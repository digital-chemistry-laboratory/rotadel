from pathlib import Path
import sqlite3
import json
from morfeus import read_xyz
from spyrmsd.rmsd import rmsd
import mdtraj as md
import numpy as np
import pandas as pd
from typing_extensions import deprecated

from aa_descriptors_library.sql import get_xyz_from_sql
from aa_descriptors_library.constants import (
    NDRD_PATH,
    NUMBER_OF_CHI_ANGLES,
    SQL_PATH,
    THREE_TO_ONE_AA,
)


def query_closest(
    pdb_file: str | Path,
    target_res_nb: int,
    charge: int,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> dict[str, float | str | dict | None]:
    """Query the rotamer in the library with minimum side-chain angles distance to the given structure.
    Args:
        pdb_file: pdb file containning the target residue
        target_res_nb: sequence number of the target residue in the pdb file
        charge: charge of the target residue
        tautomer: tautomer of the target residue if applicable ("D" or "E" for histidine)
        sql_path: path to the SQL database file
    Returns:
        Dictionary with ID, angles distance, side-chain chi angles,
        and descriptors of the closest rotamer found in the library
    """
    if sql_path is None:
        sql_path = SQL_PATH

    traj = md.load(pdb_file)
    target_name = traj.topology.residue(target_res_nb - 1).name

    def compute_chi_from_res(traj, res_nb, which_chi):
        # MDTraj functions compute all chi_i present in the pdb
        compute_chi_fct = {
            "2": md.compute_chi2,
            "3": md.compute_chi3,
            "4": md.compute_chi4,
        }
        all_atom_idxs, all_angles_rad = compute_chi_fct[str(which_chi)](traj)
        # Get the indices of all residues containing chi_i
        res_in_computed_chis = np.array(
            [
                traj.topology.atom(atom_idxs[1]).residue.index
                for atom_idxs in all_atom_idxs
            ]
        )
        # Find which computed angle corresponds to the target residue number
        try:
            corresponding_idx = np.where(res_in_computed_chis == res_nb - 1)[0][0]
        except IndexError:
            raise ValueError(
                f"Residue number {res_nb} (1-based) does not have a chi{which_chi} angle."
            )
        chi = np.degrees(all_angles_rad[0, corresponding_idx])
        return chi

    nb_chis = NUMBER_OF_CHI_ANGLES[THREE_TO_ONE_AA[target_name]]
    target_chi2 = (
        compute_chi_from_res(traj, target_res_nb, "2") if nb_chis >= 2 else None
    )
    target_chi3 = (
        compute_chi_from_res(traj, target_res_nb, "3") if nb_chis >= 3 else None
    )
    target_chi4 = (
        compute_chi_from_res(traj, target_res_nb, "4") if nb_chis >= 4 else None
    )

    closest_rot = {
        "rotamer_id": None,
        "chis_distance": float("inf"),
        "chis": {"chi2": None, "chi3": None, "chi4": None},
        "descriptors": None,
    }
    already_calculated = set()

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT rotamer_id, chi2, chi3, chi4, descriptors
            FROM rotamers_data where res = ? AND charge = ? AND tautomer IS ?
            """,
            (target_name, charge, tautomer),
        )
        for rotamer_id, chi2, chi3, chi4, descriptors in cur.fetchall():
            if (chi2, chi3, chi4) in already_calculated:
                continue
            already_calculated.add((chi2, chi3, chi4))
            chis_dist = distance_angles(
                [target_chi2, target_chi3, target_chi4], [chi2, chi3, chi4]
            )
            if chis_dist < closest_rot["chis_distance"]:
                closest_rot = {
                    "rotamer_id": rotamer_id,
                    "chis_distance": chis_dist,
                    "chis": {"chi2": chi2, "chi3": chi3, "chi4": chi4},
                    "descriptors": json.loads(descriptors),
                }

        if closest_rot["rotamer_id"] is None:
            raise ValueError(
                "No matching rotamer found. Check the specified residue number, charge, and tautomer."
            )
        closest_rot["descriptors"] = round_dict(closest_rot["descriptors"])

    return closest_rot


def distance_angles(angles1: list[float], angles2: list[float]) -> float:
    """Calculate the distance between two sets of angles.
    Args:
        angles1: first set of angles in degrees
        angles2: second set of angles in degrees
    Returns:
        Distance between the two sets of angles
    """
    if len(angles1) != len(angles2):
        raise ValueError("Angles lists must have the same length.")
    dist = 0.0
    for a1, a2 in zip(angles1, angles2):
        if a1 is None or a2 is None:
            if a1 is not None or a2 is not None:
                raise ValueError(
                    "Each pair of angles must either both have a value or both be None."
                )
            continue
        diff = abs(((a1 - a2 + 180.0) % 360.0) - 180.0)
        dist += diff**2
    return dist**0.5


def round_dict(dictionary: dict, decimals: int = 6) -> dict:
    """Recursively round up all float values in a nested dictionary.
    Args:
        d: dictionary to round up
        decimals: number of decimal places to round to
    Returns:
        The updated dictionary with rounded float values
    """
    for key, value in dictionary.items():
        if isinstance(value, dict):
            round_dict(value, decimals)
        elif isinstance(value, float):
            dictionary[key] = round(value, decimals)
        else:
            continue
    return dictionary


def query_average(
    residue: str,
    left_neighbour: str | None = None,
    right_neighbour: str | None = None,
    charge: int | None = None,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
    ndrd_path: str | Path | None = None,
) -> dict[str, float | dict]:
    """Calculate the weighted averaged descriptors for a given residue, charge, and tautomer.
    Args:
        residue: three letter code of the amino acid type
        left_neighbour: three letter code of the amino acid left from `residue`
        right_neighbour: three letter code of the amino acid right from `residue`
        charge: charge of the residue
        tautomer: tautomer of the residue if applicable ("D" or "E" for histidine)
        sql_path: path to the SQL database file
        ndrd_path: path to the csv file with NDRD data
    Returns:
        Dictionary with averaged descriptors weighted by the rotamers probability
    """
    if sql_path is None:
        sql_path = SQL_PATH
    if ndrd_path is None:
        ndrd_path = NDRD_PATH

    backbone_probs = parse_ndrd(ndrd_path, residue, left_neighbour, right_neighbour)

    def calc_weighted_desc(avg_desc, weight, desc):
        for key, value in desc.items():
            if isinstance(value, dict):
                sub_desc = avg_desc.setdefault(key, {})
                calc_weighted_desc(sub_desc, weight, value)
            else:
                avg_desc[key] = float(avg_desc.get(key, 0) + value * weight)

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            """SELECT prob, phi, psi, descriptors FROM rotamers_data WHERE res = ?
            AND (? IS NULL OR charge = ?)
            AND (? IS NULL OR tautomer IS ?)""",
            (residue, charge, charge, tautomer, tautomer),
        )

        avg_descriptors = {}
        for prob_sidechain, phi, psi, desc_str in cur.fetchall():

            # Fetch normalised backbone probability for the given phi & psi
            prob_backbone = backbone_probs.loc[
                (backbone_probs["phi"] == phi) & (backbone_probs["psi"] == psi),
                "norm_prob",
            ].values[0]

            prob_rotamer = prob_sidechain * prob_backbone
            descriptors = json.loads(desc_str)
            calc_weighted_desc(avg_descriptors, prob_rotamer, descriptors)

        return round_dict(avg_descriptors)


def parse_ndrd(
    ndrd_csv: str | Path,
    residue: str,
    left_neighbour: str | None = None,
    right_neighbour: str | None = None,
) -> pd.DataFrame:
    """Parse the NDRD csv file to get the backbone probabilities for a given residue with or without neighbours.
    Args:
        ndrd_csv: path to the csv file with NDRD data
        residue: three letter code of the amino acid type
        left_neighbour: three letter code of the amino acid left from `residue`
        right_neighbour: three letter code of the amino acid right from `residue`
    Returns:
        Normalised backbone probabilities for each phi,psi of the given residue
            - phi,psi incremented by 10° from -180° to 170°
            - probability either independent of neighbours or given one or both neighbours
            - treat trans and cis prolines together
    """
    ndrd = pd.read_csv(ndrd_csv)

    # Keep only the phi,psi incremented by 10° (increment of 10° in rotamer library while 5° in NDRD)
    ndrd_step_10 = ndrd.loc[(ndrd["phi"] % 10 == 0) & (ndrd["psi"] % 10 == 0)].copy()

    # Treat trans (PRO) and cis (CPR) prolines together
    residue = residue.upper()
    target_res = [residue] if residue != "PRO" else [residue, "CPR"]

    # Get backbone probabilities from the NDRD library for a given residue...

    phi_vals = ndrd_step_10["phi"].unique()
    psi_vals = ndrd_step_10["psi"].unique()
    phis = np.repeat(phi_vals, len(psi_vals))
    psis = np.tile(psi_vals, len(phi_vals))

    def select_rows(res, pos, neigh):
        sub_df = ndrd_step_10.loc[
            (ndrd_step_10["res"].isin(res))
            & (ndrd_step_10["neighbour_pos"] == pos)
            & (ndrd_step_10["neighbour_res"] == neigh)
        ].copy()

        # For proline, sum CPR and PRO probabilities
        if "CPR" in res:
            prob_pro = (
                sub_df[sub_df["res"] == "PRO"]
                .pivot_table(index="phi", columns="psi", values="prob", aggfunc="first")
                .reindex(index=phi_vals, columns=psi_vals)
                .to_numpy()
            )
            prob_cpr = (
                sub_df[sub_df["res"] == "CPR"]
                .pivot_table(index="phi", columns="psi", values="prob", aggfunc="first")
                .reindex(index=phi_vals, columns=psi_vals)
                .to_numpy()
            )
            prob_pro_cpr = prob_pro + prob_cpr
            log_pro_cpr = -np.log(prob_pro_cpr)
            sub_df = pd.DataFrame(
                {
                    "res": "PRO",
                    "neighbour_pos": pos,
                    "neighbour_res": neigh,
                    "phi": phis,
                    "psi": psis,
                    "prob": prob_pro_cpr.ravel(),
                    "log_prob": log_pro_cpr.ravel(),
                }
            )

        return sub_df

    if left_neighbour is None or right_neighbour is None:
        # ... independent of neighbours
        if left_neighbour is None and right_neighbour is None:
            position = "right"
            neighbour = "ALL"
        # ... given only right neighbour
        elif left_neighbour is None:
            position = "right"
            neighbour = right_neighbour
        # ... given only left neighbour
        else:
            position = "left"
            neighbour = left_neighbour
        backbone_probs = select_rows(target_res, position, neighbour)

        # Normalise
        backbone_probs["norm_prob"] = (
            backbone_probs["prob"] / backbone_probs["prob"].sum()
        )

    # ... given both neighbours
    else:
        # Calculate probabilities for triplet (Centre, Left, Right)
        # with log p(C,L,R) = log p(C,L) + log p(C,R) - log p(C,right=ALL)
        sub_CL = select_rows(target_res, "left", left_neighbour)
        sub_CR = select_rows(target_res, "right", right_neighbour)
        sub_CALL = select_rows(target_res, "right", "ALL")
        log_CL = (
            -sub_CL.pivot_table(
                index="phi", columns="psi", values="log_prob", aggfunc="first"
            )
            .reindex(index=phi_vals, columns=psi_vals)
            .to_numpy()
        )
        log_CR = (
            -sub_CR.pivot_table(
                index="phi", columns="psi", values="log_prob", aggfunc="first"
            )
            .reindex(index=phi_vals, columns=psi_vals)
            .to_numpy()
        )
        log_CRA = (
            -sub_CALL.pivot_table(
                index="phi", columns="psi", values="log_prob", aggfunc="first"
            )
            .reindex(index=phi_vals, columns=psi_vals)
            .to_numpy()
        )
        log_p_CLR = log_CL + log_CR - log_CRA
        p_CLR = np.exp(log_p_CLR)

        # Normalise
        p_norm = p_CLR / float(p_CLR.sum())

        # Save in dataframe
        backbone_probs = pd.DataFrame(
            {"phi": phis, "psi": psis, "norm_prob": p_norm.ravel()}
        )

    return backbone_probs


@deprecated("Use query_closest instead.")
def query_closest_rmsd(
    target_rotamer_xyz: str | Path,
    residue: str,
    charge: int,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> dict[str, float | str | dict | None]:
    """DEPRECATED: does not account for chi2 info because ignore hydrogens + takes time.
    Use `query_closest`, based on chi angles distance, instead.

    Query the rotamer in the library with minimum side-chain RMSD (ignoring hydrogens) to the given structure.
    Args:
        target_rotamer_xyz: xyz file of the target rotamer
        residue: three-letter code of the amino acid type of the target rotamer
        charge: charge of the target rotamer
        tautomer: tautomer of the target rotamer if applicable
        sql_path: path to the SQL database file
    Returns:
        Dictionary with ID, RMSD, side-chain chi angles, and descriptors of the closest rotamer found in the library
    """
    if sql_path is None:
        sql_path = SQL_PATH

    el_target, coords_target = read_xyz(target_rotamer_xyz)
    # Remove hydrogens
    mask_noHs = el_target != "H"
    el_target_noHs = el_target[mask_noHs]
    coords_target_noHs = coords_target[mask_noHs]

    closest_rot = {
        "rotamer_id": None,
        "rmsd": float("inf"),
        "chis": {"chi2": None, "chi3": None, "chi4": None},
        "descriptors": None,
    }
    already_calculated = set()

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT rotamer_id, chi2, chi3, chi4, descriptors FROM rotamers_data
            where res = ? AND charge = ? AND tautomer IS ?
            """,
            (residue, charge, tautomer),
        )
        for rotamer_id, chi2, chi3, chi4, descriptors in cur.fetchall():
            if (chi2, chi3, chi4) in already_calculated:
                continue
            already_calculated.add((chi2, chi3, chi4))
            el_rotamer, coords_rotamer = get_xyz_from_sql(
                cur, "sidechainH_xyz", rotamer_id
            )

            # Remove hydrogens
            el_rotamer_noHs = el_rotamer[mask_noHs]
            coords_rotamer_noHs = coords_rotamer[mask_noHs]

            rmsd_val = rmsd(
                coords_target_noHs,
                coords_rotamer_noHs,
                el_target_noHs,
                el_rotamer_noHs,
                center=True,
                minimize=True,
            )
            if rmsd_val < closest_rot["rmsd"]:
                closest_rot = {
                    "rotamer_id": rotamer_id,
                    "rmsd": rmsd_val,
                    "chis": {"chi2": chi2, "chi3": chi3, "chi4": chi4},
                    "descriptors": json.loads(descriptors),
                }

    if closest_rot["rotamer_id"] is None:
        raise ValueError(
            "No matching rotamer found. Check the specified residue, charge, and tautomer."
        )
    closest_rot["descriptors"] = round_dict(closest_rot["descriptors"])

    return closest_rot
