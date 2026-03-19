from multiprocessing import Pool
from functools import lru_cache
import json
from pathlib import Path
import sqlite3
from typing_extensions import deprecated

import mdtraj as md
import numpy as np
import pandas as pd
from spyrmsd.rmsd import rmsd

from aa_descriptors_library.io_atoms import read_geo
from aa_descriptors_library.sql import get_xyz_from_sql
from aa_descriptors_library.config import NDRD_PATH, SQL_PATH
from aa_descriptors_library.constants import (
    NUMBER_OF_CHI_ANGLES,
    ONE_TO_THREE_AA,
    PHYSIO_SPECIES,
    THREE_TO_ONE_AA,
)


@lru_cache()
def _read_ndrd(ndrd_csv: Path | str) -> pd.DataFrame:
    """Load and cache the data from NDRD csv file."""
    return pd.read_csv(ndrd_csv)


@lru_cache(maxsize=64)
def _get_rotamers_cached(
    sql_path: str, res: str, charge: int, tautomer: str | None
) -> tuple:
    """Fetch and cache rotamers for a given residue/charge/tautomer."""
    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            """SELECT rotamer_id, chi1, chi2, chi3, chi4, descriptors
            FROM rotamers_data WHERE res = ? AND charge = ? AND tautomer IS ?""",
            (res, charge, tautomer),
        )
        rows = cur.fetchall()

    if not rows:
        return (), np.empty((0, 4)), ()

    rotamer_ids, chis_list, descriptors_list = [], [], []
    seen: set[tuple] = set()
    for rotamer_id, chi1, chi2, chi3, chi4, desc in rows:
        # Ignore chi3 from Dunbrack library for proline as normally defined with 2 chis
        if res == "PRO":
            chi3 = None
        key = (chi1, chi2, chi3, chi4)
        if key in seen:
            continue
        seen.add(key)
        rotamer_ids.append(rotamer_id)
        chis_list.append(
            [
                chi1 if chi1 is not None else float("nan"),
                chi2 if chi2 is not None else float("nan"),
                chi3 if chi3 is not None else float("nan"),
                chi4 if chi4 is not None else float("nan"),
            ]
        )
        descriptors_list.append(desc)

    return tuple(rotamer_ids), np.array(chis_list, dtype=float), tuple(descriptors_list)


def _compute_chis(
    traj: "md.Trajectory", pdb_res_num: int, nb_chis: int
) -> list[float | None]:
    """Compute the chi angles for a residue in a trajectory."""

    # MDTraj functions compute all chi_i present in the PDB
    fcts = [md.compute_chi1, md.compute_chi2, md.compute_chi3, md.compute_chi4]
    chis: list[float | None] = []
    for i in range(nb_chis):
        all_atom_idxs, all_angles_rad = fcts[i](traj)
        # Get the PDB residue sequence numbers of all residues that have a chi_i angle
        all_pdb_res_nums = np.array(
            [traj.topology.atom(a[1]).residue.resSeq for a in all_atom_idxs]
        )
        # Find which computed angle corresponds to the target residue
        match = np.where(all_pdb_res_nums == pdb_res_num)[0]
        if match.size == 0:
            raise ValueError(
                f"Residue with PDB residue sequence number {pdb_res_num} does not have a chi{i + 1} angle."
            )
        chis.append(float(np.degrees(all_angles_rad[0, match[0]])))

    return chis + [None] * (4 - nb_chis)


def query_closest(
    pdb_file: str | Path,
    pdb_res_num: int | None = None,
    res_position: int | None = None,
    start_pdb_res_num: int = 1,
    charge: int,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> dict[str, float | str | dict | None]:
    """Query the rotamer in the library with minimum side-chain angles distance to the given structure.
    Give target residue either by its PDB residue sequence number (`pdb_res_num`) or by its position in the sequence
    (`res_position`, with `start_pdb_res_num` if the sequence does not start at 1 in the PDB file).
    Args:
        pdb_file: PDB file containning the target residue
        pdb_res_num: residue sequence number of the target residue in the PDB file
        res_position: position of the target residue in the amino acid sequence
        start_pdb_res_num: residue sequence number of the first residue of the target chain in the PDB file
            Only used when res_position is given
        charge: charge of the target residue
        tautomer: tautomer of the target residue if applicable ("D" or "E" for histidine)
        sql_path: path to the SQL database file
    Returns:
        Dictionary with ID, angles distance, side-chain chi angles,
        and descriptors of the closest rotamer found in the library
    """
    if sql_path is None:
        sql_path = SQL_PATH

    if (res_position is None and pdb_res_num is None) or (
        res_position is not None and pdb_res_num is not None
    ):
        raise ValueError(
            "Residue must be given by exactly one of: res_position or pdb_res_num."
        )
    if res_position is not None:
        pdb_res_num = start_pdb_res_num + res_position - 1

    traj = md.load(pdb_file)
    target_res = next(
        (r for r in traj.topology.residues if r.resSeq == pdb_res_num), None
    )
    if target_res is None:
        raise ValueError(
            f"Residue with PDB residue sequence number {pdb_res_num} not found in {pdb_file}."
        )
    target_name = target_res.name

    # Only one possibility for Ala and Gly as they do not have chi angles
    if target_name in ["ALA", "GLY"]:
        rotamer_id = "Aa0a0r0000c0" if target_name == "ALA" else "Ga0a0r0000c0"
        with sqlite3.connect(sql_path) as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT descriptors FROM rotamers_data WHERE rotamer_id = ?;",
                (rotamer_id,),
            )
            descriptors = cur.fetchone()[0]
        return {
            "rotamer_id": rotamer_id,
            "chis_distance": 0.0,
            "chis": {"chi1": None, "chi2": None, "chi3": None, "chi4": None},
            "descriptors": json.loads(descriptors),
        }

    nb_chis = NUMBER_OF_CHI_ANGLES[THREE_TO_ONE_AA[target_name]]
    target_chis = _compute_chis(traj, pdb_res_num, nb_chis)

    rotamer_ids, chi_array, descriptors_list = _get_rotamers_cached(
        str(sql_path), target_name, charge, tautomer
    )
    if not rotamer_ids:
        raise ValueError(
            "No matching rotamer found. Check the specified residue number, charge, and tautomer."
        )

    # Vectorized angular distance over all rotamers at once (NaN positions = absent chi, contributes 0)
    target_arr = np.array(
        [t if t is not None else float("nan") for t in target_chis], dtype=float
    )
    diff = np.abs(((chi_array - target_arr + 180.0) % 360.0) - 180.0)
    dists = np.sqrt(np.nansum(diff**2, axis=1))

    best_idx = int(np.argmin(dists))
    best_chis = chi_array[best_idx]

    return {
        "rotamer_id": rotamer_ids[best_idx],
        "chis_distance": float(dists[best_idx]),
        "chis": {
            "chi1": None if np.isnan(best_chis[0]) else float(best_chis[0]),
            "chi2": None if np.isnan(best_chis[1]) else float(best_chis[1]),
            "chi3": None if np.isnan(best_chis[2]) else float(best_chis[2]),
            "chi4": None if np.isnan(best_chis[3]) else float(best_chis[3]),
        },
        "descriptors": round_dict(json.loads(descriptors_list[best_idx])),
    }


def query_closest_batch(
    queries: list[dict],
    sql_path: str | Path | None = None,
    num_workers: int = 1,
) -> list[dict[str, float | str | dict | None]]:
    """Run query_closest for a list of PDB files, optionally in parallel.
    Args:
        queries: list of dicts, each containing kwargs for query_closest (except sql_path).
            Each dict must include "pdb_file" and either "pdb_res_num" or "res_position"
            (+ optional "start_pdb_res_num"), and "charge" and "tautomer"
        sql_path: path to the SQL database file
        num_workers: number of parallel worker processes
    Returns:
        List of query results in the same order as the given queries
    """
    if sql_path is None:
        sql_path = SQL_PATH
    sql_path_str = str(sql_path)

    if num_workers == 1:
        return [query_closest(**{**q, "sql_path": sql_path_str}) for q in queries]

    with Pool(processes=num_workers) as pool:
        jobs = [
            pool.apply_async(query_closest, kwds={**q, "sql_path": sql_path_str})
            for q in queries
        ]
        return [job.get() for job in jobs]


def results_closest_into_dataframes(
    queries_output: list[dict[str, float | str | dict | None]],
    pdb_ids: list[int | str],
    res_positions: list[int | str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert query_closest_batch outputs into two DataFrames.

    Args:
        queries_output: flat list of query_closest results, ordered as
            [row0_res0, row0_res1, ..., row1_res0, row1_res1, ...]
            i.e. the direct output of query_closest_batch when queries are built
            by iterating rows then residues.
        pdb_ids: one label per row (e.g. PDB IDs)
        res_positions: one label per residue position (e.g. residue position numbers)
    Returns:
        rotamers_df: indexed by pdb_ids, columns res{label}_rotamer,
            res{label}_chis_dist, res{label}_chi1..chi4 for each residue
        descriptors_df: indexed by pdb_ids, columns res{label}_{desc} for molecular
            descriptors and res{label}_{desc}_{atom} for atomic (nested) descriptors.
            If residues are not the same accross PDBs, the union of atomic descriptors
            will be saved with NaN for atoms missing in some PDBs.
    """
    n_res = len(res_positions)
    rotamer_rows = []
    descriptors_rows = []

    for pdb_id in range(len(pdb_ids)):
        rotamer_row: dict = {}
        descriptor_row: dict = {}
        for i, res_position in enumerate(res_positions):
            result = queries_output[pdb_id * n_res + i]
            prefix = f"res{res_position}"

            rotamer_row[f"{prefix}_rotamer"] = result["rotamer_id"]
            rotamer_row[f"{prefix}_chis_dist"] = result["chis_distance"]
            for chi_name, chi_val in result["chis"].items():
                rotamer_row[f"{prefix}_{chi_name}"] = chi_val

            for desc_key, desc_val in result["descriptors"].items():
                if isinstance(desc_val, dict):
                    # Atomic descriptor: flatten one column per atom
                    for atom, atom_val in desc_val.items():
                        descriptor_row[f"{prefix}_{desc_key}_{atom}"] = atom_val
                else:
                    descriptor_row[f"{prefix}_{desc_key}"] = desc_val

        rotamer_rows.append(rotamer_row)
        descriptors_rows.append(descriptor_row)

    rotamers_df = pd.DataFrame(rotamer_rows, index=pdb_ids)
    descriptors_df = pd.DataFrame(descriptors_rows, index=pdb_ids)
    rotamers_df.index.name = "variant_id"
    descriptors_df.index.name = "variant_id"

    # Group columns by residue without changing within-residue order
    res_cols: dict[int, list[str]] = {}
    seen: set[str] = set()
    for row in descriptors_rows:
        for col in row:
            if col not in seen:
                seen.add(col)
                res_position = int(col.split("_", 1)[0][3:])
                res_cols.setdefault(res_position, []).append(col)
    ordered_cols = [
        col for res_position in sorted(res_cols) for col in res_cols[res_position]
    ]
    descriptors_df = descriptors_df[ordered_cols]

    return rotamers_df, descriptors_df


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
        residue: one or three letter(s) code of the amino acid type
        left_neighbour: one or three letter(s) code of the amino acid left from `residue`
        right_neighbour: one or three letter(s) code of the amino acid right from `residue`
        charge: charge of the residue
        tautomer: tautomer of the residue if applicable ("D" or "E" for histidine)
        sql_path: path to the SQL database file
        ndrd_path: path to the csv file with NDRD data
    Returns:
        Dictionary with averaged descriptors weighted by the rotamers probability
    """
    if len(residue) == 1:
        residue = ONE_TO_THREE_AA[residue.upper()]
    if left_neighbour is not None and len(left_neighbour) == 1:
        left_neighbour = ONE_TO_THREE_AA[left_neighbour.upper()]
    if right_neighbour is not None and len(right_neighbour) == 1:
        right_neighbour = ONE_TO_THREE_AA[right_neighbour.upper()]

    if sql_path is None:
        sql_path = SQL_PATH
    if ndrd_path is None:
        ndrd_path = NDRD_PATH

    if charge is None:
        charge = PHYSIO_SPECIES[THREE_TO_ONE_AA[residue]]["charge"]
    if residue == "HIS":
        if tautomer is None:
            tautomer = PHYSIO_SPECIES[THREE_TO_ONE_AA[residue]]["tautomer"]
        if tautomer is not None and charge != 0:
            raise ValueError(
                "Histidine tautomers only available for charge 0, "
                f"but charge {charge} and tautomer {tautomer} were specified."
            )
    else:
        if tautomer is not None:
            raise ValueError("Tautomers can only be specified for histidine.")

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
            """SELECT prob, phi, psi, descriptors FROM rotamers_data
            WHERE res = ? AND charge = ? AND tautomer IS ?""",
            (residue, charge, tautomer),
        )

        avg_descriptors = {}
        all_rotamer_probs = []
        for prob_sidechain, phi, psi, desc_str in cur.fetchall():

            # Fetch normalised backbone probability for the given phi & psi
            if residue in ["ALA", "GLY"]:
                prob_backbone = 1.0
            else:
                prob_backbone = backbone_probs.loc[
                    (backbone_probs["phi"] == phi) & (backbone_probs["psi"] == psi),
                    "norm_prob",
                ].values[0]

            prob_rotamer = prob_sidechain * prob_backbone
            all_rotamer_probs.append(prob_rotamer)
            descriptors = json.loads(desc_str)
            calc_weighted_desc(avg_descriptors, prob_rotamer, descriptors)

        prob_sum = sum(all_rotamer_probs)
        if not np.isclose(prob_sum, 1.0, atol=0.01):
            raise ValueError(
                f"Sum of all rotamer probabilities is {prob_sum}, expected approximately 1.0"
            )

        return round_dict(avg_descriptors)


def parse_ndrd(
    ndrd_csv: str | Path,
    residue: str,
    left_neighbour: str | None = None,
    right_neighbour: str | None = None,
) -> pd.DataFrame:
    """Parse the NDRD csv file to get the backbone probabilities for a given residue with or without neighbours.
    Args:
        ndrd_csv: path to the csv file with NDRD data (with phi/psi incremented by 10°)
        residue: three letter code of the amino acid type
        left_neighbour: three letter code of the amino acid left from `residue`
        right_neighbour: three letter code of the amino acid right from `residue`
    Returns:
        Normalised backbone probabilities for each phi,psi of the given residue
            - phi,psi incremented by 10° from -180° to 170°
            - probability either independent of neighbours or given one or both neighbours
            - treat trans and cis prolines together
    """
    ndrd_step_10 = _read_ndrd(ndrd_csv)

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

    el_target, coords_target = read_geo(target_rotamer_xyz)
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
