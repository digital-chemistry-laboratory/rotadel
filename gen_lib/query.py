from pathlib import Path
import sqlite3
import json
from morfeus import read_xyz
from spyrmsd.rmsd import rmsd
import mdtraj as md
import numpy as np

from gen_lib.sql import get_xyz_from_sql
from gen_lib.constants import NUMBER_OF_CHI_ANGLES, THREE_TO_ONE_AA

SQL_PATH = (
    "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/gen_lib/"
    "data/merged_whole_filtered_library_sql.db"
)


def query_closest_rmsd(
    target_rotamer_xyz: str | Path,
    residue: str,
    charge: int,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> dict[str, float | str | dict | None]:
    """Query the rotamer in the library with minimum side-chain RMSD (ignoring hydrogens) to the given structure
    Args:
        target_rotamer_xyz: xyz file of the target rotamer
        residue: amino acid type of the target rotamer
        charge: charge of the target rotamer
        tautomer: tautomer of the target rotamer if applicable
        sql_path: path to the SQL database file
    Returns:
        Dictionary with ID, RMSD, and descriptors of the closest rotamer found in the library
    """
    if sql_path is None:
        sql_path = SQL_PATH

    el_target, coords_target = read_xyz(target_rotamer_xyz)
    # Remove hydrogens
    mask_noHs = el_target != "H"
    el_target_noHs = el_target[mask_noHs]
    coords_target_noHs = coords_target[mask_noHs]

    closest_rot = {"rotamer_id": None, "rmsd": float("inf"), "descriptors": None}

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT rotamer_id, descriptors FROM rotamers_data where res = ? AND charge = ? AND tautomer IS ?",
            (residue, charge, tautomer),
        )
        for rotamer_id, descriptors in cur.fetchall():
            rotamer_xyz = get_xyz_from_sql(cur, "sidechainH_xyz", rotamer_id)
            el_rotamer, coords_rotamer = (
                rotamer_xyz["elements"],
                rotamer_xyz["coordinates"],
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
                    "descriptors": json.loads(descriptors),
                }

    if closest_rot["rotamer_id"] is None:
        raise ValueError(
            "No matching rotamer found. Check the specified residue, charge, and tautomer."
        )

    return closest_rot


def query_closest_angles(
    pdb_file: str | Path,
    res_number: int,
    charge: int,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> dict[str, float | str | dict | None]:
    """Query the rotamer in the library with minimum side-chain angles distance to the given structure.
    Args:
        pdb_file: pdb file containning the target residue
        res_number: sequence number of the target residue in the pdb file
        charge: charge of the target residue
        tautomer: tautomer of the target residue if applicable
        sql_path: path to the SQL database file
    Returns:
        Dictionary with ID, angles distance, and descriptors of the closest rotamer found in the library
    """
    if sql_path is None:
        sql_path = SQL_PATH

    traj = md.load(pdb_file)
    target_name = traj.topology.residue(res_number - 1).name

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
    target_chi2 = compute_chi_from_res(traj, res_number, "2") if nb_chis >= 2 else None
    target_chi3 = compute_chi_from_res(traj, res_number, "3") if nb_chis >= 3 else None
    target_chi4 = compute_chi_from_res(traj, res_number, "4") if nb_chis >= 4 else None

    closest_rot = {
        "rotamer_id": None,
        "chis_distance": float("inf"),
        "descriptors": None,
    }

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
            chis_dist = distance_angles(
                [target_chi2, target_chi3, target_chi4], [chi2, chi3, chi4]
            )
            if chis_dist < closest_rot["chis_distance"]:
                closest_rot = {
                    "rotamer_id": rotamer_id,
                    "chis_distance": chis_dist,
                    "descriptors": json.loads(descriptors),
                }

        if closest_rot["rotamer_id"] is None:
            raise ValueError(
                "No matching rotamer found. Check the specified residue number, charge, and tautomer."
            )

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
    print(dist)
    return dist**0.5


def query_average(
    residue: str,
    charge: int | None = None,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
):
    """Calculate the weighted averaged descriptors for a given residue, charge, and tautomer.
    Args:
        residue: amino acid type
        charge: charge of the residue
        tautomer: tautomer of the residue
        sql_path: path to the SQL database file
    Returns:
        Dictionary with averaged descriptors weighted by the rotamers probability
    """
    if sql_path is None:
        sql_path = SQL_PATH

    def calc_weighted_desc(avg_desc, weight, desc):
        for key, value in desc.items():
            if isinstance(value, dict):
                sub_desc = avg_desc.setdefault(key, {})
                calc_weighted_desc(sub_desc, weight, value)
            else:
                avg_desc[key] = avg_desc.get(key, 0) + value * weight

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            """SELECT prob, descriptors FROM rotamers_data WHERE res = ?
            AND (? IS NULL OR charge = ?)
            AND (? IS NULL OR tautomer IS ?)""",
            (residue, charge, charge, tautomer, tautomer),
        )

        avg_descriptors = {}
        for prob_sidechain, desc_str in cur.fetchall():
            prob_backbone = (
                1  # replace by the P(phi,psi) once access to the NDRD library
            )
            prob_rotamer = prob_sidechain * prob_backbone
            descriptors = json.loads(desc_str)
            calc_weighted_desc(avg_descriptors, prob_rotamer, descriptors)

        return avg_descriptors
