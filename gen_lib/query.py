from pathlib import Path
import sqlite3
import json
from morfeus import read_xyz
from spyrmsd.rmsd import rmsd

from .sql import get_xyz_from_sql

SQL_PATH = (
    "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/gen_lib/"
    "data/merged_whole_filtered_library_sql.db"
)


def query_target(
    target_rotamer_xyz: str | Path,
    residue: str,
    charge: int,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> dict[str, float | str | dict | None]:
    """Query the rotamer in the library whose side-chain is the closest (ignoring hydrogens) to the given structure
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
