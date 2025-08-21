from pathlib import Path
import sqlite3
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
):
    if sql_path is None:
        sql_path = SQL_PATH

    el_target, coords_target = read_xyz(target_rotamer_xyz)
    closest = {"rotamer_id": None, "rmsd": float("inf")}

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT rotamer_id FROM rotamers_data where res = ? AND charge = ? AND tautomer IS ?",
            (residue, charge, tautomer),
        )
        for (rotamer_id,) in cur.fetchall():
            rotamer_xyz = get_xyz_from_sql(cur, "sidechainH_xyz", rotamer_id)
            el_rotamer, coords_rotamer = (
                rotamer_xyz["elements"],
                rotamer_xyz["coordinates"],
            )
            rmsd_val = rmsd(
                coords_target,
                coords_rotamer,
                el_target,
                el_rotamer,
                center=True,
                minimize=True,
            )
            if rmsd_val < closest["rmsd"]:
                closest = {"rotamer_id": rotamer_id, "rmsd": rmsd_val}

    return closest
