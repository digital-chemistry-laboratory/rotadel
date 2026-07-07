from pathlib import Path
import sqlite3
import json
import numpy as np
import pandas as pd

from rotadel.common.config import SQL_PATH


def init_sql_db(db_path: str | Path) -> None:
    """Create a SQLite database with the necessary tables for rotamer data.
    Args:
        db_path: path to the SQLite database file
    """
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS rotamers_data (
                rotamer_id TEXT PRIMARY KEY,
                res TEXT, letter TEXT, phi REAL, psi REAL,
                chi1 REAL, chi2 REAL, chi3 REAL, chi4 REAL,
                prob REAL, charge INTEGER, tautomer TEXT,
                descriptors TEXT
            )
        """
        )
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS rotamer_xyz (
                rotamer_id TEXT, atom_idx INTEGER, element TEXT,
                x REAL, y REAL, z REAL,
                PRIMARY KEY (rotamer_id, atom_idx),
                FOREIGN KEY (rotamer_id) REFERENCES rotamers_data(rotamer_id)
            )
        """
        )
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS sidechainH_xyz (
                rotamer_id TEXT, atom_idx INTEGER, element TEXT,
                x REAL, y REAL, z REAL,
                PRIMARY KEY (rotamer_id, atom_idx),
                FOREIGN KEY (rotamer_id) REFERENCES rotamers_data(rotamer_id)
            )
        """
        )
    conn.commit()


def display_rotamer_data(rotamer_id: str, sql_path: str | Path | None = None) -> None:
    """Display rotamer data from the SQL database.
    Args:
        rotamer_id: ID of the rotamer to query
        sql_path: path to the SQL database file
    """
    if sql_path is None:
        sql_path = SQL_PATH

    with sqlite3.connect(sql_path) as conn:
        cursor = conn.cursor()
    cursor.execute("SELECT * FROM rotamers_data WHERE rotamer_id = ?;", (rotamer_id,))
    columns_names = [desc[0] for desc in cursor.description]
    for row in cursor.fetchall():
        rotamer_data = dict(zip(columns_names, row))
        for col, val in rotamer_data.items():
            if col == "rotamer_id":
                print(f"{col}: {val}")
            elif col == "descriptors":
                print(f"\t{col}:")
                formatted_descriptors = json.dumps(json.loads(val), indent=4)
                for line in formatted_descriptors.splitlines():
                    print(f"\t{line}")
            else:
                print(f"\t{col}: {val}")
        print_xyz_from_sql(cursor, "rotamer_xyz", rotamer_id)
        print_xyz_from_sql(cursor, "sidechainH_xyz", rotamer_id)
        print("\n")


def print_xyz_from_sql(
    cursor: sqlite3.Cursor, table_name: str, rotamer_id: str
) -> None:
    """Print elements and xyz coordinates from table in SQL database.
    Args:
        cursor: SQLite cursor object
        table_name: name of the table containing xyz data
        rotamer_id: ID of the rotamer to query
    """
    print(f"\t{table_name}:")
    cursor.execute(
        f"SELECT element, x, y, z FROM {table_name} WHERE rotamer_id = ? ORDER BY atom_idx;",
        (rotamer_id,),
    )
    for row in cursor.fetchall():
        element, x, y, z = row
        print(f"\t\t{element:<2} ({x:.3f}, {y:.3f}, {z:.3f})")


def get_xyz_from_sql(
    table_name: str,
    rotamer_id: str,
    sql_path: str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Get XYZ coordinates of specified rotamer from SQL database.
    Args:
        table_name: name of the table containing xyz data
        rotamer_id: ID of the rotamer to query
        sql_path: path to the SQL database file
    Returns:
        Elements and xyz coordinates of the rotamer
    """
    if sql_path is None:
        sql_path = SQL_PATH

    with sqlite3.connect(sql_path) as conn:
        cursor = conn.cursor()
    cursor.execute(
        f"SELECT element, x, y, z FROM {table_name} WHERE rotamer_id = ? ORDER BY atom_idx;",
        (rotamer_id,),
    )
    rows = cursor.fetchall()
    elements = np.array([row[0] for row in rows])
    coords = np.array([[row[1], row[2], row[3]] for row in rows], dtype=float)
    return elements, coords


def get_descriptors(
    res: str | None = None,
    charge: int | None = None,
    tautomer: str | None = None,
    sql_path: str | Path | None = None,
) -> pd.DataFrame:
    """Get all descriptors from the SQL database, optionally filtered by residue, charge and tautomer.
    Args:
        residue: three letter code of the amino acid type
        charge: charge of the target residue
        tautomer: tautomer of the target residue if applicable ("D" or "E" for histidine)
        sql_path: path to the SQL database file
    Returns:
        DataFrame with rotamer IDs, residue, charge, tautomer and descriptors
    """
    if sql_path is None:
        sql_path = SQL_PATH

    query = """SELECT rotamer_id, res, charge, tautomer, descriptors FROM rotamers_data
            WHERE (? IS NULL OR res = ?)
            AND (? IS NULL OR charge = ?)
            AND (? IS NULL OR tautomer IS ?)"""
    params = (res, res, charge, charge, tautomer, tautomer)
    with sqlite3.connect(sql_path) as conn:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        base = pd.read_sql_query(query, conn, params=params)

    descriptors = [
        json.loads(d) if isinstance(d, (str, bytes)) and d else {}
        for d in base["descriptors"]
    ]
    df = pd.DataFrame(descriptors, index=base["rotamer_id"])
    df.index.name = "rotamer_id"
    df.insert(0, "residue", base["res"].values)
    df.insert(1, "charge", base["charge"].values)
    df.insert(2, "tautomer", base["tautomer"].values)

    return df


def get_descriptors_names(
    sql_path: str | Path | None = None, only_mol: bool = False
) -> list[str]:
    """Get the names of descriptors available in the SQL database.
    Args:
        sql_path: path to the SQL database file
        only_mol: if True, return only the names of molecular descriptors
    Returns:
        Names of all the descriptors or only the molecular ones
    """
    if sql_path is None:
        sql_path = SQL_PATH

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute("SELECT descriptors FROM rotamers_data LIMIT 1")
        row = cur.fetchone()
        descriptors_names = list(json.loads(row[0]).keys())

    if only_mol:
        descriptors_names = [
            descriptor
            for descriptor in descriptors_names
            if not any(
                substring in descriptor
                for substring in ("local", "fukui", "partial", "bond")
            )
        ]

    return descriptors_names


def count_rotamers_in_sql(sql_path: Path | str | None = None) -> int:
    """Count the number of rotamers in the SQL database."""
    if sql_path is None:
        sql_path = SQL_PATH

    with sqlite3.connect(sql_path) as conn:
        cur = conn.cursor()
        cur.execute("SELECT rotamer_id FROM rotamers_data;")
        rotamer_ids = [row[0] for row in cur.fetchall()]
    return len(rotamer_ids)
