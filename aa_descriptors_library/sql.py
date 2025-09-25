from pathlib import Path
import sqlite3
import json
import numpy as np


def init_sql_db(db_path: str | Path) -> None:
    """Create a SQLite database with the necessary tables for rotamer data.
    Args:
        db_path: path to the SQLite database file
    """
    db_path = Path(db_path)
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


def display_rotamer_data(cursor: sqlite3.Cursor, rotamer_id: str) -> None:
    """Display rotamer data from the SQL database.
    Args:
        cursor: SQLite cursor object
        rotamer_id: ID of the rotamer to query
    """
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
    cursor: sqlite3.Cursor, table_name: str, rotamer_id: str
) -> tuple[np.ndarray, np.ndarray]:
    """Get XYZ coordinates of specified rotamer from SQL database.
    Args:
        cursor: SQLite cursor object
        table_name: name of the table containing xyz data
        rotamer_id: ID of the rotamer to query
    Returns:
        Elements and xyz coordinates of the rotamer
    """
    cursor.execute(
        f"SELECT element, x, y, z FROM {table_name} WHERE rotamer_id = ? ORDER BY atom_idx;",
        (rotamer_id,),
    )
    rows = cursor.fetchall()
    elements = np.array([row[0] for row in rows])
    coords = np.array([[row[1], row[2], row[3]] for row in rows], dtype=float)
    return elements, coords
