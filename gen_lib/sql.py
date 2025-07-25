from pathlib import Path
import sqlite3


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
