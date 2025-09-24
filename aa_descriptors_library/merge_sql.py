import argparse
import sqlite3
from typing import Iterable
from pathlib import Path
from shutil import copyfile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge multiple SQLite databases into a single one."
    )
    parser.add_argument(
        "db_folder",
        type=Path,
        help="Path to the folder containing the SQLite .db files to merge.",
    )
    parser.add_argument(
        "merged_db_path",
        type=str,
        default=None,
        help="Path for the created merged .db file.",
    )
    args = parser.parse_args()

    return args.db_folder, args.merged_db_path


def main(db_folder: Iterable[Path], merged_db_path: Path) -> None:
    """Merge multiple SQLite databases into a single one."""
    merged_db_path = Path(merged_db_path)
    if merged_db_path.exists():
        raise FileExistsError(f"{merged_db_path} already exists.")

    db_files = sorted(Path(db_folder).glob("*.db"))

    # Create merged DB from first file for the correct schema
    copyfile(db_files[0], merged_db_path)

    # Merge all other dbs
    with sqlite3.connect(merged_db_path) as merged_conn:
        for db_file in db_files[1:]:
            with sqlite3.connect(db_file) as conn:
                for table in ["rotamers_data", "rotamer_xyz", "sidechainH_xyz"]:
                    rows = conn.execute(f"SELECT * FROM {table}").fetchall()
                    placeholders = ",".join(["?"] * len(rows[0]))
                    # try:
                    merged_conn.executemany(
                        f"INSERT INTO {table} VALUES ({placeholders})", rows
                    )
        merged_conn.commit()


if __name__ == "__main__":
    db_folder, merged_db_path = parse_args()
    main(db_folder, merged_db_path)
