import argparse
from pathlib import Path
from shutil import copyfile
import sqlite3
import time
from typing import Iterable

from aa_descriptors_library.config import SQL_PATH
from aa_descriptors_library.sql import count_rotamers_in_sql


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
        "-o",
        "--output_merged_db_path",
        type=str,
        default=None,
        help="Path for the created merged .db file.",
    )
    args = parser.parse_args()

    return args.db_folder, args.output_merged_db_path


def main(db_folder: Iterable[Path], merged_db_path: Path | str | None = None) -> None:
    """Merge multiple SQLite databases into a single one."""
    start_time = time.perf_counter()

    if merged_db_path is None:
        merged_db_path = SQL_PATH
    merged_db_path = Path(merged_db_path)
    if merged_db_path.exists():
        raise FileExistsError(f"{merged_db_path} already exists.")

    db_files = sorted(Path(db_folder).glob("*.db"))

    # Create merged db from first file for the correct schema
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

    # Print execution time and count in merged db
    elapsed_seconds = time.perf_counter() - start_time
    elapsed_minutes, remaining_seconds = divmod(elapsed_seconds, 60)
    print("Executed in: " f"{int(elapsed_minutes)} min {remaining_seconds:05.2f} sec")
    rotamer_count = count_rotamers_in_sql(merged_db_path)
    print(f"Number of rotamers in merged database: {rotamer_count}")


if __name__ == "__main__":
    db_folder, merged_db_path = parse_args()
    main(db_folder, merged_db_path)
