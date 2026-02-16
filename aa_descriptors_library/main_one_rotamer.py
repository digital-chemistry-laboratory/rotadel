import argparse
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import nullcontext
import sys
import traceback
import json
import sqlite3

from aa_descriptors_library.io import lock
from aa_descriptors_library.rotamer import Rotamer
from aa_descriptors_library.sql import init_sql_db


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        "Optimise rotamer geometry and calculate descriptors on it"
    )
    parser.add_argument(
        "rotamer_id",
        type=str,
        help="ID of the rotamer with species information",
    )
    parser.add_argument(
        "angles_json",
        type=Path,
        help="json file with the angles extracted from the Dunbrack library",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=Path,
        help="optional path for output database file (.db or .json) or folder. \
            If folder, individual rotamer_id.json will be created inside. \
            If not specified, output will be printed to stdout.",
        default=None,
    )
    parser.add_argument(
        "-r",
        "--rerun_failed",
        action="store_true",
        default=False,
        help="rerun geometry optimisations on failure",
    )
    args = parser.parse_args()

    return (
        args.rotamer_id,
        args.angles_json,
        args.output_path,
        args.rerun_failed,
    )


def main(
    rotamer_id: str,
    angles_json: PathLike | str,
    run_path: PathLike | str | None = None,
    output_path: PathLike | str | None = None,
    rerun_failed: bool = False,
    sys_exit: bool = False,
) -> None:
    """Optimise rotamer geometry and calculate descriptors on it.
    Args:
        rotamer_id: ID of the rotamer with species information
        angles_json: json file with the angles extracted from the Dunbrack library
        run_path: path of folder in which to perform the full optimisation process
        output_path: path for output database file or folder
            If .db file, write rotamer output to SQL database
            If .json file, write rotamer output to it.
            If folder, create individual rotamer_id.json inside.
            If None, print output to stdout.
        rerun_failed: whether to rerun geometry optimisations on failure
        sys_exit: whether to exit with sys.exit(1) (True) or raise exception (False) when error
    Returns:
        None, writes data to database or stdout
    """

    if output_path is not None:
        output_path = Path(output_path)
        if output_path.suffix == ".db" and not output_path.exists():
            init_sql_db(output_path)

    with open(angles_json, "r") as f:
        dunbrack_all_rotamers = json.load(f)
    rotamer_id_dunbrack = rotamer_id.split("c")[0]
    species_info = rotamer_id.split("c")[-1]
    if "t" in species_info:
        charge = int(species_info.split("t")[0])
        tautomer = species_info.split("t")[-1]
    else:
        charge = int(species_info)
        tautomer = None
    dunbrack_data = dunbrack_all_rotamers[rotamer_id_dunbrack]

    try:
        temp_dir_context = TemporaryDirectory() if run_path is None else nullcontext()
        with temp_dir_context as temp_dir:
            run_folder = Path(temp_dir) if run_path is None else Path(run_path)

            rotamer = Rotamer(
                rotamer_id, dunbrack_data, charge, tautomer, run_path=run_folder
            )

            # Check if the rotamer ID or sidechain already exists in the database
            rotamer_already_exists = False
            sidechain_already_calculated = False
            if output_path and output_path.exists() and output_path.is_file():
                if output_path.suffix == ".db":
                    with lock(output_path):
                        with sqlite3.connect(output_path) as conn:
                            rotamer_already_exists = (
                                conn.execute(
                                    "SELECT 1 FROM rotamers_data WHERE rotamer_id = ?",
                                    (rotamer_id,),
                                ).fetchone()
                                is not None
                            )

                            if not rotamer_already_exists:
                                sidechain_already_calculated = (
                                    rotamer.load_existing_sidechain_sql(conn)
                                )

                elif output_path.suffix == ".json":
                    with lock(output_path):
                        with open(output_path, "r") as f:
                            content = json.load(f)
                        rotamer_already_exists = rotamer_id in content

                    if not rotamer_already_exists:
                        sidechain_already_calculated = (
                            rotamer.load_existing_sidechain_json(output_path)
                        )

            if not rotamer_already_exists and not sidechain_already_calculated:
                if not run_folder.exists():
                    run_folder.mkdir(parents=True)

                # Optimise both whole rotamer and sidechain-H
                rotamer.opt_geometries(rerun_failed=rerun_failed)

                # Calculate descriptors on sidechain-H
                rotamer.calc_descriptors()

            # Save in database or print out
            if not rotamer_already_exists:
                if output_path is None:
                    sys.stdout.write(json.dumps(rotamer.to_dict(), indent=4))
                elif output_path.suffix == ".db":
                    with lock(output_path):
                        with sqlite3.connect(output_path) as conn:
                            rotamer.to_sql(conn)
                elif output_path.suffix == ".json":
                    with lock(output_path):
                        rotamer.to_json(output_path)
                elif output_path.suffix == "":
                    if not output_path.exists():
                        output_path.mkdir(parents=True)
                    individual_output_json = output_path / f"{rotamer_id}.json"
                    if individual_output_json.exists():
                        individual_output_json.unlink()
                    rotamer.to_json(individual_output_json)
                else:
                    raise ValueError(
                        f"The given argument `output_path` must be either: a .db file, a .json file, or be a folder. "
                        f"You gave: {output_path}"
                    )

    except Exception:
        sys.stderr.write(f"*** Rotamer {rotamer_id}\n{traceback.format_exc()}\n")
        if sys_exit:
            sys.exit(1)
        else:
            raise


if __name__ == "__main__":
    rotamer_id, angles_json, output_path, rerun_failed = parse_args()
    main(
        rotamer_id,
        angles_json,
        output_path=output_path,
        rerun_failed=rerun_failed,
        sys_exit=True,
    )
