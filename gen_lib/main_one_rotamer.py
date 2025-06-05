import argparse
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import nullcontext
import sys
import traceback
import json

from gen_lib.utils import lock
from gen_lib.rotamer import Rotamer


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
        "--output_json",
        type=Path,
        help="optional path for output JSON file. If a folder, 'data.json' will be created inside. \
            If not specified, './data.json' will be used.",
        default=None,
    )
    args = parser.parse_args()

    return args.rotamer_id, args.angles_json, args.output_json


def main(
    rotamer_id: str,
    angles_json: PathLike | str,
    run_path: PathLike | str | None = None,
    output_json: PathLike | str | None = None,
) -> None:
    """Optimise rotamer geometry and calculate descriptors on it.
    Args:
        rotamer_id: ID of the rotamer with species information
        angles_json: json file with the angles extracted from the Dunbrack library
        run_path: path of folder in which to perform the full optimisation process
        output_json: path for output JSON file (if not provided, prints to stdout)
    Returns:
        None, writes data to prodided json file or stdout
    """

    if output_json is not None:
        output_json = Path(output_json)
        if output_json.is_dir():
            raise IsADirectoryError(
                f"Provided path {output_json} must be a file not a directory."
            )

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

            already_calculated = False
            if output_json:
                with lock(output_json):
                    already_calculated = rotamer.load_existing_sidechain(output_json)

            if not already_calculated:
                if not run_folder.exists():
                    run_folder.mkdir(parents=True)

                # Optimise both whole rotamer and sidechain-H
                rotamer.opt_geometries()

                # Calculate descriptors on sidechain-H
                rotamer.calc_descriptors()

                # Save as json
                if output_json:
                    with lock(output_json):
                        rotamer.to_json(output_json)
                else:
                    sys.stdout.write(json.dumps(rotamer.to_dict(), indent=4))

    except Exception:
        # raise # Uncomment to see the error in the terminal as usual
        sys.stderr.write(f"*** Rotamer {rotamer_id}\n{traceback.format_exc()}\n")


if __name__ == "__main__":
    rotamer_id, angles_json, output_json = parse_args()
    main(rotamer_id, angles_json, output_json=output_json)
