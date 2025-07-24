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
        help="optional path for output JSON file. If folder, individual rotamer_id.json will be created inside. \
            If not specified, output will be printed to stdout.",
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
        output_json: path for output JSON file
            If json file, write rotamer output to it.
            If folder, create individual rotamer_id.json inside.
            If None, print output to stdout.
    Returns:
        None, writes data to json file or stdout
    """

    if output_json is not None:
        output_json = Path(output_json)

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

            rotamer_already_exists = False
            sidechain_already_calculated = False
            if output_json and output_json.exists() and output_json.is_file():
                with lock(output_json):
                    with open(output_json, "r") as f:
                        content = json.load(f)
                    rotamer_already_exists = rotamer_id in content

                    if not rotamer_already_exists:
                        sidechain_already_calculated = rotamer.load_existing_sidechain(
                            output_json
                        )

            if not rotamer_already_exists:
                if not sidechain_already_calculated:
                    if not run_folder.exists():
                        run_folder.mkdir(parents=True)

                    # Optimise both whole rotamer and sidechain-H
                    rotamer.opt_geometries()

                    # Calculate descriptors on sidechain-H
                    rotamer.calc_descriptors()

                # Save as json or print out
                if output_json is None:
                    sys.stdout.write(json.dumps(rotamer.to_dict(), indent=4))
                elif output_json.suffix == ".json":
                    with lock(output_json):
                        rotamer.to_json(output_json)
                elif output_json.suffix == "":
                    if not output_json.exists():
                        output_json.mkdir(parents=True)
                    individual_output_json = output_json / f"{rotamer_id}.json"
                    if individual_output_json.exists():
                        individual_output_json.unlink()
                    rotamer.to_json(individual_output_json)
                else:
                    raise ValueError(
                        f"The given argument `output_json` must be a file with .json suffix or be a folder. "
                        f"You gave: {output_json}"
                    )

    except Exception:
        # raise # Uncomment to see the error in the terminal as usual
        sys.stderr.write(f"*** Rotamer {rotamer_id}\n{traceback.format_exc()}\n")
        sys.exit(1)


if __name__ == "__main__":
    rotamer_id, angles_json, output_json = parse_args()
    main(rotamer_id, angles_json, output_json=output_json)
