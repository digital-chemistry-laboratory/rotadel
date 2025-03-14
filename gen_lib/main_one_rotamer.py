import argparse
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import nullcontext

from gen_lib.utils import extract_file_from_zip, lock
from gen_lib.rotamer import Rotamer
from gen_lib.optimisation import get_opt_structures
from gen_lib.descriptors import calc_descriptors


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        "Optimise rotamer geometry and calculate descriptors on it"
    )
    parser.add_argument(
        "rotamer_id",
        type=str,
        help="ID of the rotamer",
    )
    parser.add_argument(
        "angles_json",
        type=Path,
        help="json file with the angles extracted from the Dunbrack library",
    )
    parser.add_argument(
        "xyz_zip",
        type=Path,
        help="path to the zipped xyz folder",
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

    return args.rotamer_id, args.angles_json, args.xyz_zip, args.output_json


def main(
    rotamer_id: str,
    angles_json: PathLike | str,
    xyz_zip: PathLike | str,
    run_path: PathLike | str | None = None,
    output_json: PathLike | str | None = None,
):  # TODO: add type hint for return
    """Optimise rotamer geometry and calculate descriptors on it.
    Args:
        rotamer_id: ID of the rotamer
        angles_json: json file with the angles extracted from the Dunbrack library
        xyz_zip: path to the zipped xyz folder
        run_path: path of folder in which to perform the full optimisation process
        output_json: optional path for output JSON file.
            If a folder, 'data.json' will be created inside.
            If not specified, './data.json' will be used.
    """

    if output_json is None:
        output_json = Path("./data.json")
        if output_json.exists():
            output_json.unlink()
    else:
        output_json = Path(output_json)
        if output_json.is_dir():
            output_json = output_json / "data.json"
            if output_json.exists():
                output_json.unlink()
        else:
            if output_json.exists():
                raise FileExistsError(
                    f"JSON file for output already exists. Exit to avoid overwriting.\nProvided path: {output_json}"
                )

    rotamer = Rotamer(rotamer_id)

    temp_dir_context = TemporaryDirectory() if run_path is None else nullcontext()
    with temp_dir_context as temp_dir:
        run_folder = Path(temp_dir) if run_path is None else Path(run_path)
        if not run_folder.exists():
            run_folder.mkdir(parents=True)

        # Extract non-optimised rotamer xyz file
        name_folder = Path(xyz_zip).stem
        pre_opt_xyz = run_folder / f"{rotamer_id}_pre_opt.xyz"
        extract_file_from_zip(
            zip_archive=xyz_zip,
            file_in_zip=f"{name_folder}/{rotamer_id}.xyz",
            output_path=pre_opt_xyz,
        )

        # Optimise both whole rotamer and sidechain-H
        el_whole, coord_whole, el_sidechain, coord_sidechain = get_opt_structures(
            rotamer_key=rotamer_id,
            run_folder=run_folder,
            rotamer_xyz=pre_opt_xyz,
            angles_json=angles_json,
        )
        rotamer.set_opt_geometries(el_whole, coord_whole, el_sidechain, coord_sidechain)

    # Calculate descriptors
    descriptors = calc_descriptors(el_sidechain, coord_sidechain)
    rotamer.descriptors = descriptors

    # Save as json
    with lock(output_json):
        rotamer.to_json(output_json)


if __name__ == "__main__":
    rotamer_id, angles_json, xyz_zip, output_json = parse_args()
    main(rotamer_id, angles_json, xyz_zip, output_json)
