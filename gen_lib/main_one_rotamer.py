import argparse
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import nullcontext
from .utils import extract_file_from_zip
from .rotamer import Rotamer
from .optimisation import get_opt_structures


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
    args = parser.parse_args()
    return args.rotamer_id


def main(
    rotamer_id: str,
    angles_json: PathLike | str,
    xyz_zip: PathLike | str,
    run_path: PathLike | str | None = None,
):
    """Optimise rotamer geometry and calculate descriptors on it.
    Args:
        rotamer_id: ID of the rotamer
        angles_json: json file with the angles extracted from the Dunbrack library
        xyz_zip: path to the zipped xyz folder
        run_path: path of folder in which to perform the full optimisation process
    """
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
        rotamer = get_opt_structures(
            rotamer=rotamer,
            run_folder=run_folder,
            rotamer_xyz=pre_opt_xyz,
            angles_json=angles_json,
        )

    print(rotamer.to_dict())
    return rotamer


if __name__ == "__main__":
    rotamer_id = parse_args()
    main(rotamer_id)
