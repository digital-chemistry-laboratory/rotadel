import contextlib
import fcntl
import json
from os import PathLike
from pathlib import Path
import shutil
import zipfile

import numpy as np


@contextlib.contextmanager
def lock(file: str | PathLike):
    """Context manager for locking a file to prevent concurrent access.
    Args:
        file: name of the file to be locked (with or without suffix)
    """
    # Lock file
    file = Path(file)
    lock_path = file.with_suffix(".lock")
    lock_file = open(lock_path, "w")
    fcntl.lockf(lock_file, fcntl.LOCK_EX)
    try:
        yield
    finally:
        # Unlock file
        fcntl.lockf(lock_file, fcntl.LOCK_UN)
        lock_file.close()


def merge_json_files(json_files: list[Path | str], merged_output: Path | str):
    """Merge given JSON files into a single JSON file."""

    merged_output = Path(merged_output)
    with open(merged_output, "w") as merged_f:
        merged_f.write("{\n")
        first = True
        for file in json_files:
            with open(file, "r") as batch_f:
                data = json.load(batch_f)
            for key, value in data.items():
                if not first:
                    merged_f.write(",\n")
                json.dump(key, merged_f, indent=4)
                merged_f.write(": ")
                json.dump(value, merged_f, indent=4)
                first = False
        merged_f.write("\n}\n")


def extract_file_from_zip(
    zip_archive: str | PathLike,
    file_in_zip: str | PathLike,
    output_path: str | PathLike | None = None,
) -> None:
    """Thread-safe extraction of a specific file from a zip archive.
    Args:
        zip_archive: path to the zip archive
        file_in_zip: path of the file whithin the zip to extract
        output_path: path of the directory or file where to save the extracted file
            - if None, extract in the current directory with same file name
            - if a directory, extract in this directory with same file name
            - if a file, extract in this exact path
    """
    zip_archive = Path(zip_archive)
    with lock(zip_archive):
        with zipfile.ZipFile(zip_archive, "r") as zip_archive:
            if file_in_zip not in zip_archive.namelist():
                raise FileNotFoundError(f"'{file_in_zip}' not found in {zip_archive}")

            if output_path is None:
                output_path = Path.cwd() / Path(file_in_zip).name
            else:
                output_path = Path(output_path)
                if output_path.is_dir():
                    output_path = output_path / Path(file_in_zip).name

            with zip_archive.open(str(file_in_zip)) as source, open(
                output_path, "wb"
            ) as target:
                shutil.copyfileobj(source, target)


def update_zip_archive(
    zip_archive: str | PathLike,
    file: str | PathLike,
    file_in_zip: str | PathLike = None,
) -> None:
    """Thread-safe update of a zip archive by adding or modifying a file.
    Args:
        zip_archive: path to the zip archive
        file: path to the file to add/update in the archive
        file_in_zip: path of the file within the zip archive (if None, use the same path as provided file)
            - if exists in archive, will be overwritten by `file`
            - if doesn't exist, `file` will be added new in the archive
    """
    zip_archive = Path(zip_archive)
    file = Path(file)
    if file_in_zip is None:
        file_in_zip = file.name

    with lock(zip_archive):
        temp_zip = zip_archive.with_suffix(".tmp.zip")
        # Copy all contents except the file to be updated
        with zipfile.ZipFile(zip_archive, "r") as src_zip, zipfile.ZipFile(
            temp_zip, "w"
        ) as dst_zip:
            for item in src_zip.infolist():
                if item.filename != file_in_zip:
                    data = src_zip.read(item.filename)
                    dst_zip.writestr(item, data)

            # Add the new file
            dst_zip.write(file, file_in_zip)

        # Replace the original zip with the updated one
        shutil.move(temp_zip, zip_archive)


def read_file_from_zip(zip_archive: str | PathLike, file_in_zip: str | PathLike) -> str:
    """Read a file from a zip archive.
    Args:
        zip_archive: path to the zip archive
        file_in_zip: path of the file within the zip archive to read
    Returns:
        content of the file
    """
    zip_archive = Path(zip_archive)
    with lock(zip_archive):
        with zipfile.ZipFile(zip_archive, "r") as zip_archive:
            with zip_archive.open(str(file_in_zip)) as file:
                return file.read().decode()


def read_xyz_zip(
    xyz_zip: str | PathLike, rotamer_id: str
) -> tuple[np.ndarray, np.ndarray]:
    """Extract geometry from zipped xyz.
    Args:
        xyz_zip: path to the zipped xyz folder
        rotamer_id: ID of the rotamer to read
    Returns:
        elements and coordinates of the rotamer xyz file
    """
    file_in_zip = Path(xyz_zip).stem / Path(f"{rotamer_id}.xyz")
    content = read_file_from_zip(xyz_zip, file_in_zip)
    lines = content.strip().split("\n")
    elements = []
    coordinates = []
    for line in lines[2:]:
        element, *xyz = line.split()
        elements.append(element)
        coordinates.append(xyz)
    return np.array(elements), np.array(coordinates)
