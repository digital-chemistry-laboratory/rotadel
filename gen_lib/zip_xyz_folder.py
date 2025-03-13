# Script to zip (+ remove subfolder structure + rewrite comment line) the rotamer xyz files folder Nicolas created

import time
import zipfile
from pathlib import Path
import argparse


# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Zip xyz files (without subfolders) with modified comment line."
)
parser.add_argument("zip_file", type=str, help="Path and name of created zip file")
parser.add_argument(
    "parent_dir", type=Path, help="Directory containing subfolders/file.xyz"
)
args = parser.parse_args()
parent_dir = args.parent_dir
zip_file = args.zip_file

# Check correct arguments are given
if not zip_file.endswith(".zip") or not parent_dir.is_dir():
    print("Error: Invalid arguments")
    parser.print_usage()
    exit(1)

# Zip all xyz files in one zipped folder
start = time.perf_counter()
with zipfile.ZipFile(zip_file, "w", compression=zipfile.ZIP_DEFLATED) as zipf:
    zip_folder = Path(zip_file).stem
    for subfolder in parent_dir.iterdir():
        xyz_file = subfolder / f"{subfolder.name}.xyz"
        lines = xyz_file.read_text().splitlines(keepends=True)
        # Modify comment line of the xyz files
        lines[1] = (
            f"Non-optimised rotamer {subfolder.name} generated from Dunbrack library\n"
        )
        zipf.writestr(f"{zip_folder}/{subfolder.name}.xyz", "".join(lines))

end = time.perf_counter()
print(f"Execution time: {end - start:.6f} seconds")
