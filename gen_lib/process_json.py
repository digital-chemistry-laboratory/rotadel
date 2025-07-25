import argparse
from pathlib import Path
import json
import re
import sys
import traceback


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        "Process the output JSON files containing rotamer descriptors"
    )
    parser.add_argument(
        "json_folder",
        type=Path,
        help="folder containing the JSON files",
    )
    parser.add_argument(
        "-t",
        "--task",
        choices=["fix", "merge"],
        help="task to perform on the JSON files: 'fix' to fix incomplete JSON files, 'merge' to fix and merge them",
    )
    parser.add_argument(
        "-o",
        "--merged_output",
        type=Path,
        help="optional path for merged output JSON file",
    )
    args = parser.parse_args()

    return args.json_folder, args.task, args.merged_output


def merge_json_files(json_files: list[Path | str], merged_output: Path | str):
    """Merge given JSON files into a single JSON file."""

    merged_output = Path(merged_output)
    with open(merged_output, "w") as merged_f:
        merged_f.write("{\n")
        first = True
        for file in json_files:
            try:
                with open(file, "r") as batch_f:
                    data = json.load(batch_f)
            except json.JSONDecodeError:
                fix_json_file(file)
                with open(file, "r") as f_fixed:
                    data = json.load(f_fixed)
            for key, value in data.items():
                if not first:
                    merged_f.write(",\n")
                json.dump(key, merged_f, indent=4)
                merged_f.write(": ")
                json.dump(value, merged_f, indent=4)
                first = False
        merged_f.write("\n}\n")


def fix_json_file(json_file: Path | str) -> None:
    """Remove the last incomplete rotamer entry if JSON file is truncated."""

    json_file = Path(json_file)

    rotamer_key_pattern = r'"[A-Z]a-?\d{1,3}a-?\d{1,3}r\d{4,8}c[+-]?\d(?:t[DE])?"'
    rotamer_entry_pattern = (
        rf"({rotamer_key_pattern}:\s*{{.*?}})(?=,\s*{rotamer_key_pattern}|,\s*$|\s*}}$)"
    )

    with open(json_file, "r") as f:
        content = f.read().strip()
    try:
        json.loads(content)
        return
    except json.JSONDecodeError:
        pass

    # Find all valid rotamer entries
    valid_entries = [
        match.group(1)
        for match in re.finditer(rotamer_entry_pattern, content, re.DOTALL)
    ]
    updated_content = "{\n" + ",\n".join(valid_entries) + "\n}\n"

    all_keys = re.findall(rotamer_key_pattern, content)
    complete_keys = re.findall(rotamer_key_pattern, updated_content)
    dropped_key = next(key for key in all_keys if key not in complete_keys)
    print("Dropped incomplete rotamer:", dropped_key.strip('"'))

    # Write the updated content to a temporary file
    tmp_json = json_file.parent / f"{json_file.stem}_tmp.json"
    with open(tmp_json, "w") as f:
        f.write(updated_content)

    # Check that temporary file is valid
    try:
        with open(tmp_json, "r") as f:
            tmp_content = json.load(f)
        if len(tmp_content) != len(valid_entries):
            raise ValueError(
                f"Entry count mismatch in fixed JSON: expected {len(valid_entries)}, got {len(tmp_content)}"
            )
    except Exception:
        sys.stderr.write(f"Error fixing {json_file}: {traceback.format_exc()}\n")
        sys.exit(1)

    # Replace the original file with the fixed one
    json_file.unlink()
    tmp_json.rename(json_file)


def main(json_folder: Path | str, task: str, merged_output: Path | str | None = None):
    """Process the output JSON files containing rotamer descriptors.

    Args:
        json_folder: folder containing the output JSON files
        task: task to perform on the JSON files: 'fix' to fix incomplete JSON files, 'merge' to fix and then merge them
        merged_output: path to save the merged output JSON file
        merge: whether to merge the output JSON files into a single file.
            - False: do not merge
            - True: merge into 'merged_rotamer_descriptors.json' in the current directory
            - path given: merge into the specified file path
    """
    json_folder = Path(json_folder)
    json_files = list(json_folder.glob("*.json"))

    if task == "fix":
        for json_file in json_files:
            fix_json_file(json_file)

    elif task == "merge":
        if merged_output is None:
            merged_output = json_folder / "merged_rotamer_descriptors.json"
        merge_json_files(json_files, merged_output)


if __name__ == "__main__":
    json_folder, task, merged_output = parse_args()
    main(json_folder, task, merged_output)
