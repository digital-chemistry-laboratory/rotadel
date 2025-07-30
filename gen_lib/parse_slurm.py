import argparse
import re
from pathlib import Path
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse SLURM output files for information about failed runs."
    )
    parser.add_argument(
        "slurm_folder",
        type=Path,
        help="Path to the folder containing the SLURM output files.",
    )
    parser.add_argument(
        "-n",
        "--run_name",
        type=str,
        default=None,
        help="Name of the run to be used for the output files.",
    )
    args = parser.parse_args()

    return args.slurm_folder, args.run_name


def main(slurm_folder, run_name=None):
    # Error patterns to check
    error_patterns = {
        "DetermineConnectivity": re.compile(
            r"rdkit\.Chem\.rdDetermineBonds\.DetermineConnectivity\(NoneType\)"
        ),
        "xtbopt_NotFound": re.compile(r"FileNotFoundError"),
        ">_None": re.compile(
            r"TypeError: '>' not supported between instances of 'NoneType' and 'int'"
        ),
        "xtb_failed": re.compile(r"RuntimeError: xtb calculation failed"),
        "ConeAngle": re.compile(r"ConeAngle"),
    }

    # Parse slurm files
    parsed_info = {}
    failed_ids = []
    for slurm_file in slurm_folder.glob("slurm-*"):
        errors = Counter()
        with open(slurm_file, "r") as f:
            lines = f.readlines()
            first_line = lines[0]
            last_line = lines[-1]
            content = "".join(lines)

            # Extract batch number
            match = re.search(r"batch_\d+", first_line)
            batch = match.group(0) if match else "batch_???"

            # Extract run time
            if "DUE TO TIME LIMIT" in last_line:
                runtime = "timeout"
            else:
                match_time = re.search(r"Total execution time:\s+(.*)", last_line)
                runtime = match_time.group(1).strip() if match_time else "unknown"

        error_blocks = re.split(r"\*\*\* Rotamer ", content)[1:]
        for block in error_blocks:
            rotamer_id = block.splitlines()[0]
            failed_ids.append(rotamer_id)
            matched = False
            for name, pattern in error_patterns.items():
                if pattern.search(block):
                    errors[name] += 1
                    matched = True
                    break
            if not matched:
                errors["Other"] += 1
        parsed_info[slurm_file.name] = {
            "batch": batch,
            "error_counts": errors,
            "runtime": runtime,
        }

    error_types = [
        "DetermineConnectivity",
        "xtbopt_NotFound",
        ">_None",
        "xtb_failed",
        "ConeAngle",
        "Other",
    ]

    # Format results
    field_widths = {key: len(f"'{key}': ") + 3 for key in error_types}
    errors_record = []
    nb_tot_errors = 0
    for file, info in sorted(parsed_info.items()):
        batch = info["batch"]
        errors = info["error_counts"]
        runtime = info["runtime"]
        line = f"{file:<20} {batch:<15}"
        for error_type in error_types:
            count = errors.get(error_type, 0)
            nb_tot_errors += count
            count_str = f"{count}" if count else ""
            field = f"{error_type}: {count_str}".ljust(field_widths[error_type])
            line += field
        line += f" run: {runtime}"
        errors_record.append(line.rstrip())

    # Output error counts per file
    errors_output = "slurm_errors.out"
    if run_name:
        errors_output = f"{run_name}_{errors_output}"
    with open(errors_output, "w") as f:
        f.write(f"Total number of errors: {nb_tot_errors}\n")
        f.write("\n".join(errors_record))

    # Output IDs of failed rotamers
    failed_id_output = "failed_keys.csv"
    if run_name:
        failed_id_output = f"{run_name}_{failed_id_output}"
    with open(failed_id_output, "w") as f:
        f.write("\n".join(failed_ids))


if __name__ == "__main__":
    slurm_folder, run_name = parse_args()
    main(slurm_folder, run_name)
