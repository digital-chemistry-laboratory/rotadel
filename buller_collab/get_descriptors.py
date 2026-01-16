import argparse
import json
from pathlib import Path
import sys

import numpy as np

from aa_descriptors_library.query import query_closest
from aa_descriptors_library.utils import lock


class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder to handle NumPy types."""

    def default(self, obj):
        if isinstance(obj, (np.floating, np.integer)):
            return float(obj) if isinstance(obj, np.floating) else int(obj)
        return super().default(obj)


def main(pdb_file: str | Path, output_jsonl: str | Path | None = None) -> None:
    pdb_file = Path(pdb_file)
    variant_id = pdb_file.stem.split("_")[1]

    res_indices = [44, 45, 49, 50, 69, 84, 87, 127, 237]

    variant_data: dict = {}
    for res_idx in res_indices:
        charge = -1 if res_idx == 127 else 0
        result = query_closest(pdb_file, target_res_nb=res_idx - 3, charge=charge)
        variant_data[str(res_idx)] = result

    if output_jsonl is None:
        sys.stdout.write(
            json.dumps({str(variant_id): variant_data}, indent=4, cls=NumpyEncoder)
        )
    else:
        output_jsonl = Path(output_jsonl)
        with lock(output_jsonl):
            with output_jsonl.open("a") as f:
                f.write(
                    f"{json.dumps({str(variant_id): variant_data}, cls=NumpyEncoder)}\n"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query amino acid descriptors for PDB file."
    )
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument(
        "-o", "--output", help="Path to output JSONL file", default=None
    )
    args = parser.parse_args()
    main(args.pdb_file, args.output)
