import argparse
import pandas as pd
from pathlib import Path

from aa_descriptors_library.query import query_closest
from aa_descriptors_library.utils import lock


def flatten_dict(res_idx: int, data_dict: dict) -> dict:
    """Flatten nested dictionary."""
    out = {}
    for k, v in data_dict.items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                kk_format = str(kk).replace(", ", "-")
                out[f"{res_idx}_{k}_{kk_format}"] = vv
        else:
            out[f"{res_idx}_{k}"] = v
    return out


def clean_csv(csv_file: str) -> None:
    """Dropping NaN columns and sorting by variant_id"""
    df = pd.read_csv(csv_file, index_col=["variant_id"])
    df = df.dropna(axis="columns", how="all")
    df = df.sort_index(
        key=lambda x: x.map(lambda i: (0, int(i)) if str(i).isdigit() else (1, str(i)))
    )
    df.to_csv(csv_file)


def main(pdb_file: str | Path):
    pdb_file = Path(pdb_file)
    res_indices = [44, 45, 49, 50, 69, 84, 87, 127, 237]
    chi_csv = Path("rotamer_chis.csv")
    desc_csv = Path("descriptors.csv")

    variant_id = pdb_file.stem.split("_")[1]
    # Check if already processed
    if chi_csv.exists():
        done = set(
            pd.read_csv(chi_csv, usecols=["variant_id"])["variant_id"].astype(str)
        )
        if variant_id in done:
            return

    chis_row = {"variant_id": variant_id}
    desc_row = {"variant_id": variant_id}

    for res_idx in res_indices:
        charge = -1 if res_idx == 127 else 0
        result = query_closest(pdb_file, target_res_nb=res_idx - 3, charge=charge)

        chis_row[f"{res_idx}_rotamer_id"] = result["rotamer_id"]
        chis_row.update(flatten_dict(res_idx, result["chis"]))
        desc_row.update(flatten_dict(res_idx, result["descriptors"]))

    df_chi_row = pd.DataFrame([chis_row])
    df_desc_row = pd.DataFrame([desc_row])

    with lock(chi_csv):
        df_chi_row.to_csv(chi_csv, mode="a", header=not chi_csv.exists(), index=False)
    with lock(desc_csv):
        df_desc_row.to_csv(
            desc_csv, mode="a", header=not desc_csv.exists(), index=False
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query amino acid descriptors for PDB file."
    )
    parser.add_argument("pdb_file", help="Path to the PDB file")
    args = parser.parse_args()
    main(args.pdb_file)
