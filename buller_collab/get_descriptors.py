import argparse
import json
from pathlib import Path
import sys

from Bio.PDB import PDBParser
import numpy as np

from aa_descriptors_library.query import query_closest
from aa_descriptors_library.io import lock


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


def dist_atom_ligand_res_pdb(pdb_file, res_num, res_atom_name, lig_atom_name):
    """Calculate the distance between a residue atom and a ligand atom in a PDB file.
    Args:
        pdb_file: Path to the PDB file
        res_num: Number of the amino acid residue
        res_atom_name: Name of the atom in the residue
        lig_atom_name: Name of the atom in the ligand
    Returns:
        Distance between the specified residue atom and ligand atom
    """
    structure = PDBParser().get_structure("struct", pdb_file)[0]
    # Assuming residue in chain A
    res = structure["A"][(" ", res_num, " ")]
    # Assuming ligand in chain B and only one molecule
    ligand = structure["B"][(" ", 1, " ")]
    res_atom = res[res_atom_name]
    lig_atom = ligand[lig_atom_name]
    dist = lig_atom - res_atom
    return dist


def get_closest_D127_O(pdb_file):
    """Return atom index from rotamer descriptor library of the O in D127 closest to the ligand N atom.
    (6 for OD1, 7 for OD2).
    """
    dist_OD1 = dist_atom_ligand_res_pdb(pdb_file, 127, "OD1", "NAG")
    dist_OD2 = dist_atom_ligand_res_pdb(pdb_file, 127, "OD2", "NAG")
    if dist_OD1 < dist_OD2:
        return "OD1"
    else:
        return "OD2"
