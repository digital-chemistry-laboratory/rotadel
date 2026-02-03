from os import PathLike
from pathlib import Path

from Bio.PDB import PDBParser
from morfeus import read_xyz
from morfeus.typing import Array1DStr, Array2DFloat
import numpy as np
from openbabel import pybel


def read_geo(geo_file: PathLike | str) -> tuple[Array1DStr, Array2DFloat]:
    """Read elements and coordinates [Å] from a molecular structure file."""
    geo_file = Path(geo_file)
    if geo_file.suffix == ".xyz":
        elements, coordinates = read_xyz(geo_file)
    elif geo_file.suffix == ".pdb":
        structure = PDBParser(QUIET=True).get_structure("x", geo_file)
        atoms = list(structure.get_atoms())
        elements = np.array([a.element.strip() for a in atoms])
        coordinates = np.array([a.coord for a in atoms], dtype=np.float64)
    else:
        raise ValueError(f"Unsupported file format: {geo_file.suffix}")

    return elements, coordinates


def convert_file(
    input_file: PathLike | str,
    output_file: PathLike | str | None = None,
    output_format: str | None = None,
    delete_input: bool = True,
) -> None:
    """Convert a molecular structure file from one format to another."""
    pybel.ob.obErrorLog.SetOutputLevel(
        0
    )  # suppress the "Failed to kekulize aromatic bonds" warning
    structure = next(pybel.readfile(Path(input_file).suffix[1:], str(input_file)))
    if output_file is None and output_format is not None:
        output_file = Path(input_file).with_suffix(f".{output_format}")
    elif output_file is not None and output_format is None:
        output_format = Path(output_file).suffix[1:]
    else:
        raise ValueError("Either output_file or output_format must be provided.")
    structure.write(output_format, str(output_file), overwrite=True)
    if delete_input:
        Path(input_file).unlink()


def map_atom_indices(
    pdb_file: PathLike | str, zero_indexed=False, by_serial=False
) -> dict[str, int]:
    """Map atom names to their indices from a PDB file.
    Args:
        pdb_file: path to the PDB file
        zero_indexed: whether to use 0-based indexing (default: 1-based indexing)
        by_serial: whether to return the PDB atom serial numbers (default: indices following line order)
    Returns:
        Dictionary mapping atom names to their indices"""
    mapping = {}
    with open(pdb_file) as f:
        lines = f.readlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]

    for idx, line in enumerate(atom_lines):
        if line.startswith("ATOM"):
            parts = line.split()
            atom_name = parts[2]
            if by_serial:
                atom_serial = int(parts[1])
                mapping[atom_name] = atom_serial
            else:
                mapping[atom_name] = idx if zero_indexed else idx + 1

    return mapping


def idx_atoms_at_position(atom_mapping: dict[str, int], position: str) -> list[int]:
    """Get the indices of the heavy atoms at a given position in an amino acid.
    Args:
        atom_mapping: mapping of atom names to their indices from a PDB file
        position: position identifier to search for (e.g., "B", "G", "D", etc.)
    Returns:
        List of indices of the heavy atoms at the specified position
        (ascending order of their name, e.g. "G1" before "G2")
    """
    matches = [
        (name, index)
        for name, index in atom_mapping.items()
        if (position in name) and not name.startswith("H")
    ]

    if len(matches) > 1:
        matches.sort(key=lambda x: int(x[0][-1]))

    return [index for _, index in matches]


def get_atom_count(xyz_file: PathLike | str) -> int:
    """Get the number of atoms in an xyz file"""
    with open(xyz_file, "r") as file:
        first_line = file.readline().strip()
        return int(first_line)


def xyz_string(elements: Array1DStr, coordinates: Array2DFloat) -> str:
    """Make a XYZ string from elements and coordinates.
    Args:
        elements: elements symbols
        coordinates: coordinates [Å]
    Returns:
        XYZ string suitable for RDKit MolFromXYZBlock function
    """
    num_atoms = len(elements)
    xyz_string = f"{num_atoms}\n\n"
    for element, coords in zip(elements, coordinates):
        xyz_string += f"{element} {coords[0]} {coords[1]} {coords[2]}\n"
    return xyz_string


def rows_to_pdb_format(atom_rows: list[list[str]], add_end=True) -> str:
    """Format an atom row for a PDB file to ensure proper spacing."""
    # Renumber atom serials
    for i, row in enumerate(atom_rows, 1):
        row[1] = str(i)

    def format_atom_row(atom_row: list[str]) -> str:
        rec_type = atom_row[0]
        serial = int(atom_row[1])
        name = atom_row[2]
        resName = atom_row[3]
        chainID = atom_row[4]
        resSeq = int(atom_row[5])
        x, y, z = float(atom_row[6]), float(atom_row[7]), float(atom_row[8])
        occ = float(atom_row[9]) if len(atom_row) > 9 else 1.00
        tf = float(atom_row[10]) if len(atom_row) > 10 else 0.00
        element = atom_row[-1]
        return (
            f"{rec_type:<6}{serial:>5} {name:<5}{resName:<4}{chainID:<1}{resSeq:>4}"
            f"{x:>12.3f}{y:>8.3f}{z:>8.3f}"
            f"{occ:>6.2f}{tf:>6.2f}{element:>12}"
        )

    out_lines = [format_atom_row(row) for row in atom_rows]
    if add_end:
        out_lines.append("END")

    return "\n".join(out_lines) + "\n"


def replace_backbone(
    input_file: PathLike | str,
    output_file: PathLike | str,
    is_pro: bool = False,
) -> None:
    format = Path(input_file).suffix[1:]
    if format.lower() == "pdb":
        _replace_backbone_pdb(input_file, output_file, is_pro)
    elif format.lower() == "xyz":
        _replace_backbone_xyz(input_file, output_file, is_pro)
    else:
        raise ValueError(f"Unsupported format: {format}")


def _replace_backbone_pdb(
    input_file: PathLike | str, output_file: PathLike | str, is_pro: bool = False
) -> None:
    """Replace the backbone of the rotamer by a H in a PDB file.
    Args:
        input_file: path to the input PDB file
        output_file: path to the output PDB file
        is_pro: whether the rotamer is a proline
            If True, both the backbone N and C alpha are replaced by Hs
    Returns:
        None, writes the modified structure to the output file
    """

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(input_file, "r") as file:
        lines = file.readlines()
    atom_rows = [line.split() for line in lines if line.startswith("ATOM")]

    if is_pro:
        atoms_to_delete = {"H", "H1", "H2", "H3", "HA", "C", "O", "OXT"}
        atoms_to_replace_by_H = {"N", "CA"}
    else:
        atoms_to_delete = {"N", "H", "H1", "H2", "H3", "HA", "C", "O", "OXT"}
        atoms_to_replace_by_H = {"CA"}

    new_atom_rows = []
    for row in atom_rows:
        atom_name = row[2]
        if atom_name in atoms_to_delete:
            continue
        elif atom_name in atoms_to_replace_by_H:
            new_row = row.copy()
            new_row[2] = f"H{atom_name}"
            new_row[-1] = "H"
            if atom_name == "N":
                pro_N_row = new_row
            else:
                new_atom_rows.append(new_row)
        else:
            new_atom_rows.append(row)
    if is_pro:
        new_atom_rows.append(pro_N_row)

    pdb_str = rows_to_pdb_format(new_atom_rows)

    with open(output_file, "w") as output_file:
        output_file.write(pdb_str)


def _replace_backbone_xyz(
    input_file: PathLike | str, output_file: PathLike | str, is_pro: bool = False
) -> None:
    """Replace the backbone of the rotamer by a H in an XYZ file.
    Args:
        input_file: path to the input XYZ file
        output_file: path to the output XYZ file
        is_pro: whether the rotamer is a proline
            If True, both the backbone N and C alpha are replaced by Hs
    Returns:
        None, writes the modified structure to the output file
    """

    # Ensure that the output directory exists
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Read in original structure
    with open(input_file, "r") as file:
        lines = file.readlines()

    # Replace the comment line
    lines[1] = "backbone replaced by H\n"

    if is_pro:
        # Replace backbone N with H
        atom_line = lines[2].split()
        atom_line[0] = "H"
        lines[2] = " ".join(atom_line) + "\n"

        # Replace C alpha with H
        atom_line = lines[5].split()
        atom_line[0] = "H"
        lines[5] = " ".join(atom_line) + "\n"

        # Delete rest of backbone atoms
        lines.pop(3)
        lines.pop(3)
        lines.pop(4)
        lines.pop(4)
        lines.pop(4)
        lines.pop(-1)

        # Move the H replacing N to the end of the file
        lines.append(lines.pop(2))

        nb_atoms_deleted = 6

    else:
        # Delete backbone N and the 3 bonded Hs
        lines.pop(2)
        lines.pop(2)
        lines.pop(2)
        lines.pop(2)

        # Replace C alpha with H
        atom_line = lines[2].split()
        atom_line[0] = "H"
        lines[2] = " ".join(atom_line) + "\n"

        # Delete backbone H alpha and COO
        lines.pop(3)
        lines.pop(3)
        lines.pop(3)
        lines.pop(-1)

        nb_atoms_deleted = 8

    # Update the atom count in the first line
    new_atom_count = int(lines[0].strip()) - nb_atoms_deleted
    lines[0] = f"{new_atom_count}\n"

    # Write the modified structure to the output file
    with open(output_file, "w") as output_file:
        output_file.writelines(lines)
