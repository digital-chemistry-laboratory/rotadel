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
        coordinates = np.array([a.coord for a in atoms])
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
        matches.sort(key=lambda x: x[0])

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


def replace_backbone(
    input_file: PathLike | str, output_file: PathLike | str, is_pro: bool = False
) -> None:
    """Replace the backbone of the rotamer by a H.
    Args:
        input_file: path to the input xyz file
        output_file: path to the output xyz file
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
