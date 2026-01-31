from os import PathLike
from pathlib import Path

from morfeus.typing import Array1DStr, Array2DFloat


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
