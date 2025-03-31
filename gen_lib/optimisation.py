from os import PathLike
from pathlib import Path
import subprocess
import shutil
from morfeus import read_xyz
from morfeus.typing import Array1DStr, Array2DFloat
from typing import Any

from .utils import get_atom_count

# from .data import CHARGES_AA


def gen_dihedral_constraints(
    dunbrack_data: dict[str, Any],
) -> list[tuple[list[int], float]]:
    """Generate the dihedral contraints for the specified rotamer.
    Args:
        dunbrack_data: data with angles extracted from the Dunbrack library
    Returns:
        list of tuples for dihedral constraints (1-based atom indices)
    """
    res = dunbrack_data["res"]
    chi1, chi2, chi3, chi4 = (
        dunbrack_data["chi1"],
        dunbrack_data["chi2"],
        dunbrack_data["chi3"],
        dunbrack_data["chi4"],
    )

    if res in ["ARG", "LYS"]:
        dihedral_constraints = [
            ([1, 2, 5, 6], chi1),
            ([2, 5, 6, 7], chi2),
            ([5, 6, 7, 8], chi3),
            ([6, 7, 8, 9], chi4),
        ]
    elif res in ["GLN", "MET"]:
        dihedral_constraints = [
            ([1, 2, 5, 6], chi1),
            ([2, 5, 6, 7], chi2),
            ([5, 6, 7, 8], chi3),
        ]
    elif res == "GLU":
        dihedral_constraints = [
            ([1, 2, 5, 6], chi1),
            ([2, 5, 6, 7], chi2),
            ([5, 6, 7, 9], chi3),
        ]
    elif res in ["CYS", "PRO", "SER", "THR", "VAL"]:
        dihedral_constraints = [([1, 2, 5, 6], chi1)]
    elif res in ["ASN", "HIS", "LEU", "PHE", "TRP"]:
        dihedral_constraints = [([1, 2, 5, 6], chi1), ([2, 5, 6, 7], chi2)]
    elif res in ["ASP", "ILE"]:
        dihedral_constraints = [([1, 2, 5, 6], chi1), ([2, 5, 6, 8], chi2)]
    elif res == "TYR":
        dihedral_constraints = [([1, 2, 5, 6], chi1), ([2, 5, 6, 9], chi2)]

    return dihedral_constraints


def write_xcontrol(
    file: PathLike | str,
    type: str,
    fixed_atoms: list[int] = None,
    dihedral_constraints: list[tuple[list[int], float]] = None,
    fc: float = 0.5,
) -> None:
    """Write input file for xTB optimisation with constraining or fixing
    Args:
        file: path to the input file to create
        type: "dihedral" for constrained dihedral angles or "fix" for fixed atoms
        fixed_atoms: list of atom indices to fix in space (only when `type` is "fix")
        dihedral_constraints: list of tuples for dihedral constraints (only when `type` is "dihedral").
            Each tuple should contain:
            - A list of exactly four atom indices (1-based) involved in the dihedral angle
            - A float specifying the desired dihedral angle in degrees
        fc: force constant for constraints (only when `type` is "dihedral")
    Returns:
        None, write input file
    """

    # Block for dihedral constraint
    if type == "dihedral":
        if dihedral_constraints is None or not all(
            isinstance(item, tuple) and len(item) == 2 for item in dihedral_constraints
        ):
            raise ValueError(
                "For 'dihedral', provide a list of tuples containing the angle and a list of atom indices."
            )
        if fixed_atoms is not None:
            raise ValueError("For 'dihedral', do not provide 'fixed_atoms'.")

        # Build input for dihedral constraints
        input = "$constrain\n"
        input += f"   force constant={fc}\n"
        for atoms, angle in dihedral_constraints:
            if len(atoms) != 4:
                raise ValueError(
                    "Each dihedral constraint must have exactly 4 atom indices."
                )
            input += f"   dihedral: {', '.join(map(str, atoms))}, {angle}\n"
        input += "$end\n"

    # Block for fixed atoms
    elif type == "fix":
        if fixed_atoms is None or not all(
            isinstance(atom, int) for atom in fixed_atoms
        ):
            raise ValueError("For 'fix', provide a list of fixed atom indices.")
        if dihedral_constraints is not None:
            raise ValueError("For 'fix', do not provide 'dihedral_constraints'.")

        # Build input for fix constraints
        input = "$fix\n"
        input += f"   atoms: {','.join(map(str, fixed_atoms))}\n"
        input += "$end\n"

    else:
        raise ValueError("type must be either 'dihedral' or 'fix'.")

    # Write the input to the specified file
    with open(file, "w") as f:
        f.write(input)


def run_constraint_xtb(
    xyz_file: PathLike | str,
    path_run: PathLike | str,
    type_constraint: str,
    fixed_atoms: list[int] = None,
    dihedral_constraints: list[tuple[list[int], float]] = None,
    fc: float = 0.5,
    charge: int = 0,
    solvent: str = "ether",
) -> None:
    """Run constrained optimisation in implicit solvent with GFN2-xTB
    Args:
        xyz_file: xyz file for the starting structure
        path_run: folder to run the xTB optimisation (created if does not exist)
        type: "dihedral" for constrained dihedral angles or "fix" for fixed atoms
        fixed_atoms: list of atom indices to fix in space (only when `type` is "fix")
        dihedral_constraints: list of tuples for dihedral constraints (only when `type` is "dihedral").
            Each tuple should contain:
            - A list of exactly four atom indices (1-based) involved in the dihedral angle
            - A float specifying the desired dihedral angle in degrees
        fc: force constant for constraints (only when `type` is "dihedral")
        charge: charge of the molecule
        solvent: implicit solvent for the optimisation (default: ether)
    Returns:
        None, runs xTB optimisation
    """
    path_run = Path(path_run)
    # Create xcontrol file
    write_xcontrol(
        path_run / "xcontrol",
        type=type_constraint,
        fixed_atoms=fixed_atoms,
        dihedral_constraints=dihedral_constraints,
        fc=fc,
    )

    # Get absolute path of starting structure
    xyz_file = Path(xyz_file).resolve()

    # Run xtb from command line
    command = f"xtb {xyz_file} --opt --chrg {charge} --alpb {solvent} -I xcontrol"
    with open(path_run / "xtb.out", "w") as stdout, open(
        f"{path_run}/xtb.err", "w"
    ) as stderr:
        subprocess.run(
            command.split(),
            cwd=path_run,
            stdout=stdout,
            stderr=stderr,
        )


def replace_backbone(input_file: PathLike | str, output_file: PathLike | str) -> None:
    """Replace the backbone of the rotamer by a H.
    Args:
        input_file: path to the input xyz file
        output_file: path to the output xyz file
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

    # Delete the N backbone atom
    lines.pop(2)

    # Replace C alpha with H
    atom_line = lines[2].split()
    atom_line[0] = "H"
    lines[2] = " ".join(atom_line) + "\n"

    # Delete C and O backbone
    lines.pop(3)
    lines.pop(3)

    # Find the last O atom, delete it and the following four H atoms
    last_o_index = len(lines) - 1
    while last_o_index >= 2:
        if lines[last_o_index].startswith("O"):
            break
        last_o_index -= 1
    for _ in range(5):
        lines.pop(last_o_index)

    # Update the atom count in the first line
    new_atom_count = int(lines[0].strip()) - 8
    lines[0] = f"{new_atom_count}\n"

    # Write the modified structure to the output file
    with open(output_file, "w") as output_file:
        output_file.writelines(lines)


def get_opt_structures(
    dunbrack_data: dict[str, Any],
    run_folder: PathLike | str,
    rotamer_xyz: PathLike | str,
) -> tuple[
    Array1DStr,
    Array2DFloat,
    Array1DStr,
    Array2DFloat,
]:
    """Get the optimised geometries of side chain (with the backbone being replaced by an H) and rotamer
    Args:
        dunbrack_data: data with angles extracted from the Dunbrack library
        run_folder: path of folder in which to perform the full optimisation process
        rotamer_xyz: xyz file of the whole rotamer to optimise
    Returns:
        el_whole, coord_whole, el_sidechain, coord_sidechain: elements and coordinates of the optimised structures
    """

    run_folder = Path(run_folder)
    if not run_folder.exists():
        run_folder.mkdir(parents=True)

    # charge = CHARGES_AA[rotamer_key[0]] # TODO: deal with a.a. generated by PeptideBuilder which should be charged

    # Optimise the whole rotamer with constrained dihedral angles
    whole_folder = Path(run_folder) / "whole_rotamer"
    # Empty run folder (necessary for xtb)
    if whole_folder.exists():
        shutil.rmtree(whole_folder)
    whole_folder.mkdir(parents=True)

    dihedral_constraints = gen_dihedral_constraints(dunbrack_data)

    run_constraint_xtb(
        xyz_file=rotamer_xyz,
        path_run=whole_folder,
        type_constraint="dihedral",
        dihedral_constraints=dihedral_constraints,
    )

    # Replace backbone by H in optimised rotamer
    sidechain_folder = Path(run_folder) / "sidechain-H"
    # Empty run folder (necessary for xtb)
    if sidechain_folder.exists():
        shutil.rmtree(sidechain_folder)
    sidechain_folder.mkdir(parents=True)
    sidechain_start_xyz = sidechain_folder / "sidechain-H_start.xyz"
    replace_backbone(whole_folder / "xtbopt.xyz", sidechain_start_xyz)

    # Optimise the H in the sidechain-H structure with fixed sidechain atoms
    last_atom = get_atom_count(sidechain_start_xyz)
    fixed_atoms = list(range(2, last_atom + 1))  # atom indices are 1-based

    run_constraint_xtb(
        xyz_file=sidechain_start_xyz,
        path_run=sidechain_folder,
        type_constraint="fix",
        fixed_atoms=fixed_atoms,
    )

    # Save the optimised sidechain-H and whole rotamer in dictionary
    el_whole, coord_whole = read_xyz(whole_folder / "xtbopt.xyz")
    el_sidechain, coord_sidechain = read_xyz(sidechain_folder / "xtbopt.xyz")

    return (
        el_whole,
        coord_whole,
        el_sidechain,
        coord_sidechain,
    )
