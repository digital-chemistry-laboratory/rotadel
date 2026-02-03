from Bio.PDB import PDBIO
from morfeus.typing import Array1DStr, Array2DFloat
import numpy as np
from openmm.app import PDBFile, Modeller
import os
from os import PathLike
from pathlib import Path
from PeptideBuilder import Geometry
import PeptideBuilder
from rdkit import Chem
from rdkit.Chem import DetectChemistryProblems, MolFromPDBFile, MolFromXYZFile
from rdkit.Chem.rdDetermineBonds import DetermineBonds
from rdkit.Chem.rdmolops import GetMolFrags
import shutil
import subprocess
from typing import Any

from aa_descriptors_library.io_atoms import (
    convert_file,
    idx_atoms_at_position,
    map_atom_indices,
    read_geo,
    replace_backbone,
)
from aa_descriptors_library.constants import (
    BACKBONE_SMARTS,
    NUMBER_OF_CHI_ANGLES,
    SIDECHAIN_SMARTS,
)
from aa_descriptors_library import config


def gen_dihedral_constraints(
    dunbrack_data: dict[str, Any],
    atoms_indices: dict[str, int],
    constrain_N_CA_C_O_diangle: bool = True,
) -> list[tuple[list[int], float]]:
    """Generate the dihedral contraints for the specified whole rotamer.
    Args:
        dunbrack_data: data with angles extracted from the Dunbrack library
        atom_indices: mapping of PDB atom names to their 1-based row indices
        constrain_N_CA_C_O_diangle: whether to constrain the N-CA-C=O dihedral angle according to psi value
    Returns:
        list of tuples for dihedral constraints (1-based atom indices)
    """
    letter = dunbrack_data["letter"]
    nb_chis = NUMBER_OF_CHI_ANGLES[letter]

    dihedral_constraints = []
    N = atoms_indices["N"]
    CA = atoms_indices["CA"]
    CB = idx_atoms_at_position(atoms_indices, "B")[0]
    G = idx_atoms_at_position(atoms_indices, "G")[0]
    dihedral_constraints.append(([N, CA, CB, G], dunbrack_data["chi1"]))
    if nb_chis >= 2:
        D = idx_atoms_at_position(atoms_indices, "D")[0]
        dihedral_constraints.append(([CA, CB, G, D], dunbrack_data["chi2"]))
    if nb_chis >= 3:
        E = idx_atoms_at_position(atoms_indices, "E")[0]
        dihedral_constraints.append(([CB, G, D, E], dunbrack_data["chi3"]))
    if nb_chis >= 4:
        Z = idx_atoms_at_position(atoms_indices, "Z")[0]
        dihedral_constraints.append(([G, D, E, Z], dunbrack_data["chi4"]))

    # N-CA-C=O dihedral deduced from psi, considering the =O as being at 180° from N(i+1)
    if constrain_N_CA_C_O_diangle:
        C = atoms_indices["C"]
        O_double = atoms_indices["O"]
        N_CA_C_O_constraint = ([N, CA, C, O_double], dunbrack_data["psi"] - 180.0)
        dihedral_constraints.append(N_CA_C_O_constraint)

    return dihedral_constraints


def gen_NH3_dist_constraints(
    atoms_indices: dict[str, int],
    dist: float = 1.02,
) -> list[tuple[list[int], float]]:
    """Generate the distance constraints tuples between N and Hs of backbone NH3."""
    NH3_dist_constraints = []
    for H_name in ["H", "H2", "H3"]:
        NH3_dist_constraints.append(([atoms_indices["N"], atoms_indices[H_name]], dist))
    return NH3_dist_constraints


def gen_sidechain_constraints(
    dunbrack_data: dict[str, Any],
    atoms_indices: dict[str, int],
) -> tuple[list[tuple[list[int], float]] | None, list[int]]:
    """Generate the dihedral contraints and fixed atoms for the specified sidechain-H.
    Args:
        dunbrack_data: data with angles extracted from the Dunbrack library
        atom_indices: mapping of PDB atom names to their 1-based row indices
    Returns:
        Constraints for the sidechain-H optimisation (1-based atom indices):
            - chi2 (and chi3 for proline) dihedral constraints tuples
            - list of fixed atom indices
    """
    letter = dunbrack_data["letter"]

    # Dihedral constraint(s) for the H(s) to optimise
    if letter in ["C", "S", "T", "V"]:
        dihedral_constraints = None
    else:
        HCA = atoms_indices["HCA"]
        CB = idx_atoms_at_position(atoms_indices, "B")[0]
        G = idx_atoms_at_position(atoms_indices, "G")[0]
        D = idx_atoms_at_position(atoms_indices, "D")[0]
        dihedral_constraints = [([HCA, CB, G, D], dunbrack_data["chi2"])]
        if letter == "P":
            HN = atoms_indices["HN"]
            dihedral_constraints.append(([CB, G, D, HN], dunbrack_data["chi3"]))

    # Sidechain atoms are fixed except the H(s) on the same C as the H replacing the backbone
    atoms_to_opt = {"HCA", "HB", "HB2", "HB3"}
    if letter == "P":
        atoms_to_opt.update({"HN", "HD2", "HD3"})
    fixed_atoms = [idx for key, idx in atoms_indices.items() if key not in atoms_to_opt]

    return dihedral_constraints, fixed_atoms


def write_xcontrol(
    file: PathLike | str,
    fixed_atoms: list[int] = None,
    dihedral_constraints: list[tuple[list[int], float]] = None,
    distance_constraints: list[tuple[list[int], float]] = None,
    fc: float = 0.5,
    opt_engine: str | None = None,
) -> None:
    """Write input file for xTB optimisation with constraining or fixing
    Args:
        file: path to the input file to create
        fixed_atoms: list of atom indices to fix in space
        dihedral_constraints: list of tuples for dihedral constraints
            Each tuple should contain:
            - A list of exactly four atom indices (1-based) involved in the dihedral angle
            - A float specifying the desired dihedral angle in degrees
        distance_constraints: list of tuples for distance constraints
            Each tuple should contain:
            - A list of exactly two atom indices (1-based) involved in the distance constraint
            - A float specifying the desired distance in angstroms
        fc: force constant for constraints (only used if constraints are provided)
        opt_engine: optimisation engine to use
            - None uses xtb default engine: Approximate Normal Coordinate Rational Function optimizer (ANCopt)
            - "inertial": Fast Inertial Relaxation Engine (FIRE) (for cartesian coordinates)
    Returns:
        None, write input file
    """
    input = ""

    # Input block for constraints
    if dihedral_constraints or distance_constraints:
        for constraints_list, expected_atoms_num in zip(
            [dihedral_constraints, distance_constraints], [4, 2]
        ):
            if constraints_list and not all(
                isinstance(item, tuple)
                and len(item) == 2
                and len(item[0]) == expected_atoms_num
                for item in constraints_list
            ):
                docstring = write_xcontrol.__doc__
                raise ValueError(
                    "The provided constraints do not match the expected format. "
                    f"See docstring:\n\n{docstring}"
                )

        input += "$constrain\n"
        input += f"   force constant={fc}\n"
        for atoms, angle in dihedral_constraints if dihedral_constraints else []:
            input += f"   dihedral: {', '.join(map(str, atoms))}, {angle}\n"
        for atoms, distance in distance_constraints if distance_constraints else []:
            input += f"   distance: {', '.join(map(str, atoms))}, {distance}\n"
        input += "$end\n"

    # Input block for fixed atoms
    if fixed_atoms:
        if not all(isinstance(atom, int) for atom in fixed_atoms):
            raise ValueError("For fixing atoms, provide a list of fixed atom indices.")

        input += "$fix\n"
        input += f"   atoms: {','.join(map(str, fixed_atoms))}\n"
        input += "$end\n"

    # Input block for optimisation engine
    if opt_engine:
        input += "$opt\n"
        input += f"   engine={opt_engine}\n"
        input += "$end\n"

    # Write the input to the specified file
    with open(file, "w") as f:
        f.write(input)


def run_constraint_xtb(
    geo_file: PathLike | str,
    path_run: PathLike | str,
    fixed_atoms: list[int] = None,
    dihedral_constraints: list[tuple[list[int], float]] = None,
    distance_constraints: list[tuple[list[int], float]] = None,
    fc: float = 0.5,
    opt_engine: str | None = None,
    charge: int = 0,
    solvent: str = "ether",
) -> None:
    """Run constrained optimisation in implicit solvent with GFN2-xTB
    Args:
        geo_file: file for the starting structure
        path_run: folder to run the xTB optimisation (created if does not exist)
        fixed_atoms: list of atom indices to fix in space
        dihedral_constraints: list of tuples for dihedral constraints
            Each tuple should contain:
            - A list of exactly four atom indices (1-based) involved in the dihedral angle
            - A float specifying the desired dihedral angle in degrees
        distance_constraints: list of tuples for distance constraints
            Each tuple should contain:
            - A list of exactly two atom indices (1-based) involved in the distance constraint
            - A float specifying the desired distance in angstroms
        fc: force constant for constraints (only used if constraints are provided)
        opt_engine: optimisation engine to use
            - None uses xtb default engine: Approximate Normal Coordinate Rational Function optimizer (ANCopt)
            - "inertial": Fast Inertial Relaxation Engine (FIRE) (for cartesian coordinates)
        charge: charge of the molecule
        solvent: implicit solvent for the optimisation
    Returns:
        None, runs xTB optimisation
    """
    path_run = Path(path_run)
    # Create xcontrol file
    write_xcontrol(
        path_run / "xcontrol",
        fixed_atoms=fixed_atoms,
        dihedral_constraints=dihedral_constraints,
        distance_constraints=distance_constraints,
        fc=fc,
        opt_engine=opt_engine,
    )

    # Get absolute path of starting structure
    geo_file = Path(geo_file).resolve()

    # Run xtb from command line
    command = f"xtb {geo_file} --opt --chrg {charge} --alpb {solvent} -I xcontrol"
    with open(path_run / "xtb.out", "w") as stdout, open(
        f"{path_run}/xtb.err", "w"
    ) as stderr:
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = f"{config.OMP_NUM_THREADS},1"
        env["MKL_NUM_THREADS"] = f"{config.OMP_NUM_THREADS}"
        env["OMP_STACKSIZE"] = config.OMP_STACKSIZE
        env["OMP_MAX_ACTIVE_LEVELS"] = str(config.OMP_MAX_ACTIVE_LEVELS)
        subprocess.run(
            command.split(),
            cwd=path_run,
            stdout=stdout,
            stderr=stderr,
            env=env,
        )


def start_rotamer_geo(
    dunbrack_data: dict[str, Any],
    output_file: PathLike | str,
    charge: int,
    tautomer: str = None,
    set_N_CA_C_O_diangle: bool = True,
) -> None:
    """Build initial rotamer geometry from Dunbrack dihedral angles.
    Args:
        dunbrack_data: data with angles extracted from the Dunbrack library
        output_file: path to file to create, either PDB or XYZ
        charge: charge of the rotamer
        tautomer: tautomer of the rotamer (only for histidine, either 'D' or 'E')
        set_N_CA_C_O_diangle: whether to set the N-CA-C=O dihedral angle according to psi value
    Returns:
        None, writes the non-optimised rotamer geometry to the output file
    """
    output_file = Path(output_file)
    if output_file.suffix == ".pdb":
        convert_to_xyz = False
    elif output_file.suffix == ".xyz":
        convert_to_xyz = True
    else:
        raise ValueError("Geometry file must be either PDB or XYZ.")

    # Build initial geometry from Dunbrack dihedral angles with PeptideBuilder
    # Results in a PDB file without hydrogens
    letter = dunbrack_data["letter"]
    geo = Geometry.geometry(letter)
    geo.phi = dunbrack_data["phi"]
    geo.psi_im1 = dunbrack_data["psi"]
    if set_N_CA_C_O_diangle:
        geo.N_CA_C_O_diangle = dunbrack_data["psi"] - 180.0
    if letter == "R":
        geo.N_CA_CB_CG_diangle = dunbrack_data["chi1"]
        geo.CA_CB_CG_CD_diangle = dunbrack_data["chi2"]
        geo.CB_CG_CD_NE_diangle = dunbrack_data["chi3"]
        geo.CG_CD_NE_CZ_diangle = dunbrack_data["chi4"]
    elif letter == "I":
        geo.N_CA_CB_CG1_diangle = dunbrack_data["chi1"]
        geo.N_CA_CB_CG2_diangle = dunbrack_data["chi1"] + 121.3
        geo.CA_CB_CG1_CD1_diangle = dunbrack_data["chi2"]
    elif letter == "L":
        geo.N_CA_CB_CG_diangle = dunbrack_data["chi1"]
        geo.CA_CB_CG_CD1_diangle = dunbrack_data["chi2"]
        geo.CA_CB_CG_CD2_diangle = dunbrack_data["chi2"] - 108.4
    elif letter == "P":
        geo.N_CA_CB_CG_diangle = dunbrack_data["chi1"]
        geo.CA_CB_CG_CD_diangle = dunbrack_data["chi2"]
    elif letter == "T":
        geo.N_CA_CB_OG1_diangle = dunbrack_data["chi1"]
        geo.N_CA_CB_CG2_diangle = dunbrack_data["chi1"] + 120.3
    elif letter == "V":
        geo.N_CA_CB_CG1_diangle = dunbrack_data["chi1"]
        geo.N_CA_CB_CG2_diangle = dunbrack_data["chi1"] + 119.5
    else:
        nb_chis = NUMBER_OF_CHI_ANGLES[letter]
        geo.inputRotamers([dunbrack_data[f"chi{i}"] for i in range(1, nb_chis + 1)])
    structure = PeptideBuilder.initialize_res(geo)
    PeptideBuilder.add_terminal_OXT(structure)
    outfile = PDBIO()
    outfile.set_structure(structure)
    pdb_file = output_file.with_suffix(".pdb")
    outfile.save(str(pdb_file))

    # Add hydrogens according to species (charge and tautomer) with OpenMM
    pdb_wo_Hs = PDBFile(str(pdb_file))
    modeller = Modeller(pdb_wo_Hs.topology, pdb_wo_Hs.positions)
    if letter == "H":
        if charge == 0 and tautomer == "D":
            species = "HID"
        elif charge == 0 and tautomer == "E":
            species = "HIE"
        elif charge == 1:
            species = "HIP"
        elif charge == -1:
            species = "HIN"
        else:
            raise ValueError(
                f"Unknown charge/tautomer for histidine: {charge}, {tautomer}"
            )
    elif letter == "D" and charge == 0:
        species = "ASH"
    elif letter == "D" and charge == -1:
        species = "ASP"
    elif letter == "E" and charge == 0:
        species = "GLH"
    elif letter == "E" and charge == -1:
        species = "GLU"
    elif letter == "C" and charge == 0:
        species = "CYS"
    elif letter == "C" and charge == -1:
        species = "CYX"
    elif letter == "K" and charge == 1:
        species = "LYS"
    elif letter == "K" and charge == 0:
        species = "LYN"
    else:
        species = None
    modeller.addHydrogens(variants=[species])
    PDBFile.writeFile(modeller.topology, modeller.positions, file=str(pdb_file))

    # Convert PDB to XYZ file with Open Babel
    if convert_to_xyz:
        convert_file(pdb_file, output_file=output_file)


def get_opt_structures(
    dunbrack_data: dict[str, Any],
    charge: int,
    run_folder: PathLike | str,
    start_rotamer_file: PathLike | str,
    check: bool = True,
    tautomer: str | None = None,
    rotamer_id: str = "",
    constrain_N_CA_C_O_diangle: bool = True,
) -> tuple[
    Array1DStr,
    Array2DFloat,
    Array1DStr,
    Array2DFloat,
]:
    """Get the optimised geometries of side chain (with the backbone being replaced by an H) and rotamer.
    Args:
        dunbrack_data: data with angles extracted from the Dunbrack library
        charge: charge of the rotamer
        run_folder: path of folder in which to perform the full optimisation process
        start_rotamer_file: file of the whole rotamer to optimise (PDB or XYZ)
        check: whether to check for problems in the optimised structures
        tautomer: tautomer for histidine (used for structure checking if `check` is True)
        rotamer_id: ID of the rotamer (used for error messages if `check` is True)
        constrain_N_CA_C_O_diangle: whether to constrain the N-CA-C=O dihedral angle according to psi value
    Returns:
        el_whole, coord_whole, el_sidechain, coord_sidechain: elements and coordinates of the optimised structures
    """

    run_folder = Path(run_folder)
    if not run_folder.exists():
        run_folder.mkdir(parents=True)

    start_rotamer_file = Path(start_rotamer_file)
    file_format = start_rotamer_file.suffix
    xtbopt_file = "xtbopt" + file_format

    # Optimise the whole rotamer with constrained dihedral angles
    whole_folder = Path(run_folder) / "whole_rotamer"
    # Empty run folder (necessary for xtb)
    if whole_folder.exists():
        shutil.rmtree(whole_folder)
    whole_folder.mkdir(parents=True)

    whole_atoms_indices = map_atom_indices(start_rotamer_file)
    dihedral_constraints = gen_dihedral_constraints(
        dunbrack_data,
        whole_atoms_indices,
        constrain_N_CA_C_O_diangle=constrain_N_CA_C_O_diangle,
    )
    letter = dunbrack_data["letter"]
    # Force the Hs to stay on the NH3 for some rotamers which otherwise optimise to wrong structures
    if charge == -1 or letter == "N":
        NH3_dist_constraints = gen_NH3_dist_constraints(whole_atoms_indices)
    else:
        NH3_dist_constraints = None

    run_constraint_xtb(
        geo_file=start_rotamer_file,
        path_run=whole_folder,
        dihedral_constraints=dihedral_constraints,
        distance_constraints=NH3_dist_constraints,
        charge=charge,
        fc=1.0,
    )

    if check:
        whole_problems = has_structure_problems(
            whole_folder / xtbopt_file,
            charge=charge,
            aa_letter=letter,
            tautomer=tautomer,
        )
        if whole_problems:
            raise StructureProblemError(
                f"Wrong structure for the {rotamer_id} whole rotamer: {whole_problems}."
            )

    # Replace backbone by H in optimised rotamer
    sidechain_folder = Path(run_folder) / "sidechain-H"
    # Empty run folder (necessary for xtb)
    if sidechain_folder.exists():
        shutil.rmtree(sidechain_folder)
    sidechain_folder.mkdir(parents=True)
    start_sidechain_file = sidechain_folder / ("sidechain-H_start" + file_format)
    is_pro = dunbrack_data["res"] == "PRO"
    replace_backbone(whole_folder / xtbopt_file, start_sidechain_file, is_pro=is_pro)
    sidechain_atoms_indices = map_atom_indices(start_sidechain_file)

    # Optimise the Hs in the sidechain-H structure
    sidechain_dihedral_constraints, sidechain_fixed_atoms = gen_sidechain_constraints(
        dunbrack_data,
        sidechain_atoms_indices,
    )
    run_constraint_xtb(
        geo_file=start_sidechain_file,
        path_run=sidechain_folder,
        fixed_atoms=sidechain_fixed_atoms,
        dihedral_constraints=sidechain_dihedral_constraints,
        charge=charge,
        fc=5.0,
    )

    if check:
        sidechain_problems = has_structure_problems(
            sidechain_folder / xtbopt_file,
            charge=charge,
            aa_letter=letter,
            tautomer=tautomer,
            is_sidechainH=True,
        )
        if sidechain_problems:
            raise StructureProblemError(
                f"Wrong structure for the {rotamer_id} sidechain-H: {sidechain_problems}."
            )

    # Save the optimised sidechain-H and whole rotamer in dictionary
    el_whole, coord_whole = read_geo(whole_folder / xtbopt_file)
    el_sidechain, coord_sidechain = read_geo(sidechain_folder / xtbopt_file)

    if file_format == ".pdb":
        coord_whole = np.round(coord_whole, decimals=3)
        coord_sidechain = np.round(coord_sidechain, decimals=3)

    return (
        el_whole,
        coord_whole,
        el_sidechain,
        coord_sidechain,
        sidechain_atoms_indices,
    )


def has_structure_problems(
    geo_file: PathLike | str,
    charge: int,
    aa_letter: str,
    tautomer: str | None = None,
    is_sidechainH: bool = False,
    expected_nb_frags: int = 1,
    check_fragments: bool = True,
    check_chemistry: bool = True,
    check_smarts: bool = True,
) -> str | None:
    """Check if a structure is wrong.
    Args:
        geo_file: path to the geometry file to check (either PDB or XYZ)
        charge: charge of the molecule
        aa_letter: one-letter code of the amino acid
        tautomer: tautomer of the amino acid (only for histidine, either 'D' or 'E')
        is_sidechainH: whether the structure is a sidechain-H (default: False, whole rotamer)
        expected_nb_frags: expected number of fragments in the structure
        check_fragments: whether to check if different number of fragments than expected
        check_chemistry: whether to check for chemistry problems detected by RDKit
        check_smarts: whether to check if structure matches the expected SMARTS
    Returns:
        A message describing the problem if any of the checks fail, otherwise None
    """
    if Path(geo_file).suffix == ".pdb":
        m = MolFromPDBFile(str(geo_file), removeHs=False, sanitize=False)
    elif Path(geo_file).suffix == ".xyz":
        m = MolFromXYZFile(str(geo_file))
    else:
        raise ValueError("Geometry file must be either PDB or XYZ.")

    try:
        DetermineBonds(m, charge=charge)
    except ValueError as e:
        return f"RDKit could not determine bonds: {str(e)}"

    # Check number of fragments
    if check_fragments:
        frags = GetMolFrags(m, sanitizeFrags=False)
        if len(frags) != expected_nb_frags:
            return f"{len(frags)} fragments instead of the expected {expected_nb_frags}"

    # Check chemical problems
    if check_chemistry:
        problems = DetectChemistryProblems(m)
        if len(problems) > 0:
            return f"RDKit detected problems: {problems}"

    if check_smarts:
        # Remove implicit Hs to check SMARTS
        for atom in m.GetAtoms():
            atom.SetNoImplicit(True)
            if atom.HasProp("_MolFileHCount"):
                atom.ClearProp("_MolFileHCount")
        flags = Chem.SanitizeFlags.SANITIZE_ALL & ~Chem.SanitizeFlags.SANITIZE_ADJUSTHS
        err = Chem.SanitizeMol(m, sanitizeOps=flags, catchErrors=True)
        if err:
            return f"RDKit sanitization error code: {err}"

        # Check SMARTS match (backbone can be optimised either in zwitterion or neutral form)
        possible_smarts = []
        if aa_letter == "H" and charge == 0:
            sidechain_smarts = SIDECHAIN_SMARTS[aa_letter][charge][tautomer]
        else:
            sidechain_smarts = SIDECHAIN_SMARTS[aa_letter][charge]
        if is_sidechainH:
            sidechainH_smarts = get_smarts_sidechainH(sidechain_smarts, aa_letter)
            possible_smarts = [Chem.MolFromSmarts(sidechainH_smarts)]
        else:
            backbone_forms = (
                ["pro_zwitterion", "pro_neutral"]
                if aa_letter == "P"
                else ["zwitterion", "neutral"]
            )
            for form in backbone_forms:
                whole_smarts = BACKBONE_SMARTS[form].format(R=sidechain_smarts)
                possible_smarts.append(Chem.MolFromSmarts(whole_smarts))
        if not any(m.HasSubstructMatch(smarts) for smarts in possible_smarts):
            return "Structure does not match the amino acid SMARTS"

    return None


class StructureProblemError(RuntimeError):
    """Raised when an optimised structure fails the checks."""

    pass


def get_smarts_sidechainH(sidechain_smarts: str, aa_letter: str) -> str:
    """Adjust the sidechain SMARTS to make the sidechain-H SMARTS."""
    if aa_letter == "P":
        smarts = sidechain_smarts.replace("[CH2]", "[CH3]", 1)[::-1].replace(
            "]2HC[", "]3HC[", 1
        )[::-1]
    elif aa_letter in ["I", "T", "V"]:
        smarts = sidechain_smarts.replace("[CH]", "[CH2]", 1)
    else:
        smarts = sidechain_smarts.replace("[CH2]", "[CH3]", 1)
    return smarts
