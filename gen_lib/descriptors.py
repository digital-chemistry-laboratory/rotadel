import numpy as np
from morfeus import BuriedVolume, SASA, Sterimol, ConeAngle, Dispersion, XTB
from morfeus.typing import Array1DStr, Array2DFloat
from typing import Any
from rdkit.Chem import MolFromXYZBlock, rdDetermineBonds

from .utils import xyz_string


def move_central_atom(
    coords: Array2DFloat,
    central_index: int,
    bonded_index: int,
    target_dist: float | None = 2.28,
) -> Array2DFloat:
    """Modify the central atom coordinates such that it is at a target distance from the bonded atom
    Args:
        coords: coordinates of the molecule
        central_index: index of the atom to move away (0-indexed)
        bonded_index: index of the atom bonded to central atom (0-indexed)
        target_dist: target distance between central and bonded atoms
    Returns:
        coords: updated coordinates of the molecule
    """
    central_coord = np.array(coords[central_index])
    bonded_coord = np.array(coords[bonded_index])
    current_dist = np.linalg.norm(central_coord - bonded_coord)
    direction_vect = (central_coord - bonded_coord) / current_dist
    displacement_vect = direction_vect * (target_dist - current_dist)
    new_central_coord = central_coord + displacement_vect
    coords[central_index] = new_central_coord
    return coords


def reindex_sidechain(morfeus_dict: dict[int, np.float64]) -> dict[int, np.float64]:
    """Rewrite the dictionary returned by Morfeus to remove the H dummy atom (first entry) and reindex in 0-base.
    Args:
        morfeus_descriptor: dictionary of the atomic descriptor calculated by Morfeus
    Returns:
        Dictionary with only the sidechain atoms and 0-indexed
    """
    return {int(idx) - 1: value for idx, value in morfeus_dict.items() if idx != "1"}


def calc_descriptors(
    sidechain_el: Array1DStr, sidechain_coords: Array2DFloat
) -> dict[str, Any]:
    """Calculate stereo-electronic descriptors on optimised rotamer and sidechain.
    Args:
        sidechain_el: elements symbols of the sidechain-H geometry
        sidechain_coords: coordinates of the optimised sidechain-H geometry [Å]
    Returns:
        Dictionary of the calculated descriptors
    """
    descriptors = {}

    # Buried volume
    bv = BuriedVolume(
        sidechain_el, sidechain_coords, metal_index=1
    )  # 1-indexed in Morfeus
    descriptors["v_bur"] = bv.fraction_buried_volume

    # Solvent accessible surface area
    sasa = SASA(sidechain_el, sidechain_coords)
    sasa_polar = 0
    sasa_unpolar = 0
    for i, atom in enumerate(sidechain_el):
        sasa_atom = sasa.atom_areas[i + 1]  # 1-indexed in Morfeus
        if atom in ["C", "H"]:
            sasa_unpolar += sasa_atom
        else:
            sasa_polar += sasa_atom
    descriptors["sasa_tot"] = sasa.area
    descriptors["sasa_polar"] = sasa_polar
    descriptors["sasa_unpolar"] = sasa_unpolar

    # Sterimol parameters
    sterimol = Sterimol(sidechain_el, sidechain_coords, dummy_index=1, attached_index=2)
    descriptors["sterimol_l"] = sterimol.L_value
    descriptors["sterimol_b1"] = sterimol.B_1_value
    descriptors["sterimol_b5"] = sterimol.B_5_value

    # Cone angle
    # Needs to push away the central atom (H replacing alpha carbon)
    # so that the van der Waals spheres are not overlapping
    adjusted_coord = move_central_atom(
        sidechain_coords, central_index=0, bonded_index=1
    )
    cone_angle = ConeAngle(sidechain_el, adjusted_coord, atom_1=1)
    descriptors["cone_angle"] = cone_angle.cone_angle

    # Dispersion descriptor
    dispersion = Dispersion(sidechain_el, sidechain_coords)
    descriptors["p_int"] = dispersion.p_int

    # xTB descriptors
    xtb = XTB(sidechain_el, sidechain_coords)
    descriptors["homo"] = xtb.get_homo()
    descriptors["lumo"] = xtb.get_lumo()
    ea = xtb.get_ea()
    ip = xtb.get_ip(corrected=True)
    descriptors["electron_affinity"] = ea
    descriptors["ionization_potential"] = ip
    descriptors["electronegativity"] = (ip + ea) / 2
    descriptors["hardness"] = ip - ea
    descriptors["global_nucleophilicity"] = xtb.get_global_descriptor(
        "nucleophilicity", corrected=True
    )
    descriptors["global_electrophilicity"] = xtb.get_global_descriptor(
        "electrophilicity", corrected=True
    )
    descriptors["local_nucleophilicity"] = reindex_sidechain(
        xtb.get_fukui("local_nucleophilicity")
    )
    descriptors["local_electrophilicity"] = reindex_sidechain(
        xtb.get_fukui("local_electrophilicity")
    )
    descriptors["fukui_minus"] = reindex_sidechain(xtb.get_fukui("nucleophilicity"))
    descriptors["fukui_plus"] = reindex_sidechain(xtb.get_fukui("electrophilicity"))

    # Partial charges
    descriptors["partial_charges"] = reindex_sidechain(xtb.get_charges())

    # Bond orders
    xyz_str = xyz_string(sidechain_el, sidechain_coords)
    mol = MolFromXYZBlock(xyz_str)
    rdDetermineBonds.DetermineConnectivity(mol)
    bond_indices = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()
    ]
    bo = {}
    for atom_1, atom_2 in bond_indices:
        if atom_1 == 0 or atom_2 == 0:
            continue
        bo[(atom_1, atom_2)] = xtb.get_bond_order(
            atom_1 + 1, atom_2 + 1
        )  # 1-indexed for Morfeus
    descriptors["bond_orders"] = bo

    return descriptors
