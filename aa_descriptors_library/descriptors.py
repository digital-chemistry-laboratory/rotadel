import numpy as np
from morfeus import BuriedVolume, SASA, Sterimol, ConeAngle, Dispersion, XTB
from morfeus.typing import Array1DStr, Array2DFloat
from typing import Any

from aa_descriptors_library.constants import (
    NON_PERMUTABLE_INDICES_SIDECHAIN,
    THREE_TO_ONE_AA,
)


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
        updated coordinates of the molecule
    """
    central_coord = np.array(coords[central_index])
    bonded_coord = np.array(coords[bonded_index])
    current_dist = np.linalg.norm(central_coord - bonded_coord)
    direction_vect = (central_coord - bonded_coord) / current_dist
    displacement_vect = direction_vect * (target_dist - current_dist)
    new_central_coord = central_coord + displacement_vect
    new_coords = coords.copy()
    new_coords[central_index] = new_central_coord
    return new_coords


def clean_atomic_descriptors(
    morfeus_dict: dict[int, np.float64],
    res: str | None = None,
    is_bond_order: bool = False,
) -> dict[int, np.float64]:
    """Rewrite the dictionary returned by Morfeus to remove the permutable hydrogens.
    Args:
        morfeus_descriptor: dictionary of the atomic descriptor calculated by Morfeus
        res: three-letter code of the amino acid residue
        is_bond_order: whether the descriptor is bond orders
    Returns:
        Dictionary with only the non-permutable sidechain atoms (1-indexed)
    """
    indices_to_keep = NON_PERMUTABLE_INDICES_SIDECHAIN[THREE_TO_ONE_AA[res.upper()]]
    if is_bond_order:
        output_dict = {}
        for key, value in morfeus_dict.items():
            atom_1, atom_2 = map(int, key.split(", "))
            if atom_1 in indices_to_keep and atom_2 in indices_to_keep:
                output_dict[key] = value
    else:
        output_dict = {
            int(idx): value
            for idx, value in morfeus_dict.items()
            if int(idx) in indices_to_keep
        }

    return output_dict


def calc_descriptors(
    sidechain_el: Array1DStr,
    sidechain_coords: Array2DFloat,
    charge: int = 0,
    solvent: str = "ether",
    res: str | None = None,
) -> dict[str, Any]:
    """Calculate stereo-electronic descriptors on optimised sidechain.
    Args:
        sidechain_el: elements symbols of the sidechain-H geometry
        sidechain_coords: coordinates of the optimised sidechain-H geometry [Å]
        charge: charge of the sidechain-H geometry
        solvent: implicit solvent for the optimisation
        res: three-letter code of the amino acid residue
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
    xtb = XTB(sidechain_el, sidechain_coords, charge=charge, solvent=solvent)
    descriptors["homo"] = xtb.get_homo()
    descriptors["lumo"] = xtb.get_lumo()
    descriptors["electron_affinity"] = xtb.get_ea()
    descriptors["ionization_potential"] = xtb.get_ip()
    descriptors["electronegativity"] = xtb.get_electronegativity()
    descriptors["hardness"] = xtb.get_hardness()
    descriptors["dipole"] = xtb.get_dipole_moment()
    descriptors["global_nucleophilicity"] = xtb.get_global_descriptor("nucleophilicity")
    descriptors["global_electrophilicity"] = xtb.get_global_descriptor(
        "electrophilicity"
    )
    descriptors["local_nucleophilicity"] = clean_atomic_descriptors(
        xtb.get_fukui("local_nucleophilicity"), res=res
    )
    descriptors["local_electrophilicity"] = clean_atomic_descriptors(
        xtb.get_fukui("local_electrophilicity"), res=res
    )
    descriptors["fukui_minus"] = clean_atomic_descriptors(
        xtb.get_fukui("nucleophilicity"), res=res
    )
    descriptors["fukui_plus"] = clean_atomic_descriptors(
        xtb.get_fukui("electrophilicity"), res=res
    )
    descriptors["partial_charges"] = clean_atomic_descriptors(
        xtb.get_charges(), res=res
    )
    bond_orders = xtb.get_bond_orders()
    descriptors["bond_orders"] = clean_atomic_descriptors(
        {f"{k[0]}, {k[1]}": v for k, v in bond_orders.items()},
        res=res,
        is_bond_order=True,
    )

    return descriptors
