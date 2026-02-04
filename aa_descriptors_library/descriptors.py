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
    sidechain_atoms_mapping: dict[str, int],
    res: str,
    decimals: int | None = None,
    is_bond_order: bool = False,
) -> dict[int, np.float64]:
    """Rewrite the dictionary returned by Morfeus to remove the permutable hydrogens.
    Args:
        morfeus_descriptor: dictionary of the atomic descriptor calculated by Morfeus
        sidechain_atoms_mapping: mapping of PDB atom names to their 1-based row indices
        res: three-letter code of the amino acid residue
        decimals: number of decimals to round the values to
        is_bond_order: whether the descriptor is bond orders
    Returns:
        Dictionary with only the non-permutable sidechain atoms (1-indexed)
    """
    indices_to_keep = [
        idx
        for name, idx in sidechain_atoms_mapping.items()
        if name in NON_PERMUTABLE_INDICES_SIDECHAIN[THREE_TO_ONE_AA[res.upper()]]
    ]
    if is_bond_order:
        output_dict = {}
        for key, value in morfeus_dict.items():
            atom_1, atom_2 = map(int, key.split(", "))
            if atom_1 in indices_to_keep and atom_2 in indices_to_keep:
                output_dict[key] = (
                    round(value, decimals) if decimals is not None else value
                )
    else:
        output_dict = {
            int(idx): round(value, decimals) if decimals is not None else value
            for idx, value in morfeus_dict.items()
            if int(idx) in indices_to_keep
        }

    return output_dict


def calc_descriptors(
    sidechain_el: Array1DStr,
    sidechain_coords: Array2DFloat,
    sidechain_atoms_mapping: dict[str, int],
    res: str,
    charge: int = 0,
    solvent: str = "ether",
) -> dict[str, Any]:
    """Calculate stereo-electronic descriptors on optimised sidechain.
    Args:
        sidechain_el: elements symbols of the sidechain-H geometry
        sidechain_coords: coordinates of the optimised sidechain-H geometry [Å]
        sidechain_atoms_mapping: mapping of PDB atom names to their 1-based row indices
        res: three-letter code of the amino acid residue
        charge: charge of the sidechain-H geometry
        solvent: implicit solvent for the optimisation
    Returns:
        Dictionary of the calculated descriptors
    """
    descriptors = {}

    # Buried volume
    bv = BuriedVolume(
        sidechain_el, sidechain_coords, metal_index=1
    )  # 1-indexed in Morfeus
    descriptors["v_bur"] = round(bv.fraction_buried_volume, 4)

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
    descriptors["sasa_tot"] = round(sasa.area, 2)
    descriptors["sasa_polar"] = round(sasa_polar, 2)
    descriptors["sasa_unpolar"] = round(sasa_unpolar, 2)

    # Sterimol parameters
    sterimol = Sterimol(sidechain_el, sidechain_coords, dummy_index=1, attached_index=2)
    descriptors["sterimol_l"] = round(sterimol.L_value, 3)
    descriptors["sterimol_b1"] = round(sterimol.B_1_value, 3)
    descriptors["sterimol_b5"] = round(sterimol.B_5_value, 3)

    # Cone angle
    # Needs to push away the central atom (H replacing alpha carbon)
    # so that the van der Waals spheres are not overlapping
    adjusted_coord = move_central_atom(
        sidechain_coords, central_index=0, bonded_index=1
    )
    cone_angle = ConeAngle(
        sidechain_el, np.asarray(adjusted_coord, dtype=np.float64), atom_1=1
    )
    descriptors["cone_angle"] = round(cone_angle.cone_angle, 2)

    # Dispersion descriptor
    dispersion = Dispersion(sidechain_el, sidechain_coords)
    descriptors["p_int"] = round(dispersion.p_int, 3)

    # xTB descriptors
    xtb = XTB(sidechain_el, sidechain_coords, charge=charge, solvent=solvent)
    descriptors["homo"] = round(xtb.get_homo(), 5)
    descriptors["lumo"] = round(xtb.get_lumo(), 5)
    descriptors["electron_affinity"] = round(xtb.get_ea(), 4)
    descriptors["ionization_potential"] = round(xtb.get_ip(), 4)
    descriptors["electronegativity"] = round(xtb.get_electronegativity(), 5)
    descriptors["hardness"] = xtb.get_hardness()
    descriptors["dipole"] = xtb.get_dipole_moment()
    descriptors["global_nucleophilicity"] = xtb.get_global_descriptor("nucleophilicity")
    descriptors["global_electrophilicity"] = round(
        xtb.get_global_descriptor("electrophilicity"), 5
    )
    descriptors["local_nucleophilicity"] = clean_atomic_descriptors(
        xtb.get_fukui("local_nucleophilicity"), sidechain_atoms_mapping, res
    )
    descriptors["local_electrophilicity"] = clean_atomic_descriptors(
        xtb.get_fukui("local_electrophilicity"),
        sidechain_atoms_mapping,
        res,
        decimals=5,
    )
    descriptors["fukui_minus"] = clean_atomic_descriptors(
        xtb.get_fukui("nucleophilicity"), sidechain_atoms_mapping, res
    )
    descriptors["fukui_plus"] = clean_atomic_descriptors(
        xtb.get_fukui("electrophilicity"), sidechain_atoms_mapping, res
    )
    descriptors["partial_charges"] = clean_atomic_descriptors(
        xtb.get_charges(), sidechain_atoms_mapping, res, decimals=5
    )
    bond_orders = xtb.get_bond_orders()
    descriptors["bond_orders"] = clean_atomic_descriptors(
        {f"{k[0]}, {k[1]}": v for k, v in bond_orders.items()},
        sidechain_atoms_mapping,
        res,
        decimals=5,
        is_bond_order=True,
    )

    return descriptors
