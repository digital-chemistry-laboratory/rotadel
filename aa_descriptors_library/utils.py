from pathlib import Path

from morfeus import read_xyz
from morfeus.typing import Array1DStr, Array2DFloat
from morfeus.utils import convert_elements
from spyrmsd import rmsd, graph


def rmsd_symmetric(
    xyz_file1: str | Path | None = None,
    xyz_file2: str | Path | None = None,
    el1: Array1DStr | None = None,
    coords1: Array2DFloat | None = None,
    el2: Array1DStr | None = None,
    coords2: Array2DFloat | None = None,
    diff_order: bool = True,
) -> float:
    """Calculate RMSD between two molecular structures, handling different atom ordering.
    Args:
        xyz_file1, xyz_file2: paths to XYZ files
        el1, el2: elements
        coords1, coords2: coordinates [Å]
        diff_order: whether atoms are in different order
    Returns:
        RMSD value [Å]
    """

    if xyz_file1 is not None:
        el1, coords1 = read_xyz(xyz_file1)
    if xyz_file2 is not None:
        el2, coords2 = read_xyz(xyz_file2)

    if diff_order:
        el_nums1 = convert_elements(el1)
        el_nums2 = convert_elements(el2)
        adj1 = graph.adjacency_matrix_from_atomic_coordinates(el_nums1, coords1)
        adj2 = graph.adjacency_matrix_from_atomic_coordinates(el_nums2, coords2)
        rmsd_value = rmsd.symmrmsd(
            coords1, coords2, el_nums1, el_nums2, adj1, adj2, center=True, minimize=True
        )
    else:
        rmsd_value = rmsd.rmsd(coords1, coords2, el1, el2, center=True, minimize=True)

    return rmsd_value
