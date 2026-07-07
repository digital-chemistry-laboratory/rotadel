from pathlib import Path

from morfeus.typing import Array1DStr, Array2DFloat
from morfeus.utils import convert_elements
from spyrmsd import rmsd, graph

from aa_descriptors_library.io_atoms import read_geo


def rmsd_symmetric(
    file1: str | Path | None = None,
    file2: str | Path | None = None,
    el1: Array1DStr | None = None,
    coords1: Array2DFloat | None = None,
    el2: Array1DStr | None = None,
    coords2: Array2DFloat | None = None,
    diff_order: bool = True,
) -> float:
    """Calculate RMSD between two molecular structures, handling different atom ordering.
    Args:
        file1, file2: paths to coordinate files
        el1, el2: elements
        coords1, coords2: coordinates [Å]
        diff_order: whether atoms are in different order
    Returns:
        RMSD value [Å]
    """

    if file1 is not None:
        el1, coords1 = read_geo(file1)
    if file2 is not None:
        el2, coords2 = read_geo(file2)

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


def format_duration(seconds: float) -> str:
    """Format a duration time from seconds to h:m:s string."""
    seconds_int = int(round(seconds))
    hours = seconds_int // 3600
    minutes = (seconds_int % 3600) // 60
    secs = seconds_int % 60
    return f"{hours}h {minutes}m {secs}s"
