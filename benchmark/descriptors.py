import numpy as np

from aa_descriptors_library.query import query_average
from aa_descriptors_library.sql import get_descriptors_names


def get_descriptors_vectors_peptide(
    peptide: str,
    descriptor_names: list[str] | None = None,
    padding_len: int | None = None,
    cache: dict[tuple[str | None, str, str | None], list[float]] | None = None,
) -> np.ndarray:
    """Get the rotamer feature vectors for a given peptide sequence.
    Args:
        peptide: peptide sequence
        descriptor_names: list of descriptor names to include (if None, include all)
        padding_len: length to which L should be padded
        cache: dictionary to store descriptors of previously queried left_neighbour-residue-right_neighbour combinations
    Returns:
        L x D array,
            with L the number of residues in the peptide + optional padding, and D the number of rotamer descriptors
    """
    res_letters = list(peptide)
    if padding_len is not None and padding_len < len(res_letters):
        raise ValueError("Padding length cannot be smaller than the peptide length.")

    if descriptor_names is None:
        descriptor_names = get_descriptors_names()

        # For now, remove the atomic descriptors
        descriptor_names = [
            descriptor
            for descriptor in descriptor_names
            if not any(
                substring in descriptor
                for substring in ("local", "fukui", "partial", "bond")
            )
        ]

    vectors = []
    # Cache to avoid redundant queries
    if cache is None:
        cache = {}
    for i, res in enumerate(res_letters):
        left_neighbour = res_letters[i - 1] if i > 0 else None
        right_neighbour = res_letters[i + 1] if i < len(res_letters) - 1 else None
        cache_key = (left_neighbour, res, right_neighbour)

        if cache_key in cache:
            res_vector = cache[cache_key]
        else:
            res_descriptors = query_average(
                res,
                left_neighbour=left_neighbour,
                right_neighbour=right_neighbour,
            )
            res_vector = [
                res_descriptors[descriptor] for descriptor in descriptor_names
            ]
            cache[cache_key] = res_vector
        vectors.append(res_vector)

    # Pad with vectors of -1 values if length of peptide is smaller than padding_len
    vectors = np.asarray(vectors, dtype=float)
    if padding_len is not None and padding_len > len(res_letters):
        padding_rows = padding_len - len(res_letters)
        padding = np.full((padding_rows, len(descriptor_names)), -1.0, dtype=float)
        vectors = np.vstack([vectors, padding])

    return vectors
