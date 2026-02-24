import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from aa_descriptors_library.io import lock
from aa_descriptors_library.query import query_average
from aa_descriptors_library.sql import get_descriptors_names
from aa_descriptors_library.utils import format_duration


class PeptideDescriptorsPipeline:
    _mem_per_cpu_gb = 2.0
    _data_dir_root = Path(__file__).resolve().parent / "data"
    _splits_names = ("train", "val", "test")
    _setups_indices = [0, 1, 2, 3, 4]

    def __init__(
        self,
        dataset_name: str | None = None,
        peptide: str | None = None,
        workers: int = 1,
        batch_size: int | None = None,
        rebuild_store: bool = False,
    ) -> None:
        if dataset_name is not None and peptide is not None:
            raise ValueError("Provide either a dataset name or a peptide sequence.")

        self._dataset_name = dataset_name
        self._peptide = peptide
        self._workers = max(1, workers)
        self._batch_size = None if batch_size is None else max(1, batch_size)
        self._rebuild_store = rebuild_store

        # Cache to store descriptors of already queried left_neighbour-residue-right_neighbour combinations
        self._cache: dict[tuple[str | None, str, str | None], list[float]] = {}

        self._output_dir = (
            self._data_dir_root / dataset_name if dataset_name is not None else None
        )
        self._npz_store_path = (
            self._output_dir / f"whole_{dataset_name}_descriptors.npz"
            if self._output_dir is not None
            else None
        )
        self._dataset_splits_folder = (
            Path(__file__).resolve().parents[1]
            / "peptidy_zenodo"
            / "data"
            / dataset_name
            if dataset_name is not None
            else None
        )

    def get_peptide_descriptors_vectors(
        self,
        peptide: str | None = None,
        descriptor_names: list[str] | None = None,
        padding_len: int | None = None,
        cache_only: bool = False,
    ) -> np.ndarray:
        """Get the rotamer feature vectors for a given peptide sequence.
        Args:
            peptide: peptide sequence
            descriptor_names: list of descriptor names to include (if None, include all)
            padding_len: length to which L should be padded
            cache_only: if True, only use cache and raise error if any triplet combination is missing from cache
        Returns:
            L x D array (L = peptide length with optional padding, D = number of descriptors)
        """
        sequence = peptide if peptide is not None else self._peptide
        if sequence is None:
            raise ValueError(
                "No peptide provided. Set peptide on instance or pass peptide argument."
            )

        res_letters = list(sequence)
        if padding_len is not None and padding_len < len(res_letters):
            raise ValueError(
                "Padding length cannot be smaller than the peptide length."
            )

        if descriptor_names is None:
            descriptor_names = get_descriptors_names()

        vectors: list[list[float]] = []

        for i, res in enumerate(res_letters):
            left_neighbour = res_letters[i - 1] if i > 0 else None
            right_neighbour = res_letters[i + 1] if i < len(res_letters) - 1 else None
            cache_triplet = (left_neighbour, res, right_neighbour)

            if cache_triplet in self._cache:
                res_vector = self._cache[cache_triplet]
            else:
                if cache_only:
                    raise KeyError(
                        f"Triplet {cache_triplet} is missing from cache while building peptide '{sequence}'."
                    )
                res_descriptors = query_average(
                    res,
                    left_neighbour=left_neighbour,
                    right_neighbour=right_neighbour,
                )
                res_vector = [
                    res_descriptors[descriptor] for descriptor in descriptor_names
                ]
                self._cache[cache_triplet] = res_vector

            vectors.append(res_vector)

        vectors_arr = np.asarray(vectors, dtype=np.float32)
        if padding_len is not None and padding_len > len(res_letters):
            padding_rows = padding_len - len(res_letters)
            padding = np.full(
                (padding_rows, len(descriptor_names)), -1.0, dtype=np.float32
            )
            vectors_arr = np.vstack([vectors_arr, padding])

        return vectors_arr

    def _compute_descriptors_vectors_for_triplets(
        self,
        res_triplets: list[tuple[str | None, str, str | None]],
        descriptor_names: list[str],
    ) -> list[tuple[tuple[str | None, str, str | None], list[float]]]:
        """Compute descriptors vectors for a list of residues with neighbours triplets"""
        descriptor_names_tuple = tuple(descriptor_names)
        if self._workers <= 1:
            return [
                _compute_res_descriptor_vector_worker(triplet, descriptor_names_tuple)
                for triplet in res_triplets
            ]

        # Run in parallel across CPU cores
        with ProcessPoolExecutor(max_workers=self._workers) as executor:
            return list(
                executor.map(
                    _compute_res_descriptor_vector_worker,
                    res_triplets,
                    [descriptor_names_tuple] * len(res_triplets),
                    chunksize=64,
                )
            )

    def _gen_res_triplets(
        self, peptide: str
    ) -> set[tuple[str | None, str, str | None]]:
        """Generate the set of (left, residue, right) triplets for a given peptide sequence."""
        res_triplets: set[tuple[str | None, str, str | None]] = set()
        res_letters = list(peptide)
        for i, res in enumerate(res_letters):
            left = res_letters[i - 1] if i > 0 else None
            right = res_letters[i + 1] if i < len(res_letters) - 1 else None
            res_triplets.add((left, res, right))
        return res_triplets

    def _collect_peptides_data(
        self,
    ) -> tuple[dict[tuple[int, str], tuple[np.ndarray, np.ndarray]], list[str], int]:
        """Collect from dataset the splits data, the unique peptide sequences, and the max peptide length."""
        split_data: dict[tuple[int, str], tuple[np.ndarray, np.ndarray]] = {}
        unique_peptides: list[str] = []
        seen_peptides: set[str] = set()
        max_peptide_len = 0

        for setup_idx in self._setups_indices:
            for split in self._splits_names:
                dataset_split_path = (
                    self._dataset_splits_folder / f"setup-{setup_idx}" / f"{split}.csv"
                )
                df_split = pd.read_csv(dataset_split_path)
                peptides = df_split["sequence"].astype(str).to_numpy()
                y = df_split["label"].to_numpy(dtype=np.float32)
                split_data[(setup_idx, split)] = (peptides, y)
                for peptide in peptides:
                    if peptide not in seen_peptides:
                        seen_peptides.add(peptide)
                        unique_peptides.append(peptide)
                        max_peptide_len = max(max_peptide_len, len(peptide))

        return split_data, unique_peptides, max_peptide_len

    def _compute_peptides_descriptors_matrix(
        self,
        peptides: list[str],
        descriptor_names: list[str],
        padding_len: int,
    ) -> np.ndarray:
        """Compute descriptors for a list of peptides."""
        X_batches: list[np.ndarray] = []
        nb_total_peptides = len(peptides)

        if nb_total_peptides == 0:
            return np.empty((0, padding_len, len(descriptor_names)), dtype=np.float32)

        # Compute in batches and with cache for efficiency
        effective_batch_size = (
            nb_total_peptides if self._batch_size is None else self._batch_size
        )
        for start in range(0, nb_total_peptides, effective_batch_size):
            batch_start = time.time()
            end = min(start + effective_batch_size, nb_total_peptides)
            batch_peptides = peptides[start:end]

            batch_triplets: set[tuple[str | None, str, str | None]] = set()
            for peptide in batch_peptides:
                batch_triplets.update(self._gen_res_triplets(peptide))

            missing_triplets = [
                triplet for triplet in batch_triplets if triplet not in self._cache
            ]
            if missing_triplets:
                print(
                    f"Precomputing triplets for batch {start}:{end} -> {len(missing_triplets)} new triplets",
                    flush=True,
                )
                for (
                    triplet,
                    descriptors_vector,
                ) in self._compute_descriptors_vectors_for_triplets(
                    res_triplets=missing_triplets,
                    descriptor_names=descriptor_names,
                ):
                    self._cache[triplet] = descriptors_vector

            X_batch = np.asarray(
                [
                    self.get_peptide_descriptors_vectors(
                        peptide=peptide,
                        descriptor_names=descriptor_names,
                        padding_len=padding_len,
                        cache_only=True,
                    )
                    for peptide in batch_peptides
                ],
                dtype=np.float32,
            )
            X_batches.append(X_batch)
            print(
                f"Computed batch {start}:{end} | shape={X_batch.shape} | cache_size={len(self._cache)} "
                f"| elapsed={format_duration(time.time() - batch_start)}",
                flush=True,
            )

        if len(X_batches) == 1:
            return X_batches[0]
        return np.concatenate(X_batches, axis=0)

    def _save_peptides_store(
        self,
        peptides: list[str],
        X: np.ndarray,
        descriptor_names: list[str],
        padding_len: int,
    ) -> None:
        """Save the peptide sequences, descriptor matrix and names, and padding length to the storing file."""
        self._npz_store_path.parent.mkdir(parents=True, exist_ok=True)
        with lock(self._npz_store_path):
            np.savez_compressed(
                self._npz_store_path,
                sequences=np.asarray(peptides),
                descriptor_names=np.asarray(descriptor_names),
                padding_len=np.int32(padding_len),
                X=X,
            )

    def load_peptides_store(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """Return the peptide sequences, descriptor matrix and names, and padding length from the storing file."""
        data = np.load(self._npz_store_path, allow_pickle=False)
        return (
            data["sequences"],
            data["X"],
            data["descriptor_names"],
            int(data["padding_len"]),
        )

    def _gen_setup_splits_files(
        self,
        split_data: dict[tuple[int, str], tuple[np.ndarray, np.ndarray]],
        stored_peptides: np.ndarray,
        stored_X: np.ndarray,
        descriptor_names: list[str],
    ) -> None:
        """Generate the splits files for each setup fetching the stored dataset descriptors."""
        self._output_dir.mkdir(parents=True, exist_ok=True)
        peptide_to_idx = {
            peptide: idx for idx, peptide in enumerate(stored_peptides.tolist())
        }

        for setup_idx in self._setups_indices:
            for split in self._splits_names:
                peptides, y = split_data[(setup_idx, split)]
                indices = np.asarray(
                    [peptide_to_idx[peptide] for peptide in peptides], dtype=np.int64
                )
                X = stored_X[indices].astype(np.float32, copy=False)
                output_path = (
                    self._output_dir / f"descriptors_setup-{setup_idx}_{split}.npz"
                )
                np.savez_compressed(
                    output_path,
                    sequences=peptides,
                    descriptor_names=np.asarray(descriptor_names),
                    X=X,
                    y=y,
                )
                print(
                    f"Saved split file: {output_path} | X {X.shape} y {y.shape}",
                    flush=True,
                )

    def _check_store_structure(
        self,
        stored_descriptor_names: np.ndarray,
        expected_descriptor_names: list[str],
        stored_padding_len: int,
        expected_padding_len: int,
    ) -> None:
        if tuple(stored_descriptor_names.tolist()) != tuple(expected_descriptor_names):
            raise ValueError(
                "Descriptor store descriptor_names do not match requested descriptor_names. "
                "Use --rebuild-store to regenerate."
            )
        if stored_padding_len != expected_padding_len:
            raise ValueError(
                f"Descriptor store padding={stored_padding_len} "
                f"does not match requested padding={expected_padding_len}. "
                "Use --rebuild-store to regenerate."
            )

    def _estimate_run_ram_gb(
        self,
        nb_peptides: int,
        padding_len: int,
        nb_descriptors: int,
    ) -> tuple[float, float]:
        """Calculate estimated RAM usage in GB for the full descriptors matrix and for a single batch."""
        bytes_per_value = np.dtype(np.float32).itemsize
        full_matrix_bytes = nb_peptides * padding_len * nb_descriptors * bytes_per_value
        effective_batch_size = (
            nb_peptides
            if self._batch_size is None
            else min(self._batch_size, nb_peptides)
        )
        batch_matrix_bytes = (
            effective_batch_size * padding_len * nb_descriptors * bytes_per_value
        )
        return full_matrix_bytes / (1024**3), batch_matrix_bytes / (1024**3)

    def run_dataset_pipeline(self) -> None:
        if (
            self._dataset_name is None
            or self._output_dir is None
            or self._npz_store_path is None
        ):
            raise ValueError("dataset_name is required for dataset workflow.")

        # Get only molecular descriptors
        descriptor_names = get_descriptors_names()
        descriptor_names = [
            descriptor
            for descriptor in descriptor_names
            if not any(
                substring in descriptor
                for substring in ("local", "fukui", "partial", "bond")
            )
        ]

        total_start = time.time()
        print(
            "Starting descriptor pipeline: "
            f"dataset={self._dataset_name}, workers={self._workers}, batch_size={self._batch_size}",
            flush=True,
        )

        split_data, unique_peptides, padding_len = self._collect_peptides_data()
        estimated_full_gb, estimated_batch_gb = self._estimate_run_ram_gb(
            nb_peptides=len(unique_peptides),
            padding_len=padding_len,
            nb_descriptors=len(descriptor_names),
        )
        print(
            f"Collected splits for dataset {self._dataset_name} with {len(unique_peptides)} unique peptides.",
            flush=True,
        )
        print(
            f"Estimated descriptor RAM: full_matrix~{estimated_full_gb:.2f} GB, "
            f"per_batch~{estimated_batch_gb:.2f} GB.",
            flush=True,
        )
        if estimated_full_gb > self._mem_per_cpu_gb:
            raise MemoryError(
                f"Estimated full matrix RAM {estimated_full_gb:.2f} GB exceeds "
                f"mem-per-cpu limit={self._mem_per_cpu_gb:.2f} GB. "
                "Use a smaller batch size or increase --mem-per-cpu in Slurm script."
            )

        # Query and store descriptors for all unique peptides
        if self._npz_store_path.exists() and not self._rebuild_store:
            print(
                f"Loading existing peptide descriptors store: {self._npz_store_path}",
                flush=True,
            )
            (
                stored_peptides,
                stored_X,
                stored_descriptor_names,
                stored_padding_len,
            ) = self.load_peptides_store()
            self._check_store_structure(
                stored_descriptor_names=stored_descriptor_names,
                expected_descriptor_names=descriptor_names,
                stored_padding_len=stored_padding_len,
                expected_padding_len=padding_len,
            )

            existing_set = set(stored_peptides.tolist())
            new_peptides = [
                peptide for peptide in unique_peptides if peptide not in existing_set
            ]

            if new_peptides:
                print(
                    f"Existing store missing {len(new_peptides)} peptide(s), appending.",
                    flush=True,
                )
                new_X = self._compute_peptides_descriptors_matrix(
                    peptides=new_peptides,
                    descriptor_names=descriptor_names,
                    padding_len=padding_len,
                )
                stored_peptides = np.concatenate(
                    [stored_peptides, np.asarray(new_peptides)], axis=0
                )
                stored_X = np.concatenate([stored_X, new_X], axis=0).astype(np.float32)
                self._save_peptides_store(
                    peptides=stored_peptides.tolist(),
                    X=stored_X,
                    descriptor_names=descriptor_names,
                    padding_len=padding_len,
                )
                print(
                    f"Updated peptide store written: {self._npz_store_path}",
                    flush=True,
                )
            else:
                print(
                    "Existing store already covers all required peptides.", flush=True
                )

        else:
            if self._npz_store_path.exists() and self._rebuild_store:
                print(
                    f"Rebuilding peptide store from scratch: {self._npz_store_path}",
                    flush=True,
                )
            else:
                print(f"Building peptide store: {self._npz_store_path}", flush=True)

            stored_X = self._compute_peptides_descriptors_matrix(
                peptides=unique_peptides,
                descriptor_names=descriptor_names,
                padding_len=padding_len,
            )
            self._save_peptides_store(
                peptides=unique_peptides,
                X=stored_X,
                descriptor_names=descriptor_names,
                padding_len=padding_len,
            )
            stored_peptides = np.asarray(unique_peptides)
            print(f"Peptide store written: {self._npz_store_path}", flush=True)

        # Generate setup splits files from descriptors store
        self._gen_setup_splits_files(
            split_data=split_data,
            stored_peptides=stored_peptides,
            stored_X=stored_X,
            descriptor_names=descriptor_names,
        )

        print(
            f"Done in {format_duration(time.time() - total_start)}.",
            flush=True,
        )


def _compute_res_descriptor_vector_worker(
    res_triplet: tuple[str | None, str, str | None],
    descriptor_names: tuple[str, ...],
) -> tuple[tuple[str | None, str, str | None], list[float]]:
    """Compute the descriptors vector for a residue with neighbours."""
    left_neighbour, res, right_neighbour = res_triplet
    res_descriptors = query_average(
        res,
        left_neighbour=left_neighbour,
        right_neighbour=right_neighbour,
    )
    vector = [res_descriptors[descriptor] for descriptor in descriptor_names]
    return res_triplet, vector


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build peptide-level descriptors store for entire dataset, then "
            "generate setup train/val/test npz files from that store."
        )
    )
    parser.add_argument("dataset", help="Dataset name: 'toxicity' or 'amp'")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes for residue-triplet precomputation",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=0,
        help="Number of unique peptides per compute batch (0 to disable batching)",
    )
    parser.add_argument(
        "--rebuild-store",
        action="store_true",
        help="Force full recomputation of peptide-level store",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pipeline = PeptideDescriptorsPipeline(
        dataset_name=args.dataset,
        workers=args.workers,
        batch_size=None if args.batch_size == 0 else args.batch_size,
        rebuild_store=args.rebuild_store,
    )
    pipeline.run_dataset_pipeline()


if __name__ == "__main__":
    main()
