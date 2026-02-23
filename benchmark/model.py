import datetime
import json
import argparse
import gc
import functools
from pathlib import Path
import time

import numpy as np
import pandas as pd
import tensorflow as tf

from aa_descriptors_library.sql import get_descriptors_names
from benchmark.descriptors import get_descriptors_vectors_peptide
from peptidy_zenodo.library.models import get_predictor
from peptidy_zenodo.library.training import train_predictor
from peptidy_zenodo.library.training.hpspaces import get_sampled_hp_space
from peptidy_zenodo.library.training.training import evaluate_predictor


@functools.lru_cache(maxsize=1)
def get_filtered_descriptor_names() -> tuple[str, ...]:
    # For now, remove the atomic descriptors
    descriptor_names = get_descriptors_names()
    filtered = [
        descriptor
        for descriptor in descriptor_names
        if not any(
            substring in descriptor
            for substring in ("local", "fukui", "partial", "bond")
        )
    ]
    return tuple(filtered)


def configure_tensorflow_runtime(require_gpu: bool = True) -> None:
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print(f"Using {len(gpus)} GPU(s): {[gpu.name for gpu in gpus]}")
    else:
        message = "No GPU detected."
        if require_gpu:
            raise RuntimeError(f"{message} Exiting because GPU is required.")
        print(f"{message} Falling back to CPU.")


class ModelRunner:
    def __init__(
        self,
        model_name: str,
        dataset_name: str,
        n_combinations: int | None = None,
        start_combination_idx: int = 0,
        results_dir: str = "results",
    ) -> None:
        self.model_name = model_name
        self.dataset_name = dataset_name
        self.n_combinations = n_combinations
        self.start_combination_idx = start_combination_idx

        self.base_saving_dir = Path(results_dir)

    def run(self):
        hp_combinations_all = get_sampled_hp_space(self.model_name)
        if self.n_combinations is not None:
            hp_combinations = hp_combinations_all[
                self.start_combination_idx : self.start_combination_idx
                + self.n_combinations
            ]
        else:
            hp_combinations = hp_combinations_all[self.start_combination_idx :]

        print(
            f"Loaded {len(hp_combinations)} combinations "
            f"(start={self.start_combination_idx}, n_combinations={self.n_combinations})."
        )

        if self.n_combinations is not None:
            end_combination_idx = self.start_combination_idx + self.n_combinations - 1
            results_path = (
                self.base_saving_dir
                / f"{self.dataset_name}_batch_{self.start_combination_idx}_{end_combination_idx}.json"
            )
        else:
            results_path = self.base_saving_dir / f"{self.dataset_name}.json"
        print(f"Writing results to: {results_path}")
        if results_path.exists():
            with open(results_path, "r") as results_file:
                done_combinations = json.load(results_file)
        else:
            results_path.parent.mkdir(parents=True, exist_ok=True)
            done_combinations = dict()

        descriptor_names = list(get_filtered_descriptor_names())
        descriptor_cache: dict[tuple[str | None, str, str | None], list[float]] = {}

        start = time.time()
        for local_idx, hp_combination in enumerate(hp_combinations):
            combination_idx = hp_combination.get(
                "combination_idx",
                self.start_combination_idx + local_idx,
            )
            combination_key = f"combination-{combination_idx}"

            if combination_key in done_combinations:
                print(
                    f"Combination {combination_idx} already exists, skipping",
                )
                continue

            val_scores = list()
            test_scores = list()
            for setup_idx in range(5):
                print(
                    f"Running setup {setup_idx} of combination {combination_idx} / {len(hp_combinations)}"
                )
                tf.keras.backend.clear_session()
                tf.random.set_seed(42)

                X_train, y_train = get_X_y(
                    self.dataset_name,
                    setup_idx,
                    "train",
                    descriptor_names,
                    descriptor_cache,
                )
                X_val, y_val = get_X_y(
                    self.dataset_name,
                    setup_idx,
                    "val",
                    descriptor_names,
                    descriptor_cache,
                )
                X_test, y_test = get_X_y(
                    self.dataset_name,
                    setup_idx,
                    "test",
                    descriptor_names,
                    descriptor_cache,
                )
                add_embedding_layer = False
                predictor = get_predictor(
                    self.model_name,
                    hp_combination,
                    add_embedding_layer=add_embedding_layer,
                )
                history = train_predictor(
                    predictor,
                    X_train,
                    y_train,
                    X_val,
                    y_val,
                    hp_combination["learning_rate"],
                    hp_combination["batch_size"],
                )
                test_scores.append(evaluate_predictor(predictor, X_test, y_test))
                val_scores.append(np.min(history["val_loss"]))
                del predictor
                del X_train, y_train, X_val, y_val, X_test, y_test
                gc.collect()

            if (local_idx + 1) % 5 == 0:
                elapsed_time = time.time() - start
                print(
                    "Combination idx:",
                    combination_idx,
                    " Total time elapsed:",
                    str(datetime.timedelta(seconds=(elapsed_time))),
                )

            combination_dump = dict()
            combination_dump["hps"] = hp_combination
            combination_dump["misc"] = dict()
            combination_dump["misc"]["val_loss_mean"] = np.mean(val_scores)
            combination_dump["misc"]["val_loss_std"] = np.std(val_scores)
            for metric in test_scores[0].keys():
                combination_dump["misc"][f"test_{metric}_mean"] = np.mean(
                    [test_score[metric] for test_score in test_scores]
                )
                combination_dump["misc"][f"test_{metric}_std"] = np.std(
                    [test_score[metric] for test_score in test_scores]
                )

            combination_dump["misc"]["test_scores"] = {
                "setup-{}".format(setup_idx): test_score
                for setup_idx, test_score in enumerate(test_scores)
            }

            done_combinations[combination_key] = combination_dump

            with open(results_path, "w") as results_file:
                json.dump(done_combinations, results_file, indent=4)


def get_X_y(
    dataset_name: str,
    setup_idx: int,
    set: str,
    descriptor_names: list[str],
    descriptor_cache: dict[tuple[str | None, str, str | None], list[float]],
):
    data_path = (
        Path(__file__).resolve().parents[1]
        / "peptidy_zenodo"
        / "data"
        / dataset_name
        / f"setup-{setup_idx}"
        / f"{set}.csv"
    )
    df = pd.read_csv(data_path)
    sequences = df["sequence"].tolist()
    y = np.asarray(df["label"].tolist(), dtype=np.float32)
    X = np.asarray(
        [
            get_descriptors_vectors_peptide(
                sequence,
                descriptor_names=descriptor_names,
                padding_len=50,
                cache=descriptor_cache,
            )
            for sequence in sequences
        ],
        dtype=np.float32,
    )
    return X, y


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="cnn")
    parser.add_argument("--dataset", default="toxicity")
    parser.add_argument("--n-combinations", type=int, default=None)
    parser.add_argument("--start-combination-idx", type=int, default=0)
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--allow-cpu", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    configure_tensorflow_runtime(require_gpu=not args.allow_cpu)
    model = ModelRunner(
        model_name=args.model,
        dataset_name=args.dataset,
        n_combinations=args.n_combinations,
        start_combination_idx=args.start_combination_idx,
        results_dir=args.results_dir,
    )
    model.run()
