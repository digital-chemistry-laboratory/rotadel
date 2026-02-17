import datetime
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import tensorflow as tf

from benchmark.descriptors import get_descriptors_vectors_peptide
from peptidy_zenodo.library.models import get_predictor
from peptidy_zenodo.library.training import train_predictor
from peptidy_zenodo.library.training.hpspaces import get_sampled_hp_space
from peptidy_zenodo.library.training.training import evaluate_predictor


class ModelRunner:
    def __init__(
        self,
        # encoding_name: str,  # ref: one from("label" "onehot" "blosum62" "aa_descriptor")
        model_name: str,  # ref: "cnn"
        dataset_name: str,  # ref: one from ("amp", "toxicity")
        # n_trials: int = None,
        # start_trial_idx: int = None,
    ) -> None:
        # self.encoding_name = encoding_name
        self.model_name = model_name
        self.dataset_name = dataset_name
        # self.n_trials = n_trials
        # self.start_trial_idx = start_trial_idx

        self.base_saving_dir = Path("results")

    def run(self):
        hp_combinations = get_sampled_hp_space(self.model_name)
        # if self.n_trials is not None and self.start_trial_idx is not None:
        #     hp_combinations = hp_combinations[
        #         self.start_trial_idx : self.start_trial_idx + self.n_trials
        #     ]
        results_path = self.base_saving_dir / f"{self.dataset_name}.json"
        if results_path.exists():
            with open(results_path, "r") as results_file:
                done_combinations = json.load(results_file)
        else:
            results_path.parent.mkdir(parents=True, exist_ok=True)
            done_combinations = dict()
        start = time.time()
        for combination_idx, hp_combination in enumerate(hp_combinations):
            if f"combination-{combination_idx}" in done_combinations:
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
                    # self.encoding_name,
                    setup_idx,
                    "train",
                )
                X_val, y_val = get_X_y(
                    self.dataset_name,
                    # self.encoding_name,
                    setup_idx,
                    "val",
                )
                X_test, y_test = get_X_y(
                    self.dataset_name,
                    # self.encoding_name,
                    setup_idx,
                    "test",
                )
                add_embedding_layer = False
                # if self.encoding_name == "label":
                #     add_embedding_layer = True
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

            if combination_idx + 1 % 5 == 0:
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

            done_combinations[f"combination-{combination_idx}"] = combination_dump

            with open(results_path, "w") as results_file:
                json.dump(done_combinations, results_file, indent=4)


def get_X_y(
    dataset_name: str,
    setup_idx: int,
    set: str,
):
    df = pd.read_csv(
        f"../peptidy_zenodo/data/{dataset_name}/setup-{setup_idx}/{set}.csv"
    )
    sequences = df["sequence"].tolist()
    y = np.array(df["label"].tolist())
    X = [
        get_descriptors_vectors_peptide(sequence, padding_len=50)
        for sequence in sequences
    ]
    return X, y


if __name__ == "__main__":
    dataset_names = ["amp", "toxicity"]
    for dataset_name in dataset_names:
        model = ModelRunner(
            model_name="cnn",
            dataset_name=dataset_name,
        )
        model.run()
