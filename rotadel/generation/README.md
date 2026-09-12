# Generating the 𝚁𝚘𝚝𝙰𝙳𝚎𝙻 database

This describes how to rerun the pipeline to regenerate `database/rotadel.db` from scratch.

## 1. Obtain the Dunbrack libraries

Two libraries from the Dunbrack lab are required and must obtained on [their website](https://dunbrack.fccc.edu/lab/).

### Backbone-dependent rotamer library

- Source: <http://dunbrack.fccc.edu/bbdep2010>
- Place the rotamer library file in the `rotadel` repository:
```shell
$ mkdir -p data/Dunbrack_libraries/rotamer
$ cp <the obtained file> data/Dunbrack_libraries/rotamer/ALL_rotamers.lib
```

### Neighbor-dependent ramachandran distributions (NDRD)

- Source: <http://dunbrack.fccc.edu/ndrd>
- Place the `TCBIG` distribution file in the `rotadel` repository:
```shell
$ mkdir -p data/Dunbrack_libraries/ndrd
$ cp <the obtained TCBIG file> data/Dunbrack_libraries/ndrd/NDRD_TCBIG.txt
```

## 2. Regenerate the intermediate data files

Each step below depends on the output(s) of the previous ones and must be performed in order.

1. Extract angles and probabilities from Dunbrack rotamer library and generate ID for each rotamer:
  ```shell
  $ python rotadel/generation/dunbrack_to_json.py
  ```
   Reads `data/Dunbrack_libraries/rotamer/ALL_rotamers.lib` → writes `data/angles_Dunbrack.json`.

2. Get probability from NDRD library needed for `query_average`:
  ```shell
  $ python rotadel/generation/ndrd_to_csv.py
  ```
  Reads `data/Dunbrack_libraries/ndrd/NDRD_TCBIG.txt` → writes `database/ndrd_step10.csv`.

3. Create the rotamer IDs containing the charge and tautomer information:
  ```shell
  $ python rotadel/generation/ids_species.py
  ```
  Reads `data/angles_Dunbrack.json` → writes `data/keys_species_initial.csv`.

4. Add alanine/glycine entries (no side-chain dihedral angles and therefore absent from the Dunbrack rotamer
   library):
  ```shell
  $ python rotadel/generation/add_single_AAs.py
  ```
  Reads `data/angles_Dunbrack.json` → writes `data/angles_single_AAs.json` and `data/angles_combined.json`.

5. Curate the rotamer data and prepare batches to submit the descriptor calculation jobs:
   - Run `prepare_keys.ipynb` top to bottom.
   
   Reads `data/angles_Dunbrack.json` and `data/keys_species_initial.csv` → write `data/keys_species_cleaned.csv`
   and, in the last cell, splits it into `data/batches_keys_species_6490_cleaned/batch_*`.

## 3. Run the per-rotamer descriptor calculation

Submit the batches produced above as SLURM jobs from `jobs/`:

- `sbatch jobs/submit_batches.sh` submits one `jobs/submit_series.sh` job per batch file. Each one runs
  `rotadel.generation.main_one_rotamer` on every rotamer ID in its batch, using
  `data/angles_combined.json` as `--angles_json`, and writes one `.db` file per
  batch under `data/batches_sql_output/`.

## 4. Merge into the final database

```shell
$ python rotadel/generation/merge_sql.py data/batches_sql_output -o database/rotadel.db
```
