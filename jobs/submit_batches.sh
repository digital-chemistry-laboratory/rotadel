#!/bin/bash

# Script to submit the submitting script for each dataset batch

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

for batch_file in ../data/batches_keys_species_6490_cleaned/batch_*; do
  sbatch --export=BATCH_FILE=${batch_file} submit_series.sh
done