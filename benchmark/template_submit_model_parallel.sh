#!/bin/bash
set -euo pipefail

# Workflow variables
DATASET="toxicity"
MODEL="cnn"
NB_COMBINATIONS_PER_BATCH=12
TOTAL_NB_COMBINATIONS=108 # number of hyperparameters combinations defined in peptidy_zenodo/library/training/spaces/cnn.json
CWD="/cluster/project/jorner/lajacot/projects/aa-descriptors-library/benchmark"
RESULTS_DIR="${CWD}/results"
BATCHES_RESULTS_DIR="${RESULTS_DIR}/batches_combinations_scores_${DATASET}"

cd "${CWD}"

# Create results folder, exits if exists
if [[ -e "${BATCHES_RESULTS_DIR}" ]]; then
    echo "Error: ${BATCHES_RESULTS_DIR} already exists. Exiting to avoid overwriting."
    exit 1
fi
mkdir -p "${BATCHES_RESULTS_DIR}"

echo "Submitting batched jobs:"
echo "  dataset=${DATASET}"
echo "  model=${MODEL}"
echo "  nb_combinations_per_batch=${NB_COMBINATIONS_PER_BATCH}"
echo "  total_nb_combinations=${TOTAL_NB_COMBINATIONS}"
echo "  batch_results_dir=${BATCHES_RESULTS_DIR}"

job_ids=()

# Submit the model in parallel
# grouping into batches containing NB_COMBINATIONS_PER_BATCH of hyperparameters combinations
for ((start=0; start<TOTAL_NB_COMBINATIONS; start+=NB_COMBINATIONS_PER_BATCH)); do
    n_combinations=${NB_COMBINATIONS_PER_BATCH}
    end=$((start + n_combinations - 1))
    echo "Submitting batch from combination ${start} to ${end}"
    sbatch_output=$(sbatch "${CWD}/submit_model.sh" \
        --dataset "${DATASET}" \
        --model "${MODEL}" \
        --start-combination-idx "${start}" \
        --n-combinations "${n_combinations}" \
        --results-dir "${BATCHES_RESULTS_DIR}")
    job_id=$(echo "${sbatch_output}" | awk '{print $4}')
    job_ids+=("${job_id}")
    echo "  -> job_id=${job_id}"
done

# Make sure the previous jobs ran succesfully...
dependency=$(IFS=: ; echo "${job_ids[*]}")

# ... then merge the batches output json files into a single json file
merge_json_cmd="cd ${CWD} && \
source /cluster/project/jorner/lajacot/miniforge3/etc/profile.d/conda.sh && \
conda activate /cluster/project/jorner/lajacot/miniforge3/envs/aa && \
export PYTHONPATH=/cluster/project/jorner/lajacot/projects/aa-descriptors-library:\$PYTHONPATH && \
python -c \"import glob; from aa_descriptors_library.process_json import merge_json_files; files=sorted(glob.glob('${BATCHES_RESULTS_DIR}/${DATASET}_batch_*.json')); assert files, 'No batch files found'; merge_json_files(files, '${RESULTS_DIR}/${DATASET}_scores.json'); print('Merged', len(files), 'files into ${RESULTS_DIR}/${DATASET}_scores.json')\""

merge_output=$(sbatch \
    --dependency="afterok:${dependency}" \
    -A es_jorner \
    -t 00:30:00 \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=2GB \
    -J "merge_${DATASET}" \
    --wrap "${merge_json_cmd}")
merge_job_id=$(echo "${merge_output}" | awk '{print $4}')

echo "All batch jobs submitted."
echo "Merge job submitted with dependency afterok:${dependency}"
echo "Merge job_id=${merge_job_id}"
