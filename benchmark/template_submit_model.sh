#!/bin/bash
#SBATCH --mem-per-cpu=16GB
#SBATCH --gpus=1
#SBATCH -A es_jorner
#SBATCH -t 24:00:00
#SBATCH -J benchmark

set -euo pipefail

start_time=$(date +%s)

# Set up environment
source /cluster/project/jorner/lajacot/miniforge3/etc/profile.d/conda.sh
conda activate /cluster/project/jorner/lajacot/miniforge3/envs/aa
export PYTHONPATH=/cluster/project/jorner/lajacot/projects/aa-descriptors-library:$PYTHONPATH

# Run model
python -u model.py "$@"

# Calculate total run time
end_time=$(date +%s)
duration=$((end_time - start_time))
echo -e "\nTotal execution time: $((duration/3600))h $(((duration%3600)/60))m $((duration%60))s"
