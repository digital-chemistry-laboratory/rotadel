#!/bin/bash
#SBATCH --mem-per-cpu=2GB
#SBATCH -A es_jorner
#SBATCH -t 48:00:00
#SBATCH -J series_rot

start_time=$(date +%s)

# Set up environment
source activate /cluster/project/jorner/lajacot/miniforge3/envs/aa
export MY_MODULEPATH_ROOT=/cluster/project/jorner/modules_dcl/
module use $MY_MODULEPATH_ROOT/Core
module load xtb/bleed
export PYTHONPATH=/cluster/project/jorner/lajacot/miniforge3/envs/aa/lib/python3.11/site-packages:$PYTHONPATH
export PYTHONPATH=/cluster/project/jorner/lajacot/projects/aa-descriptors-library:$PYTHONPATH

# Files paths
angles_json="../data/angles_combined.json"
input_file="${BATCH_FILE}"
batch_number=$(basename "${input_file}" | sed 's/\.[^.]*$//' | awk -F'_' '{print $NF}')
output_path="../data/batches_sql_output/rotamers_sql_${batch_number}.db"

# Run
echo "--- Run for ${input_file} ---"
while IFS= read -r line || [ -n "$line" ]; do
    python -m aa_descriptors_library.main_one_rotamer $line -a ${angles_json} -o ${output_path} --rerun_failed
done < ${input_file}

# Calculate total run time
end_time=$(date +%s)
duration=$((end_time - start_time))
echo -e "\nTotal execution time: $((duration/3600))h $(((duration%3600)/60))m $((duration%60))s"