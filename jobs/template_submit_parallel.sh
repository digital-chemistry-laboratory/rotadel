#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH -A es_jorner
#SBATCH -t 24:00:00
#SBATCH -J run_rot

# Set above the SLURM variables

start_time=$(date +%s)

# Load the required modules and environment
mamba activate aa
export MY_MODULEPATH_ROOT=/cluster/project/jorner/modules_dcl/
module use $MY_MODULEPATH_ROOT/Core
module load xtb/bleed

# Set the required environment variables
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export PYTHONPATH=/cluster/project/jorner/lajacot/miniforge3/envs/aa/lib/python3.11/site-packages:$PYTHONPATH
export PYTHONPATH=/cluster/project/jorner/lajacot/projects/aa-descriptors-library:$PYTHONPATH
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Finds path for env_parallel command
. `which env_parallel.bash`
# paste -d '/' <(perl -pe 's/(\d+)\(x(\d+)\)/substr("$1,"x$2,0,-1)/ge' <<<$SLURM_TASKS_PER_NODE | tr ',' '\n') \
#              <(scontrol show hostnames) > ./node_list_${SLURM_JOB_ID}

# Variables for the job
input_file="${BATCH_FILE}" # Exported from the submit_batches.sh script
batch_number=$(basename "${input_file}" | sed 's/\.[^.]*$//' | awk -F'_' '{print $NF}')
angles_json="../data/angles_Dunbrack.json"
output_path="../data/batches_sql_output/rotamer_descriptors_${batch_number}.db"
# Prevent output file already existing from a previous run
if [ -f "${output_path}" ]; then
    echo "ERROR: File ${output_path} already exists. Provide a different output path." >&2
    exit 1
fi
# Run the parallel job
# Specify with --env which environment variables to export
env_parallel \
-a ${input_file} \
--env PATH \
--env PYTHONPATH \
--env LD_LIBRARY_PATH \
--env OMP_NUM_THREADS \
--env OPENBLAS_NUM_THREADS \
--env MKL_NUM_THREADS \
--env NUMEXPR_NUM_THREADS \
--env VECLIB_MAXIMUM_THREADS \
--joblog parallel_${batch_number}.log \
--wd $PWD \
--jobs ${SLURM_NTASKS} \
--timeout 300 \
"python -m aa_descriptors_library.main_one_rotamer {} ${angles_json} -o ${output_path}"
# --resume-failed \

# Calculate run times
end_time=$(date +%s)
duration=$((end_time - start_time))
echo -e "\nTotal execution time: $((duration/3600))h $(((duration%3600)/60))m $((duration%60))s"