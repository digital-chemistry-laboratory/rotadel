#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH -A es_jorner
#SBATCH -t 4:00:00
#SBATCH -J run_rot
#SBATCH -o "run_rot.out"
#SBATCH -e "run_rot.error"

# Set above the SLURM variables

start_time=$(date +%s)

# Load the required modules and environment
mamba activate aa

# Set the required environment variables
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export PYTHONPATH=/cluster/project/jorner/lajacot/miniforge3/envs/aa/lib/python3.13/site-packages:$PYTHONPATH
export PYTHONPATH=/cluster/project/jorner/lajacot/projects/aa-descriptors-library:$PYTHONPATH

# Do not change the two following lines
. `which env_parallel.bash`
paste -d '/' <(perl -pe 's/(\d+)\(x(\d+)\)/substr("$1,"x$2,0,-1)/ge' <<<$SLURM_TASKS_PER_NODE | tr ',' '\n') \
             <(scontrol show hostnames) > ./node_list_${SLURM_JOB_ID}

# Run the parallel job
input_file="data/rotamer_names.txt"
zip_file="data/xyz_folder.zip"
angles_json="data/angles_Dunbrack.json"
# Specify with --env which environment variables to export
env_parallel \
-a ${input_file} \
--env PATH \
--env PYTHONPATH \
--env LD_LIBRARY_PATH \
--joblog parallel.log \
--wd $PWD \
--jobs ${SLURM_NTASKS} \
"python main_one_rotamer.py {} ${angles_json} ${zip_file}"
#--resume-failed \

# Calculate run times
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total execution time: $duration seconds"