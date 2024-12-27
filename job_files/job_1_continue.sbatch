#!/bin/bash
#SBATCH --account=group-annelohfink
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=102
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=xrt-simulation-job1-continue
#SBATCH --output=/home/r77m975/xrb-population/out_logs/%x/%A_%a.out
#SBATCH --error=/home/r77m975/xrb-population/out_logs/%x/%A_%a.err
#SBATCH --mail-user=youssefabdulghani@montana.edu
#SBATCH --mail-type=ALL
#SBATCH --array=243-288  # Specify job array range for combinations

# Create a unique folder for the job logs
log_dir="/home/r77m975/xrb-population/out_logs/${SLURM_JOB_NAME}"
mkdir -p "$log_dir"

echo "Running on node: $(hostname)"

# Define the combinations as parameters
combinations=(
    "1.7 0.7 0.0 10.0 80.0 0.9 5000 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.2 400.0 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.2 1000 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.2 5000 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.5 400.0 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.5 1000 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.5 5000 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.9 400.0 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.9 1000 xrt"
    "1.7 0.7 0.998 6.0 0.0 0.9 5000 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.2 400.0 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.2 1000 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.2 5000 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.5 400.0 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.5 1000 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.5 5000 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.9 400.0 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.9 1000 xrt"
    "1.7 0.7 0.998 6.0 60.0 0.9 5000 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.2 400.0 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.2 1000 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.2 5000 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.5 400.0 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.5 1000 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.5 5000 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.9 400.0 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.9 1000 xrt"
    "1.7 0.7 0.998 6.0 80.0 0.9 5000 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.2 400.0 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.2 1000 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.2 5000 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.5 400.0 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.5 1000 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.5 5000 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.9 400.0 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.9 1000 xrt"
    "1.7 0.7 0.998 8.0 0.0 0.9 5000 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.2 400.0 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.2 1000 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.2 5000 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.5 400.0 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.5 1000 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.5 5000 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.9 400.0 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.9 1000 xrt"
    "1.7 0.7 0.998 8.0 60.0 0.9 5000 xrt"
)

# Adjust array index to parameter index
index=$((SLURM_ARRAY_TASK_ID - 243))

# Fetch parameters for this job
params="${combinations[index]}"

# Progress tracking
progress_file="$log_dir/progress.log"
echo "Job ${SLURM_JOB_ID}, Task ${SLURM_ARRAY_TASK_ID}: Running combination $SLURM_ARRAY_TASK_ID: $params" >> $progress_file

# Execute the simulation
apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py $params