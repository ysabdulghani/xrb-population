#!/bin/bash
#SBATCH --account=group-annelohfink
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=102
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=xrt-simulation-job-3-continue
#SBATCH --output=/home/r77m975/xrb-population/out_logs/%x-%j.out
#SBATCH --error=/home/r77m975/xrb-population/out_logs/%x-%j.err
#SBATCH --mail-user=youssefabdulghani@montana.edu
#SBATCH --mail-type=ALL

echo "Running on node: $(hostname)"

# Progress tracking setup
progress_file="/home/r77m975/xrb-population/progress-${SLURM_JOB_NAME}.log"
touch "$progress_file"

# Simulation combinations
echo 'Running combination 278 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 60.0 0.9 1000 xrt
echo 'Running combination 279 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 60.0 0.9 5000 xrt
echo 'Running combination 280 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.2 400.0 xrt
echo 'Running combination 281 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.2 1000 xrt
echo 'Running combination 282 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.2 5000 xrt
echo 'Running combination 283 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.5 400.0 xrt
echo 'Running combination 284 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.5 1000 xrt
echo 'Running combination 285 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.5 5000 xrt
echo 'Running combination 286 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.9 400.0 xrt
echo 'Running combination 287 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.9 1000 xrt
echo 'Running combination 288 of 288 for job 3' >> "$progress_file" && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 2.0 1.0 0.0 8.0 80.0 0.9 5000 xrt