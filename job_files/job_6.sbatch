#!/bin/bash
#SBATCH --account=group-annelohfink
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=102
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=xrt-simulation-job-6
#SBATCH --output=/home/r77m975/xrb-population/out_logs/job-6-%j.out
#SBATCH --error=/home/r77m975/xrb-population/out_logs/job-6-%j.err
#SBATCH --mail-user=youssefabdulghani@montana.edu
#SBATCH --mail-type=ALL

echo "Running on node: $(hostname)"

# Progress tracking setup
progress_file="/home/r77m975/xrb-population/progress.log"
touch $progress_file

echo 'Running combination 1 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.2 400.0 xrt
echo 'Running combination 2 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.2 1000 xrt
echo 'Running combination 3 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.2 5000 xrt
echo 'Running combination 4 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.5 400.0 xrt
echo 'Running combination 5 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.5 1000 xrt
echo 'Running combination 6 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.5 5000 xrt
echo 'Running combination 7 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.9 400.0 xrt
echo 'Running combination 8 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.9 1000 xrt
echo 'Running combination 9 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 60.0 0.9 5000 xrt
echo 'Running combination 10 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.2 400.0 xrt
echo 'Running combination 11 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.2 1000 xrt
echo 'Running combination 12 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.2 5000 xrt
echo 'Running combination 13 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.5 400.0 xrt
echo 'Running combination 14 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.5 1000 xrt
echo 'Running combination 15 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.5 5000 xrt
echo 'Running combination 16 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.9 400.0 xrt
echo 'Running combination 17 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.9 1000 xrt
echo 'Running combination 18 of 18 for job 6' >> /home/r77m975/xrb-population/progress.log && apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py 3.0 1.0 0.998 10.0 80.0 0.9 5000 xrt
