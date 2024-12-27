#!/bin/bash
#SBATCH --account=group-annelohfink
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=102
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=xrt-simulation-job-6
#SBATCH --output=/home/r77m975/xrb-population/out_logs/job_6/job-6-%A_%a.out
#SBATCH --error=/home/r77m975/xrb-population/out_logs/job_6/job-6-%A_%a.err
#SBATCH --mail-user=youssefabdulghani@montana.edu
#SBATCH --mail-type=ALL
#SBATCH --array=0-17

echo "Running on node: $(hostname)"

# Progress tracking setup
progress_file="/home/r77m975/xrb-population/out_logs/job_6/progress.log"
touch $progress_file

# Execute commands

    param_index=$SLURM_ARRAY_TASK_ID
    g=$(echo "3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0" | cut -d',' -f$((param_index+1)))
    T=$(echo "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0" | cut -d',' -f$((param_index+1)))
    a=$(echo "0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998,0.998" | cut -d',' -f$((param_index+1)))
    m=$(echo "10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0" | cut -d',' -f$((param_index+1)))
    i=$(echo "60.0,60.0,60.0,60.0,60.0,60.0,60.0,60.0,60.0,80.0,80.0,80.0,80.0,80.0,80.0,80.0,80.0,80.0" | cut -d',' -f$((param_index+1)))
    r=$(echo "0.2,0.2,0.2,0.5,0.5,0.5,0.9,0.9,0.9,0.2,0.2,0.2,0.5,0.5,0.5,0.9,0.9,0.9" | cut -d',' -f$((param_index+1)))
    e=$(echo "400.0,1000,5000,400.0,1000,5000,400.0,1000,5000,400.0,1000,5000,400.0,1000,5000,400.0,1000,5000" | cut -d',' -f$((param_index+1)))

    # Log progress to the unique progress file
    params="g=$g, T=$T, a=$a, m=$m, i=$i, r=$r, e=$e"
    echo "Job $SLURM_JOB_ID, Task $SLURM_ARRAY_TASK_ID: Running combination $param_index: $params" >> $progress_file

    # Execute the command
    apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python /home/r77m975/xrb-population/observational_effects.py $g $T $a $m $i $r $e xrt
    
