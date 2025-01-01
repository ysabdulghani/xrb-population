import os
import math
from itertools import product
import argparse

parser = argparse.ArgumentParser(
                    prog='sbatch_generation_arrays',
                    description='Generate job files for submission to SLURM', 
)

parser.add_argument('insturment')
parser.add_argument('partition') 
args = parser.parse_args()

insturment = args.insturment
partition = args.partition

if partition == 'priority':
    account = 'priority-annelohfink'
else:
    account = 'group-annelohfink'

# Parameters
preset_values = {
    'g': [1.7, 2.0, 3.0],
    'T': [0.5, 0.7, 1.0],
    'a': [0.0, 0.998],
    'm': [6.0, 8.0, 10.0],
    'i': [0.0, 60.0, 80.0],
    'r': [0.2, 0.5, 0.9],
    'e': [400.0, 1500.0, 10000.0]
}

# Configuration
email = "youssef.abdulghani@student.montana.edu"
cpus_per_task = 96
script_path = "/home/r77m975/xrb-population/observational_effects.py"
output_dir = f"/home/r77m975/xrb-population/job_files_{insturment}"
outlog_base_dir = f"/home/r77m975/xrb-population/out_logs_{insturment}"
time_per_combination = 600  # seconds
max_job_duration = 2 * 24 * 60 * 60  # 2 days in seconds

# Calculate combinations and jobs
combinations = math.prod(len(values) for values in preset_values.values())
combinations_per_job = max_job_duration // time_per_combination
num_jobs = math.ceil(combinations / combinations_per_job)

# Generate all parameter combinations
parameters = list(preset_values.keys())
values_list = [preset_values[key] for key in parameters]
all_combinations = list(product(*values_list))

# Divide combinations into jobs
combinations_list = [
    all_combinations[i * combinations_per_job: (i + 1) * combinations_per_job]
    for i in range(num_jobs)
]

# sbatch template
sbatch_template = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --job-name={insturment}-simulation-job-{job_id}
#SBATCH --output={outlog_path}/job-{job_id}-%A_%a-%j.out
#SBATCH --error={outlog_path}/job-{job_id}-%A_%a-%j.err
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
#SBATCH --array=0-{array_size}

echo "Running on node: $(hostname)"

# Progress tracking setup
progress_file="{progress_log_path}"
touch $progress_file

# Execute commands
{commands}
"""

# Create output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(outlog_base_dir, exist_ok=True)

for job_id, job_combinations in enumerate(combinations_list, start=1):
    # Output and progress log setup
    job_outlog_path = os.path.join(outlog_base_dir, f"job_{job_id}")
    os.makedirs(job_outlog_path, exist_ok=True)
    job_progress_log_path = os.path.join(job_outlog_path, "progress.log")

    array_size = len(job_combinations) - 1

    commands = f"""
    param_index=$SLURM_ARRAY_TASK_ID
    g=$(echo "{','.join(str(g) for g, *_ in job_combinations)}" | cut -d',' -f$((param_index+1)))
    T=$(echo "{','.join(str(T) for _, T, *_ in job_combinations)}" | cut -d',' -f$((param_index+1)))
    a=$(echo "{','.join(str(a) for _, _, a, *_ in job_combinations)}" | cut -d',' -f$((param_index+1)))
    m=$(echo "{','.join(str(m) for _, _, _, m, *_ in job_combinations)}" | cut -d',' -f$((param_index+1)))
    i=$(echo "{','.join(str(i) for _, _, _, _, i, *_ in job_combinations)}" | cut -d',' -f$((param_index+1)))
    r=$(echo "{','.join(str(r) for _, _, _, _, _, r, *_ in job_combinations)}" | cut -d',' -f$((param_index+1)))
    e=$(echo "{','.join(str(e) for _, _, _, _, _, _, e in job_combinations)}" | cut -d',' -f$((param_index+1)))

    # Log progress to the unique progress file
    params="g=$g, T=$T, a=$a, m=$m, i=$i, r=$r, e=$e"
    echo "$(date), Job $SLURM_JOB_ID, Task $SLURM_ARRAY_TASK_ID: Running combination $param_index: $params" >> $progress_file

    # Execute the command
    apptainer exec -c --mount 'type=bind,source=/home/r77m975/xrb-population/,destination=/home/r77m975/xrb-population' --pwd /home/r77m975/xrb-population/ /home/r77m975/heasoft_docker/heasoft-v6.34.sif python {script_path} $g $T $a $m $i $r $e {insturment}
    """

    sbatch_content = sbatch_template.format(
        account=account,
        partition=partition,
        cpus_per_task=cpus_per_task,
        job_id=job_id,
        outlog_path=job_outlog_path,
        email=email,
        progress_log_path=job_progress_log_path,
        array_size=array_size,
        insturment=insturment,
        commands=commands
    )

    # Save to file
    filename = os.path.join(output_dir, f"job_{job_id}_arrays_{insturment}.sbatch")
    with open(filename, "w") as file:
        file.write(sbatch_content)

print(f"Generated {num_jobs} sbatch files in {output_dir}")
