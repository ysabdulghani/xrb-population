import os
import math
from itertools import product

# Parameters
preset_values = {
    'g': [1.7, 2.0, 3.0],
    'T': [0.5, 0.7, 1.0],
    'a': [0.0, 0.998],
    'm': [6.0, 8.0, 10.0],
    'i': [0.0, 60.0, 80.0],
    'r': [0.2, 0.5, 0.9],
    'e': [400.0, 1000, 5000]
}

# Configuration
email = "youssefabdulghani@montana.edu"
cpus_per_task = 102
script_path = "/home/r77m975/xrb-population/observational_effects.py"
output_dir = "/home/r77m975/xrb-population/job_files"
outlog_path = "/home/r77m975/xrb-population/out_logs"
progress_log_path = "/home/r77m975/xrb-population/progress.log"
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
sbatch_template = f"""#!/bin/bash
#SBATCH --account=group-annelohfink
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=xrt-simulation-job-{{job_id}}
#SBATCH --output={outlog_path}/job-{{job_id}}-%j.out
#SBATCH --error={outlog_path}/job-{{job_id}}-%j.err
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL

echo "Running on node: $(hostname)"

# Progress tracking setup
progress_file="{progress_log_path}"
touch $progress_file

{{commands}}
"""

# Generate sbatch files
os.makedirs(output_dir, exist_ok=True)

for job_id, job_combinations in enumerate(combinations_list, start=1):
    commands = "\n".join(
        f"echo 'Running combination {index + 1} of {len(job_combinations)} for job {job_id}' >> {progress_log_path} && "
        f"apptainer exec /home/r77m975/heasoft_docker/heasoft-v6.34.sif python {script_path} "
        f"{g} {T} {a} {m} {i} {r} {e} xrt"
        for index, (g, T, a, m, i, r, e) in enumerate(job_combinations)
    )
    sbatch_content = sbatch_template.format(job_id=job_id, commands=commands)
    
    # Save to file
    filename = os.path.join(output_dir, f"job_{job_id}.sbatch")
    with open(filename, "w") as file:
        file.write(sbatch_content)

print(f"Generated {num_jobs} sbatch files in {output_dir}")
