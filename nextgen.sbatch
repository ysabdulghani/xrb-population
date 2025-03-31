#!/bin/bash
##
## Priority partition sbatch example.
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.

## Specify the account to use; you must use your priority account.
#SBATCH --account=group-annelohfink

## Use the priority partition.
#SBATCH --partition=nextgen

## Number of CPUs to allocate. Note that in the slurm, CPUs correspond to
## logical processors, not physical cores.
#SBATCH --cpus-per-task=384

#SBATCH --mem=64G

## Maximum job run time. Format is days-HH:MM:SS.
#SBATCH --time=0-04:00:00

## Set job name.
#SBATCH --job-name=test-xrt-simulation

## Set output log and error log filenames. 
## %j is the job number and %x is the job name set above 
#SBATCH --output=out_logs/%x-%j.out
#SBATCH --error=out_logs/%x-%j.err

## Remove one # and fill out the following lines if you want to receive emails.
##SBATCH --mail-user=      # enter your email to recieve email notifications
##SBATCH --mail-type=ALL                 # specify conditions on which to recieve emails

## Run 'man sbatch' for more information on the options above.
## Below, enter commands needed to execute your workload

# Print the name of the node
echo "Running on node: $(hostname)"

apptainer exec -c --mount 'type=bind,source=/home/r77m975/xrb-population/,destination=/home/r77m975/xrb-population' --pwd /home/r77m975/xrb-population/ /home/r77m975/heasoft_docker/heasoft-v6.33.2.sif python observational_effects.py 1.7 1.0 0.0 8.0 0.0 0.2 5000 xrt