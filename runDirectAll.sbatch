#!/bin/bash
#SBATCH --job-name runDirect      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition medium     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 2-00:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=20000     # Memory limit for each tasks (in MB)
#SBATCH -o runDirect&_%A_%a.out    # File to which STDOUT will be written
#SBATCH -e runDirect_%A_%a.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=me@uw.edu # Email to which notifications will be sent
module load R
Rscript --verbose runDirectAll.R ${SLURM_ARRAY_TASK_ID} > runDirectAll_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout