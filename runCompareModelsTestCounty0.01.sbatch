#!/bin/bash
#SBATCH --job-name runCompareModelsTestCounty0.01      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition short     # Slurm partition to use
#SBATCH --ntasks 1            # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-08:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=16000    # Memory limit for each tasks (in MB)
#SBATCH -o runCompareModelsTestCounty0.01_%j.out  # File to which STDOUT will be written
#SBATCH -e runCompareModelsTestCounty0.01_%j.err  # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=paigejo@gmail.com # Email to which notifications will be sent
 
module load R
R CMD BATCH runCompareModelsTestCounty0.01.R