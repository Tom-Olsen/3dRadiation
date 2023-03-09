#!/bin/bash
#SBATCH --job-name=myjob           # Job name
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --time=00:30:00            # maximum run time

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the module environment suitable for the job
module load foss/2019a

# Use '&' to start the first job in the background
srun -n 1 ./job1 &
srun -n 1 ./job2 

# Use 'wait' as a barrier to collect both executables when they are done. If not the batch job will finish when the job2.batch program finishes and kill job1.batch if it is still running.
wait