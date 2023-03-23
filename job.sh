#!/bin/bash
#SBATCH --job-name=myjob    # Job name
#SBATCH --nodes=1           # Number of nodes
#SBATCH --ntasks=1          # Number of MPI tasks
#SBATCH --cpus-per-task=64  # Number of cores per task
#SBATCH --mem=250G          # Memory per node
#SBATCH --time=150:00:00    # maximum run time

srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./test 1 #> output.txt

# usefull commands:
# sinfo                 shows all available nodes and their usage
# squeue -u tolsen      shows all your running jobs
# scancel <number>      cancels given job