#!/bin/bash
#SBATCH --job-name=compile              # Job name
#SBATCH --output=compile_output.log     # Output log file
#SBATCH --error=compile_error.log       # Error log file
#SBATCH --time=00:30:00                 # Time limit
#SBATCH --partition=calea               # Partition (replace with your partition)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=64              # Number of CPUs per task (adjust as needed)

# Load necessary environment and modules:
source /home/relastro-shared/scripts/relastro_env_setup.sh
module load gcc
cmake ..

# Run make to compile
make main.out