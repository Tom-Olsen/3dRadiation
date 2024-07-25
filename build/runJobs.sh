#!/bin/bash
# Main Calea:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 1"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 2"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 3"
sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 4"
sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 5"
sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 6"
sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 7"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 8"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 9"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 10"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 11"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 12"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 13"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 14"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 15"

# Main Iboga:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 1"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 2"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 3"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 4"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 5"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 6"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 7"

# Test:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./test.out"

# Benchmark:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 1"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 2"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 3"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 4"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 5"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 6"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./benchmark.out 7"

# Benchmark openMP scaling:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=1 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=8 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=16 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=24 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=32 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=40 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=48 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=56 ./benchmark.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true env OMP_NUM_THREADS=64 ./benchmark.out 0"