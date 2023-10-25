#!/bin/bash
# Calea:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 1"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 2"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 3"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 4"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 5"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 6"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 7"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 8"

# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./test.out"

# Iboga:
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 0"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 1"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 2"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 3"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 4"
# sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=32G --time=024:00:00 --partition=iboga --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./main.out 5"


# Stencil sweap:
sbatch --job-name=myjob --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=250G --time=024:00:00 --partition=calea --wrap="srun -n 1 env OMP_PLACES=threads env OMP_PROC_BIND=true ./stencilTest.out 0"