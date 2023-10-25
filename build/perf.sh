#OMP_PLACES=threads env OMP_PROC_BIND=true perf record  ./main.out
#perf report -i perf.data > report.txt
OMP_PLACES=threads env OMP_PROC_BIND=true perf stat -e cache-references,cache-misses ./main.out