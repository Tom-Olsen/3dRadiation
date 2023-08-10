OMP_PLACES=threads env OMP_PROC_BIND=true perf record  ./main.out
perf report -i perf.data > report.txt