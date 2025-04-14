#!/bin/bash
# get cpu model:
# lscpu | grep "Model name"

#OMP_PLACES=threads env OMP_PROC_BIND=true perf record ./main.out
#perf report -i perf.data > report.txt

# measure cache efficiency:
# OMP_PLACES=threads env OMP_PROC_BIND=true perf stat -e cache-references,cache-misses ./benchmark.out



# Measure memory bandwidth (of entire program):
# Adjust this based on CPU model — check via `perf list`
EVENTS="uncore_imc/cas_count_read/,uncore_imc/cas_count_write/"
OMP_PLACES=threads OMP_PROC_BIND=true perf stat -e $EVENTS ./main.out

# Record memory bandwidth (breakdown in functions): DOES NOT WORK!
# EVENTS="uncore_imc/cas_count_read/,uncore_imc/cas_count_write/"
# OMP_PLACES=threads OMP_PROC_BIND=true perf record -a -e $EVENTS -o perf.data ./main.out
# perf report -i perf.data