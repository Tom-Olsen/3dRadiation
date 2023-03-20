# shortcuts to build cmake stuff from parent directory, because I'm lazy.

all: paper test

FORCE: ;
paper: FORCE
	make paper -sCbuild -j
test: FORCE
	make test -sCbuild -j

.Phony: clean
clean:
	make clean -sCbuild

.Phony: outputClean
outputClean:
	rm -r -f output/*
	rm -r -f /mnt/ceph/tolsen/output/*

.Phony: runTest
runTest:
	env OMP_PLACES=threads env OMP_PROC_BIND=true ./test 3
# 	nohup env OMP_PLACES=threads env OMP_PROC_BIND=true ./test 3 > output.txt &
#	OMP_PLACES=threads OMP_PROC_BIND=true ./test 3

.Phony: runTestPerf
runTestPerf:
	sudo env OMP_PLACES=threads env OMP_PROC_BIND=true perf record  ./test 3
	sudo chmod -R 777 .
# For some reason 'sudo perf record <programm>' changes the ownership of files.

.Phony: perfReport
perfReport:
	sudo perf report -i perf.data > report.txt