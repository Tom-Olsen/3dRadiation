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
	rm -r output/*