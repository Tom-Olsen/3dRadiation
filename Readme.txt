Do not forget to download submodules aswell:
-git submodule init
-git submodule update

The following libraries are necessary to compile the code:

openmp:     -sudo apt-get install libomp-dev
glm:        -sudo apt-get install libglm-dev



To use the code simply incluse the 'src/Includes.hh' in your .cpp project.
This will include all needed header files.



To build the code add your .cpp file in the 'CMakeLists.txt' as an executable,
and add the c++ flags to it (analougus to 'test.cpp'):
    -add_executable(mycode mycode.cpp ${srcs})
    -target_compile_options(mycode PUBLIC -O3 -ffast-math)

Then run cmake inside the build folder:
    -cd build
    -cmake ..

Finally to make compilation easier you can also add your .cpp in the 'Makefile' in the root folder
(analougus to the 'test'):
    mycode: FORCE
	    make mycode -sCbuild -j

Finally just write 'make mycode' in the root directory to build your executable.