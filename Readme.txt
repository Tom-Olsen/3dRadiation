Do not forget to download submodules aswell:
-git submodule init
-git submodule update

The following libraries are necessary to compile the code:
openmp:     -sudo apt-get install libomp-dev

To use the code simply incluse the 'src/Radiation.h' in your .cpp project.
This will include all needed header files.
Examples for code usage can be found in exe/main.cpp.

To build the code add your .cpp file in the 'CMakeLists.txt' as an executable,
and add the c++ flags to it (analougus to 'main.cpp'):
-add_executable(myCode.out exe/main.cpp ${srcs})
-target_compile_options(myCode.out PUBLIC -O3 -ffast-math)
-target_link_libraries(myCode.out OpenMP::OpenMP_CXX)

Then run cmake inside the build folder, and afterwards the makefile:
-cd build
-cmake ..
-make -j

The src folder contains all the source code, including the two submodules eigen and glm.
The exe folder contains all .cpp files that get compiled into executable .out files.
Here you can find simple programs that contain tests for the individual modules of the code.
Have a look at these to understand the general code structure.