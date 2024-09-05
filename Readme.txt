This repository contains the 3D GRLMBRT code developed in the context of my PhD thesis, "General Relativistic Lattice Boltzmann Method for Radiation Transport", and my paper of the same title.
It is a standalone radiation transport code to simulate the behavior of radiation (photons and neutrinos) in the presence of a matter fluid in curved spacetime.
The code is not coupled to a matter fluid code yet, making all interactions between matter and radiation static.
The code is parallelized via openMP.



How To use:
The code was developed on Linux. I have not tested it on other systems and do not guarantee it will work there.

Download with submodules:
-git clone --recurse-submodules <repository-url>

Or clone first and then init submodels:
-git clone <repository-url>
-git submodule init
-git submodule update

Dont forget to create the output dir and change the global output path in ControlFlow.hh

The following libraries are necessary to compile the code:
openmp:     -sudo apt-get install libomp-dev

To use the code simply incluse the 'src/Radiation.h' in your .cpp project.
This will include all needed header files.
Examples for code usage can be found in exe/main.cpp.

To build the code add your .cpp file in the 'CMakeLists.txt' as an executable,
and add the c++ flags to it (analougus to 'main.cpp'):
-add_executable(myCode.out exe/myCode.cpp ${srcs})
-target_compile_options(myCode.out PUBLIC -O3 -ffast-math)
-target_link_libraries(myCode.out OpenMP::OpenMP_CXX)

Then run cmake inside the build folder, and afterwards the makefile:
-cd build
-cmake ..
-make <executable name>

The src folder contains all the source code, including the two submodules eigen and glm.
The exe folder contains all .cpp files that get compiled into executable .out files.
Here you can find simple programs that contain tests for the individual modules of the code.
Have a look at these to understand the general code structure.



How to perf (sudo needed if system not properly set up):
-env OMP_PLACES=threads env OMP_PROC_BIND=true perf record ./executable.out 3
-perf report -i perf.data > report.txt