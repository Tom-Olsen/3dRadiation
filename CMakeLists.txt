project(Radiation3d)
cmake_minimum_required(VERSION 3.16)

# Set output path:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Enable Debug Symbols:
set(CMAKE_BUILD_TYPE Debug) # has had no impact on performance yet

# Set cpp version:
set(CMAKE_CXX_STANDARD 17)

# Find packages:
find_package(OpenMP REQUIRED)

# Find all src files in src folder:
file(GLOB srcs src/*.cpp)

# Executables to generate and their dependencies:
add_executable(test test.cpp ${srcs})

# c++ flags:
# First option should work. Second is needed for hlr cluster.
target_compile_options(test PUBLIC -O3 -ffast-math)
# target_compile_options(test PUBLIC -O3 -ffast-math -I/home/compeng/olsen/spack/opt/spack/linux-scientific7-haswell/gcc-4.8.5/glm-0.9.9.8-zjqp4iqmf7er6yz2njayjgenkucohzy2/include/)

# libraries to link:
target_link_libraries(test OpenMP::OpenMP_CXX)