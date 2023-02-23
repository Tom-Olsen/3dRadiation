project(Radiation3d)
cmake_minimum_required(VERSION 3.16)

# Set output path:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Set cpp version:
set(CMAKE_CXX_STANDARD 17)

# Find all src files in src folder:
file(GLOB srcs src/*.cpp)

# Executables to generate and their dependencies:
add_executable(test test.cpp ${srcs})

# c++ flags:
target_compile_options(test PUBLIC -O3 -ffast-math)

# libraries to link:
find_package(OpenMP REQUIRED)
find_package(glm REQUIRED)
target_link_libraries(test jsoncpp OpenMP::OpenMP_CXX glm::glm)