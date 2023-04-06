project(Radiation3d)
cmake_minimum_required(VERSION 3.16)

# Set output path:
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

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
add_executable(gridTest tests/GridTest.cpp ${srcs})
add_executable(tensorTypesTest tests/TensorTypesTest.cpp ${srcs})
add_executable(interpolationTest tests/InterpolationTest.cpp ${srcs})

# c++ flags:
# target_compile_options(test PUBLIC)
target_compile_options(test PUBLIC -O3 -ffast-math)
target_compile_options(gridTest PUBLIC -O3 -ffast-math)

# libraries to link:
target_link_libraries(test OpenMP::OpenMP_CXX)
target_link_libraries(gridTest OpenMP::OpenMP_CXX)
target_link_libraries(tensorTypesTest OpenMP::OpenMP_CXX)
target_link_libraries(interpolationTest OpenMP::OpenMP_CXX)