project(Radiation3d)
cmake_minimum_required(VERSION 3.16)

# Set output path:
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Enable Debug Symbols:
set(CMAKE_BUILD_TYPE Debug) # has had no impact on performance yet

# Set cpp version:
set(CMAKE_CXX_STANDARD 20)

# Find packages:
find_package(OpenMP REQUIRED)

# Find all src files in src folder:
file(GLOB srcs src/*.cpp)

# Executables to generate and their dependencies:
add_executable(main.out tests/main.cpp ${srcs})
add_executable(test.out tests/test.cpp ${srcs})
add_executable(gridTest.out tests/gridTest.cpp ${srcs})
add_executable(tensorTypesTest.out tests/tensorTypesTest.cpp ${srcs})
add_executable(interpolationTest.out tests/interpolationTest.cpp ${srcs})
add_executable(metricTest.out tests/metricTest.cpp ${srcs})
add_executable(geoEqSolverTest.out tests/geodesicEquationSolverTest.cpp ${srcs})
add_executable(stencilTest.out tests/stencilTest.cpp ${srcs})
add_executable(myMathTest.out tests/myMathTest.cpp ${srcs})
add_executable(bigFunctionsTest.out tests/bigFunctionsTest.cpp ${srcs})
add_executable(meshTest.out tests/meshTest.cpp ${srcs})
add_executable(quaternionTest.out tests/quaternionTest.cpp ${srcs})

# c++ flags:
# target_compile_options(main PUBLIC)
target_compile_options(main.out PUBLIC -O3 -ffast-math)
target_compile_options(test.out PUBLIC -O3 -ffast-math)
target_compile_options(gridTest.out PUBLIC -O3 -ffast-math)
target_compile_options(tensorTypesTest.out PUBLIC -O3 -ffast-math)
target_compile_options(interpolationTest.out PUBLIC -O3 -ffast-math)
target_compile_options(metricTest.out PUBLIC -O3 -ffast-math)
target_compile_options(geoEqSolverTest.out PUBLIC -O3 -ffast-math)
target_compile_options(stencilTest.out PUBLIC -O3 -ffast-math)
target_compile_options(myMathTest.out PUBLIC -O3 -ffast-math)
target_compile_options(bigFunctionsTest.out PUBLIC -O3 -ffast-math)
target_compile_options(meshTest.out PUBLIC -O3 -ffast-math)
target_compile_options(quaternionTest.out PUBLIC -O3 -ffast-math)

# libraries to link:
target_link_libraries(main.out OpenMP::OpenMP_CXX)
target_link_libraries(test.out OpenMP::OpenMP_CXX)
target_link_libraries(gridTest.out OpenMP::OpenMP_CXX)
target_link_libraries(tensorTypesTest.out OpenMP::OpenMP_CXX)
target_link_libraries(interpolationTest.out OpenMP::OpenMP_CXX)
target_link_libraries(metricTest.out OpenMP::OpenMP_CXX)
target_link_libraries(geoEqSolverTest.out OpenMP::OpenMP_CXX)
target_link_libraries(stencilTest.out OpenMP::OpenMP_CXX)
target_link_libraries(myMathTest.out OpenMP::OpenMP_CXX)
target_link_libraries(bigFunctionsTest.out OpenMP::OpenMP_CXX)
target_link_libraries(meshTest.out OpenMP::OpenMP_CXX)
target_link_libraries(quaternionTest.out OpenMP::OpenMP_CXX)