cmake_minimum_required(VERSION 3.9)
project(pc-exp)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(src/mpi)

add_subdirectory(src/openmp)