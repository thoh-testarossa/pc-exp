cmake_minimum_required(VERSION 3.9)
project(pc-exp)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_COMPILER "gcc-7")
set(CMAKE_CXX_COMPILER "g++-7")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

SET(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/openmp)

add_executable(pc_exp_prime_openmp prime_openmp.cpp)

add_executable(pc_exp_pi_openmp pi_openmp.cpp)