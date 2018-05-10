cmake_minimum_required(VERSION 3.9)
project(pc-exp)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(pc_exp_prime_openmp src/prime_openmp.cpp)

add_executable(pc_exp_pi_openmp src/pi_openmp.cpp)