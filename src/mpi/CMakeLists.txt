cmake_minimum_required(VERSION 3.9)
project(pc-exp)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

set(OMPI_CC "gcc-7")
set(OMPI_CXX "g++-7")

SET(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/mpi)

add_executable(pc_exp_prime_mpi prime_mpi.cpp)

add_executable(pc_exp_pi_mpi pi_mpi.cpp)

add_executable(pc_exp_montCarlo_mpi montCarlo_mpi.cpp)

add_executable(pc_exp_nBody_mpi nBody_mpi.cpp)

add_executable(pc_exp_PSRS_mpi PSRS_mpi.cpp)