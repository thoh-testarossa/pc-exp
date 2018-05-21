export OMPI_CC=gcc-7
export OMPI_CXX=g++-7

mpicxx -std=c++11 -o pc_exp_pi_mpi pi_mpi.cpp
mpicxx -std=c++11 -o pc_exp_prime_mpi prime_mpi.cpp

unset OMPI_CC
unset OMPI_CXX
