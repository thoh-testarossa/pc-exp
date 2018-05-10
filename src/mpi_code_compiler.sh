export OMPI_CC=gcc-7
export OMPI_CXX=g++-7

mpicxx -o pc_exp_pi_mpi pi_mpi.cpp
mpicxx -o pc_exp_prime_mpi prime_mpi.cpp

unset OMPI_CC
unset OMPI_CXX