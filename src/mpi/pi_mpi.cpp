//
// Created by Thoh Testarossa on 2018/5/10.
//

#include <iostream>
#include <cmath>
#include "mpi.h"

#define STEP_NUM 200000000
#define STEP_LENGTH 1/STEP_NUM
#define MAXTHREADNUM 32

int main(int argv, char *argc[])
{
    int i;
    double pi = 0;
    double pi_mypart = 0;
    int thread_sum = 4;
    if(argv > 1) thread_sum = atoi(argc[1]);
    if(thread_sum > MAXTHREADNUM) thread_sum = MAXTHREADNUM;

    int numprocs, myid;

    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    for(int i = myid; i < STEP_NUM; i += numprocs)
    {
        double x = (i + 0.5) * STEP_LENGTH;
        pi_mypart += 4.0 / (1.0 + x * x);
    }

    MPI_Reduce(&pi_mypart, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(0 == myid)
        std::cout << pi / STEP_NUM << std::endl;

    MPI_Finalize();

    return 0;
}