//
// Created by Thoh Testarossa on 2018/5/10.
//

#include <iostream>
#include <cmath>
#include "omp.h"

#define STEP_NUM 200000000
#define STEP_LENGTH 1/STEP_NUM
#define MAXTHREADNUM 32

int main(int argv, char *argc[])
{
    int i;
    double x, pi = 0, sum[MAXTHREADNUM];
    int thread_sum = 1;
    if(argv > 1) thread_sum = atoi(argc[1]);
    if(thread_sum > MAXTHREADNUM) thread_sum = MAXTHREADNUM;
    omp_set_num_threads(thread_sum);


#pragma omp parallel for reduction(+:pi) schedule(dynamic)

    for(i = 0; i < STEP_NUM; i++)
    {
        x = (i + 0.5) * STEP_LENGTH;
        pi += 4.0 / (1.0 + x * x);
    };
    std::cout << pi * STEP_LENGTH << std::endl;

    return 0;
}