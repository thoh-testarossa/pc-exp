//
// Created by Thoh Testarossa on 2018/5/21.
//

#include <iostream>
#include <random>
#include "mpi.h"

#define TIMELIMIT 100000
#define TOTALSCALE 200000

#define INITDISTANCE 50
#define VMAX 20
#define P 0.514

void montCarlo_oneStep_v(int pid, int numproc, int *d_sor, int upper_d, int *v_sor, int *v_des, int vmax, int scale, double p)
{
    double p_generated;

    std::random_device r;
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::default_random_engine e1(r());

    for(int i = 0; i < scale; i++)
    {
        //Step 1: To determine a new speed value
        if(i != scale - 1)
        {
            if(v_sor[i] >= d_sor[i + 1] - d_sor[i]) v_des[i] = d_sor[i + 1] - d_sor[i] - 1;
            else v_des[i] = v_sor[i] + 1 <= vmax ? v_sor[i] + 1 : vmax;
        }
        else
        {
            if((pid < numproc - 1) && (v_sor[i] >= upper_d - d_sor[i])) v_des[i] = upper_d - d_sor[i] - 1;
            else v_des[i] = v_sor[i] + 1 <= vmax ? v_sor[i] + 1 : vmax;
        }

        //Step 2: Random slow down
        p_generated = uniform_dist(e1);
        if(p_generated <= p && v_des[i] > 0) v_des[i]--;
    }

}

void montCarlo_oneStep_d(int pid, int *d_sor, int *d_des, int *v_des, int scale)
{
    for(int i = 0; i < scale; i++) d_des[i] = d_sor[i] + v_des[i];
}

int main(int argv, char *argc[])
{
    int totalScale = TOTALSCALE;

    int *total_v, *total_d;

    //Initialization part
    MPI_Init(nullptr, nullptr);

    int numprocs, myid;

    int upper_d_from_next_process = 0, upper_d_to_previous_process = 0;

    MPI_Status s;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int *my_v_1, *my_v_2, *my_d_1, *my_d_2;
    int my_scale = myid < totalScale % numprocs ? totalScale / numprocs + 1 : totalScale / numprocs;

    my_v_1 = new int [my_scale];
    my_v_2 = new int [my_scale];
    my_d_1 = new int [my_scale];
    my_d_2 = new int [my_scale];

    //Init the local v&d table
    //Main process should do some preparation
    if(0 == myid)
    {
        //Fill the global v&d table
        total_v = new int [totalScale], total_d = new int [totalScale];

        std::random_device r;
        std::uniform_int_distribution<int> uniform_dest(0, INITDISTANCE);
        std::default_random_engine e1(r());
        total_d[0] = 0;
        for(int i = 1; i < totalScale; i++)
            total_d[i] = total_d[i - 1] + uniform_dest(e1);

        std::random_device r2;
        std::uniform_int_distribution<int> uniform_dest2(0, VMAX);
        std::default_random_engine e2(r2());
        for(int i = 0; i < totalScale; i++)
            total_v[i] = uniform_dest2(e2);

        //Distribute the table to other processers
        for(int i = numprocs - 1; i > 0; i--)
        {
            int startFrom = i < totalScale % numprocs ? (totalScale / numprocs) * i + i : (totalScale / numprocs) * i + (totalScale % numprocs);
            int localScale = i < totalScale % numprocs ? totalScale / numprocs + 1 : totalScale / numprocs;
            MPI_Send((void *)&total_v[startFrom], localScale, MPI_INT, i, INT_MAX - 1, MPI_COMM_WORLD);
            MPI_Send((void *)&total_d[startFrom], localScale, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for(int i = 0; i < my_scale; i++) my_d_2[i] = total_d[i], my_v_2[i] = total_v[i];
    }
    else
    {
        MPI_Recv((void *)my_v_2, my_scale, MPI_INT, 0, -1, MPI_COMM_WORLD, &s);
        MPI_Recv((void *)my_d_2, my_scale, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
    }

    //Calculation part
    int *d_sor = my_d_1, *d_des = my_d_2, *v_sor = my_v_1, *v_des = my_v_2;

    for(int time = 0; time < TIMELIMIT; time++)
    {
        if (myid > 0)
            //Prepare the upper d which will be sent to previous process
            upper_d_to_previous_process = time % 2 == 0 ? my_d_2[0] : my_d_1[0];
        if (myid > 0)
            //MPI_SEND the distance of the min_id obj to the "previous" process
            MPI_Send((void *)&upper_d_to_previous_process, 1, MPI_INT, myid - 1, time * numprocs + myid, MPI_COMM_WORLD);
        if (myid < numprocs - 1)
            //MPI_RECEIVE the distance of the min_id obj from the "next" process
            MPI_Recv((void *)&upper_d_from_next_process, 1, MPI_INT, myid + 1, time * numprocs + myid + 1, MPI_COMM_WORLD, &s);

        if(time % 2 == 0)
            d_sor = my_d_2, d_des = my_d_1, v_sor = my_v_2, v_des = my_v_1;
        else
            d_sor = my_d_1, d_des = my_d_2, v_sor = my_v_1, v_des = my_v_2;

        montCarlo_oneStep_v(myid, numprocs, d_sor, upper_d_from_next_process, v_sor, v_des, VMAX, my_scale, P);
        montCarlo_oneStep_d(myid, d_sor, d_des, v_des, my_scale);
    }

    //Collection part
    if(0 != myid)
    {
        MPI_Send((void *)v_des, my_scale, MPI_INT, 0, TIMELIMIT * numprocs + 2 * myid - 1, MPI_COMM_WORLD);
        MPI_Send((void *)d_des, my_scale, MPI_INT, 0, TIMELIMIT * numprocs + 2 * myid, MPI_COMM_WORLD);
    }
    else
    {
        for(int i = numprocs - 1; i > 0; i--)
        {
            int startFrom = i < totalScale % numprocs ? (totalScale / numprocs) * i + i : (totalScale / numprocs) * i + (totalScale % numprocs);
            int localScale = i < totalScale % numprocs ? totalScale / numprocs + 1 : totalScale / numprocs;
            MPI_Recv((void *)&total_v[startFrom], localScale, MPI_INT, i, TIMELIMIT * numprocs + 2 * i - 1, MPI_COMM_WORLD, &s);
            MPI_Recv((void *)&total_d[startFrom], localScale, MPI_INT, i, TIMELIMIT * numprocs + 2 * i, MPI_COMM_WORLD, &s);
        }
        for(int i = 0; i < my_scale; i++)
        {
            total_v[i] = v_des[i];
            total_d[i] = d_des[i];
        }

        //Output part
        for(int i = 0; i < totalScale; i += 100)
            std::cout << total_d[i] << " ";
        std::cout << std::endl;
    }

    MPI_Finalize();

    return 0;
}