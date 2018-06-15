//
// Created by Thoh Testarossa on 2018/6/15.
//

#include <iostream>
#include <fstream>
#include <random>
#include "mpi.h"

using namespace std;

#define GRAPH_SIZE 10
#define NOT_CONNECTED -1

//The max weight for each edge
#define DATA_BOUND 10000

#define FILENAME "input-graph.txt"
#define RESULT_FILENAME "output-result.txt"

#define VALID_EDGE_RATE 0.4

#define MASTER_PROC 0

//number of nodes

void generateGraphFile(const string &fileName, int size, double validEdgeRate)
{
    ofstream fout(fileName);
    if(fout.is_open())
    {
        fout << size << endl;

        std::random_device r;
        std::uniform_int_distribution<int> uniform_dest(1, DATA_BOUND);
        std::default_random_engine e1(r());

        std::random_device r2;
        std::uniform_real_distribution<double> uniform_dest_2(0, 1);
        std::default_random_engine e2(r2());

        int tmp;
        double p;

        for(int i = 0; i < size; i++)
        {
            for(int j = 0; j < size; j++)
            {
                p = uniform_dest_2(e2);
                if(i == j) tmp = 0;
                else if(p <= validEdgeRate) tmp = uniform_dest(e1);
                else tmp = NOT_CONNECTED;
                fout << tmp << " ";
            }
            fout << endl;
        }
    }
}

int *readGraphFromFile(const string &fileName, int &size)
{
    int *graph = nullptr;
    ifstream fin(fileName);
    if(fin.is_open())
    {
        fin >> size;
        graph = new int [size * size];
        for(int i = 0; i < size * size; i++) fin >> graph[i];
    }

    return graph;
}

int edge_No(int size, int sor, int tar)
{
    return sor * size + tar;
}

int main(int argc, char *argv[])
{
    //Prepare the graph data
    int *graph = nullptr;
    int graphSize = 0;
    string fileName = string(FILENAME);

    ifstream fin_test(fileName);
    if (fin_test.is_open())
        fin_test.close();
    else
    {
        fin_test.close();
        generateGraphFile(fileName, GRAPH_SIZE, VALID_EDGE_RATE);
    }

    graph = readGraphFromFile(fileName, graphSize);

    int numprocs, myid;

    MPI_Init(nullptr, nullptr);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(MASTER_PROC == myid)
        cout << "Computing the graph with the size " << graphSize << endl;

    //Floyd-Warshall
    for(int k = 0; k < graphSize; k++)
    {
        //Parallel computation for sor vertex
        int my_startFrom = (myid * graphSize) / numprocs, my_endTo = ((myid + 1) * graphSize) / numprocs;
        for(int i = my_startFrom; i < my_endTo; i++)
        {
            if(graph[edge_No(graphSize, i, k)] != NOT_CONNECTED)
            {
                for(int j = 0; j < graphSize; j++)
                {
                    if(graph[edge_No(graphSize, k, j)] != NOT_CONNECTED)
                    {
                        if(graph[edge_No(graphSize, i, j)] == NOT_CONNECTED || graph[edge_No(graphSize, i, j)] > graph[edge_No(graphSize, i, k)] + graph[edge_No(graphSize, k, j)])
                            graph[edge_No(graphSize, i, j)] = graph[edge_No(graphSize, i, k)] + graph[edge_No(graphSize, k, j)];
                    }
                }
            }
        }

        //Combine the result together for next iteration
        for(int p = 0; p < numprocs; p++)
        {
            int bcast_startFrom = ((p * graphSize) / numprocs) * graphSize;
            int bcast_endTo = (((p + 1) * graphSize) / numprocs) * graphSize;
            MPI_Bcast(&graph[bcast_startFrom], bcast_endTo - bcast_startFrom, MPI_INT, p, MPI_COMM_WORLD);
        }
    }

    if(MASTER_PROC == myid)
    {
        cout << "Computation finished" << endl;
        ofstream fout(string(RESULT_FILENAME));
        if(fout.is_open())
        {
            fout << graphSize << endl;
            for (int i = 0; i < graphSize; i++)
            {
                for (int j = 0; j < graphSize; j++)
                    fout << graph[edge_No(graphSize, i, j)] << " ";
                fout << endl;
            }
        }
        fout.close();
    }

    MPI_Finalize();

    return 0;
}