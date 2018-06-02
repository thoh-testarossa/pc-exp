//
// Created by Thoh Testarossa on 2018/5/25.
//

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "mpi.h"

#define MASTER 0
#define DIM 2                           /* Number of dimension          */
#define X 0                             /* coordinate x                 */
#define Y 1                             /* coordinate y                 */
#define DEBUG_GET_ARGUMENT 1
#define GENERATE_INPUT_FILE 1
#define GENERATE_OUTPUT_FILE 1

typedef double vector[DIM];             /* 2-dimensional axis system    */
const double G = 6.673e-11;

int my_rank;                            /* Rank of this process         */
int size;                               /* Total number of process      */

MPI_Datatype vectorMPI;

vector *velocities = NULL;              /* Array of total v table       */

void help()
{
    std::cout << "Usage:" << std::endl
         << "   mpirun -np <num_of_procs> <name_of_program> "
         << "<num_of_particles> <num_of_timesteps> "
         << "<time_internal> <output_frequency> "
         << std::endl;
}

void getInputArguments(int argc, char *argv[], int *num_particles, int *num_steps, double *delta_t, int *output_freq)
{
    if (argc != 5)
    {
        if (my_rank == MASTER)
            help();
        MPI_Finalize();
        exit(0);
    }

        *num_particles = strtol(argv[1], NULL, 10);
        *num_steps = strtol(argv[2], NULL, 10);
        *delta_t = strtod(argv[3], NULL);
        *output_freq = strtol(argv[4], NULL, 10);

    MPI_Bcast(num_particles, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(num_steps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(delta_t, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(output_freq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


    if (*num_particles <= 0 || *num_steps < 0 || *delta_t <= 0)
    {
        if (my_rank == MASTER)
            help();
        MPI_Finalize();
        exit(0);
    }

#  ifdef DEBUG_GET_ARGUMENT
    if (my_rank == 0)
    {
        printf("num_particles = %d\n", *num_particles);
        printf("num_steps = %d\n", *num_steps);
        printf("delta_t = %e\n", *delta_t);
        printf("output_freq = %d\n", *output_freq);
    }
#  endif
}

void generateInputFile(double *masses, vector *positions, vector *velocities, int num_particles)
{
    std::ofstream fout("generated_input.txt");

    if(!fout.is_open())
    {
        std::cout << "Cannot create input file generated_input.txt" << std::endl;
        exit(1);
    }

    for(int line = 0; line < num_particles; line++)
    {
        fout << masses[line] << " "
             << positions[line][X] << " "
             << positions[line][Y] << " "
             << velocities[line][X] << " "
             << velocities[line][Y] << " "
             << std::endl;
    }

    fout.close();
}

void generateInitConditions(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk)
{
    double mass = 10000;
    double gap = 0.01;

    int pos_square_border = (int)sqrt((double)num_particles);

    if (my_rank == MASTER)
    {
        for (int part = 0; part < num_particles; part++)
        {
            masses[part] = mass;
            positions[part][X] = (double)(part % pos_square_border) * gap;
            positions[part][Y] = (double)(part / pos_square_border) * gap;
            velocities[part][X] = 0.0;
            velocities[part][Y] = 0.0;
        }

#ifdef GENERATE_INPUT_FILE
        generateInputFile(masses, positions, velocities, num_particles);
#endif

    }

    MPI_Bcast(masses, num_particles, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(positions, num_particles, vectorMPI, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(velocities, chunk, vectorMPI, my_velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
}

void computeForce(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk)
{

    double m_g;
    vector f_part_k;
    double len, len_3, fact;

    int part = my_rank * chunk + my_particles;
    my_forces[my_particles][X] = my_forces[my_particles][Y] = 0.0;

    for (int k = 0; k < num_particles; k++)
    {
        if (k != part)
        {
            f_part_k[X] = positions[part][X] - positions[k][X];
            f_part_k[Y] = positions[part][Y] - positions[k][Y];


            len = sqrt(pow(f_part_k[X], 2) + pow(f_part_k[Y], 2));
            len_3 = pow(len, 3);


            m_g = G * masses[part] * masses[k];
            fact = m_g / len_3;

            f_part_k[X] *= fact;
            f_part_k[Y] *= fact;

            my_forces[my_particles][X] += f_part_k[X];
            my_forces[my_particles][Y] += f_part_k[Y];
        }
    }
}

void updateParticles(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double delta_t)
{

    int part;
    double fact;

    part = my_rank*chunk + my_particles;
    fact = delta_t/masses[part];

    my_positions[my_particles][X] += delta_t * my_velocities[my_particles][X];
    my_positions[my_particles][Y] += delta_t * my_velocities[my_particles][Y];
    my_velocities[my_particles][X] += fact * my_forces[my_particles][X];
    my_velocities[my_particles][Y] += fact * my_forces[my_particles][Y];

}

void generateOutputFile(double masses[], vector positions[], vector velocities[], int num_particles)
{

    if (my_rank == MASTER)
    {
        std::ofstream fout("generated_output.txt");

        if (!fout.is_open())
        {
            std::cout << "Cannot create input file generated_output.txt" << std::endl;
            exit(1);
        }

        for (int line = 0; line < num_particles; line++)
        {
            fout << masses[line] << " "
                 << positions[line][X] << " "
                 << positions[line][Y] << " "
                 << velocities[line][X] << " "
                 << velocities[line][Y] << " "
                 << std::endl;
        }

        fout.close();
    }
}

//main program
int main(int argc, char* argv[])
{
    //Definition part
    int num_particles;                  /* Number of the particles              */
    int chunk;                          /* Number of the particles per process  */
    int num_steps;                      /* Number of timesteps                  */
    int my_particles;                   /* Particles id for each process        */
    int output_freq;                    /* Test output frequency                */
    double delta_t;                     /* Internal of each timestep            */
    double *masses;                     /* Array of particles masses            */
    vector *my_positions;               /* Array of pos in this process         */
    vector *positions;                  /* Array of pos of all particles        */
    vector *my_velocities;              /* Array of v in this process           */
    vector *my_forces;                  /* Array of f in this process           */

    double start_time;                  /* MPI start time                       */
    double end_time;                    /* MPI end time                         */

    //Initialization part
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    getInputArguments(argc, argv, &num_particles, &num_steps, &delta_t, &output_freq);
    chunk = num_particles / size;

    masses = new double[num_particles];
    positions = new vector[num_particles];
    my_forces = new vector[chunk];
    my_positions = positions + my_rank * chunk;
    my_velocities = new vector[chunk];

    if (my_rank == MASTER)
        velocities = new vector[num_particles];

    MPI_Type_contiguous(DIM, MPI_DOUBLE, &vectorMPI);
    MPI_Type_commit(&vectorMPI);

    //Initial condition generation
    generateInitConditions(masses, positions, my_velocities, num_particles, chunk);

    start_time = MPI_Wtime();

    //Simulation part
    for (int steps = 1; steps <= num_steps; steps++)
    {
        for (my_particles = 0; my_particles < chunk; my_particles++)
            computeForce(my_particles, masses, my_forces, positions, num_particles, chunk);
        for (my_particles = 0; my_particles < chunk; my_particles++)
            updateParticles(my_particles, masses, my_forces, my_positions, my_velocities, num_particles, chunk,
                             delta_t);

        MPI_Allgather(MPI_IN_PLACE, chunk, vectorMPI, positions, chunk, vectorMPI, MPI_COMM_WORLD);
    }

    MPI_Gather(my_velocities, chunk, vectorMPI, velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);

#ifdef GENERATE_OUTPUT_FILE

    if(my_rank == MASTER)
        generateOutputFile(masses, positions, my_velocities, num_particles);
#endif

    end_time = MPI_Wtime();
    if (my_rank == MASTER)
    {
        std::cout << "Total time: " << end_time - start_time << " s" << std::endl
                  << "Number of particles: " << num_particles << std::endl
                  << "Number of processors: " << size << std::endl
                  ;
    }

    //Finishing the execution
    MPI_Type_free(&vectorMPI);
    free(masses);
    free(positions);
    free(my_forces);
    free(my_velocities);
    if (my_rank == MASTER)
        free(velocities);

    MPI_Finalize();
    return 0;
}
