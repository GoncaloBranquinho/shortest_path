#ifndef FOX_H
#define FOX_H

#include <mpi.h>

typedef struct { 
    int p;      // number of processes
    int q;      // square root of p, q x q grid of processes
    int rank;   // rank of this process in grid_comm
    int row;    // row where this process is located
    int col;    // col where this process is located
    MPI_Comm row_comm; 
    MPI_Comm col_comm;
    MPI_Comm grid_comm;
} Grid;

int create_cartesian_grid(Grid* grid);
void special_matrix_multiply(int* a, int* b, int* c, int n);
void fox(int n, Grid* grid, int* a, int* b, int* c);

#endif
