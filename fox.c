#include "fox.h"
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>

#define UPPER_LIMIT (INT_MAX / 2)


int create_cartesian_grid(Grid* grid) {

  MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
  grid->q = sqrt(grid->p);
  
  if (grid->q * grid->q != grid->p) {
    return 0; 
  }
  
  int dims[2], periods[2];
  dims[0] = dims[1] = grid->q;
  periods[0] = periods[1] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &(grid->grid_comm));
  MPI_Comm_rank(grid->grid_comm, &(grid->rank));
  
  int coords[2];
  MPI_Cart_coords(grid->grid_comm, grid->rank, 2, coords);
  grid->row = coords[0];
  grid->col = coords[1];
  dims[0] = 0;
  dims[1] = 1;
  MPI_Cart_sub(grid->grid_comm, dims, &(grid->row_comm));
  dims[0] = 1;
  dims[1] = 0;
  MPI_Cart_sub(grid->grid_comm, dims, &(grid->col_comm));

  return 1;
}

void special_matrix_multiply(int* a, int* b, int* c, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int idx = i*n+j;
      for (int k = 0; k < n; k++) {
        int aux = a[i*n+k] + b[k*n+j];
        c[idx] = (c[idx] < aux) ? c[idx] : aux;
      }
    }
  }
}

void fox(int n, Grid* grid,
       int* a,
       int* b,
       int* c) {
  
  int m = n/grid->q; // each submatrix has size m*m 
  int* buf = (int*)malloc(m*m*sizeof(int)); 
  int tag = 0;
  MPI_Status status;

  for (int i = 0; i < m * m; i++) {
    c[i] = UPPER_LIMIT;
  }
  
  for (int step = 0; step < grid->q; step++) {
    int chosen_submatrix = (grid->row + step) % grid->q;
    int* curr_a = (chosen_submatrix == grid->col) ? a : buf;
    MPI_Bcast(curr_a, m*m, MPI_INT, chosen_submatrix, grid->row_comm);
    special_matrix_multiply(curr_a, b, c, m);
    
    int aux, dest, source;
    aux = (grid->row+1) % grid->q;
    MPI_Cart_rank(grid->col_comm, &aux, &source);
    aux = (grid->row + grid->q - 1) % grid->q;
    MPI_Cart_rank(grid->col_comm, &aux, &dest); 
    
    MPI_Sendrecv_replace(b, m*m, MPI_INT, dest, tag, source, tag, grid->col_comm, &status);
  }
  free(buf);
}
