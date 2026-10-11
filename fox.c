#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#define UPPER_LIMIT (INT_MAX / 2)


struct { 
  int p;
  int q;
  int row;
  int col;
  MPI_Comm row_comm;
  MPI_Comm col_comm;
} Grid;


int create_cartesian_grid(Grid* grid) {

  MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
  grid->q = sqrt(grid->p);
  int dims[2], periods[2];
  dims[0] = dims[1] = q;
  periods[0] = periods[1] = 1;
  MPI_Comm grid_comm;
  // tentar trocar para 1
  MPI_Cart_create(MPI_COMM_WORLD, q, dims, periods, 0, &grid_comm);
  MPI_Comm_rank(grid->grid_comm, &(grid->rank));
  
  int coords[2];
  MPI_Cart_coords(grid->grid_comm, grid->rank, grid->q, coords);
  grid->row = coords[0];
  grid->col = coords[1];
  dims[0] = 0;
  dims[1] = 1;
  MPI_cart_sub(grid->grid_comm, dims, grid->row_comm);
  dims[0] = 1;
  dims[1] = 0;
  MPI_cart_sub(grid->grid_comm, dims, grid->col_comm);
  MPI_Comm_free(grid_comm);
  
}

void special_matrix_multiply(int* a, int* b, int* c, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int idx = i*n+j;
      c[idx] = UPPER_LIMIT;
      for (int k = 0; k < n; k++) {
        int aux = a[i*n+k] + b[k*n+j];
        c[idx] = (c[idx] < aux) ? c[idx] : aux;
      }
    }
  }
}

void fox(int n, INFO* grid,
       int* a,
       int* b,
       int* c) {
  
  int m = n/grid->q; // each submatrix has size m*m 
  int* buf = (int*)malloc(m*m*sizeof(int)); 
  int tag = 0;
  MPI_Status status;
  
  for (int step = 0; step < grid->q; step++) {
    int chosen_submatrix = (grid->row + step) % grid->q;
    if (chosen_submatrix == grid->col) {
      MPI_Bcast(a, m*m, MPI_INT, chosen_submatrix, grid->row_comm);
      Local_matrix_multiply(a, b, c);
    }
    else {
      MPI_Bcast(buf, m*m, MPI_INT, chosen_submatrix, grid->row_comm);
      Local_matrix_multiply(buf, b, c);
    }

    int aux, from, to;
    aux = (grid->row+1) % grid->q;
    MPI_Cart_rank(grid->col_comm, &aux, &from);
    aux = (grid->row + grid->q - 1) % grid->q;
    MPI_Cart_rank(grid->col_comm, &aux, &to); 
    
    MPI_Send(b, m*m, MPI_INT, to, tag, grid->col_comm);
    MPI_Recv(b, m*m, MPI_INT, from, tag, grid->col_comm, &status);
  }
  free(buf);
}
