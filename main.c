#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "fox.h"

#define UPPER_LIMIT (INT_MAX / 2)

void print_sol(int* sol, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int val = sol[i*n+j];
      if (val == UPPER_LIMIT) {
        val = 0;
      }
      if (j < n - 1) {
        printf("%d ", val);
      } else {
        printf("%d", val);
      }
    }
    printf("\n");
  }
}

void all_pairs_shortest_path(int* matrix, int* sol, int n, int k, Grid* grid) {
  int m = 1;
  int* a = matrix;
  int* curr = sol;
  int* b = (int*)malloc(k * k * sizeof(int));
   
  while (m < n-1) {
    memcpy(b, a, k * k * sizeof(int));
    fox(n, grid, a, b, curr);
    int* aux = a;
    a = curr;
    curr = aux;
    m = 2*m;
  }
  
  free(b);
  
  if (a != sol) {
    memcpy(sol, a, k*k*sizeof(int));
  }
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  Grid grid;
  if (!create_cartesian_grid(&grid)) {
    if (rank == 0) {
      fprintf(stderr, "Error: the number of processes (%d) must be a perfect square.\n", size);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  
  int n;
  int* matrix = NULL;
  int* sol = NULL;
  if (rank == 0) {
    scanf("%d", &n);

    if (n % grid.q != 0) {
      if (rank == 0) {
        fprintf(stderr, "Error: n (%d) must be divisible by q (%d).\n", n, grid.q);
      }
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    matrix = (int*)malloc(n*n*sizeof(int));
    sol = (int*)malloc(n*n*sizeof(int));

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        int cost;
        scanf("%d", &cost);
        matrix[i*n + j] = cost;
        if (i != j && cost == 0) {
          matrix[i*n + j] = UPPER_LIMIT;
        }
      }
    }  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  int m = n / grid.q;

  int* a = (int*)malloc(m * m * sizeof(int));
  int* c = (int*)malloc(m * m * sizeof(int));

  MPI_Datatype submatrix_type;

  MPI_Type_vector(m, m, n, MPI_INT, &submatrix_type);
  MPI_Type_commit(&submatrix_type);

  int new_root_rank;
  if (rank == 0) {
    new_root_rank = grid.rank;
  }
  MPI_Bcast(&new_root_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int* aux_buffer = NULL;

  if (rank == 0) {
    aux_buffer = (int*)malloc(n * n * sizeof(int));
    int ptr = 0;
    for (int r = 0; r < grid.q; r++) {
      for (int c = 0; c < grid.q; c++) {
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            aux_buffer[ptr++] = matrix[(r * m + i) * n + (c * m + j)];
          }
        }
      }
    }
  }
  MPI_Scatter(aux_buffer, m*m, MPI_INT, a, m*m, MPI_INT, new_root_rank, grid.grid_comm);

  all_pairs_shortest_path(a, c, n, m, &grid);
  
  MPI_Gather(c, m * m, MPI_INT, 
           aux_buffer, m*m, MPI_INT, 
           new_root_rank, grid.grid_comm);
  
  if (rank == 0) {
    int ptr = 0;
    for (int r = 0; r < grid.q; r++) {
      for (int c = 0; c < grid.q; c++) {
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            sol[(r * m + i) * n + (c * m + j)] = aux_buffer[ptr++];
          }
        }
      }
    }
    print_sol(sol, n);
    free(matrix);
    free(sol);
  }
  free(a); free(c);
  MPI_Type_free(&submatrix_type);
  MPI_Comm_free(&(grid.grid_comm));
  MPI_Comm_free(&(grid.row_comm));
  MPI_Comm_free(&(grid.col_comm));
  MPI_Finalize();
  return 0;
}
