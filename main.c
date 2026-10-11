#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fox.c"

#define UPPER_LIMIT (INT_MAX / 2)

void print_sol(int* sol, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int val = sol[i*n+j];
      if (val == UPPER_LIMIT) {
        val = 0;
      }
      printf("%d ", val);
    }
    printf("\n");
  }
}

void all_pairs_shortest_path(int* matrix, int* sol, int n, int k, Grid* grid) {
  int m = 1;
  int* prev = matrix;
  int* curr = sol;
  while (m < n-1) {
    fox(n, &grid, prev, prev, curr);
    int* aux = prev;
    prev = curr;
    curr = aux;
    m = 2*m;
  }
  if (prev != sol) {
    memcpy(sol, prev, k*k*sizeof(int));
  }
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  Grid grid;
  create_cartesian_grid(&grid);
  
  int n;
  int* matrix = NULL;
  int* sol = NULL;
  if (rank == 0) {
    scanf("%d", &n);

    if (n % grid.q != 0) {
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
    }
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  int m = n / grid.q;

  int* a = (int*)malloc(m * m * sizeof(int));
  int* c = (int*)malloc(m * m * sizeof(int));

  MPI_Datatype submatrix_type;

  MPI_Type_vector(m, m, n, MPI_INT, &submatrix_type);
  MPI_Type_commit(&submatrix_type);

  MPI_Scatter(matrix, m*m, submatrix_type, a, m*m, MPI_INT, 0, MPI_COMM_WORLD);

  all_pairs_shortest_path(a, c, n, m, &grid),
  
  MPI_Gather(c, m * m, MPI_INT, 
           sol, m * m, submatrix_type, 
           0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    print_sol(sol, n);
    free(matrix);
    free(sol);
  }
  free(a); free(c);
  MPI_Type_free(&submatrix_type);
  MPI_Comm_free(&(grid.row_comm));
  MPI_Comm_free(&(grid.col_comm));
  MPI_Finalize();
  return 0;
 }
