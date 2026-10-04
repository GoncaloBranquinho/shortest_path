#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UPPER_LIMIT (INT_MAX / 2)

void print_sol(int* sol, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (sol[i*n+j] == UPPER_LIMIT) {
        sol[i*n+j] = 0;
      }
      printf("%d ", sol[i*n+j]);
    }
    printf("\n");
  }
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

void all_pairs_shortest_path(int* matrix, int* sol, int n) {
  int m = 1;
  int* prev = matrix;
  int* curr = sol;
  while (m < n-1) {
    special_matrix_multiply(prev, prev, curr, n);
    int* aux = prev;
    prev = curr;
    curr = aux;
    m = 2*m;
  }
  if (prev != sol) {
    memcpy(sol, prev, n*n*sizeof(int));
  }
}

int main(int argc, char *argv[]) {

  int n;
  scanf("%d", &n);

  int* matrix = (int*)malloc(n*n*sizeof(int));
  int* sol = (int*)malloc(n*n*sizeof(int));

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

  all_pairs_shortest_path(matrix, sol, n);
  print_sol(sol, n);

  free(matrix);
  free(sol);
 }
