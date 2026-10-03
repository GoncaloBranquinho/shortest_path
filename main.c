#include <limits.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Missing input file");
    return 1;
  }

  FILE *f = fopen(argv[1], "r");

  int n;
  fscanf(f, "%d", &n);

  int matrix[n][n];

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int cost;
      fscanf(f, "%d", &cost);
      matrix[i][j] = cost;
      if (i != j && cost == 0) {
        matrix[i][j] = INT_MAX;
      }
    }
  }
  fclose(f);

 }
