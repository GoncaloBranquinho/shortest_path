typedef struct {
  int p;
  MPI_Comm comm;
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  int q;
  int my_row;
  int my_col;
  int my_rank;
} GIRD_INFO_TYPE;

void Setup_grid(GRID_INFO_TYPE* grid) {
  int old_rank;
  int dimensions[2];
  int periods[2];
  int coordinate[2];
  int varying_coords;
  MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);
  grid->q = (int)sqrt((double)grid->p);
  dimensions[0] = dimensions[1] = grid->q;
  periods[0] = periods[1] = 1;
  MPI_Cart_create_(MPI_COMM,WORLD, 2, dimensions, periods, 1, &(grid->comm));
  MPI_Comm_rank(grid->comm, &(grid->my_rank));
  MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinate);
  grid->my_row = coordinate[0];
  grid->my_col = coordinate[1];

  varying_coords[0] = 0;
  MPI_Cart_sub(grid->comm, varying_coods, &(grid->row_comm));
  varying_coords[0] = 1;
  varying_coords[1] = 0;
  MPI_Cart_sub(grid->comm, varying_coords, &(grid->col_comm));
}

void Fox(int n, GRID_INFO_TYPE* grid,
       LOCAL_MATRIX_TYPE* local_A,
       LOCAL_MATRIX_TYPE* local_B,
       LOCAL_MATRIX_TYPE* local_C) {
  LOCAL_MATRIX_TYPE* temp_A;
  int step;
  int bcast_root;
  int n_bar;
  int source;
  int dest;
  int tag = 43;
  MPI_Status status;

  n_bar = n/grid->q;
  Set_to_zero(local_C);

  source = (grid->my_row + 1) % grid->q;
  dest = (grid->my_row + grid->q-1) % grid->q;

  temp_A = Local_matrix_allocate(n_bar);

  for (step = 0; step < grid->q; step++) {
    bcast_root = (grid->my_row + step) % grid->q;
    if (bcast_root == grid->my_col) {
      MPI_Bcast(local_A, 1, DERIVED_LOCAL_MATRIX, bcast_root, grid->row_comm);
      Local_matrix_multiply(local_A, local_B, local_C);
    }
    else {
      MPI_Bcast(temp_A, 1, DERIVED_LOCAL_MATRIX, bcast_root, grid->row_comm);
      Local_matrix_multiply(temp_A, local_B, local_C);
    }
    MPI_Send(local_B, 1, DERIVED_LOCAL_MATRIX, dest, tag, grid->col_comm);
    MPI_Recv(local_B, 1, DERIVED_LOCAL_MATRIX, source, tag, grid->col_comm, &status);
  }
}
