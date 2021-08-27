#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {

  int periods[2] = {0, 0};

  int rank, size;
  int coordinates[2];
  int my_grid_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm gridComm;

  // Create our structured grid
  int nx_local;
  int ny_local;
  int mod_x;
  int mod_y;
  int Nx = 7;
  int Ny = 5;
  int gridSize = Nx * Ny;

  // Specify the size of each dimension, how many processes are
  // on each dimension.
  int dims[2];

  MPI_Dims_create(size, 2, dims);

  if (rank == 0) {
    std::cout << "Decomposition is:" << dims[0] << "," << dims[1] << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &gridComm);
  MPI_Comm_rank(gridComm, &my_grid_rank);
  MPI_Cart_coords(gridComm, my_grid_rank, 2, coordinates);

  nx_local = Nx / dims[0];
  mod_x = Nx % dims[0];

  ny_local = Ny / dims[1];
  mod_y = Ny % dims[1];

  nx_local += coordinates[0] < mod_x;
  ny_local += coordinates[1] < mod_y;

  int offset_x = 0;
  int offset_y = 0;
  for (int i_proc = 0; i_proc < coordinates[0]; i_proc++) {
    offset_x += Nx / dims[0] + (i_proc < mod_x);
  }

  for (int j_proc = 0; j_proc < coordinates[1]; j_proc++) {
    offset_y += Ny / dims[1] + (j_proc < mod_y);
  }
  // TODO: Revisit with MPI exclusive scan (inclusive scan).
  // Cumulative sums.

  // Use assert for verification of scan and below offset calculation.

  // For printing out info in order per processor rank
  for (int i = 0; i < size; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (i == rank) {
      std::cout << "Rank:" << rank << std::endl;
      std::cout << "x: " << coordinates[0] << std::endl;
      std::cout << "y: " << coordinates[1] << std::endl;
      std::cout << "x_local: " << nx_local << std::endl;
      std::cout << "y_local: " << ny_local << std::endl;
      std::cout << "offset_x: " << offset_x << std::endl;
      std::cout << "offset_y: " << offset_y << std::endl;
      for (int i = 0; i < nx_local; i++) {
        for (int j = 0; j < ny_local; j++) {
          int i_global = offset_x + i;
          int j_global = offset_y + j;
          int idx_local = i * ny_local + j;
          int idx_global = i_global * Ny + j_global;
          std::cout << i << "," << j << "->" << idx_local << " ";
          std::cout << "global: " << i_global << "," << j_global << "->"
                    << idx_global << std::endl;
        }
      }
    }
  }

  MPI_Finalize();
  return 0;
}