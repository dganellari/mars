#include <complex>
#include <iostream>
#include "adios2.h"
#include "cxxopts.hpp"
#include "mars_image.hpp"
#include "mars_image_adios_writer.hpp"

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
// Single responsibiliy principle. Function has one single purpose(makes the function also smaller).
// Open/Close
// Desingn by contract.

int main(int argc, char *argv[]) {
    int periods[3] = {0, 0, 0};

    int rank, size;
    int coordinates[3];
    int my_grid_rank;
    int Nx = 1200;
    int Ny = 1000;
    int Nz = 20;

    if (argc == 4) {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        Nz = atoi(argv[3]);
    }
    std::cout << argc;

    mars::Image image(Nx, Ny, Nz);

#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm gridComm;

    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO par_io = adios.DeclareIO("BPFile_SZ");
    adios2::Engine writer;

    // std::vector<double> data(Nx * Ny, 0);
    // std::vector<int> ranks(Nx * Ny, rank);

    // Specify the size of each dimension.
    // Can also specify how many proceses are on each dims, which in the general case
    // is not useful since MPI_Dims_create() takes care of that.
    int dims[3] = {0, 0, 0};

    MPI_Dims_create(size, 3, dims);

    // Print out how the grid has been decomposed.
    if (rank == 0) {
        std::cout << "Decomposition is:" << dims[0] << "," << dims[1] << "," << dims[2] << std::endl;
    }

    // We previously created a new communicator which we will now assign to this dims
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &gridComm);
    // We evaluate the rank of this processor passing the communicator as parameter.
    MPI_Comm_rank(gridComm, &my_grid_rank);
    // We evaluate the coordinates of this processor and store these coordinates in the array coordinates[].
    MPI_Cart_coords(gridComm, my_grid_rank, 3, coordinates);

    MPI_Barrier(MPI_COMM_WORLD);
    double tick = MPI_Wtime();

    image.setup(dims, coordinates);
    image.calc_function(mars::Mandelbrot());

    mars::Writer testWriter(par_io, writer, image);
    // Calculate how long it takes to create data.
    MPI_Barrier(MPI_COMM_WORLD);
    double tock = MPI_Wtime();
    double elapsed = tock - tick;
    if (rank == 0) {
        std::cout << "Creating data, Elapsed = " << elapsed << " seconds." << std::endl;
    }
    testWriter.setup_variables();
    testWriter.write();

    // image.calc_function([](double x, double y, double z) -> double { return x + y + z; });
    // image.calc_function(simple_func);

    MPI_Finalize();
    return 0;
#else
    // Serial Implementation

#endif
    return 0;
}