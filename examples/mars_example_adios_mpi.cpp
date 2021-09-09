#include <complex>
#include <iostream>
#include "adios2.h"
#include "cxxopts.hpp"
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_io.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_spacetime_ex.hpp"

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

typedef std::complex<double> complex;

double simple_func(const double &x, const double &y) { return std::sqrt(x * x + y * y); }

int f(double cx, double cy, double cz, int maxiter, double outofbounds) {
    complex c{cx, cy};
    complex z = c;
    int i = 0;
    while (i < maxiter) {
        double n = std::norm(z);
        if (n > outofbounds) {
            break;
        }

        z = z * z + c;
        i++;
    }

    return i;
}

int main(int argc, char *argv[]) {
    int periods[2] = {0, 0};

    int rank, size;
    int coordinates[2];
    int my_grid_rank;

    adios2::Variable<double> var_data;

#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm gridComm;
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO par_io = adios.DeclareIO("BPFile_SZ");
    adios2::Engine writer;

    // Create our structured grid
    int nx_local;
    int ny_local;
    int mod_x;
    int mod_y;
    int Nx = 8000;
    int Ny = 6000;

    std::vector<double> data(Nx * Ny, 0);
    std::vector<int> ranks(Nx * Ny, rank);

    // What is the offset of this processor compared to the origin. In this case
    // origin is (0,0).
    int offset_x = 0;
    int offset_y = 0;

    // Specify the size of each dimension.
    // Can also specify how many proceses are on each dims, which in the general case
    // is not useful since MPI_Dims_create() takes care of that.
    int dims[2] = {0, 0};

    MPI_Dims_create(size, 2, dims);

    // Print out how the grid has been decomposed.
    if (rank == 0) {
        std::cout << "Decomposition is:" << dims[0] << "," << dims[1] << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // We previously created a new communicator which we will now assign to this dims
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &gridComm);
    // We evaluate the rank of this processor passing the communicator as parameter.
    MPI_Comm_rank(gridComm, &my_grid_rank);
    // We evaluate the coordinates of this processor and store these coordinates in the array coordinates[].
    MPI_Cart_coords(gridComm, my_grid_rank, 2, coordinates);

    nx_local = Nx / dims[0];
    mod_x = Nx % dims[0];

    ny_local = Ny / dims[1];
    mod_y = Ny % dims[1];

    nx_local += coordinates[0] < mod_x;
    ny_local += coordinates[1] < mod_y;

    int recv_buf_x = 0;
    int recv_buf_y = 0;
    // // The only varying part of this code is the coordinates of the processor. i.e (i_proc).
    for (int i_proc = 0; i_proc < coordinates[0]; i_proc++) {
        offset_x += Nx / dims[0] + (i_proc < mod_x);
    }

    for (int j_proc = 0; j_proc < coordinates[1]; j_proc++) {
        offset_y += Ny / dims[1] + (j_proc < mod_y);
    }
    // For printing out info in order per processor rank
    // for (int i = 0; i < size; ++i) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     double x;
    //     double y;
    //     if (i == rank) {
    //         std::cout << "Rank:" << rank << std::endl;
    //         std::cout << "x: " << coordinates[0] << std::endl;
    //         std::cout << "y: " << coordinates[1] << std::endl;
    //         std::cout << "x_local: " << nx_local << std::endl;
    //         std::cout << "y_local: " << ny_local << std::endl;
    //         std::cout << "offset_x: " << offset_x << std::endl;
    //         std::cout << "offset_y: " << offset_y << std::endl;
    //         for (int i = 0; i < nx_local; i++) {
    //             for (int j = 0; j < ny_local; j++) {
    //                 int i_global = offset_x + i;
    //                 int j_global = offset_y + j;
    //                 int idx_local = i * ny_local + j;
    //                 int idx_global = i_global * Ny + j_global;
    //                 std::cout << i << "," << j << "->" << idx_local << " ";
    //                 std::cout << "global: " << i_global << "," << j_global << "->" << idx_global << std::endl;
    //             }
    //         }
    //     }
    // }

    double x;
    double y;
    double x_min = -2.5;
    double x_max = 1.5;
    double y_min = -1.5;
    double y_max = 1.5;
    double d_x = (x_max - x_min) / (Nx - 1);
    double d_y = (y_max - y_min) / (Ny - 1);

    MPI_Barrier(MPI_COMM_WORLD);
    double tick = MPI_Wtime();

    for (int i = 0; i < nx_local; i++) {
        for (int j = 0; j < ny_local; j++) {
            int i_global = offset_x + i;
            int j_global = offset_y + j;
            int idx_local = i * ny_local + j;
            int idx_global = i_global * Ny + j_global;
            // Scaling
            x = i_global * d_x;
            y = j_global * d_y;

            // Translation
            x += x_min;
            y += y_min;
            // data[idx_local] = simple_func(x, y);
            // data[idx_local] = f(x, y, 0, 30, 4) + rank * 30;
            data[idx_local] = f(x, y, 0, 80, 4) + rank * 80;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double tock = MPI_Wtime();
    double elapsed = tock - tick;
    if (rank == 0) {
        std::cout << "Elapsed = " << elapsed << " seconds." << std::endl;
    }
    var_data = par_io.DefineVariable<double>(
        "U",
        {static_cast<unsigned long>(Nx), static_cast<unsigned long>(Ny), 1UL},
        {static_cast<unsigned long>(offset_x), static_cast<unsigned long>(offset_y), 0UL},
        {static_cast<unsigned long>(nx_local), static_cast<unsigned long>(ny_local), 1UL});

    auto rank_data =
        par_io.DefineVariable<int>("Rank",
                                   {static_cast<unsigned long>(Nx), static_cast<unsigned long>(Ny), 1UL},
                                   {static_cast<unsigned long>(offset_x), static_cast<unsigned long>(offset_y), 0UL},
                                   {static_cast<unsigned long>(nx_local), static_cast<unsigned long>(ny_local), 1UL});
    // adios2::Attribute<std::string> color;
    // adios2::Variable<uint32_t> dimension;

    writer = par_io.Open("adios_mpi.bp", adios2::Mode::Write);
    writer.BeginStep();
    writer.Put<double>(var_data, data.data());
    writer.Put<int>(rank_data, ranks.data());
    writer.EndStep();
    writer.Close();
    MPI_Finalize();
    return 0;
#else
    // Serial Implementation

#endif
    return 0;
}