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

int main(int argc, char *argv[]) {
    int rank, size;
#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    adios2::ADIOS adios(MPI_COMM_WORLD);

    adios2::IO par_io = adios.DeclareIO("BPFile_SZ");

    Settings settings("exampleParallelImage.bp");
    ImageWriter writer(settings, par_io);

    // This just parallizes the Puts not the data creation. Since ther is just one Put operation...
    // Still takes a lot of time beacuse the data creation is not parallel.
    // Need to ask what to do with send, recv..

    // Create example data vector
    writer.new_data(300, 400, 4);
    // Open the writer, (which now is in write mode), with the settings found.
    writer.open(settings.output);
    // Write with the following steps
    writer.write(1);
    // Close writer
    writer.close();

    MPI_Finalize();
    return 0;
#else
    // Serial Implementation

#endif
    return 0;
}