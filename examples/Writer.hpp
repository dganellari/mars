#ifndef WRITER_HPP
#define WRITER_HPP
#include "Image.hpp"
#include "Mandelbrot.hpp"
#include "adios2.h"
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
namespace mars {
    class Writer {
    public:
        Writer();
        // Writer(mars::Image image);
        // Writer(mars::Image image);
#if ADIOS2_USE_MPI
        // adios2::ADIOS adios(MPI_COMM_WORLD);
        // adios2::IO par_io = adios.DeclareIO("BPFile_SZ");
        // adios2::Engine writer;
        // adios2::Variable<double> var_data;
        // mars::Image image = image;
        // void setup_variables() {
        //     var_data = par_io.DefineVariable<double>(
        //         "U",
        //         {static_cast<unsigned long>(image.Nx), static_cast<unsigned long>(image.Ny), 1UL},
        //         {static_cast<unsigned long>(image.offset_x), static_cast<unsigned long>(image.offset_y), 0UL},
        //         {static_cast<unsigned long>(image.nx_local), static_cast<unsigned long>(image.ny_local), 1UL});

        //     auto rank_data = par_io.DefineVariable<int>(
        //         "Rank",
        //         {static_cast<unsigned long>(image.Nx), static_cast<unsigned long>(image.Ny), 1UL},
        //         {static_cast<unsigned long>(image.offset_x), static_cast<unsigned long>(image.offset_y), 0UL},
        //         {static_cast<unsigned long>(image.nx_local), static_cast<unsigned long>(image.ny_local), 1UL});
        // };

        // void write(std::vector<double> data) {
        //     writer = par_io.Open("adios_mpi.bp", adios2::Mode::Write);
        //     writer.BeginStep();
        //     writer.Put<double>(var_data, data.data());
        //     // writer.Put<int>(rank_data, ranks.data());
        //     writer.EndStep();
        //     writer.Close();
        // };

#endif
    };
}  // namespace mars
#endif