#ifndef MARS_IMAGE_ADIOS_WRITER_HPP
#define MARS_IMAGE_ADIOS_WRITER_HPP
#include "adios2.h"
#include "mars_image.hpp"
#include "mars_mandelbrot.hpp"
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
namespace mars {
    class Writer {
    public:
        Writer();
        Writer(adios2::IO par_io, adios2::Engine writer, Image image)
            : new_image(image), new_par_io(par_io), new_writer(writer){};
#if ADIOS2_USE_MPI
        void setup_variables() {
            var_data = new_par_io.DefineVariable<double>("U",
                                                         {static_cast<unsigned long>(new_image.Nx),
                                                          static_cast<unsigned long>(new_image.Ny),
                                                          static_cast<unsigned long>(new_image.Nz)},
                                                         {static_cast<unsigned long>(new_image.offset_x),
                                                          static_cast<unsigned long>(new_image.offset_y),
                                                          static_cast<unsigned long>(new_image.offset_z)},
                                                         {static_cast<unsigned long>(new_image.nx_local),
                                                          static_cast<unsigned long>(new_image.ny_local),
                                                          static_cast<unsigned long>(new_image.nz_local)});

            rank_data = new_par_io.DefineVariable<int>("Rank",
                                                       {static_cast<unsigned long>(new_image.Nx),
                                                        static_cast<unsigned long>(new_image.Ny),
                                                        static_cast<unsigned long>(new_image.Nz)},
                                                       {static_cast<unsigned long>(new_image.offset_x),
                                                        static_cast<unsigned long>(new_image.offset_y),
                                                        static_cast<unsigned long>(new_image.offset_z)},
                                                       {static_cast<unsigned long>(new_image.nx_local),
                                                        static_cast<unsigned long>(new_image.ny_local),
                                                        static_cast<unsigned long>(new_image.nz_local)});
        };

        void write() {
            new_writer = new_par_io.Open("adios_mpi.bp", adios2::Mode::Write);
            new_writer.BeginStep();
            new_writer.Put<double>(var_data, new_image.data.data());
            new_writer.Put<int>(rank_data, new_image.ranks.data());
            new_writer.EndStep();
            new_writer.Close();
        };

#endif
    private:
        adios2::IO new_par_io;
        adios2::Engine new_writer;
        Image new_image;
        adios2::Variable<double> var_data;
        adios2::Variable<int> rank_data;
    };
}  // namespace mars
#endif