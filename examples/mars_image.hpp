#ifndef MARS_IMAGE_HPP
#define MARS_IMAGE_HPP
#include <vector>
#include "adios2.h"
#include "mars_mandelbrot.hpp"
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
namespace mars {
    class Image {
    public:
        Image();
        Image(const int Nx, const int Ny, const int Nz) {
            this->Nx = Nx;
            this->Ny = Ny;
            this->Nz = Nz;
        };
        void setup(int dims[], int coordinates[]) {
            nx_local = Nx / dims[0];
            mod_x = Nx % dims[0];

            ny_local = Ny / dims[1];
            mod_y = Ny % dims[1];

            nz_local = Nz / dims[2];
            mod_z = Nz % dims[2];

            nx_local += coordinates[0] < mod_x;
            ny_local += coordinates[1] < mod_y;
            nz_local += coordinates[2] < mod_z;

            // // The only varying part of this code is the coordinates of the processor. i.e (i_proc).
            for (int i_proc = 0; i_proc < coordinates[0]; i_proc++) {
                offset_x += Nx / dims[0] + (i_proc < mod_x);
            }

            for (int j_proc = 0; j_proc < coordinates[1]; j_proc++) {
                offset_y += Ny / dims[1] + (j_proc < mod_y);
            }

            for (int k_proc = 0; k_proc < coordinates[2]; k_proc++) {
                offset_z += Nz / dims[2] + (k_proc < mod_z);
            }
        };

        template <class Function>
        void calc_function(Function f) {
            double d_x = (x_max - x_min) / (Nx - 1);
            double d_y = (y_max - y_min) / (Ny - 1);
            double d_z = (z_max - z_min) / (Nz - 1);
            data.resize(Nx * Ny * Nz);
            ranks.resize(Nx * Ny * Nz);

            for (int i = 0; i < nx_local; i++) {
                for (int j = 0; j < ny_local; j++) {
                    for (int k = 0; k < nz_local; k++) {
                        int i_global = offset_x + i;
                        int j_global = offset_y + j;
                        int k_global = offset_z + k;

                        int idx_local = i * ny_local * nz_local + j * nz_local + k;
                        // int idx_global = i_global * Ny + j_global + k_global;
                        // Scaling
                        x = i_global * d_x;
                        y = j_global * d_y;
                        z = k_global * d_z;

                        // Translation
                        x += x_min;
                        y += y_min;
                        z += z_min;
                        // data[idx_local] = simple_func(x, y);
                        // data[idx_local] = f(x, y, 0, 30, 4) + rank * 30;
                        data[idx_local] = f(x, y, z);  //+ rank * 80;
                    }
                }
            }
        };

        // Fields:
        int Nx;
        int Ny;
        int Nz;
        int nx_local;
        int ny_local;
        int nz_local;
        int mod_x;
        int mod_y;
        int mod_z;
        int offset_x = 0;
        int offset_y = 0;
        int offset_z = 0;
        double x;
        double y;
        double z;
        // double x_min = -0.375;
        // double x_max = 0;
        // double y_min = 0.375;
        // double y_max = 0.75;
        double z_min = 0;
        double z_max = 0;

        // Normal for presentation.
        double x_min = -2.5;
        double x_max = 1.5;
        double y_min = -1.5;
        double y_max = 1.5;
        std::vector<double> data;
        std::vector<int> ranks;
    };
}  // namespace mars

#endif