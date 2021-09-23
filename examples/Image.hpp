#ifndef IMAGE_HPP
#define IMAGE_HPP
#include <vector>
#include "Mandelbrot.hpp"
#include "adios2.h"
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

            nx_local += coordinates[0] < mod_x;
            ny_local += coordinates[1] < mod_y;

            // // The only varying part of this code is the coordinates of the processor. i.e (i_proc).
            for (int i_proc = 0; i_proc < coordinates[0]; i_proc++) {
                offset_x += Nx / dims[0] + (i_proc < mod_x);
            }

            for (int j_proc = 0; j_proc < coordinates[1]; j_proc++) {
                offset_y += Ny / dims[1] + (j_proc < mod_y);
            }
        };

        template <class Function>
        void calc_function(Function f) {
            double d_x = (x_max - x_min) / (Nx - 1);
            double d_y = (y_max - y_min) / (Ny - 1);
            data.resize(Nx * Ny);

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
                    data[idx_local] = f(x, y, 0);  //+ rank * 80;
                }
            }
        };

        // Fields:
        int Nx = 7680;
        int Ny = 4320;
        int Nz;
        int nx_local;
        int ny_local;
        int mod_x;
        int mod_y;
        int offset_x = 0;
        int offset_y = 0;
        double x;
        double y;
        double x_min = -0.375;
        double x_max = 0;
        double y_min = 0.375;
        double y_max = 0.75;
        std::vector<double> data;
    };
}  // namespace mars

#endif