#include <adios2.h>
#include <err.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_prelepp_benchmark.hpp"
#include "mars_quality.hpp"
#include "mars_simplex.hpp"
#include "mars_utils.hpp"
#include "mars_vtk_writer.hpp"

#include "mars_longest_edge.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_test.hpp"

#include "mars_env.hpp"

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#endif  // WITH_KOKKOS

// #include "mars_fitzhugh_nagumo.hpp"
#include "mars_mesh_kokkos.hpp"
#include "mars_model_test.hpp"
#include "mars_serial_mesh_type.hpp"
#include "mars_spacetime_ex.hpp"
#include "mars_umesh_laplace.hpp"
#include "mars_umesh_st_heat_equation.hpp"

#include "mars_poisson.hpp"
#include "mars_poisson_operator.hpp"

using namespace std::chrono;

namespace mars {

    // double simple_func(const double &x, const double &y, const double &z) { return x * y * z; }

    // void create_data(const int &Nx, const int &Ny, const int &Nz) {
    //     int size = Nx * Ny * Nz;
    //     double H[3];
    //     H[0] = 0.5;
    //     H[1] = 0.25;
    //     H[2] = 1;
    //     double data[size];
    //     double x;
    //     double y;
    //     double z;

    //     for (int i = 0; i < Nx; ++i) {
    //         for (int j = 0; j < Ny; ++j) {
    //             for (int k = 0; k < Nz; ++k) {
    //                 x = i * H[0];
    //                 y = j * H[1];
    //                 z = k * H[2];
    //                 data[i * Ny * Nz + j * Nz + k] = simple_func(x, y, z);
    //             }
    //         }
    //     }
    // }

    // using VectorReal = mars::ViewVectorType<Real>;

    bool write(const std::string &path) {
        // Settings settings;
        // adios2::ADIOS adios(settings.adios_config, adios2::DebugON);
        // adios2::IO io_main = adios.DeclareIO("SimulationOutput");
        // ImageWriter main_writer(settings, io_main);

        // main_writer.open(settings.output);
        // main_writer.close();
    }

    class ST3Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_exact(p); }
    };

    class ST3RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_spacetime(p); }
    };

    template <class Mesh>
    class ST3BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex3_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };
}  // namespace mars

int main(int argc, char *argv[]) {
    // int Nx, Ny, Nz;
    // Nx = 3;
    // Ny = 3;
    // Nz = 3;
    // int size = Nx * Ny * Nz;
    // double H[3];
    // H[0] = 0.5;
    // H[1] = 0.25;
    // H[2] = 1;
    // double data[size];
    // double x;
    // double y;
    // double z;

    // for (int i = 0; i < Nx; ++i) {
    //     for (int j = 0; j < Ny; ++j) {
    //         for (int k = 0; k < Nz; ++k) {
    //             x = i * H[0];
    //             y = j * H[1];
    //             z = k * H[2];
    //             data[i * Ny * Nz + j * Nz + k] = mars::simple_func(x, y, z);
    //         }
    //     }
    // }

    Settings settings;
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");

    ImageWriter main_image(settings, io_main);

    main_image.new_data(2, 2, 1);
    main_image.open(settings.output);
    main_image.write(1);
    main_image.close();

    //
    //
    //
    //
    //
    //
    //
    //
    //
}