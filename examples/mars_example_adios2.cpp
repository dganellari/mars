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

    void create_data(const int &Nx, const int &Ny, const int &Nz) {
        int size = Nx * Ny * Nz;
        int X[size];

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                }
            }
        }
    }

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
    std::cout << "Writing results to disk..." << std::endl;

    using namespace mars;

    Env env(argc, argv);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    {
        using namespace cxxopts;
        Options options("./adios_example", "Run M.A.R.S. based applications.");

        options.add_options()("d,debug", "Enable debugging")                                     //
            ("l,level", "Number of levels", value<Integer>()->default_value("1"))                //
            ("x,nx", "Number of elements in x direction", value<Integer>()->default_value("6"))  //
            ("y,ny", "Number of elements in y direction", value<Integer>()->default_value("6"))  //
            ("z,nz", "Number of elements in z direction", value<Integer>()->default_value("6"))  //
            ("t,nt", "Number of elements in t direction", value<Integer>()->default_value("6"))  //
            ("a,adaptive", "Adaptivity", value<bool>()->default_value("false"))                  //
            ("o,output", "Enable output", value<bool>()->default_value("true"))                  //
            ("r,refine_level",
             "Number of refinements",
             value<Integer>()->default_value("1"))                                  //
            ("v,verbose", "Verbose output", value<bool>()->default_value("false"))  //
            ("h,help", "Print usage");

        auto args = options.parse(argc, argv);

        ModelTest<ParallelQuad4Mesh,
                  UMeshSTHeatEquation<ParallelQuad4Mesh>,
                  ST3BC<ParallelQuad4Mesh>,
                  ST3RHS,
                  ST3Analitcal>()
            .run(args);
    }
    return env.exit_code();
}