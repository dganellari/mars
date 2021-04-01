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

    // using VectorReal = mars::ViewVectorType<Real>;

    template <class PMesh>
    bool write(const std::string &path) {
        Settings settings;
        adios2::ADIOS adios(settings.adios_config, adios2::DebugON);
        adios2::IO io_main = adios.DeclareIO("SimulationOutput");
        ImageWriter main_writer(settings, io_main);

        main_writer.open(settings.output);
        main_writer.close();
        // int rank = 0;
        // adios2::ADIOS adios;
        // const int NSTEPS = 5;
        // // random size per process, a different size at each step
        // unsigned int Nelems;
        // const size_t Nglobal = 6;
        // // Application variables for output
        // // random size per process, 5..10 each
        // // v1 has different size on each process (but fixed over time)
        // const unsigned int Nx = rand() % 6 + 5;

        // std::vector<double> v0(Nglobal);
        // std::vector<double> v1(Nx);
        // // Local array, size is changing over time on each process
        // std::vector<double> v2;

        // // Local array, size is changing over time on each process
        // // Also, random number of processes will write it at each step
        // std::vector<double> &v3 = v2;

        // using SMesh = typename SerialMeshType<PMesh>::Type;
        // std::cout << "Writing results to disk..." << std::endl;

        // Integer n_nodes = mesh.n_nodes();

        // VectorReal::HostMirror x_host("x_host", n_nodes);
        // Kokkos::deep_copy(x_host, x);

        // SMesh serial_mesh;
        // convert_parallel_mesh_to_serial(serial_mesh, mesh);

        // // w.write(path, serial_mesh, x_host)

        // // Get io settings from the config file or
        // // create one with default settings here
        // adios2::IO io = adios.DeclareIO("Output");
        // io.SetEngine("BP3");
        // io.SetParameters({{"verbose", "4"}});

        // /*
        //  * Define local array: type, name, local size
        //  * Global dimension and starting offset must be an empty vector
        //  * Here the size of the local array is the same on every process
        //  */
        // adios2::Variable<double> varV0 = io.DefineVariable<double>("v0", {}, {}, {Nglobal});

        // /*
        //  * v1 is similar to v0 but on every process the local size
        //  * is a different value
        //  */
        // adios2::Variable<double> varV1 = io.DefineVariable<double>("v1", {}, {}, {Nx});

        //  * Define local array: type, name
        //  * Global dimension and starting offset must be an empty vector
        //  * but local size CANNOT be an empty vector.
        //  * We can use {adios2::UnknownDim} for this purpose or any number
        //  * actually since we will modify it before writing

        // adios2::Variable<double> varV2 = io.DefineVariable<double>("v2", {}, {}, {adios2::UnknownDim});

        // /*
        //  * v3 is just like v2
        //  */
        // adios2::Variable<double> varV3 = io.DefineVariable<double>("v3", {}, {}, {adios2::UnknownDim});

        // adios2::Engine writer = io.Open("solution.bp", adios2::Mode::Write);

        // for (int step = 0; step < NSTEPS; step++) {
        //     writer.BeginStep();
        //     // Write Step

        //     // // v0
        //     // for (size_t i = 0; i < Nglobal; i++) {
        //     //     v0[i] = rank * 1.0 + step * 0.1;
        //     // }
        //     // writer.Put<double>(varV0, v0.data());

        //     // // v1
        //     // for (size_t i = 0; i < Nx; i++) {
        //     //     v1[i] = rank * 1.0 + step * 0.1;
        //     // }
        //     // writer.Put<double>(varV1, v1.data());

        //     // v2

        //     // random size per process per step, 5..10 each
        //     std::cout << x.data() << std::endl;
        //     // Nelems = x.size();
        //     // x.reserve(Nelems);
        //     for (size_t i = 0; i < Nelems; i++) {
        //         v2[i] = rank * 1.0 + step * 0.1;
        //     }

        //     // Set the size of the array now because we did not know
        //     // the size at the time of definition
        //     varV2.SetSelection(adios2::Box<adios2::Dims>({}, {Nelems}));
        //     // writer.Put<double>(varV2, x.data());

        //     // v3

        //     // // random chance who writes it
        //     // unsigned int chance = rand() % 100;
        //     // /*if (step == 2)
        //     // {
        //     //     chance = 0;
        //     // }*/
        //     // bool doWrite = (chance > 60);
        //     // if (doWrite) {
        //     //     varV3.SetSelection(adios2::Box<adios2::Dims>({}, {Nelems}));
        //     //     writer.Put<double>(varV3, v3.data());
        //     // }

        //     writer.EndStep();
        // }
        // writer.Close();
        // return true;
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