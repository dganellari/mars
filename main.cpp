#include <algorithm>
#include <cassert>
#include <cstdlib>

#include "mars_config.hpp"
#include "mars_err.hpp"

#ifdef MARS_ENABLE_CXXOPTS
#include "cxxopts.hpp"
#endif
#include <iostream>
#include <numeric>

#include "mars.hpp"
#include "mars_env.hpp"

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include "mars_advection.hpp"
#include "mars_constant_viscosity_stokes.hpp"
#include "mars_poisson.hpp"
#ifdef MARS_ENABLE_AMR_BACKEND
#include "mars_test_kokkos.hpp"
#endif
#include "mars_test_mpi.hpp"
#include "mars_variable_viscosity_stokes.hpp"
#endif  // MARS_ENABLE_KOKKOS_KERNELS
#endif  // MARS_ENABLE_MPI
#include <chrono>

using namespace std::chrono;

#ifdef MARS_ENABLE_AMR_BACKEND
void run_benchmarks(int level, int refine_level) {
    using namespace mars;

    std::cout << "Generation level:" << level << std::endl;
    std::cout << "Refinement level:" << refine_level << std::endl;

#ifdef MARS_ENABLE_KOKKOS_KERNELS

    /* ParallelQuad4Mesh nsm;
    generate_cube(nsm, level, refine_level, 0); */
    ParallelMesh3 pMesh3;
    generate_cube(pMesh3, level, level, level);

    ParallelLeppBenchmark<ParallelMesh3> pb;
    pb.run(refine_level, pMesh3, "pb");

    /* Mesh3 sMesh3;
    convert_parallel_mesh_to_serial(sMesh3, pMesh3);

    PreLeppBenchmark<Mesh3> b3;
    b3.run(refine_level, sMesh3, "b3"); */
#endif
}
#endif

int main(int argc, char *argv[]) {
    using namespace mars;
#ifdef MARS_ENABLE_CXXOPTS
    using namespace cxxopts;

    /* MARS::init(argc, argv); */

    Env env(argc, argv);

    Options options("./mars_exec", "Run M.A.R.S. based applications.");

    options.add_options()("d,debug", "Enable debugging")                                //
        ("x,xDim", "Grid X Dim", value<int>()->default_value("4"))                      //
        ("y,yDim", "Grid Y Dim", value<int>()->default_value("4"))                      //
        ("z,zDim", "Grid Z Dim", value<int>()->default_value("4"))                      //
        ("b,block_size", "Vector Valued Block Size", value<int>()->default_value("1"))  //
        ("l,level", "Number of levels", value<int>()->default_value("1"))               //
        ("r,refine_level",
         "Number of refinements",
         value<int>()->default_value("1"))  //
        ("f,file",
         "File name",
         value<std::string>()->default_value("../data/write/tetrakis.MFEM"))    //
        ("a,app", "Application", value<std::string>()->default_value(""))       //
        ("v,verbose", "Verbose output", value<bool>()->default_value("false"))  //
        ("h,help", "Print usage");

    try {
        auto args = options.parse(argc, argv);

        if (args.count("help")) {
            std::cout << options.help() << std::endl;
        }

        int xDim = args["xDim"].as<int>();
        int yDim = args["yDim"].as<int>();
        int zDim = args["zDim"].as<int>();
        int block_size = args["block_size"].as<int>();
        int level = args["level"].as<int>();
        int refine_level = args["refine_level"].as<int>();
        std::string filename = args["file"].as<std::string>();
        std::string app = args["app"].as<std::string>();

        ///////////////////////////////////////////////////
        // FIXME create tests, benchmarks and separate apps for what is below
        std::map<std::string, std::function<void()>> apps;

#ifdef MARS_ENABLE_KOKKOS_KERNELS

        apps["mars_mesh_generation_kokkos_2D_a"] = [=]() { test_mars_mesh_generation_kokkos_2D(2, 4); };
        apps["mars_mesh_generation_kokkos_2D_b"] = [=]() { test_mars_mesh_generation_kokkos_2D(level + 4, level); };

        apps["mars_mesh_generation_kokkos_3D"] = [=]() { test_mars_mesh_generation_kokkos_3D(level, level, level); };
        apps["mars_mesh_generation_kokkos_1D"] = [=]() { test_mars_mesh_generation_kokkos_1D(level); };

#ifdef MARS_ENABLE_MPI

        apps["mars_distributed_mesh_generation_2D"] = [=]() {
            test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(xDim, yDim);
        };

        apps["mars_distributed_mesh_generation_hilbert_2D"] = [=]() {
            test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D<HilbertKey<Unsigned>>(xDim, yDim);
        };

        apps["mars_distributed_mesh_generation_3D"] = [=]() {
            test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(xDim, yDim, zDim);
        };

        apps["mars_distributed_mesh_generation_hilbert_3D"] = [=]() {
            test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D<HilbertKey<Unsigned>>(xDim, yDim, zDim);
        };

        apps["distributed_vector_valued"] = [=]() { test_mars_distributed_vector_valued(xDim, yDim, 0, block_size); };

        apps["distributed_vector_valued_hilbert"] = [=]() {
            test_mars_distributed_vector_valued<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(
                xDim, yDim, 0, block_size);
        };

        apps["distributed_vector_valued_degree2"] = [=]() {
            test_mars_distributed_vector_valued<ElementType::Quad4, 2>(xDim, yDim, 0, block_size);
        };

        apps["distributed_vector_valued_no_overlap_3D"] = [=]() {
            test_mars_distributed_vector_valued<ElementType::Hex8, 1, false>(xDim, yDim, zDim, block_size);
        };

        apps["distributed_vector_valued_3D"] = [=]() {
            test_mars_distributed_vector_valued<ElementType::Hex8>(xDim, yDim, zDim, block_size);
        };

        apps["distributed_vector_valued_3D_degree2"] = [=]() {
            test_mars_distributed_vector_valued<ElementType::Hex8, 2>(xDim, yDim, zDim, block_size);
        };

        apps["mfpoisson"] = [=]() {
            matrix_free_poisson<Example2Dirichlet, Example2RHS, Example2Analitcal, ElementType::Quad4>(xDim, yDim, 0);
        };

        // Quadrature for Hex8 not added yet. Only FEQuad4 currently available in Mars.
        /* apps["mfpoisson3D"] = [=]() {
            matrix_free_poisson<Example2Dirichlet, Example2RHS, Example2Analitcal, ElementType::Hex8>(xDim, yDim, zDim);
        }; */

        apps["cstokes"] = [=]() { staggered_constant_viscosty_stokes<ElementType::Quad4>(xDim, yDim, 0); };

        /* apps["cstokes3D"] = [=]() { staggered_constant_viscosty_stokes<ElementType::Hex8>(xDim, yDim, zDim); }; */

        apps["vstokes"] = [=]() { staggered_variable_viscosty_stokes<ElementType::Quad4>(xDim, yDim, 0); };

        /* apps["vstokes3D"] = [=]() { staggered_variable_viscosty_stokes<ElementType::Hex8>(xDim, yDim, zDim); }; */

        apps["advection"] = [=]() { advection(xDim, yDim); };

#endif
#endif  // MARS_ENABLE_KOKKOS_KERNELS

        if (!app.empty()) {
            auto it = apps.find(app);
            if (it == apps.end()) {
                std::cerr << "Could not find app: " << app << std::endl;
                std::cout << "Avaialble apps:\n";

                for (auto &a : apps) {
                    std::cout << a.first << "\n";
                }

            } else {
                it->second();
            }
        }

    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        std::cout << options.help() << std::endl;
    }

    return env.exit_code();
#endif
}
