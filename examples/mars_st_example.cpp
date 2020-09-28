#include <err.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>

#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
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
#include "mars_memory.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_test.hpp"

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#endif  // WITH_KOKKOS

#include "mars_mesh_kokkos.hpp"
#include "mars_model_test.hpp"
#include "mars_spacetime_ex.hpp"
#include "mars_umesh_laplace.hpp"
#include "mars_umesh_st_heat_equation.hpp"

namespace mars {
    class ST1Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex1_st_exact(p); }
    };

    class ST1RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex1_st_spacetime(p); }
    };

    template <class Mesh>
    class ST1BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex1_st_exact(p);
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

    /////////////////////////////////////////////////////////////////////////////////////////////

    class ST2Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_st_exact(p); }
    };

    class ST2RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_st_spacetime(p); }
    };

    template <class Mesh>
    class ST2BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex2_st_exact(p);
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

    /////////////////////////////////////////////////////////////////////////////////////////////

    class ST3Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_exact(p); }
    };

    class ST3RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_spacetime(p); }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace mars

int main(int argc, char *argv[]) {
    using namespace mars;
    Kokkos::initialize(argc, argv);

#ifdef MARS_USE_CUDA
    cudaDeviceSetLimit(cudaLimitStackSize,
                       32768);  // set stack to 32KB only for cuda since it is
                                // not yet supported in kokkos.
#endif

    {
        // ModelTest<ParallelMesh2, UMeshSTHeatEquation<ParallelMesh2>, ST1BC<ParallelMesh2>, ST1RHS,
        // ST1Analitcal>().run(
        //     argc, argv);

        ModelTest<ParallelMesh2, UMeshLaplace<ParallelMesh2>, Example2Dirichlet, Example2RHS, Example2Analitcal>().run(
            argc, argv);

        // ModelTest<ParallelMesh2, UMeshSTHeatEquation<ParallelMesh2>, ST2BC<ParallelMesh2>, ST2RHS,
        // ST2Analitcal>().run(
        //     argc, argv);

        // ModelTest<ParallelMesh2,
        //           UMeshSTHeatEquation<ParallelMesh2>,
        //           ZeroDirchletOnUnitCube<ParallelMesh2>,
        //           ST3RHS,
        //           ST3Analitcal>()
        //     .run(argc, argv);

        // ModelTest<ParallelMesh3,
        //           UMeshLaplace<ParallelMesh3>,
        //           ZeroDirchletOnUnitCube<ParallelMesh3>,
        //           One<ParallelMesh3>,
        //           One<ParallelMesh3>>()
        //     .run(argc, argv);
    }

    Kokkos::finalize();
}
