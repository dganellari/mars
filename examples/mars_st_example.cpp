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

#include "mars_matrix_free_operator.hpp"
#include "mars_poisson.hpp"

#include "mars_boundary_conditions.hpp"
#include "mars_fe_values.hpp"
#include "mars_identity_operator.hpp"
#include "mars_interpolate.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_model_test.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"
#include "vtu_writer.hpp"

namespace mars {

    template <class Mesh>
    class SpaceTimeMixed {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        MARS_INLINE_FUNCTION static void one_thread_eval_diag_add(const Real *J_inv, const Real &det_J, Real *val) {
            // Real g_ref[Dim], g[Dim];

            // for (int d = 0; d < Dim; ++d) {
            //     g_ref[d] = -1;
            // }

            // m_t_v_mult(J_inv, g_ref, g);
            // val[0] = dot(g, g) * det_J;

            // for (int d = 0; d < Dim; ++d) {
            //     g_ref[d] = 0;
            // }

            // for (int d = 0; d < Dim; ++d) {
            //     g_ref[d] = 1;
            //     m_t_v_mult(J_inv, g_ref, g);

            //     val[d + 1] = dot(g, g) * det_J;

            //     g_ref[d] = 0;
            // }
        }

        MARS_INLINE_FUNCTION static void one_thread_eval_add(const Real *J_inv,
                                                             const Real &det_J,
                                                             const Real *u,
                                                             Real *val) {
            Real g_ref[Dim], g[Dim], g_fe[Dim];

            ///////////////////////////////////////////////////////////////////
            ////////////////// Gradient with local basis function /////////////

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1 * u[0];
            }

            for (int i = 1; i < NFuns; ++i) {
                g_ref[i - 1] += u[i];
            }

            ///////////////////////////////////////////////////////////////////
            ////////////////// Transform gradient to physical coordinates //////

            m_t_v_mult(J_inv, g_ref, g);

            Real ut = g[Dim - 1];

            for (int i = 0; i < NFuns; ++i) {
                val[i] += ut * det_J * 1. / NFuns;
            }
        }
    };

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
        // ModelTest<ParallelMesh2, UMeshLaplace<ParallelMesh2>, Example2Dirichlet, Example2RHS,
        // Example2Analitcal>().run(
        //     argc, argv);

        ModelTest<ParallelMesh3,
                  UMeshLaplace<ParallelMesh3>,
                  ZeroDirchletOnUnitCube<ParallelMesh3>,
                  One<ParallelMesh3>,
                  One<ParallelMesh3>>()
            .run(argc, argv);
    }

    Kokkos::finalize();
}
