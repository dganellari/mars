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
        // using PMesh = ParallelMesh2;
        // using SMesh = Mesh2;

        using PMesh = ParallelMesh3;
        using SMesh = Mesh3;

        using Elem = typename PMesh::Elem;
        using SideElem = typename PMesh::SideElem;

        Integer nx = 6, ny = 6, nz = (PMesh::Dim > 2) ? 6 : 0;
        if (argc > 1) {
            nx = atol(argv[1]);
            ny = atol(argv[1]);
            nz = atol(argv[1]);

            if (PMesh::Dim <= 2) nz = 0;
        }

        PMesh mesh;
        generate_cube(mesh, nx, ny, nz);

        // ParallelBisection<PMesh> bisection(&mesh);

        // Integer n_marked = 3;
        // ViewVectorType<Integer> marked("marked", n_marked);

        // Kokkos::parallel_for(
        //     n_marked, MARS_LAMBDA(const Integer i) { marked(i) = i; });

        // bisection.refine(marked);

        ZeroDirchletOnUnitCube<PMesh> bc_fun;
        One<PMesh> rhs_fun;
        One<PMesh> an_fun;  // FIXME

        // ZeroDirchletOnUnitCube<PMesh> bc_fun;
        // Norm2Squared<PMesh> rhs_fun;

        // Example1Dirichlet bc_fun;
        // Example1RHS rhs_fun;

        // Example2Dirichlet bc_fun;
        // Example2RHS rhs_fun;
        // Example2Analitcal an_fun;

        const Integer n_nodes = mesh.n_nodes();

        ViewVectorType<Real> x("X", n_nodes);
        ViewVectorType<Real> rhs("rhs", n_nodes);
        ViewVectorType<Real> Ax("Ax", n_nodes);

        UMeshLaplace<PMesh> op(mesh);
        BoundaryConditions<PMesh> bc(mesh);
        op.init();

        auto id = std::make_shared<IdentityOperator>();
        id->init(mesh, bc_fun);

        op.set_identity(id);
        op.assemble_rhs(rhs, rhs_fun);

        Real nrm_rhs = KokkosBlas::nrm1(rhs);

        std::cout << "nrm_rhs : " << nrm_rhs << std::endl;

        bc.apply(rhs, bc_fun);
        bc.apply(x, bc_fun);

        auto prec = op.preconditioner();

        Integer num_iter = 0;
        bcg_stab(op, prec, rhs, rhs.extent(0), x, num_iter);

        // Compute Error
        ViewVectorType<Real> x_exact("X_exact", n_nodes);
        ViewVectorType<Real> diff("Diff", n_nodes);

        Interpolate<PMesh> interp(mesh);
        interp.apply(x_exact, an_fun);

        Kokkos::deep_copy(diff, x_exact);

        KokkosBlas::axpy(-1.0, x, diff);
        Real err = KokkosBlas::nrminf(diff);
        std::cout << "err : " << err << std::endl;

        ///////////////////////////////////////////////////////////////////////////

        ViewVectorType<Real>::HostMirror x_host("x_host", n_nodes);
        ViewVectorType<Real>::HostMirror rhs_host("rhs_host", n_nodes);
        Kokkos::deep_copy(x_host, x);
        Kokkos::deep_copy(rhs_host, rhs);

        SMesh serial_mesh;
        convert_parallel_mesh_to_serial(serial_mesh, mesh);

        std::cout << "n_active_elements: " << serial_mesh.n_active_elements() << std::endl;
        std::cout << "n_nodes:           " << serial_mesh.n_nodes() << std::endl;

        VTUMeshWriter<SMesh> w;
        w.write("mesh.vtu", serial_mesh, x_host);
        w.write("mesh_rhs.vtu", serial_mesh, rhs_host);

        Kokkos::deep_copy(x_host, x_exact);
        w.write("analitic.vtu", serial_mesh, x_host);
    }

    Kokkos::finalize();
}
