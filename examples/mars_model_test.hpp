#ifndef MARS_MODEL_TEST_HPP
#define MARS_MODEL_TEST_HPP

#include <err.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>

#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>

#include "mars_base.hpp"
#include "mars_globals.hpp"

#include "mars_boundary_conditions.hpp"
#include "mars_fe_values.hpp"
#include "mars_gradient_recovery.hpp"
#include "mars_identity_operator.hpp"
#include "mars_interpolate.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"
#include "vtu_writer.hpp"

namespace mars {

    template <class PMesh, class Op, class BC, class RHS, class AnalyticalFun>
    class ModelTest {
    public:
        static const int Dim = PMesh::Dim;
        using SMesh = mars::Mesh<Dim, PMesh::ManifoldDim>;

        void run(int argc, char *argv[]) {
            using Elem = typename PMesh::Elem;
            using SideElem = typename PMesh::SideElem;

            Integer ns[4] = {0, 0, 0, 0};

            Integer n = 6;

            if (argc > 1) {
                n = atol(argv[1]);
            }

            for (int i = 0; i < 4; ++i) {
                ns[i] = n;
            }

            PMesh mesh;

            if (Dim <= 3) {
                generate_cube(mesh, ns[0], ns[1], ns[2]);
            } else {
                std::cerr << "[Error] 4D not supported yet" << std::endl;
                return;
            }

            // ParallelBisection<PMesh> bisection(&mesh);

            // Integer n_marked = 3;
            // ViewVectorType<Integer> marked("marked", n_marked);

            // Kokkos::parallel_for(
            //     n_marked, MARS_LAMBDA(const Integer i) { marked(i) = i; });

            // bisection.refine(marked);

            BC bc_fun;
            RHS rhs_fun;
            AnalyticalFun an_fun;

            const Integer n_nodes = mesh.n_nodes();

            ViewVectorType<Real> x("X", n_nodes);
            ViewVectorType<Real> rhs("rhs", n_nodes);
            ViewVectorType<Real> Ax("Ax", n_nodes);

            Op op(mesh);
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

            auto prec_ptr = op.preconditioner();

            Integer num_iter = 0;
            bcg_stab(op, *prec_ptr, rhs, rhs.extent(0), x, num_iter);

            // Compute Error
            ViewVectorType<Real> x_exact("X_exact", n_nodes);
            ViewVectorType<Real> diff("Diff", n_nodes);

            Interpolate<PMesh> interp(mesh);
            interp.apply(x_exact, an_fun);

            Kokkos::deep_copy(diff, x_exact);

            KokkosBlas::axpy(-1.0, x, diff);
            Real err = KokkosBlas::nrminf(diff);
            std::cout << "err : " << err << std::endl;

            ViewVectorType<Real> error;
            GradientRecovery<PMesh> grad_rec;
            grad_rec.estimate(op.values(), x, error);

            ///////////////////////////////////////////////////////////////////////////

            std::cout << "Writing results to disk..." << std::endl;

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

            ViewVectorType<Real>::HostMirror error_host("error_host", mesh.n_elements());
            Kokkos::deep_copy(error_host, error);
            w.write("error.vtu", serial_mesh, error_host, true);
        }
    };
}  // namespace mars

#endif  // MARS_MODEL_TEST_HPP