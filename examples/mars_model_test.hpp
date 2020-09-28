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
#include <KokkosBlas1_sum.hpp>
#include <Kokkos_sort.hpp>

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
        using VectorReal = mars::ViewVectorType<Real>;
        using VectorInt = mars::ViewVectorType<Integer>;

        // void adaptive_refinement(FEValues<PMesh> &values, VectorReal &x) {
        //     auto mesh = values.mesh();

        //     ParallelBisection<PMesh> bisection(&mesh);

        //     VectorReal error;
        //     GradientRecovery<PMesh> grad_rec;
        //     grad_rec.estimate(values, x, error);

        //     const Real alpha = 0.3;
        //     const Real max_error = KokkosBlas::nrminf(error);
        //     const Real low_bound_err = alpha * max_error;

        //     VectorInt marked("marked", mesh.n_elements());

        //     Integer n_elements = mesh.n_elements();

        //     Kokkos::parallel_for(
        //         n_elements, MARS_LAMBDA(const Integer &i) { marked(i) = (error(i) >= low_bound_err) * i; });

        //     Integer n_marked = 0;

        //     Kokkos::parallel_reduce(
        //         n_elements, MARS_LAMBDA(const Integer &i, Integer &val) { val += marked(i) > 0; }, n_marked);

        //     VectorInt marked_list("maked_list", n_marked);

        //     Kokkos::sort(marked);

        //     // Kokkos::parallel_for(
        //     //     n_elements, MARS_LAMBDA(const Integer &i) { printf("%ld ", marked(i)); });
        //     Integer offset = n_elements - n_marked;

        //     Kokkos::parallel_for(
        //         n_marked, MARS_LAMBDA(const Integer i) { marked_list(i) = marked(offset + i); });

        //     // Kokkos::parallel_for(
        //     //     n_marked, MARS_LAMBDA(const Integer &i) { printf("%ld ", marked_list(i)); });

        //     bisection.refine(marked_list);

        //     auto children = mesh.get_view_children();

        //     std::cout << "nc: " << children.extent(0) << "," << children.extent(1) << std::endl;

        //     Kokkos::parallel_for(
        //         children.extent(0),
        //         MARS_LAMBDA(const Integer &i) { printf("(%ld, %ld)\n", children(i, 0), children(i, 1)); });

        //     VectorInt out("marked", mesh.n_nodes());
        //     Kokkos::parallel_for(
        //         mesh.n_nodes(), MARS_LAMBDA(const Integer &i) { out(i) = i; });

        //     VectorReal::HostMirror out_host("out_host", mesh.n_nodes());
        //     Kokkos::deep_copy(out_host, out);

        //     std::cout << "n_elements: " << n_elements << " < " << mesh.n_elements() << std::endl;

        //     SMesh serial_mesh;
        //     convert_parallel_mesh_to_serial(serial_mesh, mesh);

        //     std::cout << "n_active_elements: " << serial_mesh.n_active_elements() << std::endl;
        //     std::cout << "n_nodes:           " << serial_mesh.n_nodes() << std::endl;

        //     VTUMeshWriter<SMesh> w;
        //     w.write("mesh_refined.vtu", serial_mesh, out_host);
        // }

        void run(int argc, char *argv[]) {
            using Elem = typename PMesh::Elem;
            using SideElem = typename PMesh::SideElem;

            Integer ns[4] = {0, 0, 0, 0};

            Integer n = 6;
            bool write_output = true;

            if (argc > 1) {
                n = atol(argv[1]);
            }

            if (argc > 2) {
                write_output = atoi(argv[2]);
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

            BC bc_fun;
            RHS rhs_fun;
            AnalyticalFun an_fun;

            const Integer n_nodes = mesh.n_nodes();

            VectorReal x("X", n_nodes);
            VectorReal rhs("rhs", n_nodes);
            VectorReal Ax("Ax", n_nodes);

            Op op(mesh);
            BoundaryConditions<PMesh> bc(mesh);
            op.init();

            auto id = std::make_shared<IdentityOperator>();
            id->init(mesh, bc_fun);

            op.set_identity(id);
            op.assemble_rhs(rhs, rhs_fun);

            Real sum_rhs = KokkosBlas::sum(rhs);

            std::cout << "sum_rhs : " << sum_rhs << std::endl;

            bc.apply(rhs, bc_fun);
            bc.apply(x, bc_fun);

            auto prec_ptr = op.preconditioner();

            Integer num_iter = 0;
            bcg_stab(op, *prec_ptr, rhs, 10 * rhs.extent(0), x, num_iter);

            /////////////////////////////////////////////////////////////////////////////
            // Compute Error
            VectorReal x_exact("X_exact", n_nodes);
            VectorReal diff("Diff", n_nodes);

            Interpolate<PMesh> interp(mesh);
            interp.apply(x_exact, an_fun);

            Kokkos::deep_copy(diff, x_exact);

            KokkosBlas::axpy(-1.0, x, diff);
            Real err = KokkosBlas::nrm2(diff) / KokkosBlas::nrm2(x_exact);
            std::cout << "err : " << err << std::endl;

            ////////////////////////////////////////////////////////////////////////////
            // Compute op (x_exact) - rhs

            op.apply(x_exact, diff);

            sum_rhs = KokkosBlas::sum(rhs);

            Real norm_lapl_x = KokkosBlas::sum(diff);
            std::cout << "norm_lapl_x: " << norm_lapl_x << " == " << sum_rhs << std::endl;

            // Diff

            KokkosBlas::axpy(-1.0, rhs, diff);
            err = KokkosBlas::nrm2(diff) / KokkosBlas::nrm2(rhs);
            std::cout << "Op error : " << err << std::endl;

            ///////////////////////////////////////////////////////////////////////////

            VectorReal error;
            GradientRecovery<PMesh> grad_rec;
            grad_rec.estimate(op.values(), x, error);

            ///////////////////////////////////////////////////////////////////////////

            if (write_output) {
                std::cout << "Writing results to disk..." << std::endl;

                VectorReal::HostMirror x_host("x_host", n_nodes);
                VectorReal::HostMirror rhs_host("rhs_host", n_nodes);
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

                VectorReal::HostMirror error_host("error_host", mesh.n_elements());
                Kokkos::deep_copy(error_host, error);
                w.write("error.vtu", serial_mesh, error_host, true);

                VectorReal::HostMirror diff_host("lapl_x_host", mesh.n_nodes());
                Kokkos::deep_copy(diff_host, diff);
                w.write("lapl_x.vtu", serial_mesh, diff_host);
            }

            // adaptive_refinement(op.values(), x);
        }
    };
}  // namespace mars

#endif  // MARS_MODEL_TEST_HPP