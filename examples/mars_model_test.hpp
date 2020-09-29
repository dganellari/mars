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
// #include <Kokkos_sort.hpp>

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
        static const int ManifoldDim = PMesh::ManifoldDim;
        static const int NFuns = ManifoldDim + 1;

        using SMesh = mars::Mesh<Dim, ManifoldDim>;
        using VectorReal = mars::ViewVectorType<Real>;
        using VectorInt = mars::ViewVectorType<Integer>;
        using VectorBool = mars::ViewVectorType<bool>;

        void adaptive_refinement(FEValues<PMesh> &values, VectorReal &x) {
            auto mesh = values.mesh();

            ParallelBisection<PMesh> bisection(&mesh);

            VectorReal error;
            GradientRecovery<PMesh> grad_rec;
            grad_rec.estimate(values, x, error);

            const Real alpha = 0.3;
            const Real max_error = KokkosBlas::nrminf(error);
            const Real low_bound_err = alpha * max_error;

            Integer old_n_elements = mesh.n_elements();

            VectorBool marked("marked", mesh.n_elements());

            Kokkos::parallel_for(
                mesh.n_elements(), MARS_LAMBDA(const Integer &i) { marked(i) = (error(i) >= low_bound_err) * i; });

            VectorInt marked_list = mark_active(marked);

            // Kokkos::parallel_for(
            //     marked_list.extent(0), MARS_LAMBDA(const Integer &i) { printf("%ld ", marked_list(i)); });

            bisection.refine(marked_list);

            // VectorBool active = mesh.get_view_active();
            // VectorBool inactive("inactive", active.extent(0));

            // Kokkos::parallel_for(
            //     active.extent(0), MARS_LAMBDA(const Integer &i) { inactive(i) = !active(i); });

            VectorInt parents("parents", mesh.n_elements(), -1);

            Kokkos::parallel_for(
                old_n_elements, MARS_LAMBDA(const Integer &i) {
                    auto c = mesh.get_children(i);
                    if (c.is_valid()) {
                        parents(c(0)) = i;
                        parents(c(1)) = i;
                    }
                });

            Integer n_new = mesh.n_elements() - old_n_elements;

            // Kokkos::parallel_for(
            //     n_new,
            //     MARS_LAMBDA(const Integer &i) { printf("%ld/%ld\n", old_n_elements + i, old_n_elements + n_new); });

            VectorReal refined_x("refined_x", mesh.n_nodes());

            {
                auto elems = mesh.get_view_elements();
                auto points = mesh.get_view_points();
                auto J_inv = values.J_inv();

                Kokkos::parallel_for(
                    n_new, MARS_LAMBDA(const Integer &iter) {
                        const Integer elem_id = old_n_elements + iter;
                        Integer parent = parents(elem_id);

                        Real J_inv_e[Dim * Dim];
                        Real tr[Dim], p[Dim], p_ref[Dim];
                        Integer idx[NFuns];

                        for (int k = 0; k < Dim * Dim; ++k) {
                            J_inv_e[k] = J_inv(parent, k);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            idx[k] = elems(elem_id, k);
                        }

                        Integer n0 = elems(parent, 0);

                        for (int d = 0; d < Dim; ++d) {
                            tr[d] = points(n0, d);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                        }

                        // printf("%ld/%ld\n", old_n_elements + i, old_n_elements + n_new);
                    });
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////

            VectorInt out("marked", mesh.n_nodes());
            Kokkos::parallel_for(
                mesh.n_nodes(), MARS_LAMBDA(const Integer &i) { out(i) = i; });

            VectorReal::HostMirror out_host("out_host", mesh.n_nodes());
            Kokkos::deep_copy(out_host, out);

            SMesh serial_mesh;
            convert_parallel_mesh_to_serial(serial_mesh, mesh);

            std::cout << "n_active_elements: " << serial_mesh.n_active_elements() << std::endl;
            std::cout << "n_nodes:           " << serial_mesh.n_nodes() << std::endl;

            VTUMeshWriter<SMesh> w;
            w.write("mesh_refined.vtu", serial_mesh, out_host);
        }

        void run(int argc, char *argv[]) {
            using Elem = typename PMesh::Elem;
            using SideElem = typename PMesh::SideElem;

            PMesh mesh;

            Integer n = 6;
            bool write_output = true;

            if (argc > 1) {
                n = atol(argv[1]);
            }

            if (argc > 2) {
                write_output = atoi(argv[2]);
            }

            if (Dim <= 3) {
                Integer ns[4] = {0, 0, 0, 0};
                for (int i = 0; i < 4; ++i) {
                    ns[i] = n;
                }

                if (argc > 3) {
                    Integer mult = atoi(argv[3]);
                    if (mult) {
                        ns[Dim - 1] *= mult;
                    }
                }
                generate_cube(mesh, ns[0], ns[1], ns[2]);
            } else {
                SMesh smesh;
                read_mesh("../data/cube4d_24.MFEM", smesh);

                Bisection<SMesh> b(smesh);
                b.uniform_refine(n);
                b.clear();
                smesh.clean_up();
                convert_serial_mesh_to_parallel(mesh, smesh);
                write_output = false;
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
    };  // namespace mars
}  // namespace mars

#endif  // MARS_MODEL_TEST_HPP