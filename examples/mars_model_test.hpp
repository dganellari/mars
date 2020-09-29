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

        class Problem {
        public:
            Problem(PMesh &mesh) : mesh(mesh), op(mesh), bc(mesh), rhs("rhs", mesh.n_nodes()), write_output(true) {}

            bool solve(VectorReal &x) {
                x = VectorReal("x", mesh.n_nodes());
                bc.apply(x, bc_fun);

                auto prec_ptr = op.preconditioner();

                Integer num_iter = 0;
                if (bcg_stab(op, *prec_ptr, rhs, 10 * rhs.extent(0), x, num_iter)) {
                    if (write_output) {
                        return write(x);
                    } else {
                        return true;
                    }

                } else {
                    return false;
                }
            }

            bool write(VectorReal &x) {
                std::cout << "Writing results to disk..." << std::endl;

                Integer n_nodes = mesh.n_nodes();

                VectorReal::HostMirror x_host("x_host", n_nodes);
                VectorReal::HostMirror rhs_host("rhs_host", n_nodes);
                Kokkos::deep_copy(x_host, x);
                Kokkos::deep_copy(rhs_host, rhs);

                SMesh serial_mesh;
                convert_parallel_mesh_to_serial(serial_mesh, mesh);

                VTUMeshWriter<SMesh> w;

                if (!w.write("solution.vtu", serial_mesh, x_host)) {
                    return false;
                }

                if (!w.write("rhs.vtu", serial_mesh, rhs_host)) {
                    return false;
                }

                return true;
            }

            bool measure_actual_error(VectorReal &x) {
                /////////////////////////////////////////////////////////////////////////////
                // Compute Error
                Integer n_nodes = mesh.n_nodes();

                VectorReal x_exact("X_exact", n_nodes);
                VectorReal diff("Diff", n_nodes);
                VectorReal Ax("Ax", n_nodes);

                Interpolate<PMesh> interp(mesh);
                interp.apply(x_exact, an_fun);

                Kokkos::deep_copy(diff, x_exact);

                KokkosBlas::axpy(-1.0, x, diff);
                Real err = KokkosBlas::nrm2(diff) / KokkosBlas::nrm2(x_exact);

                std::cout << "=====================================" << std::endl;

                std::cout << "n_elements:        " << mesh.n_elements() << std::endl;
                std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
                std::cout << "n_nodes:           " << mesh.n_nodes() << std::endl;

                std::cout << "err : " << err << std::endl;

                ////////////////////////////////////////////////////////////////////////////
                // Compute op (x_exact) - rhs

                op.apply(x_exact, diff);

                Real sum_rhs = KokkosBlas::sum(rhs);

                Real norm_lapl_x = KokkosBlas::sum(diff);
                std::cout << "sum(Op(x)): " << norm_lapl_x << " == " << sum_rhs << std::endl;

                // Diff
                KokkosBlas::axpy(-1.0, rhs, diff);
                err = KokkosBlas::nrm2(diff) / KokkosBlas::nrm2(rhs);
                std::cout << "Op(x_exact) - rhs : " << err << std::endl;

                std::cout << "=====================================" << std::endl;

                if (write_output) {
                    VTUMeshWriter<SMesh> w;
                    VectorReal::HostMirror x_exact_host("x_exact_host", mesh.n_nodes());
                    Kokkos::deep_copy(x_exact_host, x_exact);
                    return w.write("x_exact.vtu", serial_mesh, x_exact_host);
                }

                return true;
            }

            void init() {
                op.init();
                auto id = std::make_shared<IdentityOperator>();
                id->init(mesh, bc_fun);
                op.set_identity(id);
                op.assemble_rhs(rhs, rhs_fun);
                bc.apply(rhs, bc_fun);

                if (write_output) {
                    convert_parallel_mesh_to_serial(serial_mesh, mesh);
                }
            }

            BC bc_fun;
            RHS rhs_fun;
            AnalyticalFun an_fun;

            PMesh &mesh;
            Op op;
            BoundaryConditions<PMesh> bc;
            VectorReal rhs;

            SMesh serial_mesh;
            bool write_output;
        };

        bool interpolate(PMesh &mesh, FEValues<PMesh> &values, const Integer &new_elem_offset, VectorReal &x) {
            VectorReal refined_x("refined_x", mesh.n_nodes());
            interpolate(mesh, values, new_elem_offset, x, refined_x);
            x = refined_x;
        }

        bool interpolate(PMesh &mesh,
                         FEValues<PMesh> &values,
                         const Integer &new_elem_offset,
                         const VectorReal &x,
                         VectorReal &refined_x) {
            Integer n_new = mesh.n_elements() - new_elem_offset;
            if (n_new == 0) return false;

            // auto parents = mesh.parent_map();

            VectorInt parents("parents", mesh.n_elements(), -1);

            Kokkos::parallel_for(
                new_elem_offset, MARS_LAMBDA(const Integer &i) {
                    auto c = mesh.get_children(i);
                    if (c.is_valid()) {
                        parents(c(0)) = i;
                        parents(c(1)) = i;
                    }
                });

            refined_x = VectorReal("refined_x", mesh.n_nodes());

            Kokkos::parallel_for(
                x.extent(0), MARS_LAMBDA(const Integer &i) { refined_x(i) = x(i); });

            {
                auto elems = mesh.get_view_elements();
                auto points = mesh.get_view_points();
                auto J_inv = values.J_inv();

                Kokkos::parallel_for(
                    n_new, MARS_LAMBDA(const Integer &iter) {
                        const Integer elem_id = new_elem_offset + iter;
                        assert(elem_id < parents.extent(0));

                        Integer parent = parents(elem_id);
                        assert(parent != -1);

                        Real J_inv_e[Dim * Dim];
                        Real tr[Dim], p[Dim], p_ref[Dim];
                        Integer idx[NFuns];
                        Real u_e[NFuns];

                        for (int k = 0; k < Dim * Dim; ++k) {
                            assert(parent < J_inv.extent(0));

                            J_inv_e[k] = J_inv(parent, k);
                            assert(J_inv_e[k] == J_inv_e[k]);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            assert(parent < elems.extent(0));

                            idx[k] = elems(parent, k);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            u_e[k] = x(idx[k]);
                        }

                        Integer n0 = elems(parent, 0);
                        assert(n0 < points.extent(0));

                        for (int d = 0; d < Dim; ++d) {
                            tr[d] = points(n0, d);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            assert(elem_id < elems.extent(0));

                            Integer c_n_id = elems(elem_id, k);

                            for (int d = 0; d < Dim; ++d) {
                                assert(c_n_id < points.extent(0));

                                p[d] = points(c_n_id, d) - tr[d];
                            }

                            Algebra<Dim>::mv_mult(J_inv_e, p, p_ref);

                            Real val = FESimplex<Dim>::fun(p_ref, u_e);
                            refined_x(c_n_id) = val;
                        }
                    });
            }
        }

        // bool refine(PMesh &mesh, FEValues<PMesh> &values, const Real tol, VectorReal &x, const bool write_output) {
        //     VectorReal error;
        //     GradientRecovery<PMesh> grad_rec;
        //     grad_rec.estimate(values, x, error);

        //     if (KokkosBlas::sum(error) < tol) {
        //         std::cout << "[Status] Reached desired tolerance" << std::endl;
        //         return false;
        //     }

        //     const Real alpha = 0.3;
        //     const Real max_error = KokkosBlas::nrminf(error);
        //     const Real low_bound_err = alpha * max_error;

        //     Integer old_n_elements = mesh.n_elements();

        //     VectorBool marked("marked", mesh.n_elements());

        //     Kokkos::parallel_for(
        //         mesh.n_elements(), MARS_LAMBDA(const Integer &i) { marked(i) = (error(i) >= low_bound_err) * i; });

        //     VectorInt marked_list = mark_active(marked);
        //     VectorInt::HostMirror marked_host("marked_host", marked_list.extent(0));
        //     Kokkos::deep_copy(marked_host, marked_list);

        //     std::vector<Integer> elems(marked_host.extent(0));

        //     for (Integer i = 0; i < marked_host.extent(0); ++i) {
        //         elems[i] = marked_host[i];
        //     }

        //     SMesh serial_mesh;
        //     convert_parallel_mesh_to_serial(serial_mesh, mesh);

        //     Bisection<SMesh> b(serial_mesh);
        //     b.refine(elems);

        //     serial_mesh.clean_up();
        //     convert_serial_mesh_to_parallel(mesh, serial_mesh);

        //     assert(serial_mesh.is_valid(true));

        //     x = VectorReal("x", mesh.n_nodes());

        //     // Kokkos::parallel_for(
        //     //     marked_list.extent(0), MARS_LAMBDA(const Integer &i) {
        //     //         assert(mesh.is_active(marked_list(i)));
        //     //         assert(marked_list(i) < mesh.n_elements());
        //     //     });

        //     // ParallelBisection<PMesh> bisection(&mesh);
        //     // bisection.refine(marked_list);

        //     // Integer n_new = mesh.n_elements() - old_n_elements;
        //     // if (n_new == 0) return false;

        //     // interpolate(mesh, values, old_n_elements, x);

        //     std::cout << "n_elements:        " << mesh.n_elements() << std::endl;
        //     std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
        //     std::cout << "n_nodes:           " << mesh.n_nodes() << std::endl;

        //     //
        //     //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //     // if (write_output) {
        //     //     VectorReal::HostMirror refined_x_host("refined_x_host", mesh.n_nodes());
        //     //     Kokkos::deep_copy(refined_x_host, x);

        //     //     SMesh serial_mesh;
        //     //     convert_parallel_mesh_to_serial(serial_mesh, mesh);

        //     //     VTUMeshWriter<SMesh> w;
        //     //     w.write("mesh_refined.vtu", serial_mesh, refined_x_host);
        //     // }

        //     return true;
        // }

        bool refine(PMesh &mesh, FEValues<PMesh> &values, const Real tol, VectorReal &x, const bool write_output) {
            ParallelBisection<PMesh> bisection(&mesh);

            VectorReal error;
            GradientRecovery<PMesh> grad_rec;
            grad_rec.estimate(values, x, error);

            if (KokkosBlas::sum(error) < tol) {
                std::cout << "[Status] Reached desired tolerance" << std::endl;
                return false;
            }

            const Real alpha = 0.3;
            const Real max_error = KokkosBlas::nrminf(error);
            const Real low_bound_err = alpha * max_error;

            Integer old_n_elements = mesh.n_elements();

            VectorBool marked("marked", mesh.n_elements());

            Kokkos::parallel_for(
                mesh.n_elements(), MARS_LAMBDA(const Integer &i) { marked(i) = (error(i) >= low_bound_err) * i; });

            VectorInt marked_list = mark_active(marked);

            Kokkos::parallel_for(
                marked_list.extent(0), MARS_LAMBDA(const Integer &i) {
                    assert(mesh.is_active(marked_list(i)));
                    assert(marked_list(i) < mesh.n_elements());
                });

            bisection.refine(marked_list);

            Integer n_new = mesh.n_elements() - old_n_elements;
            if (n_new == 0) return false;

            // x = VectorReal("x", mesh.n_nodes());

            // FIXME
            interpolate(mesh, values, old_n_elements, x);

            std::cout << "n_elements:        " << mesh.n_elements() << std::endl;
            std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
            std::cout << "n_nodes:           " << mesh.n_nodes() << std::endl;

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (write_output) {
                VectorReal::HostMirror refined_x_host("refined_x_host", mesh.n_nodes());
                Kokkos::deep_copy(refined_x_host, x);

                SMesh serial_mesh;
                convert_parallel_mesh_to_serial(serial_mesh, mesh);

                VTUMeshWriter<SMesh> w;
                w.write("mesh_refined.vtu", serial_mesh, refined_x_host);
            }

            return true;
        }

        bool solve(PMesh &mesh, const bool write_output) {
            const Integer n_nodes = mesh.n_nodes();
            VectorReal x("X", n_nodes);

            Real tol = 1e-6;
            Integer max_refinements = (Dim == 2) ? 10 : 5;
            Integer refs = 0;
            bool refined = false;
            // bool use_adaptive_refinement = true;
            bool use_adaptive_refinement = false;
            bool reset_x_to_zero = false;

            do {
                Problem problem(mesh);
                problem.write_output = write_output;
                problem.init();

                if (problem.solve(x)) {
                    problem.measure_actual_error(x);
                } else {
                    return false;
                }

                if (use_adaptive_refinement) {
                    refined = refine(mesh, problem.op.values(), tol, x, write_output);
                    mesh.clean_up();
                } else {
                    Integer elem_offset = mesh.n_elements();
                    ParallelBisection<PMesh>(&mesh).uniform_refine(1);
                    interpolate(mesh, problem.op.values(), elem_offset, x);
                    refined = true;
                }

                if (reset_x_to_zero) {
                    // reset x to 0
                    Kokkos::parallel_for(
                        x.extent(0), MARS_LAMBDA(const Integer &i) { x(i) = 0.0; });
                }

                if (!refined) {
                    std::cout << "Did not refine" << std::endl;
                    break;
                }

            } while (++refs <= max_refinements);
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
                convert_serial_mesh_to_parallel(mesh, smesh);

                ParallelBisection<PMesh>(&mesh).uniform_refine(3);
                write_output = false;
            }

            // ParallelBisection<PMesh> bisection(&mesh);
            // bisection.uniform_refine(5);

            solve(mesh, write_output);
        }
    };  // namespace mars
}  // namespace mars

#endif  // MARS_MODEL_TEST_HPP