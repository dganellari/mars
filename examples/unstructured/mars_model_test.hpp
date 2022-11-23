#ifndef MARS_MODEL_TEST_HPP
#define MARS_MODEL_TEST_HPP

#include "mars_base.hpp"
#include "mars_err.hpp"

#ifdef MARS_ENABLE_KOKKOS_KERNELS

#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <KokkosBlas1_sum.hpp>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>

#ifdef MARS_ENABLE_CXXOPTS
#include <cxxopts.hpp>
#endif //MARS_ENABLE_CXXOPTS
// #include <Kokkos_sort.hpp>

#ifdef MARS_ENABLE_SERIAL_BACKEND
#include "mars_globals.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_mesh.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_serial_utils.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"

#ifdef MARS_ENABLE_VTK
#include "vtu_writer.hpp"

#endif
#endif  // MARS_ENABLE_SERIAL_BACKEND

#include "mars_boundary_conditions.hpp"
#include "mars_copy_operator.hpp"
#include "mars_fe_values.hpp"
#include "mars_gradient_recovery.hpp"
#include "mars_identity_operator.hpp"
#include "mars_interpolate.hpp"
#include "mars_mesh_kokkos.hpp"

using namespace std::chrono;

namespace mars {

#ifdef MARS_ENABLE_SERIAL_BACKEND

    template <class Mesh>
    struct SerialMeshType {};

    template <int Dim, int ManifoldDim>
    struct SerialMeshType<
        mars::Mesh<Dim, ManifoldDim, KokkosImplementation, Simplex<Dim, ManifoldDim, KokkosImplementation>>> {
        using Type = mars::SimplicialMesh<Dim, ManifoldDim>;
    };

    template <int Dim, int ManifoldDim, int ElemType>
    struct SerialMeshType<
        mars::Mesh<Dim, ManifoldDim, KokkosImplementation, NonSimplex<ElemType, KokkosImplementation>>> {
        using Type = mars::Mesh<Dim, ManifoldDim, DefaultImplementation, NonSimplex<ElemType>>;
    };

    template <>
    struct SerialMeshType<mars::ParallelQuad4Mesh> {
        using Type = mars::Mesh<2, 2, DefaultImplementation, NonSimplex<4>>;
    };

    template <>
    struct SerialMeshType<mars::ParallelMesh2> {
        using Type = mars::SimplicialMesh<2, 2>;
    };

    template <>
    struct SerialMeshType<mars::ParallelMesh3> {
        using Type = mars::SimplicialMesh<3, 3>;
    };

    template <>
    struct SerialMeshType<mars::ParallelMesh4> {
        using Type = mars::SimplicialMesh<4, 4>;
    };

#endif  // MARS_ENABLE_SERIAL_BACKEND

    template <class PMesh, class Op, class BC, class RHS, class AnalyticalFun>
    class ModelTest {
    public:
        static const int Dim = PMesh::Dim;
        static const int ManifoldDim = PMesh::ManifoldDim;
        static const int NFuns = ManifoldDim + 1;
#ifdef MARS_ENABLE_SERIAL_BACKEND
        using SMesh = typename SerialMeshType<PMesh>::Type;
#endif  // MARS_ENABLE_SERIAL_BACKEND
        using VectorReal = mars::ViewVectorType<Real>;
        using VectorInt = mars::ViewVectorType<Integer>;
        using VectorBool = mars::ViewVectorType<bool>;

        class Problem {
        public:
            Problem(PMesh &mesh) : mesh(mesh), op(mesh), bc(mesh), rhs("rhs", mesh.n_nodes()), write_output(true) {}

            bool solve(VectorReal &x) {
                x = VectorReal("x", mesh.n_nodes());
                bc.apply(x, bc_fun);

                // auto prec_ptr = op.preconditioner();
                auto prec_ptr = std::make_shared<CopyOperator>();

                Integer num_iter = 0;

                auto start = steady_clock::now();

                bool ok = bcg_stab(op, *prec_ptr, rhs, rhs.extent(0), x, num_iter);

                auto end = steady_clock::now();
                auto diff = end - start;

                std::cout << "bcg_stab/time: " << std::chrono::duration<double>(diff).count() << " s" << std::endl;

                // FILE *file_time = fopen("timing.csv", "a");
                // fprintf(file_time, "%ld, %.4fs\n", mesh.n_nodes(), std::chrono::duration<double>(diff).count());
                // fclose(file_time);

                if (ok) {
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

#ifdef MARS_ENABLE_SERIAL_BACKEND
#ifdef MARS_ENABLE_VTK

                SMesh serial_mesh;
                convert_parallel_mesh_to_serial(serial_mesh, mesh);

                VTUMeshWriter<SMesh> w;

                if (!w.write("solution.vtu", serial_mesh, x_host)) {
                    return false;
                }

                if (!w.write("rhs.vtu", serial_mesh, rhs_host)) {
                    return false;
                }
#endif  // MARS_ENABLE_VTK
#endif  // MARS_ENABLE_SERIAL_BACKEND
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

#ifdef MARS_ENABLE_SERIAL_BACKEND
#ifdef MARS_ENABLE_VTK
                if (write_output) {
                    VTUMeshWriter<SMesh> w;
                    VectorReal::HostMirror x_exact_host("x_exact_host", mesh.n_nodes());
                    Kokkos::deep_copy(x_exact_host, x_exact);
                    return w.write("x_exact.vtu", serial_mesh, x_exact_host);
                }
#endif  // MARS_ENABLE_VTK
#endif  // MARS_ENABLE_SERIAL_BACKEND

                return true;
            }

            void init() {
                op.init();
                auto id = std::make_shared<IdentityOperator>();
                id->init(mesh, bc_fun);
                op.set_identity(id);
                op.assemble_rhs(rhs, rhs_fun);
                bc.apply(rhs, bc_fun);

#ifdef MARS_ENABLE_SERIAL_BACKEND
                if (write_output) {
                    convert_parallel_mesh_to_serial(serial_mesh, mesh);
                }
#endif  // MARS_ENABLE_SERIAL_BACKEND
            }

            BC bc_fun;
            RHS rhs_fun;
            AnalyticalFun an_fun;

            PMesh &mesh;
            Op op;
            BoundaryConditions<PMesh> bc;
            VectorReal rhs;

#ifdef MARS_ENABLE_SERIAL_BACKEND
            SMesh serial_mesh;
#endif  // MARS_ENABLE_SERIAL_BACKEND

            bool write_output;
        };

        bool interpolate(PMesh &mesh, FEValues<PMesh> &values, const Integer &new_elem_offset, VectorReal &x) {
            VectorReal refined_x("refined_x", mesh.n_nodes());
            interpolate(mesh, values, new_elem_offset, x, refined_x);
            x = refined_x;
            return true;
        }

        bool interpolate(PMesh &mesh,
                         FEValues<PMesh> &values,
                         const Integer &new_elem_offset,
                         const VectorReal &x,
                         VectorReal &refined_x) {
            Integer n_elements = mesh.n_elements();
            Integer n_new = n_elements - new_elem_offset;
            if (n_new == 0) return false;

            // auto parents = mesh.parent_map();

            VectorInt parents("parents", n_elements);
            Kokkos::deep_copy(parents, -1);

            Kokkos::parallel_for(
                new_elem_offset, MARS_LAMBDA(const Integer &i) {
                    auto c = mesh.get_children(i);
                    if (c.is_valid()) {
                        assert(c(0) < n_elements);
                        assert(c(1) < n_elements);
                        assert(c(0) > 0);
                        assert(c(1) > 0);

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
#ifdef MARS_ENABLE_SERIAL_BACKEND
#ifdef MARS_ENABLE_VTK
            if (write_output) {
                VectorReal::HostMirror refined_x_host("refined_x_host", mesh.n_nodes());
                Kokkos::deep_copy(refined_x_host, x);

                SMesh serial_mesh;
                convert_parallel_mesh_to_serial(serial_mesh, mesh);

                VTUMeshWriter<SMesh> w;
                w.write("mesh_refined.vtu", serial_mesh, refined_x_host);
            }
#endif  // MARS_ENABLE_VTK
#endif  // MARS_ENABLE_SERIAL_BACKEND

            return true;
        }

        bool solve(PMesh &mesh) {
            const Integer n_nodes = mesh.n_nodes();
            VectorReal x("X", n_nodes);

            // Integer refs{0};

            // do {
            Problem problem(mesh);
            problem.write_output = write_output;
            problem.init();

            if (problem.solve(x)) {
                problem.measure_actual_error(x);
            } else {
                return false;
            }

            // if (max_refinements == 0 || !use_adaptive_refinement) return true;
            /*


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
                            } */

            // } while (++refs <= max_refinements);
            return true;
        }


        void run(
#ifdef MARS_ENABLE_CXXOPTS
            cxxopts::ParseResult &args
#endif
            )
        {


#ifdef MARS_ENABLE_CXXOPTS
            Integer ns[4] = {
                args["nx"].as<Integer>(), args["ny"].as<Integer>(), args["nz"].as<Integer>(), args["nt"].as<Integer>()};

            use_adaptive_refinement = args["adaptive"].as<bool>();
            write_output = args["output"].as<bool>();
            max_refinements = args["refine_level"].as<Integer>();
#else
            Integer ns[4] = {10, 10, 10, 10};
#endif

            PMesh mesh;

            for (int d = Dim; d < 3; ++d) {
                ns[d] = 0;
            }

            if (Dim <= 3) {
                generate_cube(mesh, ns[0], ns[1], ns[2]);
            } else {
                assert(false);
                // SMesh smesh;
                // read_mesh("../data/cube4d_24.MFEM", smesh);
                // convert_serial_mesh_to_parallel(mesh, smesh);

                // ParallelBisection<PMesh>(&mesh).uniform_refine(3);
                /* ParallelBisection<PMesh>(&mesh).uniform_refine(4); */
                /* write_output = false; */
            }

            // ParallelBisection<PMesh> bisection(&mesh);
            // bisection.uniform_refine(5);

            solve(mesh);
        }

    private:
        Real tol{1e-6};
        Integer max_refinements{0};

        bool refined{false};
        bool use_adaptive_refinement{false};
        bool reset_x_to_zero{false};
        bool write_output{true};

    };  // namespace mars
}  // namespace mars

#endif
#endif  // MARS_MODEL_TEST_HPP
