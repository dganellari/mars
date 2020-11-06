#ifndef MARS_UMESH_ST_HEAT_EQUATION_HPP
#define MARS_UMESH_ST_HEAT_EQUATION_HPP

#include <memory>

#include "mars_base.hpp"
#include "mars_fe_simplex.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_identity_operator.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_simplex_quadrature.hpp"
#include "mars_simplex_spacetime.hpp"
#include "mars_tensor_spacetime.hpp"
#include "mars_umesh_operator.hpp"

namespace mars {

    template <class Mesh>
    class SpaceTimeMixed {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        // template <class Quadrature>
        // MARS_INLINE_FUNCTION static void one_thread_eval_diag(const Real *J_inv,
        //                                                       const Real &det_J,
        //                                                       const Quadrature &q,
        //                                                       Real *val) {
        //     Real g_ref[Dim], g_fe[Dim];
        //     const Real dx_P1 = det_J / NFuns;
        //     const Real dx_P0 = det_J;

        //     for (int d = 0; d < Dim; ++d) {
        //         g_ref[d] = -1;
        //     }

        //     Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

        //     Real u_t = g_fe[Dim - 1];
        //     val[0] = u_t * dx_P1 + Algebra<Dim - 1>::dot(g_fe, g_fe) * dx_P0;

        //     for (int d = 0; d < Dim; ++d) {
        //         g_ref[d] = 0;
        //     }

        //     for (int d = 0; d < Dim; ++d) {
        //         g_ref[d] = 1;

        //         Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

        //         u_t = g_fe[Dim - 1];
        //         val[d + 1] = u_t * dx_P1 + Algebra<Dim - 1>::dot(g_fe, g_fe) * dx_P0;

        //         g_ref[d] = 0;
        //     }

        //     for (int d = 0; d < Dim; ++d) {
        //         std::cout << val[d] << std::endl;
        //     }
        // }

        // template <class Quadrature>
        // MARS_INLINE_FUNCTION static void one_thread_eval_diag(const Real *J_inv,
        //                                                       const Real &det_J,
        //                                                       const Quadrature &q,
        //                                                       Real *val) {
        //     // assert(false);

        //     Real gi[Dim], g[Dim], g_x[Dim - 1], gj[Dim];
        //     Real pk[Dim];

        //     auto &q_points = q.points;
        //     auto &q_weights = q.weights;
        //     int n_qp = 3;  // q.n_points();

        //     for (int i = 0; i < NFuns; ++i) {
        //         val[i] = 0.0;
        //     }

        //     for (int k = 0; k < n_qp; ++k) {
        //         for (int d = 0; d < Dim; ++d) {
        //             pk[d] = q_points(k, d);
        //         }
        //         // Separate time and space dimensions

        //         assert(det_J > 0.0);
        //         const Real dx = det_J * q_weights(k) * 0.5;
        //         for (int i = 0; i < NFuns; i++) {
        //             // for each dof get the local number
        //             FESimplex<Dim>::grad(i, J_inv, gi);

        //             Real g_t = gi[Dim - 1];

        //             for (int d = 0; d < Dim - 1; ++d) {
        //                 g_x[d] = gi[d];
        //             }

        //             val[i] += (g_t * FESimplex<Dim>::fun(i, pk) + Algebra<Dim - 1>::dot(g_x, g_x)) * dx;
        //         }
        //     }
        // }

        template <class Quadrature>
        MARS_INLINE_FUNCTION static void one_thread_eval(const Real *J_inv,
                                                         const Real &det_J,
                                                         const Quadrature &,
                                                         const Real *u,
                                                         Real *val) {
            Real g_ref[Dim], g_fe[Dim], u_x[Dim - 1];
            const Real dx_P1 = det_J / NFuns;
            const Real dx_P0 = det_J;

            ///////////////////////////////////////////////////////////////////
            ////////////////// Gradient with local basis function /////////////
            // First basis function

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            // Transform
            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

            // The real u_t
            Real ut = g_fe[Dim - 1] * u[0];

            for (int d = 0; d < Dim - 1; ++d) {
                u_x[d] = g_fe[d] * u[0];
            }

            ///////////////////////////////////////////////////////////////
            // Reset to 0

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            ///////////////////////////////////////////////////////////////
            // The real \nabla_x u

            for (int i = 1; i < NFuns; ++i) {
                g_ref[i - 1] = 1;

                Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

                ut += g_fe[Dim - 1] * u[i];

                for (int d = 0; d < Dim - 1; ++d) {
                    u_x[d] += g_fe[d] * u[i];
                }

                g_ref[i - 1] = 0;
            }

            ////////////////////////////////////////////////////////////////////////
            // Integrate (u_t, v) [Set]
            for (int i = 0; i < NFuns; ++i) {
                val[i] = ut * dx_P1;
            }

            ///////////////////////////////////////////////////////////////////////
            // Integrate (u_x, v_x)  [Add]

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

            val[0] += Algebra<Dim - 1>::dot(u_x, g_fe) * dx_P0;

            //////////////////////////////////////////////////

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0.0;
            }

            for (int i = 1; i < NFuns; ++i) {
                g_ref[i - 1] = 1;

                Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

                val[i] += Algebra<Dim - 1>::dot(u_x, g_fe) * dx_P0;

                g_ref[i - 1] = 0;
            }
        }
    };

    template <class Mesh>
    class UMeshSTHeatEquation final : public UMeshOperator<Mesh> {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;
        using Super = mars::UMeshOperator<Mesh>;

        static constexpr int Dim = Mesh::Dim;
        // static constexpr int NFuns = Mesh::Dim + 1;
        static constexpr int NFuns = Elem::NNodes;

        UMeshSTHeatEquation(Mesh &mesh) : Super(mesh) {}

        SimplexQuadratureLinear<Dim> quad_;

        void init() override {
            this->values().init();

            quad_ = SimplexQuadratureLinear<Dim>::make();
            auto prec = std::make_shared<JacobiPreconditioner>();
            prec->quad_ = quad_;

            this->set_precontitioner(prec);
            this->preconditioner()->init(this->values());
        }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {
            // For Kokkos-Cuda
            auto values = this->values();
            auto det_J = values.det_J();
            auto J_inv = values.J_inv();
            auto mesh = values.mesh();

            ViewMatrixType<Integer> elems = mesh.get_view_elements();
            auto active = mesh.get_view_active();

            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            auto quad = quad_;

            SpaceTime<Elem> st;

            Kokkos::parallel_for(
                "UMeshSTHeatEquation::apply", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real Au[NFuns];
                    Integer idx[NFuns];
                    Real J_inv_e[Dim * Dim];

                    if (!active(i)) return;  // ACTIVE

                    for (int k = 0; k < (Dim * Dim); ++k) {
                        J_inv_e[k] = J_inv(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        u[k] = x(idx[k]);
                    }

                    // SpaceTimeMixed<Mesh>::one_thread_eval(J_inv_e, det_J(i), quad, u, Au);
                    st.one_thread_eval(J_inv_e, det_J(i), u, Au);

                    for (Integer k = 0; k < NFuns; ++k) {
                        Kokkos::atomic_add(&op_x(idx[k]), Au[k]);
                    }
                });

            if (this->identity()) {
                this->identity()->apply(x, op_x);
            }
        }

        class JacobiPreconditioner final : public UMeshJacobiPreconditioner<Mesh> {
        public:
            void init(FEValues<Mesh> &values) override {
                auto mesh = values.mesh();
                ViewMatrixType<Integer> elems = mesh.get_view_elements();
                auto active = mesh.get_view_active();

                auto det_J = values.det_J();
                auto J_inv = values.J_inv();

                ViewVectorType<Real> inv_diag("inv_diag", mesh.n_nodes());

                auto quad = quad_;

                SpaceTime<Elem> st;

                Kokkos::parallel_for(
                    "JacobiPreconditioner::init", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                        Integer idx[NFuns];
                        Real val[NFuns];
                        Real J_inv_e[Dim * Dim];

                        if (!active(i)) return;  // ACTIVE

                        const Real det_J_e = det_J(i);

                        assert(det_J_e > 0.0);

                        for (Integer k = 0; k < (Dim * Dim); ++k) {
                            J_inv_e[k] = J_inv(i, k);
                        }

                        for (Integer k = 0; k < NFuns; ++k) {
                            idx[k] = elems(i, k);
                        }
                        // SpaceTimeMixed<Mesh>::one_thread_eval_diag(J_inv_e, det_J_e, quad, val);
                        st.one_thread_eval_diag(J_inv_e, det_J_e, val);
                        //
                        for (Integer k = 0; k < NFuns; ++k) {
                            assert(val[k] != 0.0);
                            Real inv_val = val[k];
                            assert(inv_val == inv_val);
                            // inv_diag(idx[k]) += inv_val;

                            Kokkos::atomic_add(&inv_diag(idx[k]), inv_val);
                        }
                    });

                Kokkos::parallel_for(

                    mesh.n_nodes(), MARS_LAMBDA(const Integer d) {
                        inv_diag(d) = 1. / std::sqrt(inv_diag(d));
                        // std::cout << inv_diag(d) << std::endl;
                    });

                this->inv_diag_ = inv_diag;
            }

            SimplexQuadratureLinear<Dim> quad_;
        };
    };

}  // namespace mars

#endif  // MARS_UMESH_ST_HEAT_EQUATION_HPP
