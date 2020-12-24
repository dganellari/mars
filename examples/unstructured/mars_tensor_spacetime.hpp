#ifndef MARS_TENSOR_SPACETIME_HPP
#define MARS_TENSOR_SPACETIME_HPP

#include "mars_base.hpp"
#include "mars_fe_simplex.hpp"
#include "mars_globals.hpp"

#include "mars_quad4.hpp"
#include "mars_spacetime.hpp"

namespace mars {

    template <Integer Type, class Implementation>
    class SpaceTime<NonSimplex<Type, Implementation>> {
    public:
        using Elem = mars::NonSimplex<Type, Implementation>;
        static constexpr int NFuns = Elem::NNodes;
        static constexpr int Dim = Elem::Dim;

        MARS_INLINE_FUNCTION void one_thread_eval_diag(const Real *J_inv, const Real &det_J, Real *val) const {
            // assert(false);
            // for (int i = 0; i < NFuns; ++i) {
            //     val[i] = 1.0;
            // }

            Real gi[Dim], g[Dim], g_x[Dim - 1];
            Real pk[Dim];

            auto &q_points = q.q_p;
            auto &q_weights = q.q_w;
            int n_qp = q.n_points();

            for (int i = 0; i < NFuns; ++i) {
                val[i] = 0.0;
            }

            for (int k = 0; k < n_qp; ++k) {
                for (int d = 0; d < Dim; ++d) {
                    pk[d] = q_points(k, d);
                }

                // Separate time and space dimensions

                assert(det_J > 0.0);
                const Real dx = det_J * q_weights(k);

                for (int i = 0; i < NFuns; i++) {
                    // for each dof get the local number
                    FEQuad4<Real>::Grad::affine_f(i, J_inv, pk, gi);
                    Real g_t = gi[Dim - 1];

                    for (int d = 0; d < Dim - 1; ++d) {
                        g_x[d] = gi[d];
                    }
                    val[i] += (g_t * FEQuad4<Real>::Fun::f(i, pk) + Algebra<Dim - 1>::dot(g_x, g_x)) * dx;
                }
            }
        }

        MARS_INLINE_FUNCTION void one_thread_eval(const Real *J_inv,
                                                  const Real &det_J,
                                                  const Real *u,
                                                  Real *val) const {
            Real gi[Dim], g[Dim], g_x[Dim - 1];
            Real pk[Dim];

            auto &q_points = q.q_p;
            auto &q_weights = q.q_w;
            int n_qp = q.n_points();

            for (int i = 0; i < NFuns; ++i) {
                val[i] = 0.0;
            }

            for (int k = 0; k < n_qp; ++k) {
                for (int d = 0; d < Dim; ++d) {
                    pk[d] = q_points(k, d);
                }

                ////////////////////////
                // Compute physical gradient of solution once per quadrature point
                FEQuad4<Real>::Grad::ref(pk, u, gi);
                Algebra<Dim>::m_t_v_mult(J_inv, gi, g);
                ////////////////////////
                // Separate time and space dimensions
                Real g_t = g[Dim - 1];

                for (int d = 0; d < Dim - 1; ++d) {
                    g_x[d] = g[d];
                }

                assert(det_J > 0.0);
                const Real dx = det_J * q_weights(k);

                for (int i = 0; i < NFuns; i++) {
                    // for each dof get the local number
                    FEQuad4<Real>::Grad::affine_f(i, J_inv, pk, gi);
                    val[i] += (g_t * FEQuad4<Real>::Fun::f(i, pk) / NFuns + Algebra<Dim - 1>::dot(g, gi)) * dx;
                }
            }
        }

        SpaceTime() : q(FEQuad4<Real, 2>::Quadrature::make()) {}
        FEQuad4<Real, 2>::Quadrature q;
    };

    using SpaceTimeQuad4 = SpaceTime<NonSimplex<4, KokkosImplementation>>;
    // using JacobianQuad4 = Jacobian_st<NonSimplex<4, KokkosImplementation>>;
}  // namespace mars

#endif  // MARS_TENSOR_SPACETIME_HPP
