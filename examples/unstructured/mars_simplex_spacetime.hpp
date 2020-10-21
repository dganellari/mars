#ifndef MARS_SIMPLEX_SPACETIME_HPP
#define MARS_SIMPLEX_SPACETIME_HPP

#include "mars_base.hpp"
#include "mars_fe_simplex.hpp"
#include "mars_globals.hpp"
#include "mars_simplex_quadrature.hpp"

#include "mars_quad4.hpp"
#include "mars_spacetime.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim, class Implementation>
    class SpaceTime<Simplex<Dim, ManifoldDim, Implementation> > {
    public:
        using Elem = mars::Simplex<Dim, ManifoldDim, Implementation>;
        static constexpr int NFuns = Elem::NNodes;

        MARS_INLINE_FUNCTION void one_thread_eval_diag(const Real *J_inv, const Real &det_J, Real *val) const {
            // assert(false);

            Real gi[Dim], g[Dim], g_x[Dim - 1], gj[Dim];
            Real pk[Dim];

            auto &q_points = q.points;
            auto &q_weights = q.weights;
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
                const Real dx = det_J * q_weights(k) * 0.5;
                for (int i = 0; i < NFuns; i++) {
                    // for each dof get the local number
                    FESimplex<Dim>::grad(i, J_inv, gi);

                    Real g_t = gi[Dim - 1];

                    for (int d = 0; d < Dim - 1; ++d) {
                        g_x[d] = gi[d];
                    }

                    val[i] += (g_t * FESimplex<Dim>::fun(i, pk) + Algebra<Dim - 1>::dot(g_x, g_x)) * dx;
                }
            }
        }

        MARS_INLINE_FUNCTION void one_thread_eval(const Real *J_inv,
                                                  const Real &det_J,
                                                  const Real *u,
                                                  Real *val) const {
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

        SpaceTime() : q(SimplexQuadratureLinear<Dim>::make()) {}
        SimplexQuadratureLinear<Dim> q;
    };

}  // namespace mars

#endif MARS_SIMPLEX_SPACETIME_HPP