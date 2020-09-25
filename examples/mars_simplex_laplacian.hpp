#ifndef MARS_SIMPLEX_LAPLACIAN_HPP
#define MARS_SIMPLEX_LAPLACIAN_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {

    template <class Mesh>
    class SimplexLaplacian {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        MARS_INLINE_FUNCTION static void m_t_v_mult(const Real *A, const Real *x, Real *y) {
            for (int d1 = 0; d1 < Dim; ++d1) {
                y[d1] = 0;

                for (int d2 = 0; d2 < Dim; ++d2) {
                    y[d1] += A[d1 + d2 * Dim] * x[d2];
                }
            }
        }

        MARS_INLINE_FUNCTION static Real dot(const Real *l, const Real *r) {
            Real ret = 0.0;
            for (Integer i = 0; i < Dim; ++i) {
                ret += l[i] * r[i];
            }

            return ret;
        }

        MARS_INLINE_FUNCTION static void one_thread_eval_diag(const Real *J_inv, const Real &det_J, Real *val) {
            Real g_ref[Dim], g[Dim];

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            assert(det_J > 0.0);
            assert(has_non_zero<Dim>(J_inv));

            m_t_v_mult(J_inv, g_ref, g);
            val[0] = dot(g, g) * det_J;

            assert(val[0] > 0.0);

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 1;
                m_t_v_mult(J_inv, g_ref, g);

                val[d + 1] = dot(g, g) * det_J;

                assert(val[d + 1] > 0.0);

                g_ref[d] = 0;
            }
        }

        MARS_INLINE_FUNCTION static void one_thread_eval(const Real *J_inv,
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

            ///////////////////////////////////////////////////////////////////
            ////////////////// evaluate bilinear form ////////////////////////

            ///////////// Evaluate for Phi_0(x) = 1 - x_0 - ... x_n

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            m_t_v_mult(J_inv, g_ref, g_fe);

            val[0] = dot(g, g_fe) * det_J;

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            ///////////// Evaluate for Phi_i(x) = x_{i-1}
            for (int i = 1; i < NFuns; ++i) {
                int d = i - 1;

                g_ref[d] = 1.0;

                m_t_v_mult(J_inv, g_ref, g_fe);

                val[i] = dot(g, g_fe) * det_J;

                g_ref[d] = 0.0;
            }
        }
    };

}  // namespace mars

#endif  // MARS_SIMPLEX_LAPLACIAN_HPP
