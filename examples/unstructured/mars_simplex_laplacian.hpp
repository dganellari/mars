#ifndef MARS_SIMPLEX_LAPLACIAN_HPP
#define MARS_SIMPLEX_LAPLACIAN_HPP

#include "mars_base.hpp"
#include "mars_fe_simplex.hpp"
#include "mars_globals.hpp"

#include "mars_laplacian.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim, class Implementation>
    class Laplacian<Simplex<Dim, ManifoldDim, Implementation> > {
    public:
        using Elem = mars::Simplex<Dim, ManifoldDim, Implementation>;
        static constexpr int NFuns = Elem::NNodes;

        MARS_INLINE_FUNCTION static void one_thread_eval_diag(const Real *J_inv, const Real &det_J, Real *val) {
            Real g_ref[Dim], g[Dim];

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            assert(det_J > 0.0);
            assert(has_non_zero<Dim>(J_inv));

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g);
            val[0] = Algebra<Dim>::dot(g, g) * det_J;

            assert(val[0] > 0.0);

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 1;
                Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g);

                val[d + 1] = Algebra<Dim>::dot(g, g) * det_J;

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

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g);

            ///////////////////////////////////////////////////////////////////
            ////////////////// evaluate bilinear form ////////////////////////

            ///////////// Evaluate for Phi_0(x) = 1 - x_0 - ... x_n

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

            val[0] = Algebra<Dim>::dot(g, g_fe) * det_J;

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            ///////////// Evaluate for Phi_i(x) = x_{i-1}
            for (int i = 1; i < NFuns; ++i) {
                int d = i - 1;

                g_ref[d] = 1.0;

                Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g_fe);

                val[i] = Algebra<Dim>::dot(g, g_fe) * det_J;

                g_ref[d] = 0.0;
            }
        }
    };

}  // namespace mars

#endif  // MARS_SIMPLEX_LAPLACIAN_HPP
