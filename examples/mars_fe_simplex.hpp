#ifndef MARS_FE_SIMPLEX_HPP
#define MARS_FE_SIMPLEX_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {

    template <int Dim>
    class Algebra {
    public:
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
    };

    template <int Dim>
    class FESimplex {
    public:
        MARS_INLINE_FUNCTION static Real fun(const int i, const Real *p) {
            Real ret = 0.0;

            if (i == 0) {
                ret = 1.0;
                for (int d = 0; d < Dim; ++d) {
                    ret -= p[d];
                }
            } else {
                ret = p[i - 1];
            }

            return ret;
        }

        MARS_INLINE_FUNCTION static void grad(const Real *J_inv, const Real *u, Real *g) {
            Real g_ref[Dim];

            ///////////////////////////////////////////////////////////////////
            ////////////////// Gradient with local basis function /////////////

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1 * u[0];
            }

            for (int i = 1; i < Dim + 1; ++i) {
                g_ref[i - 1] += u[i];
            }

            ///////////////////////////////////////////////////////////////////
            ////////////////// Transform gradient to physical coordinates //////

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g);
        }
    };
}  // namespace mars

#endif  // MARS_FE_SIMPLEX_HPP