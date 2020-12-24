#ifndef MARS_FE_SIMPLEX_HPP
#define MARS_FE_SIMPLEX_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

#include "mars_invert.hpp"
#include "mars_jacobian.hpp"

namespace mars {

    template <int Dim>
    MARS_INLINE_FUNCTION bool has_non_zero(const Real *J_inv) {
        bool ret = false;
        for (int k = 0; k < (Dim * Dim); ++k) {
            if (J_inv[k] * J_inv[k] > 0.0) {
                ret = true;
                break;
            }
        }

        return ret;
    }

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

        MARS_INLINE_FUNCTION static void mv_mult(const Real *A, const Real *x, Real *y) {
            for (int d1 = 0; d1 < Dim; ++d1) {
                y[d1] = 0;

                for (int d2 = 0; d2 < Dim; ++d2) {
                    y[d1] += A[d1 * Dim + d2] * x[d2];
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

    template <Integer Dim, Integer ManifoldDim, class Implementation>
    class Jacobian<Simplex<Dim, ManifoldDim, Implementation> > {
    public:
        using Elem = mars::Simplex<Dim, ManifoldDim, Implementation>;
        static const int NFuns = Elem::NNodes;

        template <class View>
        MARS_INLINE_FUNCTION void static compute(const Integer *idx,
                                                 const View &points,
                                                 Real *J,
                                                 Real *J_inv,
                                                 Real &det_J) {
            Real p0[Dim], pk[Dim];

            for (int d = 0; d < Dim; ++d) {
                p0[d] = points(idx[0], d);
            }

            for (int k = 1; k < NFuns; ++k) {
                const int km1 = k - 1;

                for (int d = 0; d < Dim; ++d) {
                    pk[d] = points(idx[k], d);
                }

                for (int d = 0; d < Dim; ++d) {
                    J[d * Dim + km1] = pk[d] - p0[d];
                }
            }

            Invert<Dim>::apply(J, J_inv, det_J);

            assert(has_non_zero<Dim>(J_inv));
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

        MARS_INLINE_FUNCTION static bool inside(const Real *p, const Real &tol) {
            Real sum_p = 0.0;
            for (int i = 0; i < Dim; ++i) {
                if (p[i] <= -tol) {
                    return false;
                }

                sum_p += p[i];
            }

            return sum_p <= 1.0 + tol;
        }

        MARS_INLINE_FUNCTION static Real fun(const Real *p, const Real *u) {
            assert(inside(p, 1e-8));

            Real ret = 1;

            for (int i = 0; i < Dim; ++i) {
                ret -= p[i];
            }

            ret *= u[0];

            for (int i = 0; i < Dim; ++i) {
                ret += p[i] * u[i + 1];
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

        MARS_INLINE_FUNCTION static void grad(const int i, const Real *J_inv, Real *g) {
            Real g_ref[Dim];

            ///////////////////////////////////////////////////////////////////
            ////////////////// Gradient with local basis function /////////////

            if (i == 0) {
                for (int d = 0; d < Dim; ++d) {
                    g_ref[d] = -1;
                }
            } else {
                for (int d = 0; d < Dim; ++d) {
                    g_ref[d] = 0.0;
                }

                g_ref[i - 1] = 1.0;
            }

            ///////////////////////////////////////////////////////////////////
            ////////////////// Transform gradient to physical coordinates //////

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g);
        }
    };
}  // namespace mars

#endif  // MARS_FE_SIMPLEX_HPP
