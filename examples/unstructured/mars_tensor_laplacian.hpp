#ifndef MARS_TENSOR_LAPLACIAN_HPP
#define MARS_TENSOR_LAPLACIAN_HPP

#include "mars_base.hpp"
#include "mars_fe_simplex.hpp"
#include "mars_globals.hpp"

#include "mars_laplacian.hpp"
#include "mars_quad4.hpp"

namespace mars {

    template <Integer Type, class Implementation>
    class Laplacian<NonSimplex<Type, Implementation> > {
    public:
        using Elem = mars::NonSimplex<Type, Implementation>;
        static constexpr int NFuns = Elem::NNodes;
        static constexpr int Dim = Elem::Dim;

        MARS_INLINE_FUNCTION void one_thread_eval_diag(const Real *J_inv, const Real &det_J, Real *val) const {
            // assert(false);
            for (int i = 0; i < NFuns; ++i) {
                val[i] = 1.0;
            }
        }

        MARS_INLINE_FUNCTION void one_thread_eval(const Real *J_inv,
                                                  const Real &det_J,
                                                  const Real *u,
                                                  Real *val) const {
            Real gi[Dim], g[Dim];
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

                assert(det_J > 0.0);
                const Real dx = det_J * q_weights(k);

                for (int i = 0; i < NFuns; i++) {
                    // for each dof get the local number
                    FEQuad4<Real>::Grad::affine_f(i, J_inv, pk, gi);

                    val[i] += Algebra<Dim>::dot(g, gi) * dx;
                }
            }
        }

        Laplacian() : q(FEQuad4<Real, 2>::Quadrature::make()) {}
        FEQuad4<Real, 2>::Quadrature q;
    };

    template <Integer Type, class Implementation>
    class Jacobian<NonSimplex<Type, Implementation> > {
    public:
        using Elem = mars::NonSimplex<Type, Implementation>;
        static const int Dim = Elem::Dim;

        template <class View>
        MARS_INLINE_FUNCTION void static compute(const Integer *idx,
                                                 const View &points,
                                                 Real *J,
                                                 Real *J_inv,
                                                 Real &det_J) {
            static const int ref_idx[4] = {0, 1, 3, 4};
            static const int NCorners = Dim + 1;

            Real p0[Dim], pk[Dim];

            for (int d = 0; d < Dim; ++d) {
                p0[d] = points(idx[0], d);
            }

            for (int k = 1; k < NCorners; ++k) {
                const int km1 = k - 1;
                const int r = idx[ref_idx[k]];

                for (int d = 0; d < Dim; ++d) {
                    pk[d] = points(r, d);
                }

                for (int d = 0; d < Dim; ++d) {
                    J[d * Dim + km1] = pk[d] - p0[d];
                }
            }

            Invert<Dim>::apply(J, J_inv, det_J);

            assert(has_non_zero<Dim>(J_inv));
        }
    };

}  // namespace mars

#endif  // MARS_TENSOR_LAPLACIAN_HPP
