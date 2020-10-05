#ifndef MARS_INVERT_HPP
#define MARS_INVERT_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {
    template <int Dim>
    class Invert {};

    template <>
    class Invert<1> {
    public:
        MARS_INLINE_FUNCTION static bool apply(const Real *m, Real *m_inv, Real &det_m) {
            m_inv[0] = 1. / m[0];
            det_m = m[0];
            return det_m != 0.0;
        }
    };

    template <>
    class Invert<2> {
    public:
        MARS_INLINE_FUNCTION static Real det(const Real *m) { return m[0] * m[3] - m[2] * m[1]; }

        MARS_INLINE_FUNCTION static bool invert(const Real *mat, Real *mat_inv, const Real &det) {
            mat_inv[0] = mat[3] / det;
            mat_inv[1] = -mat[1] / det;
            mat_inv[2] = -mat[2] / det;
            mat_inv[3] = mat[0] / det;
            return true;
        }

        MARS_INLINE_FUNCTION static bool apply(const Real *m, Real *m_inv, Real &det_m) {
            det_m = det(m);
            return invert(m, m_inv, det_m);
        }
    };

    template <>
    class Invert<3> {
    public:
        MARS_INLINE_FUNCTION static Real det(const Real *m) {
            return m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[0] * m[5] * m[7] -
                   m[1] * m[3] * m[8] - m[2] * m[4] * m[6];
        }

        MARS_INLINE_FUNCTION static bool invert(const Real *mat, Real *mat_inv, const Real &det) {
            assert(det != 0.);

            if (det == 0.) {
                return false;
            }

            mat_inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / det;
            mat_inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / det;
            mat_inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / det;
            mat_inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / det;
            mat_inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / det;
            mat_inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / det;
            mat_inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / det;
            mat_inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / det;
            mat_inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / det;
            return true;
        }

        MARS_INLINE_FUNCTION static bool apply(const Real *m, Real *m_inv, Real &det_m) {
            det_m = det(m);
            return invert(m, m_inv, det_m);
        }
    };

    template <>
    class Invert<4> {
    public:
        MARS_INLINE_FUNCTION static bool apply(const Real *m, Real *m_inv, Real &det_m) {
            m_inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] +
                       m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

            m_inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] -
                       m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

            m_inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] +
                       m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

            m_inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] -
                        m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

            m_inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] -
                       m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

            m_inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] +
                       m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

            m_inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] -
                       m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

            m_inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] +
                        m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

            m_inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] +
                       m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

            m_inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] -
                       m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

            m_inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] +
                        m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

            m_inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] -
                        m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

            m_inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] -
                       m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

            m_inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] +
                       m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

            m_inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] -
                        m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

            m_inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] +
                        m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

            det_m = m[0] * m_inv[0] + m[1] * m_inv[4] + m[2] * m_inv[8] + m[3] * m_inv[12];

            if (det_m == 0) return false;

            const Real inv_det = 1.0 / det_m;

            for (int i = 0; i < 16; i++) {
                m_inv[i] *= inv_det;
            }

            return true;
        }
    };

}  // namespace mars

#endif  // MARS_INVERT_HPP
