#ifndef MARS_LAPLACIAN_EXAMPLES
#define MARS_LAPLACIAN_EXAMPLES

#include <cmath>

#include "mars_base.hpp"
#include "mars_config.hpp"

namespace mars {

    /* Examples for -laplacian u = f */

    // Smooth solution
    MARS_INLINE_FUNCTION Real ex1_exact(const Real *points, Real &result) {
        return points[0] * points[1] * (1 - points[0]) * (1 - points[1]);
    }

    MARS_INLINE_FUNCTION Real ex1_laplacian(const Real *points, Real &result) {
        return -2.0 * points[0] * (1 - points[0]) - 2.0 * points[1] * (1 - points[1]);
    }

    /* ---------------------------------------------------------------------------- */

    // Peak solution
    MARS_INLINE_FUNCTION Real ex2_exact(const Real *points, Real &result) {
        return std::exp(-1000 * ((points[0] - 0.5) * (points[0] - 0.5) + (points[1] - 0.5) * (points[1] - 0.5)));
    }

    MARS_INLINE_FUNCTION Real ex2_laplacian(const Real *points, Real &result) {
        return -((std::exp(-1000 * (points[0] - 1 / 2) * (points[0] - 1 / 2) -
                           1000 * (points[1] - 1 / 2) * (points[1] - 1 / 2)) *
                      (2000 * points[0] - 1000) * (2000 * points[0] - 1000) -
                  2000 * std::exp(-1000 * (points[0] - 1 / 2) * (points[0] - 1 / 2) -
                                  1000 * (points[1] - 1 / 2) * (points[1] - 1 / 2))) +
                 (std::exp(-1000 * (points[0] - 1 / 2) * (points[0] - 1 / 2) -
                           1000 * (points[1] - 1 / 2) * (points[1] - 1 / 2)) *
                      (2000 * points[1] - 1000) * (2000 * points[1] - 1000) -
                  2000 * std::exp(-1000 * (points[0] - 1 / 2) * (points[0] - 1 / 2) -
                                  1000 * (points[1] - 1 / 2) * (points[1] - 1 / 2))));
    }

    /* ---------------------------------------------------------------------------- */
    // Mild wave front

    MARS_INLINE_FUNCTION Real ex3_exact(const Real *points, Real &result) {
        Real alpha = 20;
        Real x_c = -0.05;
        Real y_c = -0.05;
        Real r_0 = 0.7;

        Real r = sqrt((points[0] - x_c) * (points[0] - x_c) + (points[1] - y_c) * (points[1] - y_c));

        return atan(alpha * (r - r_0));
    }

    MARS_INLINE_FUNCTION Real ex3_laplacian(const Real *points, Real &result) { return }

}  // namespace mars

#endif  // MARS_LAPLACIAN_EXAMPLESs