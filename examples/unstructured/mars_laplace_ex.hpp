#ifndef MARS_LAPLACIAN_EXAMPLES
#define MARS_LAPLACIAN_EXAMPLES

#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include <Kokkos_ArithTraits.hpp>
#include <cmath>

#include "mars_base.hpp"
#include "mars_config.hpp"

namespace mars {

    /* Examples for -laplacian u = f */

    // Smooth solution
    MARS_INLINE_FUNCTION Real ex1_exact(const Real *points) {
        return points[0] * points[1] * (1 - points[0]) * (1 - points[1]);
    }

    MARS_INLINE_FUNCTION Real ex1_laplacian(const Real *points) {
        return -2.0 * points[0] * (1 - points[0]) - 2.0 * points[1] * (1 - points[1]);
    }

    /* ---------------------------------------------------------------------------- */

    // Peak solution
    MARS_INLINE_FUNCTION Real ex2_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return Kokkos::ArithTraits<Real>::exp(-1000 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
    }

    MARS_INLINE_FUNCTION Real ex2_laplacian(const Real *x) {
        return -(Kokkos::ArithTraits<Real>::exp(-1000.0 * (x[0] - 1 / 2.0) * (x[0] - 1 / 2.0) -
                                                1000.0 * (x[1] - 1 / 2.0) * (x[1] - 1 / 2.0)) *
                     (2000 * x[0] - 1000.0) * (2000 * x[0] - 1000.0) -
                 2000.0 * Kokkos::ArithTraits<Real>::exp(-1000.0 * (x[0] - 1 / 2.0) * (x[0] - 1 / 2.0) -
                                                         1000.0 * (x[1] - 1 / 2.0) * (x[1] - 1 / 2.0))) -
               (Kokkos::ArithTraits<Real>::exp(-1000.0 * (x[0] - 1 / 2.0) * (x[0] - 1 / 2.0) -
                                               1000.0 * (x[1] - 1 / 2.0) * (x[1] - 1 / 2.0)) *
                    (2000 * x[1] - 1000.0) * (2000 * x[1] - 1000.0) -
                2000 * Kokkos::ArithTraits<Real>::exp(-1000.0 * (x[0] - 1 / 2.0) * (x[0] - 1 / 2.0) -
                                                      1000.0 * (x[1] - 1 / 2.0) * (x[1] - 1 / 2.0)));
    }

    /* ---------------------------------------------------------------------------- */
    // Mild wave front

    MARS_INLINE_FUNCTION Real ex3_exact(const Real *points) {
        Real alpha = 20;
        Real x_c = -0.05;
        Real y_c = -0.05;
        Real r_0 = 0.7;

        Real r = Kokkos::ArithTraits<Real>::sqrt((points[0] - x_c) * (points[0] - x_c) +
                                                 (points[1] - y_c) * (points[1] - y_c));

        return Kokkos::ArithTraits<Real>::atan(alpha * (r - r_0));
    }

    MARS_INLINE_FUNCTION Real ex3_laplacian(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return -(20 / (((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                              (y + 1 / 20.0) * (y + 1 / 20.0)) -
                         14) *
                            (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                  (y + 1 / 20.0) * (y + 1 / 20.0)) -
                             14) +
                        1) *
                       Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                       (y + 1 / 20.0) * (y + 1 / 20.0))) -
                 (5 * (2 * x + 1 / 10) * (2 * x + 1 / 10)) /
                     (((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                             (y + 1 / 20.0) * (y + 1 / 20.0)) -
                        14) *
                           (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                 (y + 1 / 20.0) * (y + 1 / 20.0)) -
                            14) +
                       1) *
                      Kokkos::ArithTraits<Real>::pow(
                          ((x + 1 / 20.0) * (x + 1 / 20.0) + (y + 1 / 20.0) * (y + 1 / 20.0)), (3 / 2))) -
                 (200 * (2 * x + 1 / 10) * (2 * x + 1 / 10) *
                  (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                        (y + 1 / 20.0) * (y + 1 / 20.0)) -
                   14)) /
                     (((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                             (y + 1 / 20.0) * (y + 1 / 20.0)) -
                        14) *
                           (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                 (y + 1 / 20.0) * (y + 1 / 20.0)) -
                            14) +
                       1) *
                      ((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                             (y + 1 / 20.0) * (y + 1 / 20.0)) -
                        14) *
                           (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                 (y + 1 / 20.0) * (y + 1 / 20.0)) -
                            14) +
                       1) *
                      ((x + 1 / 20.0) * (x + 1 / 20.0) + (y + 1 / 20.0) * (y + 1 / 20.0)))) -
               (20 / (((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                             (y + 1 / 20.0) * (y + 1 / 20.0)) -
                        14) *
                           (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                 (y + 1 / 20.0) * (y + 1 / 20.0)) -
                            14) +
                       1) *
                      Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                      (y + 1 / 20.0) * (y + 1 / 20.0))) -
                (5 * (2 * y + 1 / 10) * (2 * y + 1 / 10)) /
                    (((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                            (y + 1 / 20.0) * (y + 1 / 20.0)) -
                       14) *
                          (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                (y + 1 / 20.0) * (y + 1 / 20.0)) -
                           14) +
                      1) *
                     Kokkos::ArithTraits<Real>::pow(((x + 1 / 20.0) * (x + 1 / 20.0) + (y + 1 / 20.0) * (y + 1 / 20.0)),
                                                    (3 / 2))) -
                (200 * (2 * y + 1 / 10) * (2 * y + 1 / 10) *
                 (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                       (y + 1 / 20.0) * (y + 1 / 20.0)) -
                  14)) /
                    (((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                            (y + 1 / 20.0) * (y + 1 / 20.0)) -
                       14) *
                          (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                (y + 1 / 20.0) * (y + 1 / 20.0)) -
                           14) +
                      1) *
                     ((20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                            (y + 1 / 20.0) * (y + 1 / 20.0)) -
                       14) *
                          (20 * Kokkos::ArithTraits<Real>::sqrt((x + 1 / 20.0) * (x + 1 / 20.0) +
                                                                (y + 1 / 20.0) * (y + 1 / 20.0)) -
                           14) +
                      1) *
                     ((x + 1 / 20.0) * (x + 1 / 20.0) + (y + 1 / 20.0) * (y + 1 / 20.0))));
    }

}  // namespace mars
#endif
#endif  // MARS_LAPLACIAN_EXAMPLESs