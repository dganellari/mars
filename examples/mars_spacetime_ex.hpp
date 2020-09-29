#ifndef MARS_SPACETIME_EXAMPLES
#define MARS_SPACETIME_EXAMPLES

#include <Kokkos_ArithTraits.hpp>
#include <cmath>

#include "mars_base.hpp"
#include "mars_config.hpp"

namespace mars {
    /* Examples for heat equation u_t - laplacian u = f */

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// 2D examples /////////////////////////////////////////////////

    /* ---------------------------------------------------------------------------- */

    /* BC --> zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex1_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return Kokkos::ArithTraits<Real>::sin(M_PI * x) * Kokkos::ArithTraits<Real>::sin(M_PI * y);
    }

    MARS_INLINE_FUNCTION Real ex1_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return M_PI * Kokkos::ArithTraits<Real>::cos(M_PI * y) * Kokkos::ArithTraits<Real>::sin(M_PI * x) -
               (-M_PI * M_PI * Kokkos::ArithTraits<Real>::sin(M_PI * y) * Kokkos::ArithTraits<Real>::sin(M_PI * x));
    }

    /* ---------------------------------------------------------------------------- */

    /* BC --> zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex2_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return x * Kokkos::ArithTraits<Real>::cos(3. * M_PI * x / 2.) * Kokkos::ArithTraits<Real>::sin(3. * y);
    }

    MARS_INLINE_FUNCTION Real ex2_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return M_PI * 3.0 / 4.0 * Kokkos::ArithTraits<Real>::sin(3. * y) *
                   (4 * Kokkos::ArithTraits<Real>::sin(3. * M_PI * x / 2.) +
                    3. * M_PI * x * Kokkos::ArithTraits<Real>::cos(3. * M_PI * x / 2.)) +
               3. * x * Kokkos::ArithTraits<Real>::cos(3. * M_PI * x / 2.) * Kokkos::ArithTraits<Real>::cos(3. * y);
    }

    /* ---------------------------------------------------------------------------- */

    /* BC --> zero dirichlet */

    MARS_INLINE_FUNCTION Real ex3_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return (x * x - x) * (y * y - y) * Kokkos::ArithTraits<Real>::exp(-100 * ((x - y) * (x - y)));
    }

    MARS_INLINE_FUNCTION Real ex3_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];

        return (-Kokkos::ArithTraits<Real>::exp(-100.0 * (y - x) * (y - x)) * (-x * x + x) * (2.0 * y - 1) -
                Kokkos::ArithTraits<Real>::exp(-100 * (y - x) * (y - x)) * (-y * y + y) * (200 * y - 200 * x) *
                    (-x * x + x)) -
               (Kokkos::ArithTraits<Real>::exp(-100.0 * (y - x) * (y - x)) * (-y * y + y) * (200.0 * y - 200.0 * x) *
                    (200.0 * y - 200.0 * x) * (-x * x + x) -
                200.0 * Kokkos::ArithTraits<Real>::exp(-100.0 * (y - x) * (y - x)) * (-y * y + y) * (-x * x + x) -
                2.0 * Kokkos::ArithTraits<Real>::exp(-100.0 * (y - x) * (y - x)) * (-y * y + y) *
                    (200.0 * y - 200.0 * x) * (2.0 * x - 1) -
                2.0 * Kokkos::ArithTraits<Real>::exp(-100 * (y - x) * (y - x)) * (-y * y + y));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// 3D examples /////////////////////////////////////////////////

    // Smooth solution
    /* BC -->  zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex4_3D_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];

        return (1 - x) * x * x * (1 - y) * y * y * (1 - z) * z * z;
    }

    MARS_INLINE_FUNCTION Real ex4_3D_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];

        return (-z * z * x * x * y * y * (x - 1.0) * (y - 1.0) -
                2 * z * x * x * y * y * (z - 1.0) * (x - 1.0) * (y - 1.0)) -
               ((-4.0 * z * z * x * y * y * (z - 1.0) * (y - 1.0) -
                 2 * z * z * y * y * (z - 1.0) * (x - 1.0) * (y - 1.0)) +
                (-4.0 * z * z * x * x * y * (z - 1.0) * (x - 1.0) -
                 2 * z * z * x * x * (z - 1.0) * (x - 1.0) * (y - 1.0)));
    }

    /* ---------------------------------------------------------------------------- */

    // Moving peak
    /* BC -->  zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex5_3D_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];

        return (x * x - x) * (y * y - y) * (z * z - z) *
               Kokkos::ArithTraits<Real>::exp(-100.0 * ((x - z) * (x - z) + (y - z) * (y - z)));
    }

    MARS_INLINE_FUNCTION Real ex5_3D_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];

        return (Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                    (-x * (x) + x) * (-y * (y) + y) * (2 * z - 1) -
                Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100 * (z - y) * (z - y)) * (-z * (z) + z) *
                    (-x * (x) + x) * (-y * (y) + y) * (200.0 * x - 400 * z + 200.0 * y)) -
               ((2.0 * Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (-y * (y) + y) +
                 200.0 * Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (-x * (x) + x) * (-y * (y) + y) +
                 2 * Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (200.0 * z - 200.0 * x) * (-y * (y) + y) * (2 * x - 1) -
                 Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (200.0 * z - 200.0 * x) * (200.0 * z - 200.0 * x) * (-x * (x) + x) *
                     (-y * (y) + y)) +
                (2.0 * Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (-x * (x) + x) +
                 200.0 * Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (-x * (x) + x) * (-y * (y) + y) +
                 2 * Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (200.0 * z - 200.0 * y) * (-x * (x) + x) * (2 * y - 1) -
                 Kokkos::ArithTraits<Real>::exp(-100.0 * (z - x) * (z - x) - 100.0 * (z - y) * (z - y)) *
                     (-z * (z) + z) * (200.0 * z - 200.0 * y) * (200.0 * z - 200.0 * y) * (-x * (x) + x) *
                     (-y * (y) + y)));
    }

    /* ---------------------------------------------------------------------------- */
    // Moving Gaussian peak --> Domain (-1,1)^2 x (O,5)

    /* BC -->  zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex6_3D_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];

        Real A = 1;
        Real sigma = 5e-02;
        Real r = 0.5;
        Real p = 5;

        Real gamma_1 = x - r * (Kokkos::ArithTraits<Real>::cos(2 * M_PI * z / p));
        Real gamma_2 = y - r * (Kokkos::ArithTraits<Real>::sin(2 * M_PI * z / p));

        Real num = gamma_1 * gamma_1 + gamma_2 * gamma_2;

        return A * Kokkos::ArithTraits<Real>::exp(-num / (2.0 * sigma * sigma));
    }

    MARS_INLINE_FUNCTION Real ex6_3D_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];

        return (-Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2)) *
                (80 * M_PI * Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) *
                     (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) -
                 80 * M_PI * Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) *
                     (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2))) -
               ((Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2)) *
                     (400 * x - 200 * Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5)) *
                     (400 * x - 200 * Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5)) -
                 400 * Kokkos::ArithTraits<Real>::exp(
                           -200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) *
                               (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) -
                           200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2) *
                               (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2))) +
                (Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2)) *
                     (400 * y - 200 * Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5)) *
                     (400 * y - 200 * Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5)) -
                 400 * Kokkos::ArithTraits<Real>::exp(
                           -200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) *
                               (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * z) / 5) / 2) -
                           200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2) *
                               (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * z) / 5) / 2))));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// 4D examples /////////////////////////////////////////////////

    // Moving peak solution
    /* BC -->  zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex7_4D_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];
        const Real w = points[3];

        return 100.0 * (x * x - x) * (y * y - y) * (z * z - z) * (w * w - w) *
               Kokkos::ArithTraits<Real>::exp(-100.0 * ((x - w) * (x - w) + (y - w) * (y - w) + (z - w) * (z - w)));
    }

    MARS_INLINE_FUNCTION Real ex7_4D_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];
        const Real w = points[3];

        return (Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                               100.0 * (w - z) * (w - z)) *
                    (-w * w + w) * (-y * y + y) * (-z * z + z) * (-100.0 * x * x + 100.0 * x) *
                    (200 * x - 600 * w + 200 * y + 200 * z) -
                Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                               100.0 * (w - z) * (w - z)) *
                    (-y * y + y) * (-z * z + z) * (2 * w - 1) * (-100.0 * x * x + 100.0 * x)) -
               ((Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (200.0 * w - 200.0 * x) * (200.0 * w - 200.0 * x) * (-y * y + y) * (-z * z + z) *
                     (-100.0 * x * x + 100.0 * x) -
                 200.0 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (-y * y + y) * (-z * z + z) * (-100.0 * x * x + 100.0 * x) -
                 200.0 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (-y * y + y) * (-z * z + z) -
                 2.0 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (200.0 * w - 200.0 * x) * (-y * y + y) * (-z * z + z) * (200.0 * x - 100.0)) +
                (Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (200.0 * w - 200.0 * y) * (200.0 * w - 200.0 * y) * (-y * y + y) * (-z * z + z) *
                     (-100.0 * x * x + 100.0 * x) -
                 200.0 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (-y * y + y) * (-z * z + z) * (-100.0 * x * x + 100.0 * x) -
                 2 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (200.0 * w - 200.0 * y) * (-z * z + z) * (2.0 * y - 1.0) *
                     (-100.0 * x * x + 100.0 * x) -
                 2 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (-z * z + z) * (-100.0 * x * x + 100.0 * x)) +
                (Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (200.0 * w - 200.0 * z) * (200.0 * w - 200.0 * z) * (-y * y + y) * (-z * z + z) *
                     (-100.0 * x * x + 100.0 * x) -
                 200.0 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (-y * y + y) * (-z * z + z) * (-100.0 * x * x + 100.0 * x) -
                 2 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (200.0 * w - 200.0 * z) * (-y * y + y) * (2.0 * z - 1.0) *
                     (-100.0 * x * x + 100.0 * x) -
                 2 *
                     Kokkos::ArithTraits<Real>::exp(-100.0 * (w - x) * (w - x) - 100.0 * (w - y) * (w - y) -
                                                    100.0 * (w - z) * (w - z)) *
                     (-w * w + w) * (-y * y + y) * (-100.0 * x * x + 100.0 * x)));
    }
    /* ---------------------------------------------------------------------------- */
    // Moving Gaussian peak --> Domain (-1,1)^3 x (O,5)

    /* BC -->  zero dirichlet + natural neumann on upper bound */

    MARS_INLINE_FUNCTION Real ex8_4D_st_exact(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];
        const Real w = points[3];

        Real A = 1;
        Real sigma = 5e-02;
        Real r = 0.5;
        Real p = 5;

        Real gamma_1 = x - r * (Kokkos::ArithTraits<Real>::cos(2 * M_PI * w / p));
        Real gamma_2 = y - r * (Kokkos::ArithTraits<Real>::sin(2 * M_PI * w / p));
        Real gamma_3 = z - r * (Kokkos::ArithTraits<Real>::cos(2 * M_PI * w / p));

        Real num = gamma_1 * gamma_1 + gamma_2 * gamma_2 + gamma_3 * gamma_3;

        return A * Kokkos::ArithTraits<Real>::exp(-num / (2.0 * sigma * sigma));
    }

    MARS_INLINE_FUNCTION Real ex8_4D_st_spacetime(const Real *points) {
        const Real x = points[0];
        const Real y = points[1];
        const Real z = points[2];
        const Real w = points[3];

        return (-Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2)) *
                (80 * M_PI * Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) *
                     (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                 80 * M_PI * Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) *
                     (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) +
                 80 * M_PI * Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) *
                     (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2))) -
               ((Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2)) *
                     (400 * x - 200 * Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5)) *
                     (400 * x - 200 * Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5)) -
                 400 * Kokkos::ArithTraits<Real>::exp(
                           -200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                               (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                           200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                               (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                           200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                               (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2))) +
                (Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2)) *
                     (400 * y - 200 * Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5)) *
                     (400 * y - 200 * Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5)) -
                 400 * Kokkos::ArithTraits<Real>::exp(
                           -200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                               (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                           200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                               (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                           200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                               (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2))) +
                (Kokkos::ArithTraits<Real>::exp(-200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                                                    (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                                                200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                                                    (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2)) *
                     (400 * z - 200 * Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5)) *
                     (400 * z - 200 * Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5)) -
                 400 * Kokkos::ArithTraits<Real>::exp(
                           -200 * (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                               (x - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                           200 * (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) *
                               (z - Kokkos::ArithTraits<Real>::cos((2 * M_PI * w) / 5) / 2) -
                           200 * (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2) *
                               (y - Kokkos::ArithTraits<Real>::sin((2 * M_PI * w) / 5) / 2))));
    }

}  // namespace mars

#endif  // MARS_SPACETIME_EXAMPLES