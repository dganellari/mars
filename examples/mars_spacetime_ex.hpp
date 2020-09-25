#ifndef MARS_SPACETIME_EXAMPLES
#define MARS_SPACETIME_EXAMPLES

#include <Kokkos_ArithTraits.hpp>
#include <cmath>

#include "mars_base.hpp"
#include "mars_config.hpp"

namespace mars {
    /* Examples for heat equation u_t - laplacian u = f */

    /* ---------------------------------------------------------------------------- */

    /* BC --> zero dirichlet on the whole boundary */

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

}  // namespace mars

#endif  // MARS_SPACETIME_EXAMPLES