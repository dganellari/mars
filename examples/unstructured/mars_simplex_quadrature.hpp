#ifndef MARS_SIMPLEX_QUADRATURE_HPP
#define MARS_SIMPLEX_QUADRATURE_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {

    template <typename T, int N>
    using ViewQWeights =
        Kokkos::View<T[N], Kokkos::LayoutRight, KokkosSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    template <typename T, Integer XDim_, Integer YDim_>
    using ViewQPoints =
        Kokkos::View<T[XDim_][YDim_], Kokkos::LayoutRight, KokkosSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    template <int Dim, int Order>
    class SimplexQuadrature {};

    template <int Dim>
    class SimplexQuadrature<Dim, 1> {
    public:
        static SimplexQuadrature make() {
            SimplexQuadrature ret;
            ret.init();
            return ret;
        }

        ViewQPoints<Real, Dim + 1, 2> points;
        ViewQWeights<Real, Dim + 1> weights;

        MARS_INLINE_FUNCTION static constexpr int n_points() { return 6; }
        MARS_INLINE_FUNCTION static constexpr int dim() { return 2; }

        SimplexQuadrature() : points("q_points"), weights("q_weights") {}

        void init() {
            auto points_tmp = points;
            auto weights_tmp = weights;

            Kokkos::parallel_for(
                1, MARS_LAMBDA(const int &) {
                    weights_tmp(0) = 1. / (Dim + 1);

                    for (int p = 0; p < Dim; ++p) {
                        weights_tmp(p + 1) = 1. / (Dim + 1);
                        points_tmp(p + 1, p) = 1.0;
                    }
                });
        }
    };

    template <>
    class SimplexQuadrature<2, 2> {
    public:
        static SimplexQuadrature make() {
            SimplexQuadrature ret;
            ret.init();
            return ret;
        }

        ViewQPoints<Real, 6, 2> points;
        ViewQWeights<Real, 6> weights;

        MARS_INLINE_FUNCTION static constexpr int n_points() { return 6; }
        MARS_INLINE_FUNCTION static constexpr int dim() { return 2; }

        SimplexQuadrature() : points("q_points"), weights("q_weights") {}

        void init() {
            auto points_tmp = points;
            auto weights_tmp = weights;

            Kokkos::parallel_for(
                1, MARS_LAMBDA(const int &) {
                    Real pts[6][2] = {{0.5, 0.5},
                                      {0.5, 0.0},
                                      {0.0, 0.5},
                                      {1.0 / 6.0, 1.0 / 6.0},
                                      {1.0 / 6.0, 2.0 / 3.0},
                                      {2.0 / 3.0, 1.0 / 6.0}};

                    Real w[6] = {1.0 / 30.0, 1.0 / 30.0, 1.0 / 30.0, 0.3, 0.3, 0.3};

                    for (int p = 0; p < n_points(); ++p) {
                        weights_tmp(p) = w[p];

                        for (int d = 0; d < dim(); ++d) {
                            points_tmp(p, d) = pts[p][d];
                        }
                    }
                });
        }
    };

    template <>
    class SimplexQuadrature<2, 6> {
    public:
        static SimplexQuadrature make() {
            SimplexQuadrature ret;
            ret.init();
            return ret;
        }

        ViewQPoints<Real, 13, 2> points;
        ViewQWeights<Real, 13> weights;

        MARS_INLINE_FUNCTION static constexpr int n_points() { return 13; }
        MARS_INLINE_FUNCTION static constexpr int dim() { return 2; }

        void init() {
            auto points_tmp = points;
            auto weights_tmp = weights;

            Kokkos::parallel_for(
                1, MARS_LAMBDA(const int &) {
                    Real pts[13][2] = {{0.063089014491502228340331602870819, 0.063089014491502228340331602870819},
                                       {0.063089014491502228340331602870819, 0.87382197101699554331933679425836},
                                       {0.87382197101699554331933679425836, 0.063089014491502228340331602870819},
                                       {0.24928674517091042129163855310702, 0.24928674517091042129163855310702},
                                       {0.24928674517091042129163855310702, 0.50142650965817915741672289378596},
                                       {0.50142650965817915741672289378596, 0.24928674517091042129163855310702},
                                       {0.053145049844816947353249671631398, 0.31035245103378440541660773395655},
                                       {0.053145049844816947353249671631398, 0.63650249912139864723014259441205},
                                       {0.31035245103378440541660773395655, 0.053145049844816947353249671631398},
                                       {0.31035245103378440541660773395655, 0.63650249912139864723014259441205},
                                       {0.63650249912139864723014259441205, 0.053145049844816947353249671631398},
                                       {0.63650249912139864723014259441205, 0.31035245103378440541660773395655}};

                    Real w[13] = {0.050844906370206816920936809106869,
                                  0.050844906370206816920936809106869,
                                  0.050844906370206816920936809106869,
                                  0.11678627572637936602528961138558,
                                  0.11678627572637936602528961138558,
                                  0.11678627572637936602528961138558,
                                  0.082851075618373575193553456420442,
                                  0.082851075618373575193553456420442,
                                  0.082851075618373575193553456420442,
                                  0.082851075618373575193553456420442,
                                  0.082851075618373575193553456420442,
                                  0.082851075618373575193553456420442};

                    for (int p = 0; p < n_points(); ++p) {
                        weights_tmp(p) = w[p];

                        for (int d = 0; d < dim(); ++d) {
                            points_tmp(p, d) = pts[p][d];
                        }
                    }
                });
        }
    };

    template <int Dim>
    using SimplexQuadratureLinear = mars::SimplexQuadrature<Dim, 1>;

    template <int Dim>
    using SimplexQuadratureQuartic = mars::SimplexQuadrature<Dim, 6>;

}  // namespace mars

#endif  // MARS_SIMPLEX_QUADRATURE_HPP
