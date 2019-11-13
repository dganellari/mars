#ifndef mars_gauss_points_hpp
#define mars_gauss_points_hpp

#include "Kokkos_Vector.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#include "generation/mars_utils_kokkos.hpp"

namespace mars{

    template<int Dim, int Order>
    class QGauss {};

    template<>
    class QGauss<3, 2> {
    public:
        static constexpr const Real a = 0.333333333333;
        static constexpr const Real b = 1.0;
        static constexpr const Real c = 0.0;
        static constexpr const int N = 8;

        MARS_INLINE_FUNCTION constexpr QGauss()
        : points{
            c,c,c,
            b,c,c,
            c,b,c,
            c,c,b,
            a,a,c,
            a,c,a,
            a,c,c,
            a,a,a
        },
        weights{
            0.025, 0.025, 0.025, 0.025, 0.225, 0.225, 0.225, 0.225
        }
        {}

        MARS_INLINE_FUNCTION constexpr Real sum_weights() const
        {
            return weights[0] + weights[1] + weights[2] + weights[3] +  weights[4] +  weights[5] + weights[6] + weights[7];
        }

        MARS_INLINE_FUNCTION constexpr const Real *point(const int i) const
        {
            return &points[i*3];
        }


        Real points[8 * 3];
        Real weights[8];
    };

    // class gauss_quadrature_rule{
    // public:
    //     static const int n_qG2_points = 4;
    //     static const int Dim = 3;

    //     void reserve_qp_points()
    //     {
    //         points_ = ViewMatrixType<Real>("n_qG2_pts", n_qG2_points, Dim);
    //     }

    //     ViewMatrixType<Real> points_;
    //     Kokkos::View<Real*> wts_;

    // };
}

#endif
