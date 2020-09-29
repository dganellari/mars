#ifndef MARS_FE_VALUES_HPP
#define MARS_FE_VALUES_HPP

#include "mars_invert.hpp"

#include <Kokkos_ArithTraits.hpp>

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

    template <class Mesh>
    class FEValues {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static const int Dim = Mesh::Dim;
        static const int ManifoldDim = Mesh::ManifoldDim;
        static const int NFuns = Elem::ElemType;
        static const int NQPoints = 1;
        // template <class Quadrature>
        FEValues(Mesh &mesh) : mesh_(mesh) {}

        void init() {
            ViewMatrixType<Integer> elems = mesh_.get_view_elements();
            ViewMatrixType<Real> points = mesh_.get_view_points();
            auto active = mesh_.get_view_active();

            auto ne = elems.extent(0);
            auto nen = elems.extent(1);

            det_J_ = ViewVectorType<Real>("det_J", mesh_.n_elements());
            J_inv_ = ViewMatrixType<Real>("J_inv", mesh_.n_elements(), Dim * Dim);

            auto det_J = det_J_;
            auto J_inv = J_inv_;

            const Integer n_nodes = mesh_.n_nodes();

            Kokkos::parallel_for(
                "FEValues::init", mesh_.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real J[Dim * Dim], J_inv_e[Dim * Dim];
                    Real p0[Dim], pk[Dim];

                    for (int k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        assert(idx[k] < n_nodes);
                    }

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

                    Real e_det_J = 0.0;
                    Invert<Dim>::apply(J, J_inv_e, e_det_J);

                    // Turn off inactive contributions
                    // det_J(i) = Kokkos::ArithTraits<Real>::abs(e_det_J) * active(i);
                    det_J(i) = Kokkos::ArithTraits<Real>::abs(e_det_J);

                    for (int k = 0; k < (Dim * Dim); ++k) {
                        J_inv(i, k) = J_inv_e[k];
                    }

                    assert(has_non_zero<Dim>(J_inv_e));
                    // assert(has_non_zero<Dim>(&J_inv(i, 0)));

                    assert(e_det_J == e_det_J);
                    assert(e_det_J != 0.0);
                    assert(e_det_J > 0.0);

                    if (e_det_J < 0) {
                        static const int dim = Dim;
                        printf("found element with wrong orientation\n");
                        Integer temp = elems(i, dim);
                        elems(i, dim) = elems(i, dim - 1);
                        elems(i, dim - 1) = temp;
                    }
                });

            //     Kokkos::parallel_for(
            //         "print_det_J", mesh_.n_elements(), MARS_LAMBDA(const Integer i) { printf("%g\n", det_J(i)); });

            Real measure = KokkosBlas::nrm1(det_J_);

            std::cout << "measure: " << measure << std::endl;
        }

        MARS_INLINE_FUNCTION Mesh &mesh() { return mesh_; }
        MARS_INLINE_FUNCTION const Mesh &mesh() const { return mesh_; }

        MARS_INLINE_FUNCTION ViewVectorType<Real> det_J() const { return det_J_; }
        MARS_INLINE_FUNCTION ViewMatrixType<Real> J_inv() const { return J_inv_; }

    private:
        Mesh &mesh_;

    public:
        ViewMatrixType<Real> J_inv_;
        ViewVectorType<Real> det_J_;
    };

}  // namespace mars

#endif  // MARS_FE_VALUES_HPP
