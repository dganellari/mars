#ifndef MARS_FITZHUGH_NAGUMO_HPP
#define MARS_FITZHUGH_NAGUMO_HPP

#include "mars_simplex_quadrature.hpp"

namespace mars {

    template <class Mesh>
    class FitzHughNagumo {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;
        using Super = mars::UMeshOperator<Mesh>;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        FitzHughNagumo(FEValues<Mesh> &values) : values_(values) {}

        void init() { quad_ = SimplexQuadratureQuartic<Dim>::make(); }
        void jacobian(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
            // Kokkos::parallel_for(
            //     "UMeshOperator::assemble_rhs", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
            //         Integer idx[NFuns];
            //         Real p[Dim];
            //         Real det_J_e = det_J(i);

            //         for (Integer k = 0; k < NFuns; ++k) {
            //             idx[k] = elems(i, k);
            //         }

            //         for (Integer k = 0; k < NFuns; ++k) {
            //             for (Integer d = 0; d < Dim; ++d) {
            //                 p[d] = points(idx[k], d);
            //             }

            //             const Real val = f(p);
            //             const Real scaled_val = val * det_J_e / NFuns;
            //             Kokkos::atomic_add(&rhs(idx[k]), scaled_val);
            //         }
            //     });
        }

        void fun(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
            // Kokkos::parallel_for(
            //     "UMeshOperator::assemble_rhs", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
            //         Integer idx[NFuns];
            //         Real p[Dim];
            //         Real det_J_e = det_J(i);

            //         for (Integer k = 0; k < NFuns; ++k) {
            //             idx[k] = elems(i, k);
            //         }

            //         for (Integer k = 0; k < NFuns; ++k) {
            //             for (Integer d = 0; d < Dim; ++d) {
            //                 p[d] = points(idx[k], d);
            //             }

            //             const Real val = f(p);
            //             const Real scaled_val = val * det_J_e / NFuns;
            //             Kokkos::atomic_add(&rhs(idx[k]), scaled_val);
            //         }
            //     });
        }

    private:
        FEValues<Mesh> &values_;
        SimplexQuadratureQuartic<Dim> quad_;
    };

}  // namespace mars

#endif  // MARS_FITZHUGH_NAGUMO_HPP
