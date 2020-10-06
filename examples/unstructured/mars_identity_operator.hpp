#ifndef MARS_IDENTITY_OPERATOR_HPP
#define MARS_IDENTITY_OPERATOR_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {

    class IdentityOperator {
    public:
        template <class Mesh, class BC>
        void init(Mesh &mesh, BC bc) {
            static const int Dim = Mesh::Dim;

            auto points = mesh.points();
            ViewVectorType<bool> is_boundary("is_boundary", mesh.n_nodes());

            Kokkos::parallel_for(
                "IdentityOperator::init", mesh.n_nodes(), MARS_LAMBDA(const Integer i) {
                    Real p[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points(i, d);
                    }
                    is_boundary(i) = bc.is_boundary(p);
                });

            is_boundary_ = is_boundary;
        }

        void apply(const ViewVectorType<Real> &input, ViewVectorType<Real> &x) {
            auto is_boundary = is_boundary_;

            Kokkos::parallel_for(
                "IdentityOperator::apply", is_boundary_.extent(0), MARS_LAMBDA(const Integer i) {
                    x(i) = x(i) * (!is_boundary(i)) + input(i) * is_boundary(i);
                });
        }

    private:
        ViewVectorType<bool> is_boundary_;
    };
}  // namespace mars

#endif  // MARS_IDENTITY_OPERATOR_HPP
