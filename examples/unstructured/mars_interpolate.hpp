#ifndef MARS_INTERPOLATE_HPP
#define MARS_INTERPOLATE_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {

    template <class Mesh>
    class Interpolate {
    public:
        Interpolate(Mesh &mesh) : mesh_(mesh) {}

        template <typename F>
        void apply(ViewVectorType<Real> &v, F fun) {
            static const int Dim = Mesh::Dim;

            ViewMatrixType<Real> points = mesh_.get_view_points();

            Kokkos::parallel_for(
                "Interpolate::apply", mesh_.n_nodes(), MARS_LAMBDA(const Integer i) {
                    Real p[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points(i, d);
                    }

                    v(i) = fun(p);
                });
        }

    private:
        Mesh &mesh_;
    };

}  // namespace mars

#endif  // MARS_INTERPOLATE_HPP
