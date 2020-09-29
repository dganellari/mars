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
            auto values = this->values();
            auto det_J = values.det_J();
            auto J_inv = values.J_inv();
            auto mesh = values.mesh();

            ViewMatrixType<Integer> elems = mesh.get_view_elements();

            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            auto quad = quad_;

            Kokkos::parallel_for(
                "FitzHughNagumo::apply", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real I_ion[NFuns];
                    Integer idx[NFuns];
                    Real J_inv_e[Dim * Dim];

                    for (int k = 0; k < (Dim * Dim); ++k) {
                        s J_inv_e[k] = J_inv(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        u[k] = x(idx[k]);
                    }

                    Real u_elem;

                    // for (Integer i = 0; i < quad.n_points(); ++i) {
                    //     for (Integer k = 0; k < NFuns; ++k) {
                    //         u_elem += FESimplex::fun(k, quad.points[i])
                    //     }
                    // }

                    // SpaceTimeMixed<Mesh>::one_thread_eval(J_inv_e, det_J(i), quad, u, Au);

                    // for (Integer k = 0; k < NFuns; ++k) {
                    //     Kokkos::atomic_add(&op_x(idx[k]), Au[k]);
                    // }
                });

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
