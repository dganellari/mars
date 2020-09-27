#ifndef MARS_GRADIENT_RECOVERY_HPP
#define MARS_GRADIENT_RECOVERY_HPP

#include "mars_fe_simplex.hpp"
#include "mars_fe_values.hpp"

namespace mars {

    template <class Mesh>
    class GradientRecovery {
    public:
        static const int Dim = Mesh::Dim;
        static const int ManifoldDim = Mesh::ManifoldDim;
        static const int NFuns = ManifoldDim + 1;

        // Compute gradient from solution
        static void grad(const FEValues<Mesh> &values, const ViewVectorType<Real> &u, ViewMatrixType<Real> &grad_u) {
            auto J_inv = values.J_inv();

            ViewMatrixType<Integer> elems = values.mesh().get_view_elements();

            auto mesh = values.mesh();

            if (grad_u.extent(0) != mesh.n_elements()) {
                grad_u = ViewMatrixType<Real>("grad_u", mesh.n_elements(), Dim);
            }

            Kokkos::parallel_for(
                "GradientRecovery::grad", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u_e[NFuns];
                    Real J_inv_e[Dim * Dim];
                    Real grad_u_e[Dim];

                    for (Integer k = 0; k < NFuns; ++k) {
                        u_e[k] = u(elems(i, k));
                    }

                    for (int k = 0; k < Dim * ManifoldDim; ++k) {
                        J_inv_e[k] = J_inv(i, k);
                    }

                    FESimplex<ManifoldDim>::grad(J_inv_e, u_e, grad_u_e);

                    for (int k = 0; k < Dim; ++k) {
                        grad_u(i, k) = grad_u_e[k];
                    }
                });
        }

        static void project(const FEValues<Mesh> &values,
                            const ViewMatrixType<Real> &grad_p0,
                            ViewMatrixType<Real> &grad_p1) {
            auto det_J = values.det_J();

            ViewMatrixType<Integer> elems = values.mesh().get_view_elements();

            auto mesh = values.mesh();

            ViewVectorType<Real> lumped_mass("lumped_mass", mesh.n_nodes());

            if (grad_p1.extent(0) != mesh.n_nodes()) {
                grad_p1 = ViewMatrixType<Real>("grad_p1", mesh.n_nodes(), Dim);
            }

            Kokkos::parallel_for(
                "GradientRecovery::project", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real grad_p[Dim];

                    for (int k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                    }

                    const Real w = det_J(i) / NFuns;

                    for (int d = 0; d < Dim; ++d) {
                        grad_p[d] = grad_p0(i, d) * w;
                    }

                    for (int k = 0; k < NFuns; ++k) {
                        Kokkos::atomic_add(&lumped_mass(idx[k]), w);

                        for (int d = 0; d < Dim; ++d) {
                            Kokkos::atomic_add(&grad_p1(idx[k], d), grad_p[d]);
                        }
                    }
                });

            Kokkos::parallel_for(
                "GradientRecovery::project_remove_mass", mesh.n_nodes(), MARS_LAMBDA(const Integer i) {
                    for (int d = 0; d < Dim; ++d) {
                        grad_p1(i, d) /= lumped_mass(i);
                    }
                });
        }

        static void estimate(const FEValues<Mesh> &values,
                             const ViewMatrixType<Real> &grad_p0,
                             const ViewMatrixType<Real> &grad_p1,
                             ViewVectorType<Real> &error) {
            auto det_J = values.det_J();

            ViewMatrixType<Integer> elems = values.mesh().get_view_elements();

            auto mesh = values.mesh();

            if (error.extent(0) != mesh.n_elements()) {
                error = ViewVectorType<Real>("grad_p1", mesh.n_elements());
            }

            Kokkos::parallel_for(
                "GradientRecovery::estimate", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real grad_p0_e[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        grad_p0_e[d] = grad_p0(i, d);
                    }

                    for (int k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                    }

                    Real error_e = 0.0;

                    for (int k = 0; k < NFuns; ++k) {
                        Real err_k = 0.0;
                        for (int d = 0; d < Dim; ++d) {
                            Real diff_d = grad_p1(idx[k], d) - grad_p0_e[d];
                            err_k += diff_d * diff_d;
                        }

                        error_e += err_k * det_J(i);
                    }

                    error(i) = error_e;
                });
        }

        static void estimate(const FEValues<Mesh> &values, const ViewVectorType<Real> &u, ViewVectorType<Real> &error) {
            ViewMatrixType<Real> grad_p0, grad_p1;
            grad(values, u, grad_p0);
            project(values, grad_p0, grad_p1);
            estimate(values, grad_p0, grad_p1, error);
        }
    };

}  // namespace mars

#endif  // MARS_GRADIENT_RECOVERY_HPP
