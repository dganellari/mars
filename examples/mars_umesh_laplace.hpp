#ifndef MARS_UMESH_LAPLACE_HPP
#define MARS_UMESH_LAPLACE_HPP

#include <memory>

#include "mars_base.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_identity_operator.hpp"
#include "mars_simplex_laplacian.hpp"

namespace mars {

    template <class Mesh>
    class UMeshLaplace final {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        UMeshLaplace(Mesh &mesh) : values_(mesh) {}

        class FakeComm {
        public:
            static Real sum(const Real &v) { return v; }
        };

        FakeComm comm() { return FakeComm(); }

        void init() {
            values_.init();
            preconditioner_.init(values_);
        }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
            // For Kokkos-Cuda
            auto det_J = values_.det_J();
            auto J_inv = values_.J_inv();
            auto mesh = values_.mesh();

            ViewMatrixType<Integer> elems = mesh.get_view_elements();

            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            Kokkos::parallel_for(
                "UMeshLaplace::apply", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real Au[NFuns];
                    Integer idx[NFuns];
                    Real J_inv_e[Dim * Dim];

                    for (int k = 0; k < (Dim * Dim); ++k) {
                        J_inv_e[k] = J_inv(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        u[k] = x(idx[k]);
                    }

                    SimplexLaplacian<Mesh>::one_thread_eval(J_inv_e, det_J(i), u, Au);

                    for (Integer k = 0; k < NFuns; ++k) {
                        Kokkos::atomic_add(&op_x(idx[k]), Au[k]);
                    }
                });

            id_->apply(x, op_x);
        }

        template <class F>
        void assemble_rhs(ViewVectorType<Real> &rhs, F f) {
            auto det_J = values_.det_J();

            ViewMatrixType<Integer> elems = values_.mesh().get_view_elements();
            ViewMatrixType<Real> points = values_.mesh().get_view_points();

            auto mesh = values_.mesh();

            Kokkos::parallel_for(
                "UMeshLaplace::assemble_rhs", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real p[Dim];
                    Real det_J_e = det_J(i);

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        for (Integer d = 0; d < Dim; ++d) {
                            p[d] = points(idx[k], d);
                        }

                        const Real val = f(p);
                        const Real scaled_val = val * det_J_e / NFuns;
                        Kokkos::atomic_add(&rhs(idx[k]), scaled_val);
                    }
                });
        }

        class JacobiPreconditioner {
        public:
            void init(FEValues<Mesh> &values) {
                auto mesh = values.mesh();
                ViewMatrixType<Integer> elems = values.mesh().get_view_elements();
                auto det_J = values.det_J();
                auto J_inv = values.J_inv();

                ViewVectorType<Real> inv_diag("inv_diag", mesh.n_nodes());

                Kokkos::parallel_for(
                    "JacobiPreconditioner::init", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                        Integer idx[NFuns];
                        Real val[NFuns];
                        Real J_inv_e[Dim * Dim];
                        const Real det_J_e = det_J(i);

                        assert(det_J_e > 0.0);

                        for (Integer k = 0; k < (Dim * Dim); ++k) {
                            J_inv_e[k] = J_inv(i, k);
                        }

                        for (Integer k = 0; k < NFuns; ++k) {
                            idx[k] = elems(i, k);
                        }

                        SimplexLaplacian<Mesh>::one_thread_eval_diag(J_inv_e, det_J_e, val);

                        for (Integer k = 0; k < NFuns; ++k) {
                            assert(val[k] != 0.0);
                            Real inv_val = 1. / val[k];

                            assert(inv_val == inv_val);

                            Kokkos::atomic_add(&inv_diag(idx[k]), inv_val);
                        }
                    });

                inv_diag_ = inv_diag;
            }

            void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
                auto n = inv_diag_.extent(0);
                auto inv_diag = inv_diag_;

                Kokkos::parallel_for(
                    "JacobiPreconditioner::apply", n, MARS_LAMBDA(const Integer i) { op_x(i) = inv_diag(i) * x(i); });

                // Kokkos::parallel_for(
                //     "JacobiPreconditioner::copy", n, MARS_LAMBDA(const Integer i) { op_x(i) = x(i); });

                id_->apply(x, op_x);
            }

            void set_identity(std::shared_ptr<IdentityOperator> id) { id_ = id; }

        private:
            ViewVectorType<Real> inv_diag_;
            std::shared_ptr<IdentityOperator> id_;
        };

        void set_identity(const std::shared_ptr<IdentityOperator> &id) {
            id_ = id;
            preconditioner_.set_identity(id);
        }

        inline JacobiPreconditioner &preconditioner() { return preconditioner_; }

        FEValues<Mesh> values_;
        JacobiPreconditioner preconditioner_;
        std::shared_ptr<IdentityOperator> id_;
    };

}  // namespace mars

#endif  // MARS_UMESH_LAPLACE_HPP
