#ifndef MARS_UMESH_OPERATOR_HPP
#define MARS_UMESH_OPERATOR_HPP

#include <memory>

#include "mars_base.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_identity_operator.hpp"
#include "mars_simplex_laplacian.hpp"

namespace mars {

    class FakeComm {
    public:
        static Real sum(const Real &v) { return v; }
    };

    template <class Mesh>
    class UMeshPreconditioner {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        virtual ~UMeshPreconditioner() {}
        virtual void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) = 0;
        virtual void init(FEValues<Mesh> &values) = 0;

        void set_identity(std::shared_ptr<IdentityOperator> id) { id_ = id; }
        inline const std::shared_ptr<IdentityOperator> &identity() const { return id_; }

    private:
        std::shared_ptr<IdentityOperator> id_;
    };

    template <class Mesh>
    class JacobiPreconditioner : public UMeshPreconditioner<Mesh> {
    public:
        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {
            auto n = inv_diag_.extent(0);
            auto inv_diag = inv_diag_;

            Kokkos::parallel_for(
                "JacobiPreconditioner::apply", n, MARS_LAMBDA(const Integer i) { op_x(i) = inv_diag(i) * x(i); });

            // Kokkos::parallel_for(
            //     "JacobiPreconditioner::copy", n, MARS_LAMBDA(const Integer i) { op_x(i) = x(i); });

            this->identity()->apply(x, op_x);
        }

        void set_identity(std::shared_ptr<IdentityOperator> id) { id_ = id; }

    private:
        ViewVectorType<Real> inv_diag_;
    };

    template <class Mesh>
    class UMeshOperator {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        UMeshOperator(Mesh &mesh) : values_(mesh) {}
        virtual ~UMeshOperator() {}

        FakeComm comm() { return FakeComm(); }

        virtual void init() = 0;
        virtual void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) = 0;

        template <class F>
        void assemble_rhs(ViewVectorType<Real> &rhs, F f) {
            auto det_J = values_.det_J();

            ViewMatrixType<Integer> elems = values_.mesh().get_view_elements();
            ViewMatrixType<Real> points = values_.mesh().get_view_points();

            auto mesh = values_.mesh();

            Kokkos::parallel_for(
                "UMeshOperator::assemble_rhs", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
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

        void set_identity(const std::shared_ptr<IdentityOperator> &id) {
            id_ = id;
            preconditioner_->set_identity(id);
        }

        inline UMeshPreconditioner<Mesh> &preconditioner() { return preconditioner_; }
        inline const std::shared_ptr<IdentityOperator> &identity() const { return id_; }

    protected:
        FEValues<Mesh> values_;
        std::shared_ptr<UMeshPreconditioner<Mesh>> preconditioner_;
        std::shared_ptr<IdentityOperator> id_;
    };

}  // namespace mars

#endif  // MARS_UMESH_OPERATOR_HPP
