#ifndef MARS_UMESH_OPERATOR_HPP
#define MARS_UMESH_OPERATOR_HPP

#include <memory>

#include "mars_base.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_identity_operator.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_tensor_laplacian.hpp"

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
        static constexpr int NFuns = Elem::NNodes;

        virtual ~UMeshPreconditioner() {}
        virtual void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) = 0;
        virtual void precondition_left(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) = 0;
        virtual void precondition_right(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) = 0;
        virtual void init(FEValues<Mesh> &values) = 0;

        void set_identity(std::shared_ptr<IdentityOperator> id) { id_ = id; }
        inline const std::shared_ptr<IdentityOperator> &identity() const { return id_; }

    private:
        std::shared_ptr<IdentityOperator> id_;
    };

    template <class Mesh>
    class UMeshJacobiPreconditioner : public UMeshPreconditioner<Mesh> {
    public:
        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {
            auto n = inv_diag_.extent(0);
            auto inv_diag = inv_diag_;

            Kokkos::parallel_for(
                "UMeshJacobiPreconditioner::apply", n, MARS_LAMBDA(const Integer i) { op_x(i) = inv_diag(i) * x(i); });

            // Kokkos::parallel_for(
            //     "UMeshJacobiPreconditioner::copy", n, MARS_LAMBDA(const Integer i) { op_x(i) = x(i); });

            this->identity()->apply(x, op_x);
        }

        void precondition_left(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {
            auto n = inv_diag_.extent(0);
            auto inv_diag = inv_diag_;

            Kokkos::parallel_for(
                "UMeshJacobiPreconditioner::precondition_left", n, MARS_LAMBDA(const Integer i) {
                    // const Real inv_diag_i = Kokkos::ArithTraits<Real>::abs(inv_diag(i));
                    const Real inv_diag_i = inv_diag(i);
                    assert(inv_diag_i > 0.0);
                    op_x(i) = Kokkos::ArithTraits<Real>::sqrt(inv_diag_i) * x(i);
                });

            this->identity()->apply(x, op_x);
        }

        void precondition_right(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {
            precondition_left(x, op_x);
        }

    protected:
        ViewVectorType<Real> inv_diag_;
    };

    template <class Mesh>
    class UMeshOperator {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Elem::NNodes;

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
            auto active = mesh.get_view_active();

            Kokkos::parallel_for(
                "UMeshOperator::assemble_rhs", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real p[Dim];

                    if (!active(i)) return;  // ACTIVE

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

        inline void set_identity(const std::shared_ptr<IdentityOperator> &id) {
            id_ = id;
            preconditioner_->set_identity(id);
        }

        inline void set_precontitioner(const std::shared_ptr<UMeshPreconditioner<Mesh>> &p) { preconditioner_ = p; }

        inline std::shared_ptr<UMeshPreconditioner<Mesh>> &preconditioner() { return preconditioner_; }
        inline const std::shared_ptr<IdentityOperator> &identity() const { return id_; }

        inline FEValues<Mesh> &values() { return values_; }
        inline const FEValues<Mesh> &values() const { return values_; }

    private:
        FEValues<Mesh> values_;
        std::shared_ptr<UMeshPreconditioner<Mesh>> preconditioner_;
        std::shared_ptr<IdentityOperator> id_;
    };

}  // namespace mars

#endif  // MARS_UMESH_OPERATOR_HPP
