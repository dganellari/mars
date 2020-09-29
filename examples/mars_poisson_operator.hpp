#ifndef MARS_POISSON_OPERATOR_HPP
#define MARS_POISSON_OPERATOR_HPP

#include "mars_matrix_free_operator.hpp"
#include "mars_poisson.hpp"

namespace mars {

    class PoissonOperator final : public Operator {
    public:
        using DMType = mars::DMQ2;

        PoissonOperator(DMType &dm, const context &c) : Operator(c) {}

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override;
    };

    class CopyOperator final {
    public:
        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) { Kokkos::deep_copy(op_x, x); }
    };

}  // namespace mars

#endif