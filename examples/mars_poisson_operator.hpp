#ifndef MARS_POISSON_OPERATOR_HPP
#define MARS_POISSON_OPERATOR_HPP

#include "mars_matrix_free_operator.hpp"
#include "mars_poisson.hpp"

namespace mars {

    template<class DM, Integer INPUT, Integer OUTPUT>
    class PoissonOperator final : public Operator<DM> {
    public:
        PoissonOperator(const context &c, DM &dm) : Super(c, dm) {}

        void init() override {
            this.values().init();
        }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override;

    };

    class CopyOperator final {
    public:
        static void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) { Kokkos::deep_copy(op_x, x); }
    };

}  // namespace mars

#endif
