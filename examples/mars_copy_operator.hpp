#ifndef MARS_COPY_OPERATOR_HPP
#define MARS_COPY_OPERATOR_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {

    class CopyOperator final {
    public:
        static void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) { Kokkos::deep_copy(op_x, x); }
        static void precondition_left(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) { apply(x, op_x); }
        static void precondition_right(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) { apply(x, op_x); }
    };

}  // namespace mars

#endif  // MARS_COPY_OPERATOR_HPP
