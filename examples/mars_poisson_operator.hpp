#ifndef MARS_POISSON_OPERATOR_HPP
#define MARS_POISSON_OPERATOR_HPP

#include "mars_matrix_free_operator.hpp"
#include "mars_poisson.hpp"

namespace mars {

class PoissonOperator final : public Operator {
public:
  using DMType = mars::DMQ2;

  PoissonOperator(DMType &dm, const context &c) : Operator(c) {}

  void apply(const ViewVectorType<Real> &x,
             ViewVectorType<Real> &op_x) override;
};

} // namespace mars

#endif