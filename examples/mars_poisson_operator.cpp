#include "mars_poisson_operator.hpp"
#include "mars_precon_conjugate_grad.hpp"

namespace mars {

    void PoissonOperator::apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
        // gather ghost data
        // form_operator(dm, ...);
        // scatter ghost data
        // apply identity operator for BC (constrained nodes)
    }
}  // namespace mars