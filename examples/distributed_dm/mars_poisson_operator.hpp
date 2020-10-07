#ifndef MARS_POISSON_OPERATOR_HPP
#define MARS_POISSON_OPERATOR_HPP

#include "mars_copy_operator.hpp"
#include "mars_matrix_free_operator.hpp"

namespace mars {

    template <class DM, Integer INPUT, Integer OUTPUT, Integer RHS>
    class PoissonOperator final : public Operator<DM> {
    public:
        template <Integer idx>
        using DMDataType = typename DM::template UserDataType<idx>;

        PoissonOperator(const context &c, DM &dm) : Operator<DM>(c, dm) {}

        void init() override { this->values().init(); }

        using InputVector = ViewVectorType<DMDataType<INPUT>>;
        using OutputVector = ViewVectorType<DMDataType<OUTPUT>>;
        using RhsVector = ViewVectorType<DMDataType<RHS>>;

        void apply(const InputVector &x, OutputVector &op_x) override {
            /* auto context = this->ctx(); */
            auto fe_values = this->values();
            auto dm = fe_values.dm();

            dm.template set_locally_owned_data<INPUT>(x);

            // specify the tuple indices of the tuplelements that are needed to gather.
            // if no index specified it gathers all views of the tuple. All data.
            dm.template gather_ghost_data<INPUT>(this->ctx());

            fe_values.template form_operator<INPUT, OUTPUT>();

            scatter_add_ghost_data<DM, OUTPUT>(dm, this->ctx());

            // cleanup
            /* dm.boundary_dof_iterate<INPUT>(MARS_LAMBDA(const Integer local_dof, INPUTDT &value) { */
            dm.template boundary_dof_iterate<INPUT>(MARS_LAMBDA(const Integer local_dof, DMDataType<INPUT> &value) {
                dm.template get_dof_data<OUTPUT>(local_dof) = value;
            });

            dm.template get_locally_owned_data<OUTPUT>(op_x);
        }

        template <class F>
        void assemble_rhs(F f, RhsVector &rhs) {
            auto fe_values = this->values();
            auto dm = fe_values.dm();

            fe_values.template assemble_local_rhs<F, RHS>(f);

            scatter_add_ghost_data<DM, RHS>(dm, this->ctx());

            dm.template get_locally_owned_data<RHS>(rhs);
        }
    };

}  // namespace mars

#endif
