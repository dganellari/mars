#ifndef MARS_POISSON_OPERATOR_HPP
#define MARS_POISSON_OPERATOR_HPP

#ifdef WITH_KOKKOS_KERNELS
#include "mars_copy_operator.hpp"
#include "mars_matrix_free_operator.hpp"

namespace mars {

    template <Integer INPUT, Integer OUTPUT, Integer RHS, class DM, class FEM>
    class PoissonOperator /*final*/ : public Operator<DM, FEM> {
    public:
        template <Integer idx>
        using DMDataType = typename DM::template UserDataType<idx>;

        PoissonOperator(const context &c, DM &dm, FEM &fe) : Operator<DM, FEM>(c, dm, fe) {}

        void init() override { this->values().init(); }

        using InputVector = ViewVectorType<DMDataType<INPUT>>;
        using OutputVector = ViewVectorType<DMDataType<OUTPUT>>;
        using RhsVector = ViewVectorType<DMDataType<RHS>>;

        void apply(const InputVector &x, OutputVector &op_x) override {
            /* auto context = this->ctx(); */
            auto fe_values = this->values();
            auto dm = fe_values.dm();
            auto dof_handler = dm.get_dof_handler();

            dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) { dm.template get_dof_data<OUTPUT>(i) = 0.0; });

            dm.template set_locally_owned_data<INPUT>(x);

            // specify the tuple indices of the tuplelements that are needed to gather.
            // if no index specified it gathers all views of the tuple. All data.
            dm.template gather_ghost_data<INPUT>();

            fe_values.template form_operator<INPUT, OUTPUT>();

            scatter_add_ghost_data<DM, OUTPUT>(dm);

            // cleanup
            /* dm.boundary_dof_iterate<INPUT>(MARS_LAMBDA(const Integer local_dof, INPUTDT &value) { */
            dof_handler.boundary_dof_iterate(MARS_LAMBDA(const Integer local_dof) {
                dm.template get_dof_data<OUTPUT>(local_dof) = dm.template get_dof_data<INPUT>(local_dof);
            });

            dm.template get_locally_owned_data<OUTPUT>(op_x);
        }

        template <class F>
        void assemble_rhs(F f, RhsVector &rhs) {
            auto fe_values = this->values();
            auto dm = fe_values.dm();

            fe_values.template assemble_local_rhs<F, RHS>(f);

            scatter_add_ghost_data<DM, RHS>(dm);

            dm.template get_locally_owned_data<RHS>(rhs);
        }
    };

}  // namespace mars

#endif
#endif