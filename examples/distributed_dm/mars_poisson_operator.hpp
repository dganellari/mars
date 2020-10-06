#ifndef MARS_POISSON_OPERATOR_HPP
#define MARS_POISSON_OPERATOR_HPP

#include "mars_matrix_free_operator.hpp"

namespace mars {

    template<class DM, Integer INPUT, Integer OUTPUT>
    class PoissonOperator final : public Operator<DM> {
    public:

        template <Integer idx>
        using DMDataType = typename DM::template UserDataType<idx>;

        PoissonOperator(const context &c, DM &dm) : Operator<DM>(c, dm) {}

        void init() override {
            this->values().init();
        }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {

            /* auto context = this->ctx(); */
            auto fe_values = this->values();
            auto dm = fe_values.dm();

            dm.template set_locally_owned_data<INPUT>(x);

            /* Kokkos::parallel_for(
                "printglobaldatavalues", locall_owned_dof_size,
                MARS_LAMBDA(const Integer i) {
                  printf("i: %li, gdata: %lf - rank: %i\n", i, x(i), proc_num);
                }); */

            // specify the tuple indices of the tuplelements that are needed to gather.
            // if no index specified it gathers all views of the tuple. All data.
            dm.template gather_ghost_data<INPUT>(this->ctx());

            fe_values.template form_operator<INPUT, OUTPUT>();

            // iterate through the local dofs and print the local number and the data
            /* dm.dof_iterate(
                MARS_LAMBDA(const Integer i) {
                    printf("lid: %li, u: %lf, v: %lf, rank: %i\n", i,
                           dm.get_dof_data<u>(i), dm.get_dof_data<v>(i), proc_num);
                }); */
            /* using dm_tuple = typename DM::user_tuple;
            dm_tuple boundary_data = dm.template scatter_ghost_data<OUTPUT>(this->ctx());

            dm.template scatter_add<OUTPUT>(boundary_data); */

            scatter_add_ghost_data<DM, OUTPUT>(dm, this->ctx());

            // cleanup
            /* dm.boundary_dof_iterate<INPUT>(MARS_LAMBDA(const Integer local_dof, INPUTDT &value) { */
            dm.template boundary_dof_iterate<INPUT>(
                    MARS_LAMBDA(const Integer local_dof, DMDataType<INPUT> &value) {
                dm.template get_dof_data<OUTPUT>(local_dof) = value;
            });

            dm.template get_locally_owned_data<OUTPUT>(op_x);
        }


        template <class F, typename T>
        void assemble_rhs(ViewVectorType<T> &rhs, F f) {
            auto fe_values = this->values();
            fe_values.template assemble_rhs<F>(rhs, f);
        }

    };

    class CopyOperator final {
    public:
        static void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) { Kokkos::deep_copy(op_x, x); }
    };

}  // namespace mars

#endif
