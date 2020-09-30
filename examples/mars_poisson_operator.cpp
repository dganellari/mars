/* #include "mars_poisson_operator.hpp"
#include "mars_precon_conjugate_grad.hpp"

namespace mars {

    template<class DM, Integer INPUT, Integer OUTPUT>
    void PoissonOperator<DM, INPUT, OUTPUT>::apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
        // from x to DM<INPUT> (global to local)
        // gather ghost data
        // form_operator(dm, ...);
        // scatter ghost data
        // apply identity operator for BC (constrained nodes)
        // from DM<OUTPUT> to op_x (local to global)

        [>using INPUTDT = typename DM::UserDataType<INPUT>;<]
        [>using OUTPUTDT = typename DM::UserDataType<OUTPUT>;<]

        auto context = this->ctx();
        auto fe_values = this->values();
        auto dm = fe_values.dm();

        dm.set_locally_owned_data<INPUT>(x);

        [>Kokkos::parallel_for(
            "printglobaldatavalues", locall_owned_dof_size,
            MARS_LAMBDA(const Integer i) {
              printf("i: %li, gdata: %lf - rank: %i\n", i, x(i), proc_num);
            });<]

        // specify the tuple indices of the tuplelements that are needed to gather.
        // if no index specified it gathers all views of the tuple. All data.
        dm.gather_ghost_data<INPUT>(context);

        fe_values.form_operator(dm);

        // iterate through the local dofs and print the local number and the data
        [>dm.dof_iterate(
            MARS_LAMBDA(const Integer i) {
                printf("lid: %li, u: %lf, v: %lf, rank: %i\n", i,
                       dm.get_dof_data<u>(i), dm.get_dof_data<v>(i), proc_num);
            });<]

        [>scatter_add_ghost_data<DM, OUTPUT>(dm, context);<]

        // cleanup
        [>dm.boundary_dof_iterate<INPUT>(MARS_LAMBDA(const Integer local_dof, INPUTDT &value) {<]
        dm.boundary_dof_iterate<INPUT>(MARS_LAMBDA(const Integer local_dof, double &value) {
            dm.get_elem_data<OUTPUT>(local_dof) = value;
        });

        dm.get_locally_owned_data<OUTPUT>(op_x);
    }
}  // namespace mars */
