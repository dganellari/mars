#ifndef GENERATION_MARS_DISTRIBUTED_DM_HPP_
#define GENERATION_MARS_DISTRIBUTED_DM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_base_data_management.hpp"
#include "mars_distributed_finite_element.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <typename DHandler, typename... T>
    class DM : public BDM<T...> {
    public:
        static constexpr Integer Degree = DHandler::Degree;

        using DofHandler = DHandler;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = BDM<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        MARS_INLINE_FUNCTION
        DM(DofHandler d) : dof_handler(d) {
            SuperDM::template reserve_user_data(
                user_data, "user_data tuple", get_dof_handler().get_local_dof_enum().get_elem_size());
        }

        MARS_INLINE_FUNCTION
        const user_tuple &get_user_data() const { return user_data; }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_dof_data() const {
            return std::get<idx>(user_data);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_dof_data(const Integer i) const {
            return std::get<idx>(user_data)(i);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_data(const Integer i) const {
            return std::get<idx>(user_data)(i);
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx>
        void gather_ghost_data() {
            // gather operation: fill the data from the received ghost data
            SuperDM::template gather_ghost_data<dataidx...>(get_dof_handler(), user_data);
        }

        template <Integer... dataidx>
        user_tuple scatter_ghost_data() {
            return SuperDM::template scatter_ghost_data<dataidx...>(get_dof_handler(), user_data);
        }

        template <Integer... dataidx>
        void scatter_add_ghost_data() {
            // scatter the data to the procs and keep them in a boundary data tuple
            // again if no template argument is specified all the data is scattered.
            // if not all of them then be careful since the tuple is not initialized on
            // the others example: dm_tuple boundary_data =
            // dm.scatter_ghost_data<1>();
            auto boundary_data = scatter_ghost_data<dataidx...>();

            // use the scattered data "boundary_data" to do ops like max, add or min in
            // the dof contributions. Otherwise you can use predifined features like
            // scatter_add as following. careful to use the same template argument for
            // specifing the data as in the scatter_ghost_data since otherwise you might
            // try to access uninitialized tuplelement and get seg faults. example::
            // dm.scatter_add<1>(boundary_data); If: dm.scatter_add<0>(boundary_data) then
            // seg faults.
            SuperDM::template scatter_add<dataidx...>(get_dof_handler(), boundary_data, user_data);
            /*dm.scatter_max<u>(boundary_data);*/
            /*dm.scatter_min<u>(boundary_data);*/
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void get_locally_owned_data(const ViewVectorType<H> &x) {
            using namespace Kokkos;

            assert(get_dof_handler().get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = get_dof_handler().get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_dof_handler().get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_dof_handler().get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);
                    x(i) = dof_data(local);
                });
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void set_locally_owned_data(const ViewVectorType<H> &x) {
            using namespace Kokkos;

            assert(get_dof_handler().get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = get_dof_handler().get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_dof_handler().get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_dof_handler().get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);
                    dof_data(local) = x(i);
                });
        }

        /* building the stencil is the responsibility of the specialized DM. */
        template <typename ST, bool Orient = false>
        ST build_stencil() {
            return mars::build_stencil<ST, Orient>(get_dof_handler());
        }

        /* building the FE dof map*/
        auto build_fe_dof_map() { return mars::build_fe_dof_map(get_dof_handler()); }

        MARS_INLINE_FUNCTION
        const DofHandler &get_dof_handler() const { return dof_handler; }

        template <Integer idx, typename F>
        void owned_data_iterate(F f) const {
            SuperDM::template owned_data_iterate<idx>(*this, f);
        }

        template <Integer idx, typename F>
        void data_iterate(F f) const {
            SuperDM::template data_iterate<idx>(*this, f);
        }

    private:
        DofHandler dof_handler;
        // data associated to the dof data.
        user_tuple user_data;
    };

}  // namespace mars

#endif
#endif

#endif
