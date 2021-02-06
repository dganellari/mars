#ifndef GENERATION_MARS_DISTRIBUTED_SDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_SDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_base_data_management.hpp"
#include "mars_distributed_staggered_dof_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <typename DHandler, typename... T>
    class SDM : public BDM<T...> {
    public:
        static constexpr Integer Degree = DHandler::Degree;

        using SDofHandler = DHandler;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = BDM<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;


        MARS_INLINE_FUNCTION
        SDM(DofHandler<typename SDofHandler::Mesh, SDofHandler::Degree> d) : dof_handler(SDofHandler(d)) {
            SuperDM::template reserve_user_data(
                vdata, "separated_user_data tuple", dof_handler.get_local_dofs().extent(0));
        }

        MARS_INLINE_FUNCTION
        SDM(SDofHandler d) : dof_handler(d) {
            SuperDM::template reserve_user_data(
                vdata, "separated_user_data tuple", dof_handler.get_local_dofs().extent(0));
        }

        MARS_INLINE_FUNCTION
        const user_tuple &get_user_data() const { return vdata; }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_dof_data() const {
            return std::get<idx>(vdata);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_dof_data(const Integer local_dof) const {
            return std::get<idx>(vdata)(get_dof_handler().get_dof_index(local_dof));
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_data(const Integer i) const {
            return std::get<idx>(vdata)(i);
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx>
        void gather_ghost_data() {
            // gather operation: fill the data from the received ghost data
            SuperDM::template gather_ghost_data<dataidx...>(get_dof_handler(), vdata);
        }

        template <Integer... dataidx>
        user_tuple scatter_ghost_data() {
            return SuperDM::template scatter_ghost_data<dataidx...>(get_dof_handler(), vdata);
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
            SuperDM::template scatter_add<dataidx...>(get_dof_handler(), boundary_data, vdata);
            /*dm.scatter_max<u>(boundary_data);*/
            /*dm.scatter_min<u>(boundary_data);*/
        }

        /* Face numbering on the stencil => ordering in the stencil stencil[1,0,3,2]
                       ----3----
                       |       |
                       0   x   1
                       |       |
                       ----2---- */
        // building the stencil is the responsibility of the specialized DM.
        template <typename ST, bool Orient = false>
        ST build_stencil() {
            return mars::build_stencil<ST, Orient>(*this);
        }

        /* building the FE dof map*/
        auto build_fe_dof_map() { return mars::build_fe_dof_map(*this); }

        MARS_INLINE_FUNCTION
        const SDofHandler &get_dof_handler() const { return dof_handler; }

        template <Integer idx, typename F>
        void owned_data_iterate(F f) const {
            SuperDM::template owned_data_iterate<idx>(*this, f);
        }

        template <Integer idx, typename F>
        void data_iterate(F f) const {
            SuperDM::template data_iterate<idx>(*this, f);
        }

    private:
        SDofHandler dof_handler;
        // data assigned to each separated local dof
        user_tuple vdata;
    };

    template <class DofHandler, typename... T>
    using VDM = SDM<VolumeDofHandler<DofHandler>, T...>;

    template <class DofHandler, typename... T>
    using FDM = SDM<FaceDofHandler<DofHandler>, T...>;

    template <class DofHandler, typename... T>
    using CDM = SDM<CornerDofHandler<DofHandler>, T...>;

}  // namespace mars

#endif
#endif

#endif
