#ifndef GENERATION_MARS_DISTRIBUTED_DM_HPP_
#define GENERATION_MARS_DISTRIBUTED_DM_HPP_

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
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

        DM(const DofHandler& d) : dof_handler(d), halo_plan(d) {
            SuperDM::template reserve_user_data<>(
                user_data, "user_data tuple", get_dof_handler().get_dof_size());
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
            halo_plan.template start_gather<dataidx...>(user_data);
            halo_plan.wait();
        }

        template <Integer... dataidx>
        user_tuple scatter_ghost_data() {
            halo_plan.template start_scatter<dataidx...>(user_data);
            return halo_plan.wait();
        }

        template <Integer... dataidx>
        void scatter_add_ghost_data() {
            auto boundary_data = scatter_ghost_data<dataidx...>();
            SuperDM::template scatter_add<dataidx...>(get_dof_handler(), boundary_data, user_data);
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void get_locally_owned_data(const ViewVectorType<H> &x) {
            SuperDM::template get_locally_owned_data<idx>(get_dof_handler(), x, user_data);
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void set_locally_owned_data(const ViewVectorType<H> &x) {
            SuperDM::template set_locally_owned_data<idx>(get_dof_handler(), user_data, x);
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

        template<Integer... fields> void start_gather() { halo_plan.template start_gather<fields...>(user_data); }
        template<Integer... fields> void start_scatter() { halo_plan.template start_scatter<fields...>(user_data); }
        user_tuple wait() { return halo_plan.wait(); }
        const HaloExchangePlan<T...>& exchange_plan() const { return halo_plan; }

        struct DeviceView {
            typename DofHandler::DeviceView dofs;
            user_tuple fields;
            MARS_INLINE_FUNCTION auto get_dof_handler() const { return dofs; }
            template<std::size_t i> MARS_INLINE_FUNCTION auto get_dof_data() const { return std::get<i>(fields); }
            template<std::size_t i> MARS_INLINE_FUNCTION auto& get_dof_data(Integer j) const { return std::get<i>(fields)(j); }
            template<std::size_t i> MARS_INLINE_FUNCTION auto& get_data(Integer j) const { return std::get<i>(fields)(j); }
        };
        DeviceView device_view() const { return {dof_handler.device_view(), user_data}; }

    private:
        DofHandler dof_handler;
        HaloExchangePlan<T...> halo_plan;
        // data associated to the dof data.
        user_tuple user_data;
    };

}  // namespace mars

#endif
#endif

#endif
