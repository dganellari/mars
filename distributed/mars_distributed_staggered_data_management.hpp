#ifndef GENERATION_MARS_DISTRIBUTED_SDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_SDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_base_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <Integer Label, class Mesh, Integer degree, typename... T>
    class SDM : public BDM<T...> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = BDM<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        static constexpr Integer dofLabel = Label;
        static constexpr Integer Degree = degree;

        MARS_INLINE_FUNCTION
        SDM(SDofHandler<Label, Mesh, degree> d) : dof_handler(d) {
            SuperDM::template reserve_user_data(vdata, "separated_user_data tuple", dof_handler.get_local_dofs().extent(0));
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_dof_data() const {
            return std::get<idx>(vdata);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_dof_data(const Integer local_dof) const {
            return std::get<idx>(vdata)(get_dof_index(local_dof));
        }

        MARS_INLINE_FUNCTION
        const user_tuple &get_data() const { return vdata; }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_data(const Integer i) const {
            return std::get<idx>(vdata)(i);
        }

        template <Integer idx, typename F>
        void owned_data_iterate(F f) const {
            Kokkos::parallel_for(
                "owned_separated_dof_iter", get_owned_dofs().extent(0), MARS_LAMBDA(const Integer i) {
                    const Integer local_dof = get_owned_dof(i);
                    f(get_dof_data<idx>(local_dof));
                });
        }

        template <Integer idx, typename F>
        void data_iterate(F f) const {
            Kokkos::parallel_for(
                "owned_separated_dof_iter", get_dof_size(), MARS_LAMBDA(const Integer i) {
                    f(get_data<idx>(i));
                });
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
        auto build_fe_dof_map() {
            return mars::build_fe_dof_map(*this);
        }

       template <Integer... dataidx>
        void scatter_add(user_tuple &boundary_user_data) {
            SuperDM::template fill_buffer_data<1, dataidx...>(
                vdata, boundary_user_data, get_boundary_dofs(), get_local_dof_map());
        }

        template <Integer... dataidx>
        void scatter_max(user_tuple &boundary_user_data) {
            SuperDM::template fill_buffer_data<2, dataidx...>(
                vdata, boundary_user_data, get_boundary_dofs(), get_local_dof_map());
        }

        template <Integer... dataidx>
        void scatter_min(user_tuple &boundary_user_data) {
            SuperDM::template fill_buffer_data<3, dataidx...>(
                vdata, boundary_user_data, get_boundary_dofs(), get_local_dof_map());
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx>
        void gather_ghost_data() {
            const context &context = get_dof_handler().get_context();

            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = get_view_scan_recv_mirror()(size);
            user_tuple ghost_user_data;
            SuperDM::template reserve_user_data<dataidx...>(ghost_user_data, "ghost_user_data", ghost_size);

            // prepare the buffer to send the boundary data
            const Integer buffer_size = get_boundary_dofs().extent(0);
            user_tuple buffer_data;
            SuperDM::template reserve_user_data<dataidx...>(buffer_data, "buffer_data", buffer_size);

            SuperDM::template fill_buffer_data<0, dataidx...>(
                vdata, buffer_data, get_boundary_dofs(), get_local_dof_map());

            SuperDM::template exchange_ghost_dofs_data<dataidx...>(
                context, ghost_user_data, buffer_data, get_view_scan_recv_mirror().data(), get_view_scan_send_mirror().data());

            SuperDM::template fill_user_data<0, dataidx...>(
                vdata, ghost_user_data, get_ghost_dofs(), get_local_dof_map());
        }

        template <Integer... dataidx>
        user_tuple scatter_ghost_data() {
            const context &context = get_dof_handler().get_context();

            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = get_view_scan_recv_mirror()(size);
            user_tuple ghost_buffer_data;
            SuperDM::template reserve_user_data<dataidx...>(ghost_buffer_data, "ghost_user_data", ghost_size);

            SuperDM::template fill_user_data<1, dataidx...>(
                vdata, ghost_buffer_data, get_ghost_dofs(), get_local_dof_map());

            const Integer boundary_size = get_boundary_dofs().extent(0);
            user_tuple boundary_user_data;
            SuperDM::template reserve_user_data<dataidx...>(boundary_user_data, "boundary_user_data", boundary_size);

            // prepare the buffer to send the boundary data
            SuperDM::template exchange_ghost_dofs_data<dataidx...>(context,
                                                                   boundary_user_data,
                                                                   ghost_buffer_data,
                                                                   get_view_scan_send_mirror().data(),
                                                                   get_view_scan_recv_mirror().data());

            return boundary_user_data;
        }

        MARS_INLINE_FUNCTION
        virtual Integer get_dof_size() const override { return get_dof_handler().get_local_dofs().extent(0); }

        MARS_INLINE_FUNCTION
        virtual Integer get_owned_dof_size() const override { return get_dof_handler().get_owned_dofs().extent(0); }

        MARS_INLINE_FUNCTION
        const SDofHandler<Label, Mesh, degree> &get_dof_handler() const { return dof_handler; }

    private:
        SDofHandler<Label, Mesh, degree> dof_handler;
        // data assigned to each separated local dof
        user_tuple vdata;
    };

    template <class Mesh, Integer degree, typename... T>
    using VDM = SDM<DofLabel::lVolume, Mesh, degree, T...>;

    template <class Mesh, Integer degree, typename... T>
    using FDM = SDM<DofLabel::lFace, Mesh, degree, T...>;

    template <class Mesh, Integer degree, typename... T>
    using CDM = SDM<DofLabel::lCorner, Mesh, degree, T...>;

}  // namespace mars

#endif
#endif

#endif
