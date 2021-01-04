#ifndef GENERATION_MARS_DISTRIBUTED_DM_HPP_
#define GENERATION_MARS_DISTRIBUTED_DM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_base_data_management.hpp"
#include "mars_distributed_stencil.hpp"
#include "mars_distributed_finite_element.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class DM : public BDM<Mesh, degree, T...> {
    public:
        static constexpr Integer Degree = degree;

        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = BDM<Mesh, degree, T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        MARS_INLINE_FUNCTION
        DM(DofHandler<Mesh, degree> d) : SuperDM(d) {
            SuperDM::template reserve_user_data(
                user_data, "user_data tuple", SuperDM::get_dof_handler().get_local_dof_enum().get_elem_size());
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

        template <Integer... dataidx>
        void scatter_add(user_tuple &boundary_user_data) {
            SuperDM::template fill_buffer_data<1, dataidx...>(
                user_data,
                boundary_user_data,
                SuperDM::get_dof_handler().get_boundary_dofs(),
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        template <Integer... dataidx>
        void scatter_max(user_tuple &boundary_user_data) {
            SuperDM::template fill_buffer_data<2, dataidx...>(
                user_data,
                boundary_user_data,
                SuperDM::get_dof_handler().get_boundary_dofs(),
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        template <Integer... dataidx>
        void scatter_min(user_tuple &boundary_user_data) {
            SuperDM::template fill_buffer_data<3, dataidx...>(
                user_data,
                boundary_user_data,
                SuperDM::get_dof_handler().get_boundary_dofs(),
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx>
        void gather_ghost_data() {
            const context &context = SuperDM::get_dof_handler().get_context();
            // exchange the ghost dofs first since it will be used to find the address
            // of the userdata based on the sfc code.
            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = SuperDM::get_dof_handler().get_view_scan_recv_mirror()(size);
            user_tuple ghost_user_data;
            SuperDM::template reserve_user_data<dataidx...>(ghost_user_data, "ghost_user_data", ghost_size);

            // prepare the buffer to send the boundary data
            const Integer buffer_size = SuperDM::get_dof_handler().get_boundary_dofs().extent(0);
            user_tuple buffer_data;
            SuperDM::template reserve_user_data<dataidx...>(buffer_data, "buffer_data", buffer_size);

            SuperDM::template fill_buffer_data<0, dataidx...>(
                user_data,
                buffer_data,
                SuperDM::get_dof_handler().get_boundary_dofs(),
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());

            SuperDM::template exchange_ghost_dofs_data<dataidx...>(
                context,
                ghost_user_data,
                buffer_data,
                SuperDM::get_dof_handler().get_view_scan_recv_mirror().data(),
                SuperDM::get_dof_handler().get_view_scan_send_mirror().data());

            // use the received ghost data and the sfc to put them to the unified local data
            SuperDM::template fill_user_data<0, dataidx...>(
                user_data,
                ghost_user_data,
                SuperDM::get_dof_handler().get_ghost_dofs(),
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        template <Integer... dataidx>
        user_tuple scatter_ghost_data() {
            const context &context = SuperDM::get_dof_handler().get_context();
            // exchange the ghost dofs first since it will be used to find the address
            // of the userdata based on the sfc code.
            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = SuperDM::get_dof_handler().get_view_scan_recv_mirror()(size);
            user_tuple ghost_buffer_data;
            SuperDM::template reserve_user_data<dataidx...>(ghost_buffer_data, "ghost_user_data", ghost_size);

            SuperDM::template fill_user_data<1, dataidx...>(
                user_data,
                ghost_buffer_data,
                SuperDM::get_dof_handler().get_ghost_dofs(),
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());

            const Integer boundary_size = SuperDM::get_dof_handler().get_boundary_dofs().extent(0);
            user_tuple boundary_user_data;
            SuperDM::template reserve_user_data<dataidx...>(boundary_user_data, "boundary_user_data", boundary_size);

            // prepare the buffer to send the boundary data
            SuperDM::template exchange_ghost_dofs_data<dataidx...>(
                context,
                boundary_user_data,
                ghost_buffer_data,
                SuperDM::get_dof_handler().get_view_scan_send_mirror().data(),
                SuperDM::get_dof_handler().get_view_scan_recv_mirror().data());

            return boundary_user_data;
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void get_locally_owned_data(const ViewVectorType<H> &x) {
            using namespace Kokkos;

            assert(SuperDM::get_dof_handler().get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = SuperDM::get_dof_handler().get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc =
                SuperDM::get_dof_handler().get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local =
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local();
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

            assert(SuperDM::get_dof_handler().get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = SuperDM::get_dof_handler().get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc =
                SuperDM::get_dof_handler().get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local =
                SuperDM::get_dof_handler().get_local_dof_enum().get_view_sfc_to_local();
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
            return mars::build_stencil<ST, Orient>(SuperDM::get_dof_handler());
        }

        /* building the FE dof map*/
        auto build_fe_dof_map() {
            return mars::build_fe_dof_map(SuperDM::get_dof_handler());
        }

    private:
        // data associated to the dof data.
        user_tuple user_data;
    };

}  // namespace mars

#endif
#endif

#endif
