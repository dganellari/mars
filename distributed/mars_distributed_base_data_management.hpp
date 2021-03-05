#ifndef GENERATION_MARS_DISTRIBUTED_BDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_BDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    template <typename... T>
    class BDM {
    public:
        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        template <Integer... dataidx>
        MARS_INLINE_FUNCTION static void reserve_user_data(user_tuple &tuple,
                                                           std::string view_desc,
                                                           const Integer size) {
            expand_tuple<resize_view_functor, user_tuple, dataidx...>(resize_view_functor(view_desc, size), tuple);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void parallel_for_data(const Integer size, H f) {
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename ElementType, Integer Op = 0>
        struct FillBufferData {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                buffer_data(i) = user_data(sfc_index);
            }
        };

        template <typename ElementType>
        struct FillBufferData<ElementType, 1> {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                Kokkos::atomic_add(&user_data(sfc_index), buffer_data(i));
            }
        };

        template <typename ElementType>
        struct FillBufferData<ElementType, 2> {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                Kokkos::atomic_fetch_max(&user_data(sfc_index), buffer_data(i));
            }
        };

        template <typename ElementType>
        struct FillBufferData<ElementType, 3> {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                Kokkos::atomic_fetch_min(&user_data(sfc_index), buffer_data(i));
            }
        };

        template <Integer Op>
        struct FillBufferDataFunctor {
            FillBufferDataFunctor(std::string d, size_t s, ViewVectorType<Integer> b, ViewVectorType<Integer> stl)
                : desc(d), size(s), boundary_sfc(b), sfc_to_local(stl) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                Kokkos::parallel_for(
                    desc, size, FillBufferData<ElementType, Op>(el_1, el_2, boundary_sfc, sfc_to_local));
            }

            std::string desc;
            size_t size;

            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;
        };

        template <Integer Op, Integer... dataidx>
        MARS_INLINE_FUNCTION static void fill_buffer_data(user_tuple &udata,
                                                          user_tuple &buffer_data,
                                                          const ViewVectorType<Integer> &boundary,
                                                          const ViewVectorType<Integer> &map) {
            const Integer size = boundary.extent(0);
            expand_tuple<FillBufferDataFunctor<Op>, user_tuple, dataidx...>(
                FillBufferDataFunctor<Op>("fill_buffer_data", size, boundary, map), buffer_data, udata);
        }

        struct ExchangeGhostDofsData {
            ExchangeGhostDofsData(const context &c, Integer *sr, Integer *ss)
                : con(c), sc_rcv_mirror(sr), sc_snd_mirror(ss) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                con->distributed->i_send_recv_view(el_1, sc_rcv_mirror, el_2, sc_snd_mirror);
            }

            Integer *sc_rcv_mirror;
            Integer *sc_snd_mirror;
            const context &con;
        };

        template <Integer... dataidx>
        MARS_INLINE_FUNCTION static void exchange_ghost_dofs_data(const context &c,
                                                                  user_tuple &recv_data,
                                                                  user_tuple &send_data,
                                                                  Integer *recv_mirror,
                                                                  Integer *send_mirror) {
            expand_tuple<ExchangeGhostDofsData, user_tuple, dataidx...>(
                ExchangeGhostDofsData(c, recv_mirror, send_mirror), recv_data, send_data);
        }

        // gather operation: fill the data from the received ghost data
        template <typename ElementType, bool Op = 0>
        struct FillUserData {
            ElementType ghost_data;
            ElementType user_data;
            ViewVectorType<Integer> ghost_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillUserData(ElementType gd, ElementType ud, ViewVectorType<Integer> gs, ViewVectorType<Integer> stl)
                : ghost_data(gd), user_data(ud), ghost_sfc(gs), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer local_sfc_index = sfc_to_local(ghost_sfc(i));
                user_data(local_sfc_index) = ghost_data(i);
            }
        };

        // scatter operation: fill the to be sent ghost buffer from data
        template <typename ElementType>
        struct FillUserData<ElementType, 1> {
            ElementType ghost_data;
            ElementType user_data;
            ViewVectorType<Integer> ghost_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillUserData(ElementType gd, ElementType ud, ViewVectorType<Integer> gs, ViewVectorType<Integer> stl)
                : ghost_data(gd), user_data(ud), ghost_sfc(gs), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer local_sfc_index = sfc_to_local(ghost_sfc(i));
                ghost_data(i) = user_data(local_sfc_index);
            }
        };

        template <bool Op>
        struct FillUserDataFunctor {
            FillUserDataFunctor(std::string d, size_t s, ViewVectorType<Integer> b, ViewVectorType<Integer> stl)
                : desc(d), size(s), ghost_sfc(b), sfc_to_local(stl) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                Kokkos::parallel_for(desc, size, FillUserData<ElementType, Op>(el_1, el_2, ghost_sfc, sfc_to_local));
            }

            std::string desc;
            size_t size;

            ViewVectorType<Integer> ghost_sfc;
            ViewVectorType<Integer> sfc_to_local;
        };

        /* fill the user data views from the ghost data that were filled in FillBufferData and received through mpi. */
        template <bool Op, Integer... dataidx>
        MARS_INLINE_FUNCTION static void fill_user_data(user_tuple &udata,
                                                        user_tuple &ghost_user_data,
                                                        const ViewVectorType<Integer> &ghost_sfc,
                                                        const ViewVectorType<Integer> &map) {
            const Integer size = ghost_sfc.extent(0);
            expand_tuple<FillUserDataFunctor<Op>, user_tuple, dataidx...>(
                FillUserDataFunctor<Op>("fill_user_data", size, ghost_sfc, map), ghost_user_data, udata);
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx, typename H>
        static void gather_ghost_data(const H &dof_handler, user_tuple &user_data) {
            const context &context = dof_handler.get_context();
            // exchange the ghost dofs first since it will be used to find the address
            // of the userdata based on the sfc code.
            /* int proc_num = rank(context); */
            int size = num_ranks(context);

            Integer ghost_size = dof_handler.get_view_scan_recv_mirror()(size);
            user_tuple ghost_user_data;
            reserve_user_data<dataidx...>(ghost_user_data, "ghost_user_data", ghost_size);

            /* prepare the buffer to send the boundary data */
            const Integer buffer_size = dof_handler.get_boundary_dof_size();
            user_tuple buffer_data;
            reserve_user_data<dataidx...>(buffer_data, "buffer_data", buffer_size);

            fill_buffer_data<0, dataidx...>(
                user_data, buffer_data, dof_handler.get_boundary_dofs(), dof_handler.get_local_dof_map());

            exchange_ghost_dofs_data<dataidx...>(context,
                                                 ghost_user_data,
                                                 buffer_data,
                                                 dof_handler.get_view_scan_recv_mirror().data(),
                                                 dof_handler.get_view_scan_send_mirror().data());

            /* use the received ghost data and the sfc to put them to the unified local data */
            fill_user_data<0, dataidx...>(
                user_data, ghost_user_data, dof_handler.get_ghost_dofs(), dof_handler.get_local_dof_map());
        }

        template <Integer... dataidx, typename H>
        static void scatter_add(const H &handler, user_tuple &boundary_user_data, user_tuple &user_data) {
            fill_buffer_data<1, dataidx...>(
                user_data, boundary_user_data, handler.get_boundary_dofs(), handler.get_local_dof_map());
        }

        template <Integer... dataidx, typename H>
        void scatter_max(const H &handler, user_tuple &boundary_user_data, user_tuple &user_data) {
            fill_buffer_data<2, dataidx...>(
                user_data, boundary_user_data, handler.get_boundary_dofs(), handler.get_local_dof_map());
        }

        template <Integer... dataidx, typename H>
        void scatter_min(const H &handler, user_tuple &boundary_user_data, user_tuple &user_data) {
            fill_buffer_data<3, dataidx...>(
                user_data, boundary_user_data, handler.get_boundary_dofs(), handler.get_local_dof_map());
        }

        template <Integer... dataidx, typename H>
        static user_tuple scatter_ghost_data(const H &dof_handler, user_tuple &user_data) {
            const context &context = dof_handler.get_context();
            /* int proc_num = rank(context); */
            int size = num_ranks(context);

            Integer ghost_size = dof_handler.get_view_scan_recv_mirror()(size);
            user_tuple ghost_buffer_data;
            reserve_user_data<dataidx...>(ghost_buffer_data, "ghost_user_data", ghost_size);

            fill_user_data<1, dataidx...>(
                user_data, ghost_buffer_data, dof_handler.get_ghost_dofs(), dof_handler.get_local_dof_map());

            const Integer boundary_size = dof_handler.get_boundary_dofs().extent(0);
            user_tuple boundary_user_data;
            reserve_user_data<dataidx...>(boundary_user_data, "boundary_user_data", boundary_size);

            // prepare the buffer to send the boundary data
            exchange_ghost_dofs_data<dataidx...>(context,
                                                 boundary_user_data,
                                                 ghost_buffer_data,
                                                 dof_handler.get_view_scan_send_mirror().data(),
                                                 dof_handler.get_view_scan_recv_mirror().data());

            return boundary_user_data;
        }

        template <Integer idx, typename DM, typename F>
        static void owned_data_iterate(const DM &dm, F f) {
            auto handler = dm.get_dof_handler();
            Kokkos::parallel_for(
                "owned_separated_dof_iter", handler.get_owned_dof_size(), MARS_LAMBDA(const Integer i) {
                    const Integer local_dof = handler.get_owned_dof(i);
                    f(local_dof, dm.template get_dof_data<idx>(local_dof));
                });
        }

        template <Integer idx, typename DM, typename F>
        void data_iterate(const DM &dm, F f) {
            auto handler = dm.get_dof_handler();
            Kokkos::parallel_for(
                "owned_separated_dof_iter", handler.get_dof_size(), MARS_LAMBDA(const Integer i) {
                    f(dm.get_data<idx>(i));
                });
        }
    };

    // gather operation: fill the data from the received ghost data
    template <typename H, typename T>
    void gather_ghost_data(const H &dof_handler, ViewVectorType<T> &data) {
        assert(data.extent(0) == dof_handler.get_local_dof_enum().get_elem_size());
        using SuperDM = BDM<T>;
        auto tuple = std::make_tuple(data);
        SuperDM::template gather_ghost_data<0>(dof_handler, tuple);
    }

    // scatter add but using the handler and the view instead of a dm.
    template <typename H, typename T>
    void scatter_add_ghost_data(const H &dof_handler, ViewVectorType<T> &data) {
        assert(data.extent(0) == dof_handler.get_local_dof_enum().get_elem_size());
        using SuperDM = BDM<T>;
        auto tuple = std::make_tuple(data);
        auto boundary_data = SuperDM::template scatter_ghost_data<0>(dof_handler, tuple);
        SuperDM::template scatter_add<0>(dof_handler, boundary_data, tuple);
    }

    template <class DM, Integer... dataidx>
    void scatter_add_ghost_data(DM &dm) {
        dm.template scatter_add_ghost_data<dataidx...>();
    }

}  // namespace mars

#endif
#endif

#endif
