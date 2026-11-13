#ifndef GENERATION_MARS_DISTRIBUTED_BDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_BDM_HPP_

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_dof_management.hpp"
#include "mars_distributed_halo_plan.hpp"

namespace mars {

    template <typename... T>
    class BDM {
    public:
        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        template <Integer... dataidx>
        static void reserve_user_data(user_tuple &tuple,
                                                           std::string view_desc,
                                                           const Integer size) {
            expand_tuple<resize_view_functor, user_tuple, dataidx...>(resize_view_functor(view_desc, size), tuple);
        }

        template <typename H>
        void parallel_for_data(const Integer size, H f) {
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H, typename ElementType, Integer Op = 0>
        struct FillBufferData {
            ElementType buffer_data;
            ElementType user_data;
            H dof_handler;

            FillBufferData(ElementType bf, ElementType ud, H d) : buffer_data(bf), user_data(ud), dof_handler(d) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                auto block_local = dof_handler.get_boundary_local_index(i);
                buffer_data(i) = user_data(block_local);
            }
        };

        template <typename H, typename ElementType>
        struct FillBufferData<H, ElementType, 1> {
            ElementType buffer_data;
            ElementType user_data;
            H dof_handler;

            FillBufferData(ElementType bf, ElementType ud, H d) : buffer_data(bf), user_data(ud), dof_handler(d) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                auto block_local = dof_handler.get_boundary_local_index(i);
                Kokkos::atomic_add(&user_data(block_local), buffer_data(i));
            }
        };

        template <typename H, typename ElementType>
        struct FillBufferData<H, ElementType, 2> {
            ElementType buffer_data;
            ElementType user_data;
            H dof_handler;

            FillBufferData(ElementType bf, ElementType ud, H d) : buffer_data(bf), user_data(ud), dof_handler(d) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                auto block_local = dof_handler.get_boundary_local_index(i);
                Kokkos::atomic_fetch_max(&user_data(block_local), buffer_data(i));
            }
        };

        template <typename H, typename ElementType>
        struct FillBufferData<H, ElementType, 3> {
            ElementType buffer_data;
            ElementType user_data;
            H dof_handler;

            FillBufferData(ElementType bf, ElementType ud, H d) : buffer_data(bf), user_data(ud), dof_handler(d) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                auto block_local = dof_handler.get_boundary_local_index(i);
                Kokkos::atomic_fetch_min(&user_data(block_local), buffer_data(i));
            }
        };

        template <Integer Op, typename H>
        struct FillBufferDataFunctor {
            FillBufferDataFunctor(std::string d, size_t s, H hd) : desc(d), size(s), dof_handler(hd) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                Kokkos::parallel_for(desc, size, FillBufferData<H, ElementType, Op>(el_1, el_2, dof_handler));
            }

            std::string desc;
            size_t size;
            H dof_handler;
        };

        template <typename H, Integer Op, Integer... dataidx>
        static void fill_buffer_data(user_tuple &udata, user_tuple &buffer_data, const H &dof_handler) {
            const Integer size = dof_handler.get_boundary_dof_size();
            const auto view = dof_handler.device_view();
            using DeviceView = decltype(view);
            expand_tuple<FillBufferDataFunctor<Op, DeviceView>, user_tuple, dataidx...>(
                FillBufferDataFunctor<Op, DeviceView>("fill_buffer_data", size, view), buffer_data, udata);
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx, typename H>
        static void gather_ghost_data(const H &dof_handler, user_tuple &user_data) {
            gather_ghost_data<dataidx...>(dof_handler, user_data, []() {});
        }

        // gather operation: fill the data from the received ghost data
        template <Integer... dataidx, typename H, typename F>
        static void gather_ghost_data(const H &dof_handler, user_tuple &user_data, F callback) {
            HaloExchangePlan<T...> plan(dof_handler, false);
            plan.template start_gather<dataidx...>(user_data);
            try { callback(); } catch (...) { plan.wait(); throw; }
            plan.wait();
        }

        template <Integer... dataidx, typename H>
        static void scatter_add(const H &handler, user_tuple &boundary_user_data, user_tuple &user_data) {
            fill_buffer_data<H, 1, dataidx...>(user_data, boundary_user_data, handler);
        }

        template <Integer... dataidx, typename H>
        void scatter_max(const H &handler, user_tuple &boundary_user_data, user_tuple &user_data) {
            fill_buffer_data<H, 2, dataidx...>(user_data, boundary_user_data, handler);
        }

        template <Integer... dataidx, typename H>
        void scatter_min(const H &handler, user_tuple &boundary_user_data, user_tuple &user_data) {
            fill_buffer_data<H, 3, dataidx...>(user_data, boundary_user_data, handler);
        }

        template <Integer... dataidx, typename H>
        static user_tuple scatter_ghost_data(const H &dof_handler, user_tuple &user_data) {
            return scatter_ghost_data<dataidx...>(dof_handler, user_data, []() {});
        }

        template <Integer... dataidx, typename H, typename F>
        static user_tuple scatter_ghost_data(const H &dof_handler, user_tuple &user_data, F callback) {
            HaloExchangePlan<T...> plan(dof_handler, false);
            plan.template start_scatter<dataidx...>(user_data);
            try { callback(); } catch (...) { plan.wait(); throw; }
            return plan.wait();
        }

        template <Integer idx, typename DM, typename F>
        static void owned_data_iterate(const DM &dm, F f) {
            auto handler = dm.get_dof_handler().device_view();
            auto values = dm.template get_dof_data<idx>();
            Kokkos::parallel_for(
                "owned_separated_dof_iter", handler.get_owned_dof_size(), MARS_LAMBDA(const Integer i) {
                    const Integer local_dof = handler.get_owned_dof(i);
                    f(local_dof, values(local_dof));
                });
        }

        template <Integer idx, typename DM, typename F>
        static void data_iterate(const DM &dm, F f) {
            auto handler = dm.get_dof_handler().device_view();
            auto values = dm.template get_dof_data<idx>();
            Kokkos::parallel_for(
                "owned_separated_dof_iter", handler.get_dof_size(), MARS_LAMBDA(const Integer i) {
                    f(values(i));
                });
        }

        template <Integer dataidx, typename H, typename DH>
        static void set_locally_owned_data(const DH &dhandler, user_tuple &dof_data, const ViewVectorType<H> &x) {
            assert(dhandler.get_owned_dof_size() == x.extent(0));

            const auto handler = dhandler.device_view();
            const auto values = std::get<dataidx>(dof_data);
            dhandler.owned_dof_iterate(MARS_LAMBDA(const Integer i) {
                const auto block_local = handler.get_owned_dof(i);
                values(block_local) = x(i);
            });
        }

        template <Integer dataidx, typename H, typename DH>
        static void get_locally_owned_data(const DH &dhandler, const ViewVectorType<H> &x, user_tuple &dof_data) {
            assert(dhandler.get_owned_dof_size() == x.extent(0));

            const auto handler = dhandler.device_view();
            const auto values = std::get<dataidx>(dof_data);
            dhandler.owned_dof_iterate(MARS_LAMBDA(const Integer i) {
                const auto block_local = handler.get_owned_dof(i);
                x(i) = values(block_local);
            });
        }
    };

    // gather operation: fill the data from the received ghost data
    template <typename H, typename T, typename F>
    void gather_ghost_data(const H &dof_handler, ViewVectorType<T> &data, F callback) {
        assert(data.extent(0) == dof_handler.get_dof_size());
        using SuperDM = BDM<T>;
        auto tuple = std::make_tuple(data);
        SuperDM::template gather_ghost_data<0>(dof_handler, tuple, callback);
    }

    // gather operation: fill the data from the received ghost data
    template <typename H, typename T>
    void gather_ghost_data(const H &dof_handler, ViewVectorType<T> &data) {
        gather_ghost_data(dof_handler, data, []() {});
    }

    // scatter add but using the handler and the view instead of a dm.
    template <typename H, typename T, typename F>
    void scatter_add_ghost_data(const H &dof_handler, ViewVectorType<T> &data, F callback) {
        assert(data.extent(0) == dof_handler.get_dof_size());
        using SuperDM = BDM<T>;
        auto tuple = std::make_tuple(data);
        auto boundary_data = SuperDM::template scatter_ghost_data<0>(dof_handler, tuple, callback);
        SuperDM::template scatter_add<0>(dof_handler, boundary_data, tuple);
    }

    // scatter add but using the handler and the view instead of a dm.
    template <typename H, typename T>
    void scatter_add_ghost_data(const H &dof_handler, ViewVectorType<T> &data) {
        scatter_add_ghost_data(dof_handler, data, []() {});
    }

    template <class DM, Integer... dataidx>
    void scatter_add_ghost_data(DM &dm) {
        dm.template scatter_add_ghost_data<dataidx...>();
    }

    template <typename H, typename T>
    void set_locally_owned_data(const H &dof_handler, ViewVectorType<T> &data, ViewVectorType<T> owned) {
        assert(data.extent(0) == dof_handler.get_dof_size());
        using SuperDM = BDM<T>;
        auto tuple = std::make_tuple(data);
        SuperDM::template set_locally_owned_data<0>(dof_handler, tuple, owned);
    }

    template <typename H, typename T>
    void get_locally_owned_data(const H &dof_handler, ViewVectorType<T> owned, ViewVectorType<T> &data) {
        assert(data.extent(0) == dof_handler.get_dof_size());
        using SuperDM = BDM<T>;
        auto tuple = std::make_tuple(data);
        SuperDM::template get_locally_owned_data<0>(dof_handler, owned, tuple);
    }

}  // namespace mars

#endif
#endif

#endif
