#ifndef GENERATION_MARS_DISTRIBUTED_BDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_BDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class BDM {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        BDM(DofHandler<Mesh, degree> d) : dof_handler(d) {}

        template <Integer... dataidx>
        MARS_INLINE_FUNCTION void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size) {
            expand_tuple<resize_view_functor, dataidx...>(resize_view_functor(view_desc, size), tuple);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void parallel_for_data(const Integer size, H f) {
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename F, Integer... dataidx>
        MARS_INLINE_FUNCTION void expand_tuple(const F &f, user_tuple &t) {
            if (sizeof...(dataidx) == 0) {
                apply_impl(f, t);
            } else {
                for_each_arg<F, 0, dataidx...>(f, t);
            }
        }

        template <typename F, Integer... dataidx>
        MARS_INLINE_FUNCTION void expand_tuple(const F &f, user_tuple &t, user_tuple &v) {
            if (sizeof...(dataidx) == 0) {
                apply_impl(f, t, v);
            } else {
                for_each_arg<F, 0, dataidx...>(f, t, v);
            }
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
        MARS_INLINE_FUNCTION void fill_buffer_data(user_tuple &buffer_data,
                                                   const ViewVectorType<Integer> &boundary,
                                                   const ViewVectorType<Integer> map) {
            const Integer size = boundary.extent(0);
            expand_tuple<FillBufferDataFunctor<Op>, dataidx...>(
                FillBufferDataFunctor<Op>("fill_buffer_data", size, boundary, map,
                buffer_data,
                user_data);
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
        MARS_INLINE_FUNCTION void fill_user_data(user_tuple &ghost_user_data,
                                                 const ViewVectorType<Integer> &ghost_sfc) {
            const Integer size = ghost_sfc.extent(0);
            expand_tuple<FillUserDataFunctor<Op>, dataidx...>(
                FillUserDataFunctor<Op>("fill_user_data", size, ghost_sfc, get_dof_handler().get_local_dof_enum().get_view_sfc_to_local()),
                ghost_user_data,
                user_data);
        }

        template <Integer... dataidx>
        void scatter_add(user_tuple &boundary_user_data) {
            fill_buffer_data<1, dataidx...>(boundary_user_data,
                                            get_dof_handler().get_boundary_dofs(),
                                            get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        template <Integer... dataidx>
        void scatter_max(user_tuple &boundary_user_data) {
            fill_buffer_data<2, dataidx...>(boundary_user_data,
                                            get_dof_handler().get_boundary_dofs(),
                                            get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        template <Integer... dataidx>
        void scatter_min(user_tuple &boundary_user_data) {
            fill_buffer_data<3, dataidx...>(boundary_user_data, get_dof_handler().get_boundary_dofs(), get_dof_handler().get_local_dof_enum().get_view_sfc_to_local());
        }

        template <Integer idx,
                  Integer face_nr = -1,
                  typename F,
                  typename H = typename std::tuple_element<idx, tuple>::type>
        void boundary_dof_iterate(F f) {
            using namespace Kokkos;
            constexpr Integer Type = simplex_type::ElemType;

            const Integer size = get_dof_handler().get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_dof_handler().get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_dof_handler().get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            const Integer xdim = get_dof_handler().get_local_dof_enum().get_XDim();
            const Integer ydim = get_dof_handler().get_local_dof_enum().get_YDim();
            const Integer zdim = get_dof_handler().get_local_dof_enum().get_ZDim();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);

                    if (is_boundary_sfc<Type, face_nr>(sfc, xdim, ydim, zdim)) {
                        f(local, dof_data(local));
                    }
                });
        }

        const DofHandler<Mesh, degree> &get_dof_handler() { return dof_handler; }

    private:
        DofHandler<Mesh, degree> dof_handler;
    };

}  // namespace mars

#endif
#endif

#endif
