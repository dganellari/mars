#ifndef GENERATION_MARS_DISTRIBUTED_DM_HPP_
#define GENERATION_MARS_DISTRIBUTED_DM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class DM : public DOFM<Mesh, degree> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = DOFM<Mesh, degree>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        DM(Mesh *mesh, const context &c) : SuperDM(mesh, c) {}

        template <Integer... dataidx>
        MARS_INLINE_FUNCTION void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size) {
            expand_tuple<resize_view_functor, dataidx...>(resize_view_functor(view_desc, size), tuple);
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
        MARS_INLINE_FUNCTION void fill_buffer_data(user_tuple &buffer_data, const ViewVectorType<Integer> &boundary) {
            const Integer size = boundary.extent(0);
            expand_tuple<FillBufferDataFunctor<Op>, dataidx...>(
                FillBufferDataFunctor<Op>("fill_buffer_data", size, boundary, SuperDM::get_local_dof_enum().get_view_sfc_to_local()),
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
                FillUserDataFunctor<Op>("fill_user_data", size, ghost_sfc, SuperDM::get_local_dof_enum().get_view_sfc_to_local()),
                ghost_user_data,
                user_data);
        }

        template <Integer... dataidx>
        void gather_ghost_data(const context &context) {
            using namespace Kokkos;

            Kokkos::Timer timer;

            // exchange the ghost dofs first since it will be used to find the address
            // of the userdata based on the sfc code.

            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = SuperDM::get_view_scan_recv_mirror()(size);
            user_tuple ghost_user_data;
            reserve_user_data<dataidx...>(ghost_user_data, "ghost_user_data", ghost_size);

            // prepare the buffer to send the boundary data
            const Integer buffer_size = SuperDM::get_boundary_dofs().extent(0);
            user_tuple buffer_data;
            reserve_user_data<dataidx...>(buffer_data, "buffer_data", buffer_size);

            fill_buffer_data<0, dataidx...>(buffer_data, SuperDM::get_boundary_dofs());

            expand_tuple<ExchangeGhostDofsData, dataidx...>(
                ExchangeGhostDofsData(context, SuperDM::get_view_scan_recv_mirror().data(), SuperDM::get_view_scan_send_mirror().data()),
                ghost_user_data,
                buffer_data);

            // use the received ghost data and the sfc to put them to the unified local data
            fill_user_data<0, dataidx...>(ghost_user_data, SuperDM::get_ghost_dofs());

            /* print_nth_tuple<1>(proc_num); */
        }

        template <Integer... dataidx>
        user_tuple scatter_ghost_data(const context &context) {
            using namespace Kokkos;

            Kokkos::Timer timer;

            // exchange the ghost dofs first since it will be used to find the address
            // of the userdata based on the sfc code.

            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = SuperDM::get_view_scan_recv_mirror()(size);
            user_tuple ghost_buffer_data;
            reserve_user_data<dataidx...>(ghost_buffer_data, "ghost_user_data", ghost_size);

            fill_user_data<1, dataidx...>(ghost_buffer_data, SuperDM::get_ghost_dofs());

            const Integer boundary_size = SuperDM::get_boundary_dofs().extent(0);
            user_tuple boundary_user_data;
            reserve_user_data<dataidx...>(boundary_user_data, "boundary_user_data", boundary_size);

            // prepare the buffer to send the boundary data
            expand_tuple<ExchangeGhostDofsData, dataidx...>(
                ExchangeGhostDofsData(context, SuperDM::get_view_scan_send_mirror().data(), SuperDM::get_view_scan_recv_mirror().data()),
                boundary_user_data,
                ghost_buffer_data);
            /* print_nth_tuple<1>(proc_num); */

            return boundary_user_data;
        }

        template <Integer... dataidx>
        void scatter_add(user_tuple &boundary_user_data) {
            fill_buffer_data<1, dataidx...>(boundary_user_data, SuperDM::get_boundary_dofs());
        }

        template <Integer... dataidx>
        void scatter_max(user_tuple &boundary_user_data) {
            fill_buffer_data<2, dataidx...>(boundary_user_data, SuperDM::get_boundary_dofs());
        }

        template <Integer... dataidx>
        void scatter_min(user_tuple &boundary_user_data) {
            fill_buffer_data<3, dataidx...>(boundary_user_data, SuperDM::get_boundary_dofs());
        }

        virtual void enumerate_dofs(const context &context) override {
            SuperDM::enumerate_dofs(context);
            // reserve the user data tuple with the size as the local dof size.
            reserve_user_data(user_data, "user_data tuple", SuperDM::get_local_dof_enum().get_elem_size());
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void get_locally_owned_data(const ViewVectorType<H> &x) {
            using namespace Kokkos;

            assert(SuperDM::get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = SuperDM::get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = SuperDM::get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = SuperDM::get_local_dof_enum().get_view_sfc_to_local();
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

            assert(SuperDM::get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = SuperDM::get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = SuperDM::get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = SuperDM::get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);
                    dof_data(local) = x(i);
                });
        }

        template <Integer idx,
                  Integer face_nr = -1,
                  typename F,
                  typename H = typename std::tuple_element<idx, tuple>::type>
        void boundary_dof_iterate(F f) {
            using namespace Kokkos;
            constexpr Integer Type = simplex_type::ElemType;

            const Integer size = SuperDM::get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = SuperDM::get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = SuperDM::get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            const Integer xdim = SuperDM::get_local_dof_enum().get_XDim();
            const Integer ydim = SuperDM::get_local_dof_enum().get_YDim();
            const Integer zdim = SuperDM::get_local_dof_enum().get_ZDim();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);

                    if (is_boundary_sfc<Type, face_nr>(sfc, xdim, ydim, zdim)) {
                        f(local, dof_data(local));
                    }
                });
        }

    private:
       // data associated to the dof data.
        user_tuple user_data;
    };

}  // namespace mars

#endif
#endif

#endif
