#ifndef GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_
#define GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_

#ifdef MARS_ENABLE_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_distributed_mesh_management.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    template <class Mesh, typename... T>
    class UserData {
        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using simplex_type = typename Mesh::Elem;
        /* using MM = MeshManager<Mesh>; */
        using SfcKeyType = typename Mesh::SfcKeyType;
        using KeyType = typename Mesh::KeyType;

    public:
        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        UserData(Mesh m) : mesh(m) {
            const Integer size = get_mesh().get_chunk_size();
            reserve_user_data(user_data_, "user_data", size);
        }

        MARS_INLINE_FUNCTION UserData(Mesh m, const user_tuple &data) : mesh(m), user_data_(data) {}

        void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size) {
            apply_impl(resize_view_functor(view_desc, size), tuple);
        }

        template <typename ElementType>
        struct FillBufferData {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<KeyType> boundary_lsfc_index;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<KeyType> bd)
                : buffer_data(bf), user_data(ud), boundary_lsfc_index(bd) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const auto lsfc_index = boundary_lsfc_index(i);
                buffer_data(i) = user_data(lsfc_index);
            }
        };

        struct fill_buffer_data_functor {
            fill_buffer_data_functor(std::string d, size_t s, ViewVectorType<KeyType> b)
                : desc(d), size(s), boundary_lsfc_index(b) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                Kokkos::parallel_for(desc, size, FillBufferData<ElementType>(el_1, el_2, boundary_lsfc_index));
            }

            std::string desc;
            size_t size;

            ViewVectorType<KeyType> boundary_lsfc_index;
        };

        void fill_buffer_data(user_tuple &buffer_data) {
            const Integer size = get_mesh().get_view_boundary().extent(0);

            reserve_user_data(buffer_data, "buffer_data", size);

            auto boundary_lsfc_index = get_mesh().get_view_boundary_sfc_index();
            apply_impl(
                fill_buffer_data_functor("fill_buffer_data", size, boundary_lsfc_index), buffer_data, user_data_);
        }

        struct exchange_ghost_data_functor {
            exchange_ghost_data_functor(const context &c,
                                        ViewVectorType<Integer>::HostMirror sr,
                                        ViewVectorType<Integer>::HostMirror ss)
                : con(c), sc_rcv_mirror(sr), sc_snd_mirror(ss) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                con->distributed->i_send_recv_view(el_1, sc_rcv_mirror.data(), el_2, sc_snd_mirror.data());
            }

            ViewVectorType<Integer>::HostMirror sc_rcv_mirror;
            ViewVectorType<Integer>::HostMirror sc_snd_mirror;

            const context &con;
        };

        void exchange_ghost_data(const context &context) {
            using namespace Kokkos;

            Kokkos::Timer timer;

            // exchange the ghost layer first since it will be used to find the address
            // of the userdata based on the sfc code.
            /* exchange_ghost_layer(context); */

            int proc_num = rank(context);
            int size = num_ranks(context);

            auto scan_send_mirror = get_mesh().get_view_scan_send_mirror();
            auto scan_recv_mirror = get_mesh().get_view_scan_recv_mirror();
            Integer ghost_size = scan_recv_mirror(size);
            reserve_user_data(ghost_user_data_, "ghost_user_data", ghost_size);

            user_tuple buffer_data;
            fill_buffer_data(buffer_data);

            Kokkos::fence();

            apply_impl(exchange_ghost_data_functor(context, scan_recv_mirror, scan_send_mirror),
                       ghost_user_data_,
                       buffer_data);

            /* print_nth_tuple<1>(proc_num); */
        }

        template <Integer I, typename H = typename std::tuple_element<I, tuple>::type>
        void print_nth_tuple(const int proc) {
            using namespace Kokkos;

            ViewVectorType<Integer> scan_ghost = get_view_scan_ghost();
            ViewVectorType<Integer> ghost = get_view_ghost();
            ViewVectorType<H> data = std::get<I>(ghost_user_data_);
            Integer ghost_size = data.extent(0);

            Integer xDim = get_mesh().get_XDim();
            Integer yDim = get_mesh().get_YDim();
            Integer zDim = get_mesh().get_ZDim();

            parallel_for(
                "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                    const Integer r = find_owner_processor(scan_ghost, i, 1, proc);

                    double point[3];
                    get_vertex_coordinates_from_sfc<simplex_type::ElemType, SfcKeyType>(
                        ghost(i), point, xDim, yDim, zDim);

                    Octant o = get_octant_from_sfc<simplex_type::ElemType, SfcKeyType>(ghost(i));
                    printf(
                        "ghost data: %li - %li - %li - (%lf, %lf) - data: %lf - proc: "
                        "%li - rank: %i\n",
                        i,
                        ghost(i),
                        elem_index(o.x, o.y, o.z, xDim, yDim),
                        point[0],
                        point[1],
                        data(i),
                        r,
                        proc);
                });
        }

        template <typename ElementType>
        struct InitialCondition {
            ElementType user_data;
            Integer proc;

            InitialCondition(ElementType ud, Integer p) : user_data(ud), proc(p) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const { user_data(i) = proc; }
        };

        struct InitData {
            InitData(std::string d, size_t s, tuple t) : desc(d), size(s), tup(t) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, std::size_t I) const {
                using ET = typename ElementType::value_type;
                constexpr std::size_t TypeIndex = TypeIdx<tuple, ET>::value;

                Kokkos::parallel_for(
                    desc + std::to_string(I), size, InitialCondition<ElementType>(el_1, std::get<TypeIndex>(tup)));
            }

            std::string desc;
            size_t size;
            tuple tup;
        };

        void init_user_data(T... args) {
            const Integer size = get_mesh().get_chunk_size();

            apply_impl(InitData("init_data", size, std::forward_as_tuple(args...)), user_data_);
        }

        /* template <typename H>
      struct InitCond
      {
          user_tuple tuple;
          H init_cond;

          InitCond(user_tuple t, H f) : tuple(t), init_cond(f) {}

          void operator()(int i) const
          {
              init_cond(tuple, i);
          }
      };
    */

        template <typename H>
        void parallel_for_data(const Integer size, H f) {
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        void set_init_cond(H f) {
            elem_iterate(f);
        }

        template <typename H, typename S>
        void elem_iterate_reduce(H f, S s) {
            const Integer size = get_mesh().get_chunk_size();
            Kokkos::parallel_reduce("elem_reduce", size, f, s);
        }

        template <typename H>
        void elem_iterate(H f) {
            get_mesh().elem_iterate(f);
        }

        template <typename H>
        void face_iterate(H f) const {
            get_mesh().face_iterate(f);
        }

        MARS_INLINE_FUNCTION
        Integer get_ghost_elem(const Integer i) const { return get_mesh().get_ghost_sfc(i); }

        const ViewVectorType<Integer> &get_view_boundary() const {
            return get_mesh().get_view_boundary();
        }

        const ViewVectorType<Integer> &get_view_ghost() const {
            return get_mesh().get_view_ghost();
        }

        const ViewVectorType<Integer> &get_view_scan_ghost() const {
            return get_mesh().get_view_scan_ghost();
        }

        MARS_INLINE_FUNCTION
        const user_tuple &get_user_data() const { return user_data_; }

        /* template<std::size_t idx, typename H = typename NthType<idx, T...>::type>
         */
        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_elem_data(const int i) const {
            return std::get<idx>(user_data_)(i);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION void set_elem_data(const int i, const H value) const {
            std::get<idx>(user_data_)(i) = value;
        }

        /* template<std::size_t idx, typename H = NthType<idx, T...>> */
        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_data() const {
            return std::get<idx>(user_data_);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_ghost_elem_data(const int i) const {
            return std::get<idx>(ghost_user_data_)(i);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_ghost_data() const {
            return std::get<idx>(ghost_user_data_);
        }

        // return the sfc from the sfc index from either ghost or local index
        template <bool Ghost = false>
        MARS_INLINE_FUNCTION Integer get_sfc(const Integer sfc_index) {
            if (Ghost)
                return get_ghost_elem(sfc_index);
            else
                return get_mesh().get_sfc(sfc_index);
        }

        MARS_INLINE_FUNCTION
        Mesh get_mesh() const { return mesh; }

    private:
        /* MM mesh_manager_; */
        Mesh mesh;
        // ghost data layer
        user_tuple user_data_;
        user_tuple ghost_user_data_;
    };

    template <class UserData>
    void exchange_ghost_user_data(UserData &data) {
        const context &context = data.get_mesh().get_context();
        int size = num_ranks(context);

        if (data.get_mesh().get_view_scan_recv_mirror()(size) ==
            data.get_mesh().get_view_ghost().extent(0)) {
            // std::cout << "Exchange the ghost data..." << std::endl;
            data.exchange_ghost_data(context);
        } else {
            errorx(1,
                   "Not allowed to call exchange ghost data before the ghost layer "
                   "creation. Please call create_ghost_layer method first!");
        }
    }

}  // namespace mars

#endif
#endif

#endif
