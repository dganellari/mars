#ifndef GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_
#define GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_

#ifdef MARS_ENABLE_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_sparse_exchange.hpp"
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

        UserData(Mesh m, const user_tuple &data) : mesh(m), user_data_(data) {}

        void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size) {
            apply_impl(resize_view_functor(view_desc, size), tuple);
        }

        void exchange_ghost_data(const context& context) {
            exchange_fields(context, std::index_sequence_for<T...>{});
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

        ViewVectorType<KeyType> get_view_boundary() const {
            return mesh.get_view_boundary();
        }

        ViewVectorType<KeyType> get_view_ghost() const {
            return mesh.get_view_ghost();
        }

        ViewVectorType<Integer> get_view_scan_ghost() const {
            return mesh.get_view_scan_ghost();
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

        struct DeviceView {
            typename Mesh::DeviceView mesh;
            user_tuple local, ghosts;
            MARS_INLINE_FUNCTION auto get_mesh() const { return mesh; }
            template<std::size_t i> MARS_INLINE_FUNCTION auto get_data() const { return std::get<i>(local); }
            template<std::size_t i> MARS_INLINE_FUNCTION auto get_ghost_data() const { return std::get<i>(ghosts); }
            template<std::size_t i> MARS_INLINE_FUNCTION auto& get_elem_data(Integer j) const { return std::get<i>(local)(j); }
            template<std::size_t i> MARS_INLINE_FUNCTION auto& get_ghost_elem_data(Integer j) const { return std::get<i>(ghosts)(j); }
        };
        DeviceView device_view() const { return {mesh.device_view(), user_data_, ghost_user_data_}; }

    private:
        template<std::size_t I>
        void prepare_field(const context& context) {
            auto& exchange = std::get<I>(exchanges_);
            if (!exchange) {
                const auto comm = context->distributed->mpi_comm();
                const auto ranks = num_ranks(context);
                using Scalar = typename std::tuple_element<I, tuple>::type;
                exchange = std::make_unique<mpi::SparseExchange<Scalar>>(
                    comm, I, mpi::sparse_segments(mesh.get_view_scan_send_mirror().data(), ranks, 1, comm),
                    mpi::sparse_segments(mesh.get_view_scan_recv_mirror().data(), ranks, 1, comm),
                    mesh.get_boundary_size(), mesh.get_ghost_size());
                std::get<I>(ghost_user_data_) = exchange->recv_buffer();
            }
            const auto buffer = exchange->send_buffer();
            const auto values = std::get<I>(user_data_);
            const auto indices = mesh.get_view_boundary_sfc_index();
            Kokkos::parallel_for("cell_halo_pack", indices.extent(0), KOKKOS_LAMBDA(Integer i) {
                buffer(i) = values(indices(i));
            });
        }
        template<std::size_t... I>
        void exchange_fields(const context& context, std::index_sequence<I...>) {
            (prepare_field<I>(context), ...);
            (std::get<I>(exchanges_)->start(), ...);
            (std::get<I>(exchanges_)->wait(), ...);
        }
        std::tuple<std::unique_ptr<mpi::SparseExchange<T>>...> exchanges_;

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
