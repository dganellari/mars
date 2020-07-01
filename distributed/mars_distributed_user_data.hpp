#ifndef GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_
#define GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars
{

/* template <Integer Dim, Integer ManifoldDim, Integer Type, typename ...T> */
template <class Mesh, typename ...T>
class UserData
{
    /* using Mesh = mars::Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>; */
    /* using user_tuple = std::tuple<ViewVectorType<T>... >; */
    using user_tuple = ViewsTuple<T...>;
    using tuple  = std::tuple<T...>;

    using simplex_type = typename Mesh::Elem;

    template<Integer idx>
    using type = typename std::tuple_element<idx, tuple>::type;

public:
    MARS_INLINE_FUNCTION UserData(Mesh *mesh) : host_mesh(mesh), mesh(nullptr)
    {
        const Integer size = host_mesh->get_chunk_size();
        reserve_user_data(user_data_, "user_data", size);
        copy_mesh_to_device();
    }

    MARS_INLINE_FUNCTION UserData(Mesh *mesh, const user_tuple& data) :
               host_mesh(mesh), mesh(nullptr), user_data_(data)
    {
        copy_mesh_to_device();
    }

    MARS_INLINE_FUNCTION void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size)
    {
        /* constexpr Integer data_size = std::tuple_size<std::decay<decltype(tuple)>::type>::value; */
        /* reserve_view_tuple(tuple, size, view_desc); */
        apply_impl(resize_view_functor(view_desc, size), tuple);
    }

    template <typename ElementType>
    struct FillBufferData
    {
        ElementType buffer_data;
        ElementType user_data;
        ViewVectorType<Integer> boundary_lsfc_index;

        FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd)
            : buffer_data(bf), user_data(ud), boundary_lsfc_index(bd)
        {
        }

        MARS_INLINE_FUNCTION
        void operator()(Integer i) const
        {
            const Integer lsfc_index = boundary_lsfc_index(i);
            buffer_data(i) = user_data(lsfc_index);
        }
    };

    struct fill_buffer_data_functor
    {
        fill_buffer_data_functor(std::string d, size_t s, ViewVectorType<Integer> b) : desc(d), size(s), boundary_lsfc_index(b) {}

        template <typename ElementType>
        void operator()(ElementType &el_1, ElementType &el_2) const
        {
            Kokkos::parallel_for(desc, size,
                                 FillBufferData<ElementType>(el_1, el_2, boundary_lsfc_index));

        }

        std::string desc;
        size_t size;

        ViewVectorType<Integer> boundary_lsfc_index;
    };

    MARS_INLINE_FUNCTION void
    fill_buffer_data(user_tuple &buffer_data)
    {
        const Integer size = host_mesh->get_view_boundary().extent(0);

        reserve_user_data(buffer_data, "buffer_data", size);

        ViewVectorType<Integer> boundary_lsfc_index = host_mesh->get_view_boundary_sfc_index();
        apply_impl(fill_buffer_data_functor("fill_buffer_data", size, boundary_lsfc_index), buffer_data, user_data_);
    }

    void exchange_ghost_counts(const context &context)
    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        int proc_num = rank(context);
        int size = num_ranks(context);


        std::vector<Integer> send_count(size, 0);
        std::vector<Integer> receive_count(size, 0);

        scan_send_mirror = create_mirror_view(host_mesh->get_view_scan_boundary());
        Kokkos::deep_copy(scan_send_mirror, host_mesh->get_view_scan_boundary());

        proc_count = 0;
        for (int i = 0; i < size; ++i)
        {
            Integer count = scan_send_mirror(i + 1) - scan_send_mirror(i);
            if (count > 0)
            {
                send_count[i] = count;
                ++proc_count;
                std::cout<<"****ToProc: "<<i<< " count:"<<count<< " Proc: "<<proc_num<<std::endl;
            }
        }


        context->distributed->i_send_recv_vec(send_count, receive_count, proc_count);

        for (int i = 0; i < size; ++i)
        {
            if (receive_count[i] > 0)
            {
                std::cout << "-----FromProc: " << i << " count:" << receive_count[i]<< " Proc: "<<proc_num<<std::endl;
            }
        }

        //create the scan recv mirror view from the receive count
        reserve_scan_ghost(size + 1);

        scan_recv_mirror = create_mirror_view(get_view_scan_ghost());
        make_scan_index_mirror(scan_recv_mirror,receive_count);
        Kokkos::deep_copy(get_view_scan_ghost(), scan_recv_mirror);
    }

    void exchange_ghost_layer(const context &context)
    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        int proc_num = rank(context);
        int size = num_ranks(context);

        /* auto count_mirror = create_mirror_view(host_mesh->get_view_scan_boundary());
           Kokkos::deep_copy(count_mirror, host_mesh->get_view_scan_boundary()); */

        Integer ghost_size = scan_recv_mirror(size);

        reserve_ghost(ghost_size);

        std::cout<<"Starting mpi send receive for the ghost layer"<<std::endl;
        context->distributed->i_send_recv_view(get_view_ghost(), scan_recv_mirror.data(),
                host_mesh->get_view_boundary(), scan_send_mirror.data(), proc_count);
        std::cout<<"Ending mpi send receive for the ghost layer"<<std::endl;
/*
        parallel_for(
                "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer rank = host_mesh->find_owner_processor(get_view_scan_ghost(), i, 1, proc_num);

                printf(" ghost: %i - %li - proc: %li - rank: %li\n", i, get_view_ghost()(i),
                        rank , proc_num);
                });
*/
    }

    struct exchange_ghost_data_functor
    {
        exchange_ghost_data_functor(const context &c, ViewVectorType<Integer>::HostMirror sr,
                                    ViewVectorType<Integer>::HostMirror ss, Integer p) : con(c), sc_rcv_mirror(sr), sc_snd_mirror(ss), proc(p) {}

        template <typename ElementType>
        void operator()(ElementType &el_1, ElementType &el_2) const
        {
            con->distributed->i_send_recv_view(el_1, sc_rcv_mirror.data(),
                                               el_2, sc_snd_mirror.data(), proc);
        }

        Integer proc;
        ViewVectorType<Integer>::HostMirror sc_rcv_mirror;
        ViewVectorType<Integer>::HostMirror sc_snd_mirror;

        const context &con;
    };

    void exchange_ghost_data(const context &context)
    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        //exchange the ghost layer first since it will be used to find the address of the userdata based on the sfc code.
        /* exchange_ghost_layer(context); */

        int proc_num = rank(context);
        int size = num_ranks(context);

        Integer ghost_size = scan_recv_mirror(size);
        reserve_user_data(ghost_user_data_, "ghost_user_data", ghost_size);

        user_tuple buffer_data;
        fill_buffer_data(buffer_data);

        apply_impl(exchange_ghost_data_functor(context, scan_recv_mirror, scan_send_mirror, proc_count),
                   ghost_user_data_, buffer_data);

        /* print_nth_tuple<1>(proc_num); */

    }

    template <Integer I, typename H = typename std::tuple_element<I, tuple>::type>
    void print_nth_tuple(const int proc)
    {
        using namespace Kokkos;

        ViewVectorType<Integer> scan_ghost = get_view_scan_ghost();
        ViewVectorType<Integer> ghost = get_view_ghost();
        ViewVectorType<H> data = std::get<I>(ghost_user_data_);
        Integer ghost_size = data.extent(0);

        Integer xDim = host_mesh->get_XDim();
        Integer yDim = host_mesh->get_YDim();
        Integer zDim = host_mesh->get_ZDim();

        parallel_for(
            "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer r = find_owner_processor(scan_ghost, i, 1, proc);

                double point[3];
                get_vertex_coordinates_from_sfc<simplex_type::ElemType>(ghost(i), point, xDim, yDim, zDim);

                Octant o = get_octant_from_sfc<simplex_type::ElemType>(ghost(i));
                /* printf("ghost data: %li - %li - %li - (%lf, %lf) - data: %lf - proc: %li - rank: %i\n", i, ghost(i), elem_index(o.x, o.y, o.z, xDim, yDim), point[0], point[1], data(i), r, proc); */
            });
    }

    template <typename ElementType>
    struct InitialCondition
    {
        ElementType user_data;
        Integer proc;

        InitialCondition(ElementType ud, Integer p) : user_data(ud), proc(p)
        {
        }

        MARS_INLINE_FUNCTION
        void operator()(Integer i) const
        {
            user_data(i) = proc;
        }
    };


    struct InitData
    {
        InitData(std::string d, size_t s, tuple t) : desc(d), size(s), tup(t) {}

        template <typename ElementType>
        void operator()(ElementType &el_1, std::size_t I) const
        {
            using ET = typename ElementType::value_type;
            constexpr std::size_t TypeIndex = TypeIdx<tuple, ET>::value;

            Kokkos::parallel_for(desc + std::to_string(I), size,
                                 InitialCondition<ElementType>(el_1, std::get<TypeIndex>(tup)));
        }

        std::string desc;
        size_t size;
        tuple tup;
    };

    MARS_INLINE_FUNCTION void
    init_user_data(T... args)
    {
        const Integer size = host_mesh->get_chunk_size();

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
    MARS_INLINE_FUNCTION void
    parallel_for_data(const Integer size, H f)
    {
        Kokkos::parallel_for("init_initial_cond", size, f);
    }

    template <typename H>
    MARS_INLINE_FUNCTION void
    set_init_cond(H f)
    {
        const Integer size = host_mesh->get_chunk_size();
        Kokkos::parallel_for("init_initial_cond", size, f);
    }

    template <typename H>
    MARS_INLINE_FUNCTION void
    elem_iterate(H f)
    {
        const Integer size = host_mesh->get_chunk_size();
        Kokkos::parallel_for("elem_iterate", size, f);
    }

    template <typename H>
    struct FaceIterate
    {
        FaceIterate(Mesh *m, H f, ViewVectorType<Integer> gl,
                    ViewVectorType<Integer> sg, Integer p, Integer x, Integer y, Integer z)
            : mesh(m), func(f), ghost_layer(gl), scan_ghost(sg), proc(p), xDim(x), yDim(y), zDim(z) {}

        template <Integer dir>
        MARS_INLINE_FUNCTION void iterate(const Integer i, const Octant &ref_octant) const
        {
            for (int side = 0; side < 2; ++side)
            {
                Integer face_nr;

                if (side == 0)
                    face_nr = 2 * dir + 1;
                else
                    face_nr = 2 * dir;

                Octant nbh_oc = face_nbh<simplex_type::ElemType>(ref_octant, face_nr,
                                                                 xDim, yDim, zDim);

                bool ghost = false;
                Integer index;

                if (nbh_oc.is_valid())
                {
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(nbh_oc);

                    Integer owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                    assert(owner_proc >= 0);

                    /* if the face neighbor element is ghost then do a binary search
                     * on the ghost layer to find the index */
                    if (proc != owner_proc)
                    {
                        ghost = true;

                        /* to narrow down the range of search we use the scan ghost
                        and the owner proc of the ghost. */
                        const int start_index = scan_ghost(owner_proc);
                        const int last_index = scan_ghost(owner_proc + 1) - 1;

                        /* as opposed to the whole range: */
                        /* const int start_index = 0;
                        const int last_index = ghost_layer.extent(0) -1; */

                        index = binary_search(ghost_layer, start_index, last_index, enc_oc);
                        assert(index >= 0);
                    }
                    else
                    {
                        //using the sfc (global) to local mapping of the mesh.
                        index = mesh->get_index_of_sfc_elem(enc_oc);
                        assert(index >= 0);
                    }

                    /* printf("Index: %li, o.x: %li, y: %li, elem-index: %li, owner_proc: %li, proc: %li , o.x: %li, y: %li, index: %li, ghost: %i\n", index, ref_octant.x, ref_octant.y, elem_index(ref_octant.x, ref_octant.y, ref_octant.z, xDim, yDim), owner_proc, proc, o.x, o.y, elem_index(o.x, o.y, o.z, xDim, yDim), face.get_second_side().is_ghost()); */
                }

                bool boundary = nbh_oc.shares_boundary_side<simplex_type::ElemType>(xDim, yDim, zDim);

                if ((side == 0 && nbh_oc.is_valid()) || ghost || boundary)
                {
                    Face<simplex_type::ElemType, dir> face;
                    int otherside = side ^ 1;

                    face.get_side(side).set_elem_id(i);
                    face.get_side(side).set_boundary(boundary);
                    //if it is the side element of the ref octant.
                    face.get_side(side).set_origin();

                    if (!boundary)
                    {
                        face.get_side(otherside).set_elem_id(index);
                        face.get_side(otherside).set_ghost(ghost);
                    }

                    func(face);
                }
            }
        }

        MARS_INLINE_FUNCTION
        void operator()(const Integer i) const
        {
            const Integer oc = mesh->get_view_sfc()(i);
            Octant ref_octant = get_octant_from_sfc<simplex_type::ElemType>(oc);

            iterate<0>(i, ref_octant);
            iterate<1>(i, ref_octant);
            //TODO: 3D part
        }

        Mesh *mesh;
        H func;

        ViewVectorType<Integer> ghost_layer;
        ViewVectorType<Integer> scan_ghost;

        Integer proc;
        Integer xDim;
        Integer yDim;
        Integer zDim;
        };

        template <typename H>
        void face_iterate(H f)
        {
            Integer xDim = host_mesh->get_XDim();
            Integer yDim = host_mesh->get_YDim();
            Integer zDim = host_mesh->get_ZDim();

            const Integer size = host_mesh->get_chunk_size();

            Kokkos::parallel_for("face_iterate", size, FaceIterate<H>(mesh, f, get_view_ghost(), get_view_scan_ghost(), host_mesh->get_proc(), xDim, yDim, zDim));
        }

        Mesh *get_mesh() const
        {
            return mesh;
        }

        Mesh *get_host_mesh() const
        {
            return host_mesh;
        }

        void reserve_ghost(const Integer n_elements)
        {
            ghost_ = ViewVectorType<Integer>("ghost_", n_elements);
        }

        void reserve_scan_ghost(const Integer n_elements)
        {
            scan_ghost_ = ViewVectorType<Integer>("scan_ghost_", n_elements);
        }

        MARS_INLINE_FUNCTION
        void set_view_ghost(const ViewVectorType<Integer> &b)
        {
            ghost_ = b;
        }


        MARS_INLINE_FUNCTION
        const Integer get_ghost_elem(const Integer i) const
        {
            return ghost_(i);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_ghost() const
        {
            return ghost_;
        }

        MARS_INLINE_FUNCTION
        void set_view_scan_ghost(const ViewVectorType<Integer> &b)
        {
            scan_ghost_ = b;
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_ghost() const
        {
            return scan_ghost_;
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror  &get_view_scan_recv_mirror() const
        {
            return scan_recv_mirror;
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_send_mirror() const
        {
            return scan_send_mirror;
        }

        MARS_INLINE_FUNCTION
        const user_tuple &get_user_data() const
        {
            return user_data_;
        }

        /* template<std::size_t idx, typename H = typename NthType<idx, T...>::type> */
        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION
            H &
            get_elem_data(const int i) const
        {
            return std::get<idx>(user_data_)(i);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION
        void set_elem_data(const int i, const H value) const
        {
            std::get<idx>(user_data_)(i) = value;
        }

        /* template<std::size_t idx, typename H = NthType<idx, T...>> */
        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_data() const
        {
            return std::get<idx>(user_data_);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION
            H &
            get_ghost_elem_data(const int i) const
        {
            return std::get<idx>(ghost_user_data_)(i);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_ghost_data() const
        {
            return std::get<idx>(ghost_user_data_);
        }

        //does ont perform a deep copy of the view containted in the mesh. Just the mesh object.
        void copy_mesh_to_device()
        {
            Mesh *tmp = (Mesh *)Kokkos::kokkos_malloc(sizeof(Mesh));
            Mesh mCopy = *host_mesh;
            Mesh *oldDeviceMesh = mesh;
            Kokkos::parallel_for(
                "CreateDistributedMeshObject", 1, KOKKOS_LAMBDA(const int &) {
                    // two local copies for m and tmp since this->m, this->mesh host pointers
                    new ((Mesh *)tmp) Mesh(mCopy);
                    if (oldDeviceMesh)
                        oldDeviceMesh->~Mesh();
                    // it only works on a copy since **m is still a host pointer and fails on the device.
                });

            mesh = tmp; //make the mesh pointer a device one so that this init func is not neccessary anymore
        }

    private:
        Mesh *mesh; //device mesh
        Mesh *host_mesh; //host mesh that is copied to device.

        //ghost and boundary layers
        ViewVectorType<Integer> ghost_;
        ViewVectorType<Integer> scan_ghost_;

        //ghost data layer
        user_tuple user_data_;
        user_tuple ghost_user_data_;

        //mirror view on the mesh scan boundary view used for the mpi send receive
        ViewVectorType<Integer>::HostMirror scan_send_mirror;

        //mirror view on the scan_ghost view
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;

        Integer proc_count;
    };

    template <class UserData, Integer Type>
    void create_ghost_layer(const context &context, UserData &data)
    {
        std::cout << "Building the ghost layer (boundary element set)..." << std::endl;
        data.get_host_mesh()->template build_boundary_element_sets<Type>();

        data.exchange_ghost_counts(context);
        data.exchange_ghost_layer(context);
    }

    template <class UserData>
    void exchange_ghost_user_data(const context &context, UserData &data)
    {
        int size = num_ranks(context);

        if (data.get_view_scan_recv_mirror()(size) == data.get_view_ghost().extent(0))
        {
            std::cout << "Exchange the ghost data..." << std::endl;
            data.exchange_ghost_data(context);
        }
        else
        {
           errx(1, "Not allowed to call exchange ghost data before the ghost layer creation. Please call create_ghost_layer method first!");
        }
    }

} // namespace mars

#endif
#endif

#endif
