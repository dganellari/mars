#ifndef GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_
#define GENERATION_MARS_DISTRIBUTED_USERDATA_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"
#include "mars_distributed_utils.hpp"

namespace mars
{

/* template <Integer Dim, Integer ManifoldDim, Integer Type, typename ...T> */
template <class Mesh, typename ...T>
class UserData
{
    /* using Mesh = mars::Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>; */
    using user_tuple = std::tuple<ViewVectorType<T>... >;

public:
    MARS_INLINE_FUNCTION UserData(Mesh *mesh) : mesh(mesh)
    {
        const Integer size = mesh->get_chunk_size();
        reserve_user_data(user_data_, "user_data", size);
    }

    MARS_INLINE_FUNCTION UserData(Mesh *mesh, const user_tuple& data) :
                mesh(mesh), user_data_(data)
    {
    }

    MARS_INLINE_FUNCTION void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size)
    {
        /* constexpr Integer data_size = std::tuple_size<std::decay<decltype(tuple)>::type>::value; */
        /* reserve_view_tuple(tuple, size, view_desc); */
        apply_impl(resize_view_functor(view_desc, size), tuple);
    }

    /* MARS_INLINE_FUNCTION void reserve_ghost_user_data(const Integer size)
    {
        const Integer data_size = std::tuple_size<decltype(user_data)>::value;

        for (int i = 0; i < data_size; ++i)
        {
            std::get<i>(ghost_user_data) = typename std::tuple_element<i, decltype(ghost_user_data)>::type("ghost_user_data_" + i, size);
        }
    }
 */
    MARS_INLINE_FUNCTION void reserve_ghost_count(const Integer size)
    {
        send_count.assign(size, 0);
        receive_count.assign(size, 0);
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
        const Integer size = mesh->get_view_boundary().extent(0);

        reserve_user_data(buffer_data, "buffer_data", size);

        ViewVectorType<Integer> boundary_lsfc_index = mesh->get_view_boundary_sfc_index();
        apply_impl(fill_buffer_data_functor("fill_buffer_data", size, boundary_lsfc_index), buffer_data, user_data_);
    }

    void exchange_ghost_counts(const context &context)
    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        int proc_num = rank(context);
        int size = num_ranks(context);

        reserve_ghost_count(size);
        scan_send_mirror = create_mirror_view(mesh->get_view_scan_boundary());
        Kokkos::deep_copy(scan_send_mirror, mesh->get_view_scan_boundary());

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
    }

    // returns the prefix sum of C
    template <typename C>
    void make_scan_index_mirror(const ViewVectorType<Integer>::HostMirror &out, C const &c)
    {
        static_assert(
            std::is_integral<typename C::value_type>::value,
            "make_index only applies to integral types");

        out(0) = 0;
        std::partial_sum(c.begin(), c.end(), out.data() + 1);
    }

    void exchange_ghost_layer(const context &context)
    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        int proc_num = rank(context);
        int size = num_ranks(context);

        /* auto count_mirror = create_mirror_view(mesh->get_view_scan_boundary());
           Kokkos::deep_copy(count_mirror, mesh->get_view_scan_boundary()); */

        reserve_scan_ghost(size + 1);

        scan_recv_mirror = create_mirror_view(get_view_scan_ghost());
        make_scan_index_mirror(scan_recv_mirror,receive_count);
        Kokkos::deep_copy(get_view_scan_ghost(), scan_recv_mirror);

        Integer ghost_size = scan_recv_mirror(size);

        reserve_ghost(ghost_size);

        std::cout<<"Starting mpi send receive for the ghost layer"<<std::endl;
        context->distributed->i_send_recv_view(get_view_ghost(), scan_recv_mirror.data(),
                mesh->get_view_boundary(), scan_send_mirror.data(), proc_count);
        std::cout<<"Ending mpi send receive for the ghost layer"<<std::endl;
/*
        parallel_for(
                "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer rank = mesh->find_owner_processor(get_view_scan_ghost(), i, 1, proc_num);

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
        exchange_ghost_layer(context);

        int proc_num = rank(context);
        int size = num_ranks(context);

        Integer ghost_size = scan_recv_mirror(size);
        reserve_user_data(ghost_user_data_, "ghost_user_data", ghost_size);

        user_tuple buffer_data;
        fill_buffer_data(buffer_data);

        apply_impl(exchange_ghost_data_functor(context, scan_recv_mirror, scan_send_mirror, proc_count),
                   ghost_user_data_, buffer_data);

        parallel_for(
            "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer rank = mesh->find_owner_processor(get_view_scan_ghost(), i, 1, proc_num);
                std::stringstream stream;
                stream << " ghost data: " << i << " - " << get_view_ghost()(i) << " data: " << std::get<2>(ghost_user_data_)(i)
                       << " proc: " << rank << " - rank: " << proc_num << std::endl;
                std::cout << stream.str();
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
        InitData(std::string d, size_t s, Integer p) : desc(d), size(s), proc(p) {}

        template <typename ElementType>
        void operator()(ElementType &el_1) const
        {
            Kokkos::parallel_for(desc, size,
                                 InitialCondition<ElementType>(el_1, proc));
        }

        std::string desc;
        size_t size;
        Integer proc;
    };

    /* TODO: make it for T... values */
    MARS_INLINE_FUNCTION void
    init_user_data()
    {
        const Integer size = mesh->get_chunk_size();
        const Integer proc = mesh->get_proc();

        apply_impl(InitData("init_data", size, proc), user_data_);
    }

    template <typename elem_type>
    MARS_INLINE_FUNCTION static void init_cond(elem_type v, int i);

    template <typename ElementType>
    struct InitCond
    {
        /* InitCond(user_tuple t, H ic) : tuple(t), init_cond(ic) {} */
        InitCond(ElementType t) : view(t) {}

        void operator()(int i) const
        {
            init_cond(view, i);
        }

        ElementType view;
    };

    struct InitCondData
    {
        InitCondData(std::string d, size_t s) : desc(d), size(s) {}

        template <typename ElementType>
        void operator()(ElementType &el_1) const
        {
            Kokkos::parallel_for(desc, size,
                                 InitCond<ElementType>(el_1));
        }

        std::string desc;
        size_t size;
    };

    /* template <typename H> */
    MARS_INLINE_FUNCTION void
    set_init_cond()
    {
        const Integer size = mesh->get_chunk_size();
        apply_impl(InitCondData("init_data", size), user_data_);
    }

    Mesh *
    get_mesh() const
    {
        return mesh;
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

private:
    Mesh *mesh;

     //ghost and boundary layers
    ViewVectorType<Integer> ghost_;
    ViewVectorType<Integer> scan_ghost_;

    //ghost data layer
    /* ViewVectorType<T> user_data; */
    /* ViewVectorType<T> ghost_user_data; */
    user_tuple user_data_;
    user_tuple ghost_user_data_;

    std::vector<Integer> send_count;
    //mirror view on the mesh scan boundary view used for the mpi send receive
    ViewVectorType<Integer>::HostMirror scan_send_mirror;

    std::vector<Integer> receive_count;
    //mirror view on the scan_ghost view
    ViewVectorType<Integer>::HostMirror scan_recv_mirror;

    Integer proc_count;
};

template <class UserData>
void exchange_ghost_user_data(const context &context, UserData &data)
{
    using namespace Kokkos;

    data.exchange_ghost_counts(context);
    data.exchange_ghost_data(context);
}

/* template <typename Functor>
inline void init_data(UserData data, Mesh mesh, Functor init_cond)
{
    init_cond();
} */


} // namespace mars

#endif
#endif

#endif
