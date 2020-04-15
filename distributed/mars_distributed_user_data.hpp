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

template <Integer Dim, Integer ManifoldDim, Integer Type, typename T>
class UserData
{
    using Mesh = mars::Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;

public:
    MARS_INLINE_FUNCTION UserData(Mesh *mesh) : mesh(mesh)
    {
        reserve_user_data();
    }

    MARS_INLINE_FUNCTION UserData(Mesh *mesh, const ViewVectorType<T>& data) : 
                mesh(mesh), user_data(data)
    {
    }

    MARS_INLINE_FUNCTION void reserve_user_data()
    {
        user_data = ViewVectorType<T>("user_data", mesh->get_chunk_size());
    }

    MARS_INLINE_FUNCTION void reserve_ghost_user_data(const Integer size)
    {  
        ghost_user_data = ViewVectorType<T>("ghost_user_data", size);
    }

    MARS_INLINE_FUNCTION void reserve_ghost_count(const Integer size)
    {
        send_count.assign(size, 0);
        receive_count.assign(size, 0);
    }

    struct FillBufferData
    {
        ViewVectorType<T> buffer_data;
        ViewVectorType<T> user_data;
        ViewMatrixType<Integer> boundary;

        FillBufferData(ViewVectorType<T> bf, ViewVectorType<T> ud, ViewMatrixType<Integer> bd)
            : buffer_data(bf), user_data(ud), boundary(bd)
        {
        }

        MARS_INLINE_FUNCTION 
        void operator()(Integer i) const
        {
            const Integer lsfc_index = boundary(i, 1);
            buffer_data(i) = user_data(lsfc_index);
        }
    };

    MARS_INLINE_FUNCTION void
    fill_buffer_data(ViewVectorType<T> buffer_data)
    {
        ViewMatrixType<Integer> boundary = mesh->get_view_boundary();

        Kokkos::parallel_for("fill_buffer_data", buffer_data.extent(0),
                     FillBufferData(buffer_data, user_data, boundary));
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
    void make_scan_index_mirror(const ViewVectorType<Integer>::HostMirror& out, C const &c)
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

        mesh->reserve_scan_ghost(size + 1);

        scan_recv_mirror = create_mirror_view(mesh->get_view_scan_ghost());
        make_scan_index_mirror(scan_recv_mirror,receive_count);
        Kokkos::deep_copy(mesh->get_view_scan_ghost(), scan_recv_mirror);

        Integer ghost_size = scan_recv_mirror(size);

        mesh->reserve_ghost(ghost_size);

        //context->distributed->i_send_recv_view(mesh->ghost_, mesh->scan_ghost_, mesh->boundary_, 
                   // mesh->scan_boundary_unsigned, proc_count);

/*         for (int i = 0; i < size; ++i)
        {
            if (receive_count[i] > 0)
            {
                std::cout << "-----FromProc: " << i << " count:" << receive_count[i] << " Proc: " << proc_num << std::endl;
            }
        } */
    }

    void exchange_ghost_data(const context &context, const ViewVectorType<T> &buffer_data)
    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        exchange_ghost_layer(context);
        
        int proc_num = rank(context);
        int size = num_ranks(context);

        Integer ghost_size = scan_recv_mirror(size);
        reserve_ghost_user_data(ghost_size);
        
       // context->distributed->i_send_recv_view(send_count, receive_count, proc_count);

        /*         for (int i = 0; i < size; ++i)
        {
            if (receive_count[i] > 0)
            {
                std::cout << "-----FromProc: " << i << " count:" << receive_count[i] << " Proc: " << proc_num << std::endl;
            }
        } */
    }

    Mesh *get_mesh() const
    {
        return mesh;
    }

private:
    Mesh *mesh;
    ViewVectorType<T> user_data;
    ViewVectorType<T> ghost_user_data;

    std::vector<Integer> send_count;
    ViewVectorType<Integer>::HostMirror scan_send_mirror;

    std::vector<Integer> receive_count;
    ViewVectorType<Integer>::HostMirror scan_recv_mirror;

    Integer proc_count;
};

template <Integer Dim, Integer ManifoldDim, Integer Type, typename T>
void exchange_ghost_user_data(const context &context, UserData<Dim, ManifoldDim, Type, T> &data)
{
    using namespace Kokkos;

    const Integer size = data.get_mesh()->get_view_boundary().extent(0);
    ViewVectorType<T> buffer_data("buffer_data", size);

    data.fill_buffer_data(buffer_data);
    data.exchange_ghost_counts(context);
    data.exchange_ghost_data(context, buffer_data);
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
