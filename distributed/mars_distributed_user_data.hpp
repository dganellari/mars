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

    MARS_INLINE_FUNCTION void reserve_ghost_user_data()
    {   //should be done based on the receive count total. not on the boundary
     /*    const Integer size = mesh->get_view_boundary().extent(0);
        ghost_user_data = ViewVectorType<T>("ghost_user_data", size); */
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
        auto count_mirror = create_mirror_view(mesh->get_view_scan_boundary());
        Kokkos::deep_copy(count_mirror, mesh->get_view_scan_boundary()); 

        Integer proc_count = 0;
        for (int i = 0; i < size; ++i)
        {
            Integer count = count_mirror(i + 1) - count_mirror(i);
            if (count > 0)
            {
                send_count[i] = count;
                ++proc_count;
            }
        }

        context->distributed->i_send_recv_vec(send_count, receive_count, proc_count);
    }

    Mesh* get_mesh() const
    {
        return mesh;
    }

private:
    Mesh *mesh;
    ViewVectorType<T> user_data;
    ViewVectorType<T> ghost_user_data;

    std::vector<Integer> send_count;
    std::vector<Integer> receive_count;
    std::vector<Integer> scan_receive_count;
};

template <Integer Dim, Integer ManifoldDim, Integer Type, typename T>
void exchange_ghost_user_data(const context &context, UserData<Dim, ManifoldDim, Type, T> &data)
{
    using namespace Kokkos;

   

    //reserve_ghost_user_data();
    const Integer size = data.get_mesh()->get_view_boundary().extent(0);
    ViewVectorType<T> buffer_data("buffer_data", size);

    data.fill_buffer_data(buffer_data);
    data.exchange_ghost_counts(context);
   // data.exchange_ghost_data(context, data, buffer_data);
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
