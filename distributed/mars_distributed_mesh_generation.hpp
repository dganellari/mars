#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"
namespace mars
{

template <Integer Dim, Integer ManifoldDim, Integer Type>
using DMesh = Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;

template <Integer Dim, Integer ManifoldDim, Integer Type>
void broadcast_gp_np(const context &context, DMesh<Dim, ManifoldDim, Type> &mesh, const std::vector<int> &counts, ViewVectorType<unsigned int> elems, Integer n__anchor_nodes)
{
    int proc_num = rank(context);
    int size = num_ranks(context);

    ViewVectorType<Integer> GpNp = ViewVectorType<Integer>("global_static_partition", 2 * (size + 1));

    if (proc_num == 0)
    {
        auto GpNp_host = create_mirror_view(GpNp);

        auto elem_view_host = create_mirror_view(elems);
        deep_copy(elem_view_host, elems);

        for (int i = 0; i < size; ++i)
        {
            //counts[0] is always chunk size
            GpNp_host(2 * i) = elem_view_host(i * counts[0]); // acc sum scan giving the first element per process
            GpNp_host(2 * i + 1) = counts[i];
        }
        //insert the last element of the sfc adding 1 to it (to make the last element not part of the linearization) for the binary search to work properly
        GpNp_host(2 * size) = elem_view_host(n__anchor_nodes - 1) + 1;
        deep_copy(GpNp, GpNp_host);
    }

    context->distributed->broadcast(GpNp); //broadcast to all processors.
    std::cout << "MPI broadcast ended!" << std::endl;

    mesh.set_view_gp(GpNp);

    /* parallel_for(
       "print_elem_gp:", size+1, KOKKOS_LAMBDA(const int i) {
       printf(" elch: (%li-%li) - %i - %i\n", GpNp(2*i),  GpNp(2*i+1), i, proc_num);
       }); */
}

template <Integer Dim, Integer ManifoldDim, Integer Type>
bool generate_distributed_cube(const context &context, DMesh<Dim, ManifoldDim, Type> &mesh,
                               const Integer xDim, const Integer yDim, const Integer zDim)
{
    using namespace Kokkos;

    Kokkos::Timer timer;

    int proc_num = rank(context);
    // std::cout << "rank -:    " << proc_num << std::endl;

    int size = num_ranks(context);
    // std::cout << "size - :    " << size << std::endl;

    unsigned int n__anchor_nodes = 0;

    switch (Type)
    {
    case ElementType::Quad4:
    {
        n__anchor_nodes = xDim * yDim;
        break;
    }
    case ElementType::Hex8:
    {
        n__anchor_nodes = xDim * yDim * zDim;
        break;
    }
    default:
    {
        std::cout << "Not yet implemented for other element types!" << std::endl;
        return false;
    }
    }

    //unsigned int chunk_size = (unsigned int)ceil((double)n__anchor_nodes / size);
    //Integer chunk_size = n__anchor_nodes / size + (n__anchor_nodes % size != 0);
    //Integer last_chunk_size = chunk_size - (chunk_size * size - n__anchor_nodes);
    Integer chunk_size = n__anchor_nodes / size;
    Integer last_chunk_size = chunk_size - (chunk_size * size - n__anchor_nodes);

    assert(chunk_size > 0);

    if (chunk_size <= 0)
    {
        errx(1, " Invalid number of mpi processes. Defined more mpi processes than mesh elements to be generated!");
    }

    SFC morton;

    bool root = mars::rank(context) == 0;
    if (root)
    {
        std::cout << "chunk_size - :    " << chunk_size << std::endl;
        std::cout << "n__anchor_nodes size:: - :    " << n__anchor_nodes << std::endl;
        std::cout << "last_chunk_size size:: - :    " << last_chunk_size << std::endl;

        morton.generate_sfc_elements<Type>(xDim, yDim, zDim);
    }

    std::vector<int> counts(size);

    for (int i = 0; i < size; ++i)
    {
        if (i == size - 1)
            counts[i] = last_chunk_size;
        else
            counts[i] = chunk_size;
    }

    //set the chunk size to the remainder for the last mpi processes.
    if (proc_num == size - 1)
    {
        chunk_size = last_chunk_size;
    }

    ViewVectorType<unsigned int> local = ViewVectorType<unsigned int>("local_partition_sfc", chunk_size);

    context->distributed->scatterv_gids(morton.get_view_elements(), local, counts);

    std::cout << "MPI Scatter ended!" << std::endl;

    /*  parallel_for(
        "print_elem_chunk",chunk_size, KOKKOS_LAMBDA(const int i) {
        printf(" elch: %u-%i\n", local(i), proc_num);
        }); */

    mesh.set_view_sfc(local);
    mesh.set_chunk_size(chunk_size);
    mesh.set_proc(proc_num);

    broadcast_gp_np(context, mesh, counts, morton.get_view_elements(), n__anchor_nodes);

    assert(Dim <= 3);
    assert(ManifoldDim <= Dim);

    Kokkos::Timer timer_gen;

    //the mesh construct depends on template parameters.
    bool gen_pts = mesh.template generate_points<Type>(xDim, yDim, zDim);

    bool gen_elm = mesh.template generate_elements<Type>(xDim, yDim, zDim);

    /* mesh.template build_ghost_element_sets<Type>(xDim, yDim, zDim); */

    mesh.template build_boundary_element_sets<Type>(xDim, yDim, zDim);

    double time_gen = timer_gen.seconds();
    std::cout << "Distributed Generation kokkos took: " << time_gen << " seconds. Process: " << proc_num << std::endl;

    if (!gen_pts || !gen_elm)
        std::cerr << "Not implemented for other dimensions yet" << std::endl;

    double time = timer.seconds();
    std::cout << "Total distributed generation  kokkos took: " << time << " seconds. Process: " << proc_num << std::endl;

    return (gen_pts && gen_elm);
}

} // namespace mars

#endif
#endif

#endif // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
