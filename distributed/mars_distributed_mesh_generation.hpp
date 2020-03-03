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
bool generate_distributed_cube(const context &context, DMesh<Dim, ManifoldDim, Type> &mesh,
                               const Integer xDim, const Integer yDim, const Integer zDim)
{
    using namespace Kokkos;

	Kokkos::Timer timer;

    int proc_num =  rank(context);
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

    unsigned int chunk_size = (unsigned int)ceil((double)n__anchor_nodes / size);
    unsigned int last_chunk_size  = chunk_size - (chunk_size * size - n__anchor_nodes);
    
    SFC morton;

    bool root = mars::rank(context) == 0;
    if (root)
    {
        std::cout << "chunk_size - :    " << chunk_size << std::endl;
        std::cout << "n__anchor_nodes size:: - :    " << n__anchor_nodes << std::endl;
        std::cout << "last_chunk_size size:: - :    " << last_chunk_size << std::endl;

        morton.generate_sfc_elements<Type>(xDim, yDim, zDim);

        /*  parallel_for(
            "print_elem", xDim * yDim, KOKKOS_LAMBDA(const int i) {
                printf(" el: %u-%i\n", morton.get_view_elements()(i), i);
            }); */
    }

    std::vector<int> counts(size);

    for (int i = 0; i < size; ++i)
    {
        if (i == size - 1)
            counts[i] = last_chunk_size;
        else
            counts[i] = chunk_size;

        //printf(" count: %d - ", counts[i]);
    }
    //printf(" endcount\n");

    //set the chunk size to the remainder for the last mpi processes.
    if(proc_num == size - 1)
    {
        chunk_size = last_chunk_size;
    }

    ViewVectorType<unsigned int> local = ViewVectorType<unsigned int>("local_partition_sfc", chunk_size);
    
    context->distributed->scatterv_gids(morton.get_view_elements(), local, counts);

    std::cout << "MPI Scatter ended!"<< std::endl;


  /*  parallel_for(
        "print_elem_chunk",chunk_size, KOKKOS_LAMBDA(const int i) {
            printf(" elch: %u-%i\n", local(i), proc_num);
        }); */

     mesh.set_view_sfc(local);
     mesh.set_chunk_size(chunk_size);

    assert(Dim <= 3);
    assert(ManifoldDim <= Dim);

	Kokkos::Timer timer_gen;

    //the mesh construct depends on template parameters. 
    bool gen_pts = mesh.template generate_points<Type>(xDim, yDim, zDim);

    bool gen_elm = mesh.template generate_elements<Type>(xDim, yDim, zDim);

   	double time_gen = timer_gen.seconds();
	std::cout << "Distributed Generation kokkos took: " << time_gen << " seconds. Process: "<<proc_num << std::endl;

      if (!gen_pts || !gen_elm)
        std::cerr << "Not implemented for other dimensions yet" << std::endl;

	double time = timer.seconds();
	std::cout << "Total distributed generation  kokkos took: " << time << " seconds. Process: "<<proc_num << std::endl;
    
    return (gen_pts && gen_elm);
}

} // namespace mars

#endif
#endif

#endif // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
