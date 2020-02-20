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
    std::cout << "rank -:    " << proc_num << std::endl;

    int size = num_ranks(context);
    std::cout << "size - :    " << size << std::endl;

    unsigned int chunk_size = 0;
    unsigned int last_chunk_size = 0;
    unsigned int n__anchor_nodes = 0;
    ViewVectorType<unsigned int> local;

    switch (Type)
    {
    case ElementType::Quad4:
    {
        //std::cout << "ElementType:: - :    " << ElementType::Quad4 << std::endl;

        n__anchor_nodes = xDim * yDim;
        chunk_size = (unsigned int)ceil((double)n__anchor_nodes / size);
        last_chunk_size = chunk_size - (chunk_size * size - n__anchor_nodes);
        local = ViewVectorType<unsigned int>("local_partition_sfc", chunk_size);
        break;
    }
    case ElementType::Hex8:
    {
        //std::cout << "ElementType:: - :    " << ElementType::Quad4 << std::endl;

        n__anchor_nodes = xDim * yDim * zDim;
        chunk_size = (unsigned int)ceil((double)n__anchor_nodes / size);
        last_chunk_size = chunk_size - (chunk_size * size - n__anchor_nodes);
        local = ViewVectorType<unsigned int>("local_partition_sfc", chunk_size);
        break;
    }
    default:
    {
        std::cout << "Not yet implemented for other element types!" << std::endl;
        return false;
    }
    }

    SFC morton;

    bool root = mars::rank(context) == 0;
    if (root)
    {
        std::cout << "chunk_size - :    " << chunk_size << std::endl;

        morton.generate_sfc_elements<Type>(xDim, yDim, zDim);

       /*  parallel_for(
            "print_elem", xDim * yDim, KOKKOS_LAMBDA(const int i) {
                printf(" el: %u-%i\n", morton.get_view_elements()(i), i);
            }); */
    }

    context->distributed->scatter_gids(morton.get_view_elements(), local);

    std::cout << "MPI Scatter ended!"<< std::endl;


   parallel_for(
        "print_elem_chunk",chunk_size, KOKKOS_LAMBDA(const int i) {
            printf(" elch: %u-%i\n", local(i), proc_num);
        });

     mesh.set_view_sfc(local);
    
    //set the chunk size to the remainder for the last mpi processes.
    if(proc_num == size-1)
        mesh.set_chunk_size(last_chunk_size);
    else
        mesh.set_chunk_size(chunk_size);

    assert(Dim <= 3);
    assert(ManifoldDim <= Dim);

	Kokkos::Timer timer_gen;

    //the mesh construct depends on template parameters. 
    bool gen_pts = mesh.template generate_points<Type>(xDim, yDim, zDim);

    bool gen_elm = mesh.generate_elements(xDim, yDim, zDim, Type);

   	double time_gen = timer_gen.seconds();
	std::cout << "Generation 2D distributed kokkos took: " << time_gen << " seconds. Process: "<<proc_num << std::endl;

      if (!gen_pts || !gen_elm)
        std::cerr << "Not implemented for other dimensions yet" << std::endl;

	double time = timer.seconds();
	std::cout << "Total Generation 2D distributed kokkos took: " << time << " seconds. Process: "<<proc_num << std::endl;


    ViewMatrixType<Real> poi = mesh.get_view_points();

    parallel_for(
        "print_elem_chunk1", mesh.get_view_points().extent(0), KOKKOS_LAMBDA(const int i) {
            printf(" pt: [(%f, %f) - %i]\n", poi(i, 0), poi(i, 1), proc_num);
        });

          ViewMatrixType<Integer> eeel = mesh.get_view_elements();

    parallel_for(
        "print_elem_chunk", mesh.get_view_elements().extent(0), KOKKOS_LAMBDA(const int i) {
            printf(" pt: [(%li, %li, %li, %li) - %i]\n", eeel(i, 0), eeel(i, 1), eeel(i, 2), eeel(i, 3), proc_num);
        });

    return (gen_pts && gen_elm);
    return true;
}

} // namespace mars

#endif
#endif

#endif // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
