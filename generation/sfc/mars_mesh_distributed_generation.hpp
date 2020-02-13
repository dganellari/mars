#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "mars_mesh_kokkos.hpp"
#include "mars_sfc_generation.hpp"
#include "mars_utils_kokkos.hpp"
namespace mars
{

template <Integer Dim, Integer ManifoldDim, Integer Type>
bool generate_distributed_cube(const context &context,
                               Mesh<Dim, ManifoldDim, KokkosImplementation, NonSimplex<Type, KokkosImplementation>> &mesh,
                               const Integer xDim, const Integer yDim, const Integer zDim)
{
    using namespace Kokkos;

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
        std::cout << "ElementType:: - :    " << ElementType::Quad4 << std::endl;

        n__anchor_nodes = xDim * yDim;
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

        parallel_for(
            "print_elem", xDim * yDim, KOKKOS_LAMBDA(const int i) {
                printf(" el: %u-%i\n", morton.get_view_elements()(i), i);
            });
    }

     context->distributed->scatter_gids(morton.get_view_elements(), local);

   parallel_for(
        "print_elem_chunk",chunk_size, KOKKOS_LAMBDA(const int i) {
            printf(" elch: %u-%i\n", local(i), proc_num);
        });

    using Elem = mars::NonSimplex<Type, KokkosImplementation>;

    assert(Dim <= 3);
    assert(ManifoldDim <= Dim);

    bool gen_pts = mesh.generate_points(xDim, yDim, zDim, Type);

    bool gen_elm = mesh.generate_elements(xDim, yDim, zDim, Type);

    if (!gen_pts || !gen_elm)
        std::cerr << "Not implemented for other dimensions yet" << std::endl;

    return (gen_pts && gen_elm);
}

} // namespace mars

#endif
#endif

#endif // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
