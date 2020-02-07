#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_

#include "mars_parallel_descriptor.hpp"
#ifdef WITH_KOKKOS
#include "mars_sfc_generation.hpp"

namespace mars
{

template <Integer Dim, Integer ManifoldDim, Integer Type>
bool generate_distributed_cube(
    Mesh<Dim, ManifoldDim, KokkosImplementation, NonSimplex<Type, KokkosImplementation>> &mesh,
    const Integer xDim, const Integer yDim, const Integer zDim)
{

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

#endif // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
