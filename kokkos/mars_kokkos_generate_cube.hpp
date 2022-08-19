#ifndef MARS_KOKKOS_GENERATE_CUBE_HPP
#define MARS_KOKKOS_GENERATE_CUBE_HPP

#include "mars_mesh_kokkos.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim, class KokkosImplementation>
    bool generate_cube(
        Mesh<Dim, ManifoldDim, KokkosImplementation, mars::Simplex<Dim, ManifoldDim, KokkosImplementation>>& mesh,
        const Integer xDim,
        const Integer yDim,
        const Integer zDim) {
        using Elem = mars::Simplex<Dim, ManifoldDim>;

        assert(Dim <= 3);
        assert(ManifoldDim <= Dim);

        bool gen_pts = mesh.generate_points(xDim, yDim, zDim);

        bool gen_elm = mesh.generate_elements(xDim, yDim, zDim);

        if (!gen_pts || !gen_elm) std::cerr << "Not implemented for other dimensions yet" << std::endl;

        return (gen_pts && gen_elm);
    }

    // non simplex cube generation.
    template <Integer Dim, Integer ManifoldDim, Integer Type>
    bool generate_cube(Mesh<Dim, ManifoldDim, KokkosImplementation, NonSimplex<Type, KokkosImplementation>>& mesh,
                       const Integer xDim,
                       const Integer yDim,
                       const Integer zDim) {
        using Elem = mars::NonSimplex<Type, KokkosImplementation>;

        assert(Dim <= 3);
        assert(ManifoldDim <= Dim);

        bool gen_pts = mesh.generate_points(xDim, yDim, zDim, Type);

        bool gen_elm = mesh.generate_elements(xDim, yDim, zDim, Type);

        if (!gen_pts || !gen_elm) std::cerr << "Not implemented for other dimensions yet" << std::endl;

        return (gen_pts && gen_elm);
    }

    inline bool generate_square(ParallelMesh2& mesh, const Integer xDim, const Integer yDim) {
        return generate_cube(mesh, xDim, yDim, 0);
    }

    inline bool generate_line(ParallelMesh1& mesh, const Integer xDim) { return generate_cube(mesh, xDim, 0, 0); }

}  // namespace mars

#endif