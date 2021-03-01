#ifndef MARS_SERIAL_MESH_TYPE_HPP
#define MARS_SERIAL_MESH_TYPE_HPP

namespace mars {
    template <class Mesh>
    struct SerialMeshType {};

    template <int Dim, int ManifoldDim>
    struct SerialMeshType<
        mars::Mesh<Dim, ManifoldDim, KokkosImplementation, Simplex<Dim, ManifoldDim, KokkosImplementation>>> {
        using Type = mars::Mesh<Dim, ManifoldDim>;
    };

    template <int Dim, int ManifoldDim, int ElemType>
    struct SerialMeshType<
        mars::Mesh<Dim, ManifoldDim, KokkosImplementation, NonSimplex<ElemType, KokkosImplementation>>> {
        using Type = mars::Mesh<Dim, ManifoldDim, DefaultImplementation, NonSimplex<ElemType>>;
    };

    template <>
    struct SerialMeshType<mars::ParallelQuad4Mesh> {
        using Type = mars::Mesh<2, 2, DefaultImplementation, NonSimplex<4>>;
    };

    template <>
    struct SerialMeshType<mars::ParallelMesh2> {
        using Type = mars::Mesh<2, 2>;
    };

    template <>
    struct SerialMeshType<mars::ParallelMesh3> {
        using Type = mars::Mesh<3, 3>;
    };

    template <>
    struct SerialMeshType<mars::ParallelMesh4> {
        using Type = mars::Mesh<4, 4>;
    };
}  // namespace mars

#endif  // MARS_SERIAL_MESH_TYPE_HPP
