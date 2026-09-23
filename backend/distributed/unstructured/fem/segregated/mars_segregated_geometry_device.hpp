#pragma once
#include "mars_segregated_geometry.hpp"
#include <cuda_runtime.h>

namespace mars::segregated {

// Accepts ElementDomain's existing SoA coordinates and local-node connectivity.
template<class KeyType, class RealType> struct TetMeshView {
    const KeyType* nodes[4];
    const RealType *x, *y, *z;
    int node_count, element_count;
};

template<class KeyType, class RealType>
__global__ void build_tet_geometry(TetMeshView<KeyType,RealType> mesh,
    TetGeometry<RealType>* geometry, int* error)
{
    const int e = blockIdx.x*blockDim.x+threadIdx.x;
    if (e >= mesh.element_count) return;
    RealType coordinates[12];
    for (int n = 0; n < 4; ++n) {
        const auto node = mesh.nodes[n][e];
        if (node < KeyType(0) || node >= KeyType(mesh.node_count)) { atomicExch(error,1); return; }
        coordinates[3*n] = mesh.x[node]; coordinates[3*n+1] = mesh.y[node]; coordinates[3*n+2] = mesh.z[node];
    }
    if (!tet_geometry(coordinates,geometry[e])) atomicExch(error,1);
}

// Zero sums before scattering uniquely owned elements/faces. For MPI, reverse-add
// and publish complete nodal sums before normalization; never scatter published sums again.
template<class KeyType, class RealType>
__global__ void accumulate_dual_volumes(TetMeshView<KeyType,RealType> mesh,
    const TetGeometry<RealType>* geometry, int begin, int end, RealType* volumes)
{
    const int e = begin+blockIdx.x*blockDim.x+threadIdx.x;
    if (e >= end) return;
    for (int n = 0; n < 4; ++n) atomicAdd(volumes+mesh.nodes[n][e],geometry[e].volume/RealType(4));
}

template<int Components, class KeyType, class RealType>
__global__ void accumulate_interior_gradient(TetMeshView<KeyType,RealType> mesh,
    const TetGeometry<RealType>* geometry, const RealType* field, int begin, int end,
    bool shifted, bool incremental, RealType* numerator)
{
    const int e = begin+blockIdx.x*blockDim.x+threadIdx.x;
    if (e >= end) return;
    RealType values[4*Components], local[4*Components*3];
    for (int n = 0; n < 4; ++n)
        for (int c = 0; c < Components; ++c) values[n*Components+c] = field[mesh.nodes[n][e]*Components+c];
    tet_gradient_numerator<Components>(geometry[e], values, shifted, incremental, local);
    for (int n = 0; n < 4; ++n)
        for (int c = 0; c < Components*3; ++c)
            atomicAdd(numerator+mesh.nodes[n][e]*Components*3+c,local[n*Components*3+c]);
}

template<int Components, class KeyType, class RealType>
__global__ void accumulate_boundary_gradient(TetMeshView<KeyType,RealType> mesh,
    const TetGeometry<RealType>* geometry, const TetBoundaryFace* faces, int count,
    const RealType* field, bool shifted, bool incremental, RealType* numerator, int* error)
{
    const int f = blockIdx.x*blockDim.x+threadIdx.x;
    if (f >= count) return;
    const auto face = faces[f];
    if (face.element < 0 || face.element >= mesh.element_count || face.ordinal < 0 || face.ordinal >= 4) {
        atomicExch(error,1); return;
    }
    RealType values[3*Components], local[3*Components*3], area[3];
    for (int n = 0; n < 3; ++n)
        for (int c = 0; c < Components; ++c)
            values[n*Components+c] = field[mesh.nodes[tet_face_node(face.ordinal,n)][face.element]*Components+c];
    tet_boundary_area(geometry[face.element],face.ordinal,area);
    tri_gradient_numerator<Components>(area, values, shifted, incremental, local);
    for (int n = 0; n < 3; ++n)
        for (int c = 0; c < Components*3; ++c)
            atomicAdd(numerator+mesh.nodes[tet_face_node(face.ordinal,n)][face.element]*Components*3+c,
                      local[n*Components*3+c]);
}

template<int Components, class RealType>
__global__ void normalize_gradient(const RealType* numerator, const RealType* volumes,
    int nodes, RealType relaxation, bool initialized, RealType* gradient, int* error)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i >= nodes*Components*3) return;
    const RealType previous = initialized && relaxation != RealType(1) ? gradient[i] : RealType(0);
    if (!finish_gradient(numerator[i],volumes[i/(Components*3)],previous,relaxation,initialized,gradient[i]))
        atomicExch(error,1);
}
} // namespace mars::segregated
