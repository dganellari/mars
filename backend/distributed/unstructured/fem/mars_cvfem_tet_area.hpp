#pragma once

#include "mars_cvfem_kernel.hpp"   // Tet4CVFEM
#include "mars_element_traits.hpp" // d_tetLRSCV
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// Tet SCS area-vector precompute.
//
// Counterpart of the hex precomputeAreaVectorsKernel. Fills d_areaVec_{x,y,z}
// with elementCount*6 vectors -- one per (element, SCS edge) -- so the
// matrix-free NS scatter kernels (advection, divergence, grad) can compute a
// face flux as a single dot product v . A.
//
// Geometry: median-dual construction. For each edge (L,R) of the tet, the SCS
// face is the quadrilateral [edge midpoint M, centroid of one adjacent tet
// face, tet centroid C, centroid of the other adjacent tet face]. The four
// points wind around the dual face; its area-normal is A.
//
// Orientation: A points from L toward R, matching the hex convention that the
// NS kernels rely on (flux v.A leaves L with -flux, enters R with +flux). We
// fix the sign at the end by checking A . (xR - xL) > 0.
//
// The 6 edges and the two tet faces each edge belongs to. A tet face is named
// by the node it is OPPOSITE. Edge (L,R) is shared by the two faces that are
// opposite the OTHER two nodes. Edge order matches d_tetLRSCV / Tet4CVFEM:
//   ip0:(0,1) ip1:(1,2) ip2:(0,2) ip3:(0,3) ip4:(1,3) ip5:(2,3)
// For edge (L,R) with the remaining nodes {a,b}, the two adjacent face
// centroids are centroid(L,R,a) and centroid(L,R,b).
__device__ __constant__ int d_tetEdgeOtherNodes[6][2] = {
    {2, 3},  // edge (0,1): other nodes 2,3
    {0, 3},  // edge (1,2): other nodes 0,3
    {1, 3},  // edge (0,2): other nodes 1,3
    {1, 2},  // edge (0,3): other nodes 1,2
    {0, 2},  // edge (1,3): other nodes 0,2
    {0, 1}   // edge (2,3): other nodes 0,1
};

// Area-normal of the planar quad (p0,p1,p2,p3) by splitting into two triangles
// sharing the p0-p2 diagonal. Same triangulation idea as the hex helper.
template<typename RealType>
__device__ inline void quad_area_normal_tet(const RealType p0[3],
                                            const RealType p1[3],
                                            const RealType p2[3],
                                            const RealType p3[3],
                                            RealType A[3])
{
    // Triangle 1: (p0, p1, p2); Triangle 2: (p0, p2, p3).
    // Each triangle's area-normal is 0.5 * (e1 x e2); sum the two.
    RealType a1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    RealType b1[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    RealType a2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    RealType b2[3] = {p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]};

    A[0] = RealType(0.5) * ((a1[1]*b1[2] - a1[2]*b1[1]) + (a2[1]*b2[2] - a2[2]*b2[1]));
    A[1] = RealType(0.5) * ((a1[2]*b1[0] - a1[0]*b1[2]) + (a2[2]*b2[0] - a2[0]*b2[2]));
    A[2] = RealType(0.5) * ((a1[0]*b1[1] - a1[1]*b1[0]) + (a2[0]*b2[1] - a2[1]*b2[0]));
}

template<typename KeyType, typename RealType>
__global__ void precomputeTetAreaVectorsKernel(
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    size_t numElements,
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    RealType* __restrict__ d_areaVec_x,
    RealType* __restrict__ d_areaVec_y,
    RealType* __restrict__ d_areaVec_z)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    constexpr int nnodes = 4;
    constexpr int nscs   = 6;

    KeyType nodes[nnodes] = {d_conn0[e], d_conn1[e], d_conn2[e], d_conn3[e]};
    RealType coords[nnodes][3];
    #pragma unroll
    for (int n = 0; n < nnodes; ++n) {
        coords[n][0] = d_x[nodes[n]];
        coords[n][1] = d_y[nodes[n]];
        coords[n][2] = d_z[nodes[n]];
    }

    // Tet centroid (shared by all 6 SCS faces).
    RealType C[3];
    #pragma unroll
    for (int d = 0; d < 3; ++d)
        C[d] = RealType(0.25) * (coords[0][d] + coords[1][d] + coords[2][d] + coords[3][d]);

    #pragma unroll
    for (int ip = 0; ip < nscs; ++ip) {
        int L = d_tetLRSCV[ip * 2];
        int R = d_tetLRSCV[ip * 2 + 1];
        int a = d_tetEdgeOtherNodes[ip][0];
        int b = d_tetEdgeOtherNodes[ip][1];

        // Edge midpoint.
        RealType M[3];
        // Centroids of the two adjacent tet faces (L,R,a) and (L,R,b).
        RealType Fa[3], Fb[3];
        #pragma unroll
        for (int d = 0; d < 3; ++d) {
            M[d]  = RealType(0.5) * (coords[L][d] + coords[R][d]);
            Fa[d] = (coords[L][d] + coords[R][d] + coords[a][d]) / RealType(3);
            Fb[d] = (coords[L][d] + coords[R][d] + coords[b][d]) / RealType(3);
        }

        // Dual quad winds M -> Fa -> C -> Fb. Its area-normal is the SCS area.
        RealType A[3];
        quad_area_normal_tet<RealType>(M, Fa, C, Fb, A);

        // Force the L->R orientation the NS kernels expect: A must have a
        // positive component along (xR - xL).
        RealType dirx = coords[R][0] - coords[L][0];
        RealType diry = coords[R][1] - coords[L][1];
        RealType dirz = coords[R][2] - coords[L][2];
        RealType dot  = A[0]*dirx + A[1]*diry + A[2]*dirz;
        if (dot < RealType(0)) { A[0] = -A[0]; A[1] = -A[1]; A[2] = -A[2]; }

        size_t off = e * nscs + ip;
        d_areaVec_x[off] = A[0];
        d_areaVec_y[off] = A[1];
        d_areaVec_z[off] = A[2];
    }
}

template<typename KeyType, typename RealType>
void precomputeTetAreaVectorsGpu(
    const KeyType* d_conn0, const KeyType* d_conn1,
    const KeyType* d_conn2, const KeyType* d_conn3,
    size_t numElements,
    const RealType* d_x, const RealType* d_y, const RealType* d_z,
    RealType* d_areaVec_x, RealType* d_areaVec_y, RealType* d_areaVec_z,
    cudaStream_t stream = 0)
{
    int blockSize = 256;
    int numBlocks = (numElements + blockSize - 1) / blockSize;
    precomputeTetAreaVectorsKernel<KeyType, RealType><<<numBlocks, blockSize, 0, stream>>>(
        d_conn0, d_conn1, d_conn2, d_conn3,
        numElements,
        d_x, d_y, d_z,
        d_areaVec_x, d_areaVec_y, d_areaVec_z);
}

} // namespace fem
} // namespace mars
