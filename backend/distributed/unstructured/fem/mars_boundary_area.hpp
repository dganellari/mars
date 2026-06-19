#pragma once

// Per-node boundary median-dual area vectors S_bnd for the collocated Rhie-Chow continuity.
//
// The interior continuity a^pu (assembleRhieChowDivGradKernel) is a divergence over the INTERIOR
// median-dual SCS faces only. At an OPEN boundary node the control volume's exterior facet on dOmega
// is missing, so the boundary mass flux rho*(u.n)*|S_bnd| is dropped from div(u). That term is
// identically zero on a closed domain (walls/lid: u.n=0) but nonzero at an inlet/outlet -- without it
// the inlet flux has nowhere to go and the interior velocity collapses (host-validated in
// scripts/rhie_chow_channel.py). S_bnd closes the open dual cells:
//   S_bnd_P = sum over boundary triangles touching P of (1/3) * outward-area-vector
// A boundary triangle is a tet face that belongs to exactly ONE tet. Its outward area vector is
// 0.5*(e1 x e2) flipped to point AWAY from the tet's opposite (4th) node. Each of the 3 vertices gets
// 1/3 of it (the median-dual split of a triangle is three equal-area parts).
//
// Boundary-face detection is a fresh thrust dedup (4 sorted faces/tet -> sort -> a face is boundary iff
// it is unique in the sorted list). NOT FaceTopology / getBoundaryNodes (both unreliable for tets here).
// Single-rank: every owned element is local. Multi-rank fold (owner-once + reverseExchangeNodeHaloAdd)
// is deferred -- the model-operator driver is single-rank.

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

namespace mars {
namespace fem {

// Tet face f is OPPOSITE local node f; these are its 3 face-local node indices.
__device__ __constant__ int d_tetBndFaceNodes[4][3] = {
    {1, 2, 3},  // face opposite node 0
    {0, 2, 3},  // face opposite node 1
    {0, 1, 3},  // face opposite node 2
    {0, 1, 2}   // face opposite node 3
};

// Sorted global-node triple (so shared faces collide) + the owning (element,localFace) payload.
// Assumes a CONFORMING manifold mesh: each face is in 1 (boundary) or 2 (interior) tets. A
// non-manifold face (>2 tets) would be misclassified interior by the neighbor-diff test below.
struct BndFace {
    unsigned long long k0, k1, k2;
    unsigned long long payload;     // localElemIdx*4 + localFace
};
__host__ __device__ inline bool operator<(const BndFace& a, const BndFace& b) {
    if (a.k0 != b.k0) return a.k0 < b.k0;
    if (a.k1 != b.k1) return a.k1 < b.k1;
    return a.k2 < b.k2;
}
__host__ __device__ inline bool sameKey(const BndFace& a, const BndFace& b) {
    return a.k0 == b.k0 && a.k1 == b.k1 && a.k2 == b.k2;
}

template<typename KeyType>
__global__ void buildTetFacesKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2,
                                     const KeyType* c3, size_t startElem, size_t numLocal, BndFace* faces)
{
    size_t kk = blockIdx.x * blockDim.x + threadIdx.x;
    if (kk >= numLocal) return;
    size_t e = startElem + kk;
    KeyType n[4] = {c0[e], c1[e], c2[e], c3[e]};
    for (int f = 0; f < 4; ++f) {
        unsigned long long a = (unsigned long long)n[d_tetBndFaceNodes[f][0]];
        unsigned long long b = (unsigned long long)n[d_tetBndFaceNodes[f][1]];
        unsigned long long c = (unsigned long long)n[d_tetBndFaceNodes[f][2]];
        if (a > b) { unsigned long long t = a; a = b; b = t; }   // sort the triple
        if (b > c) { unsigned long long t = b; b = c; c = t; }
        if (a > b) { unsigned long long t = a; a = b; b = t; }
        BndFace& bf = faces[kk * 4 + f];
        bf.k0 = a; bf.k1 = b; bf.k2 = c;
        bf.payload = (unsigned long long)kk * 4 + (unsigned long long)f;
    }
}

template<typename KeyType, typename RealType>
__global__ void scatterBoundaryAreaKernel(
    const BndFace* faces, size_t nFaces,
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3, size_t startElem,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    RealType* sBndX, RealType* sBndY, RealType* sBndZ)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nFaces) return;
    // boundary iff this sorted face key differs from BOTH neighbors (appears exactly once)
    bool bnd = (i == 0 || !sameKey(faces[i], faces[i - 1])) &&
               (i == nFaces - 1 || !sameKey(faces[i], faces[i + 1]));
    if (!bnd) return;

    unsigned long long p = faces[i].payload;
    size_t kk = p / 4; int f = (int)(p % 4);
    size_t e = startElem + kk;
    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n4[4]; for (int k = 0; k < 4; ++k) n4[k] = (int)cc[k][e];
    int fn[3] = {n4[d_tetBndFaceNodes[f][0]], n4[d_tetBndFaceNodes[f][1]], n4[d_tetBndFaceNodes[f][2]]};
    int opp = n4[f];                                            // face f is opposite local node f

    RealType X[3][3];
    for (int k = 0; k < 3; ++k) { X[k][0] = nodeX[fn[k]]; X[k][1] = nodeY[fn[k]]; X[k][2] = nodeZ[fn[k]]; }
    RealType e1[3] = {X[1][0]-X[0][0], X[1][1]-X[0][1], X[1][2]-X[0][2]};
    RealType e2[3] = {X[2][0]-X[0][0], X[2][1]-X[0][1], X[2][2]-X[0][2]};
    RealType A[3] = {RealType(0.5)*(e1[1]*e2[2]-e1[2]*e2[1]),
                     RealType(0.5)*(e1[2]*e2[0]-e1[0]*e2[2]),
                     RealType(0.5)*(e1[0]*e2[1]-e1[1]*e2[0])};
    // outward = away from the opposite node
    RealType fcx = (X[0][0]+X[1][0]+X[2][0]) / RealType(3);
    RealType fcy = (X[0][1]+X[1][1]+X[2][1]) / RealType(3);
    RealType fcz = (X[0][2]+X[1][2]+X[2][2]) / RealType(3);
    RealType dot = A[0]*(fcx-nodeX[opp]) + A[1]*(fcy-nodeY[opp]) + A[2]*(fcz-nodeZ[opp]);
    if (dot < RealType(0)) { A[0] = -A[0]; A[1] = -A[1]; A[2] = -A[2]; }

    RealType third = RealType(1) / RealType(3);
    for (int k = 0; k < 3; ++k) {
        atomicAdd(&sBndX[fn[k]], A[0] * third);
        atomicAdd(&sBndY[fn[k]], A[1] * third);
        atomicAdd(&sBndZ[fn[k]], A[2] * third);
    }
}

// Fill d_sBnd{X,Y,Z} (nNodes, caller need not pre-zero -- we zero here). Single rank: startElem=0.
template<typename KeyType, typename RealType>
void computeBoundaryAreaVectorsGpu(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    size_t startElem, size_t numLocal,
    const RealType* d_x, const RealType* d_y, const RealType* d_z, int nNodes,
    RealType* d_sBndX, RealType* d_sBndY, RealType* d_sBndZ, cudaStream_t stream = 0)
{
    const int blockSize = 256;
    const size_t nFaces = 4 * numLocal;
    thrust::device_vector<BndFace> faces(nFaces);
    int eBlocks = (int)((numLocal + blockSize - 1) / blockSize);
    buildTetFacesKernel<KeyType><<<eBlocks, blockSize, 0, stream>>>(
        c0, c1, c2, c3, startElem, numLocal, thrust::raw_pointer_cast(faces.data()));
    thrust::sort(thrust::cuda::par.on(stream), faces.begin(), faces.end());
    cudaMemsetAsync(d_sBndX, 0, nNodes * sizeof(RealType), stream);
    cudaMemsetAsync(d_sBndY, 0, nNodes * sizeof(RealType), stream);
    cudaMemsetAsync(d_sBndZ, 0, nNodes * sizeof(RealType), stream);
    int fBlocks = (int)((nFaces + blockSize - 1) / blockSize);
    scatterBoundaryAreaKernel<KeyType, RealType><<<fBlocks, blockSize, 0, stream>>>(
        thrust::raw_pointer_cast(faces.data()), nFaces, c0, c1, c2, c3, startElem,
        d_x, d_y, d_z, d_sBndX, d_sBndY, d_sBndZ);
}

} // namespace fem
} // namespace mars
