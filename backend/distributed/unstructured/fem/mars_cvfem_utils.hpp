#pragma once

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>

namespace mars {
namespace fem {

// GPU kernel: Build DOF mapping from ownership array
// ownership: 0=ghost, 1=owned, 2=shared
// nodeToDof: output mapping (-1 for ghosts)
template<typename KeyType>
__global__ void buildDofMappingKernel(
    const uint8_t* __restrict__ ownership,
    const int* __restrict__ prefixSum,
    int* __restrict__ nodeToDof,
    size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    if (ownership[i] != 0) {
        nodeToDof[i] = prefixSum[i];
    } else {
        nodeToDof[i] = -1;
    }
}

// Build DOF mapping on GPU
template<typename KeyType>
int buildDofMappingGpu(
    const uint8_t* d_ownership,
    int* d_nodeToDof,
    size_t numNodes,
    cudaStream_t stream = 0)
{
    // Create temporary array for prefix sum
    thrust::device_vector<int> d_isOwned(numNodes);
    thrust::device_vector<int> d_prefixSum(numNodes);

    // Mark owned nodes (ownership != 0)
    thrust::transform(
        thrust::device,
        thrust::device_pointer_cast(d_ownership),
        thrust::device_pointer_cast(d_ownership + numNodes),
        d_isOwned.begin(),
        [] __device__ (uint8_t o) { return (o != 0) ? 1 : 0; }
    );

    // Exclusive prefix sum to get DOF indices
    thrust::exclusive_scan(thrust::device, d_isOwned.begin(), d_isOwned.end(), d_prefixSum.begin());

    // Get total number of DOFs
    int lastOwned, lastPrefix;
    cudaMemcpy(&lastOwned, thrust::raw_pointer_cast(d_isOwned.data() + numNodes - 1), sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&lastPrefix, thrust::raw_pointer_cast(d_prefixSum.data() + numNodes - 1), sizeof(int), cudaMemcpyDeviceToHost);
    int numDofs = lastPrefix + lastOwned;

    // Launch kernel to build mapping
    int blockSize = 256;
    int numBlocks = (numNodes + blockSize - 1) / blockSize;
    buildDofMappingKernel<KeyType><<<numBlocks, blockSize, 0, stream>>>(
        d_ownership,
        thrust::raw_pointer_cast(d_prefixSum.data()),
        d_nodeToDof,
        numNodes
    );

    return numDofs;
}

// GPU kernel: Initialize CVFEM fields from coordinates
template<typename RealType>
__global__ void initFieldsKernel(
    const RealType* __restrict__ x,
    const RealType* __restrict__ y,
    const RealType* __restrict__ z,
    RealType* __restrict__ phi,
    RealType* __restrict__ gamma,
    RealType* __restrict__ beta,
    RealType* __restrict__ grad_phi_x,
    RealType* __restrict__ grad_phi_y,
    RealType* __restrict__ grad_phi_z,
    size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    RealType xi = x[i];
    RealType yi = y[i];
    RealType zi = z[i];

    // phi = sin(x) + 3*cos(y) + 4*sin(5*x*y*z)
    RealType xyz5 = 5.0 * xi * yi * zi;
    phi[i] = sin(xi) + 3.0 * cos(yi) + 4.0 * sin(xyz5);

    // grad_phi
    RealType cos5xyz = cos(xyz5);
    grad_phi_x[i] = cos(xi) + 20.0 * yi * zi * cos5xyz;
    grad_phi_y[i] = 20.0 * xi * zi * cos5xyz - 3.0 * sin(yi);
    grad_phi_z[i] = 20.0 * xi * yi * cos5xyz;

    // Constants
    gamma[i] = 0.1;
    beta[i] = 1.234;
}

// Initialize fields on GPU
template<typename RealType>
void initFieldsGpu(
    const RealType* d_x,
    const RealType* d_y,
    const RealType* d_z,
    RealType* d_phi,
    RealType* d_gamma,
    RealType* d_beta,
    RealType* d_grad_phi_x,
    RealType* d_grad_phi_y,
    RealType* d_grad_phi_z,
    size_t numNodes,
    cudaStream_t stream = 0)
{
    int blockSize = 256;
    int numBlocks = (numNodes + blockSize - 1) / blockSize;
    initFieldsKernel<RealType><<<numBlocks, blockSize, 0, stream>>>(
        d_x, d_y, d_z,
        d_phi, d_gamma, d_beta,
        d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
        numNodes
    );
}

// hexLRSCV constant for sparsity building
__device__ __constant__ int d_hexLRSCV[24] = {
    0, 1, 1, 2, 2, 3, 0, 3,
    4, 5, 5, 6, 6, 7, 4, 7,
    0, 4, 1, 5, 2, 6, 3, 7
};

// GPU kernel: Count NNZ per row for reduced sparsity pattern
template<typename KeyType>
__global__ void countNnzPerRowKernel(
    const KeyType* __restrict__ conn0,
    const KeyType* __restrict__ conn1,
    const KeyType* __restrict__ conn2,
    const KeyType* __restrict__ conn3,
    const KeyType* __restrict__ conn4,
    const KeyType* __restrict__ conn5,
    const KeyType* __restrict__ conn6,
    const KeyType* __restrict__ conn7,
    const int* __restrict__ nodeToDof,
    int* __restrict__ nnzPerRow,
    size_t numElements,
    int numDofs)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType nodes[8];
    nodes[0] = conn0[e];
    nodes[1] = conn1[e];
    nodes[2] = conn2[e];
    nodes[3] = conn3[e];
    nodes[4] = conn4[e];
    nodes[5] = conn5[e];
    nodes[6] = conn6[e];
    nodes[7] = conn7[e];

    // For each SCS, mark adjacency
    for (int scs = 0; scs < 12; ++scs) {
        int nodeL = d_hexLRSCV[scs * 2];
        int nodeR = d_hexLRSCV[scs * 2 + 1];

        int dofL = nodeToDof[nodes[nodeL]];
        int dofR = nodeToDof[nodes[nodeR]];

        // Atomic increment for each unique edge (simplified - may overcount)
        if (dofL >= 0 && dofR >= 0 && dofL != dofR) {
            atomicAdd(&nnzPerRow[dofL], 1);
            atomicAdd(&nnzPerRow[dofR], 1);
        }
    }
}

// =============================================================================
// Area vector pre-computation kernel
// =============================================================================
// Pre-computes area vectors for all SCS (12 per element) to avoid redundant
// geometry computation in the assembly kernel.

// SCS face quad node indices (into 27-point subdivision)
__device__ __constant__ int d_hexEdgeFacetTable[12][4] = {
    {20, 8, 12, 26},  {24, 9, 12, 26},  {10, 12, 26, 23}, {11, 25, 26, 12},
    {13, 20, 26, 17}, {17, 14, 24, 26}, {17, 15, 23, 26}, {16, 17, 26, 25},
    {19, 20, 26, 25}, {20, 18, 24, 26}, {22, 23, 26, 24}, {21, 25, 26, 23}
};

// Subdivide hex into 27 points for SCS geometry
// MUST match mars_cvfem_hex_kernel.hpp::subdivide_hex_8 exactly!
template<typename RealType>
__device__ inline void subdivide_hex_8_generic(const RealType coords[8][3], RealType coordv[27][3])
{
    // Copy 8 corner nodes
    for (int n = 0; n < 8; ++n) {
        for (int d = 0; d < 3; ++d) {
            coordv[n][d] = coords[n][d];
        }
    }

    for (int d = 0; d < 3; ++d) {
        // 12 edge midpoints (indices 8-11, 13-16, 18-19, 21-22)
        coordv[8][d] = RealType(0.5) * (coords[0][d] + coords[1][d]);   // edge 1
        coordv[9][d] = RealType(0.5) * (coords[1][d] + coords[2][d]);   // edge 2
        coordv[10][d] = RealType(0.5) * (coords[2][d] + coords[3][d]);  // edge 3
        coordv[11][d] = RealType(0.5) * (coords[3][d] + coords[0][d]);  // edge 4
        coordv[13][d] = RealType(0.5) * (coords[4][d] + coords[5][d]);  // edge 5
        coordv[14][d] = RealType(0.5) * (coords[5][d] + coords[6][d]);  // edge 6
        coordv[15][d] = RealType(0.5) * (coords[6][d] + coords[7][d]);  // edge 7
        coordv[16][d] = RealType(0.5) * (coords[7][d] + coords[4][d]);  // edge 8
        coordv[18][d] = RealType(0.5) * (coords[1][d] + coords[5][d]);  // edge 9
        coordv[19][d] = RealType(0.5) * (coords[0][d] + coords[4][d]);  // edge 10
        coordv[21][d] = RealType(0.5) * (coords[3][d] + coords[7][d]);  // edge 11
        coordv[22][d] = RealType(0.5) * (coords[2][d] + coords[6][d]);  // edge 12

        // 6 face centers (indices 12, 17, 20, 23, 24, 25)
        coordv[12][d] = RealType(0.25) * (coords[0][d] + coords[1][d] + coords[2][d] + coords[3][d]); // face 0
        coordv[17][d] = RealType(0.25) * (coords[4][d] + coords[5][d] + coords[6][d] + coords[7][d]); // face 1
        coordv[20][d] = RealType(0.25) * (coords[0][d] + coords[1][d] + coords[4][d] + coords[5][d]); // face 2
        coordv[23][d] = RealType(0.25) * (coords[2][d] + coords[3][d] + coords[6][d] + coords[7][d]); // face 3
        coordv[24][d] = RealType(0.25) * (coords[1][d] + coords[2][d] + coords[5][d] + coords[6][d]); // face 4
        coordv[25][d] = RealType(0.25) * (coords[0][d] + coords[3][d] + coords[4][d] + coords[7][d]); // face 5

        // Volume centroid (index 26)
        coordv[26][d] = RealType(0.0);
        for (int n = 0; n < 8; ++n) {
            coordv[26][d] += coords[n][d];
        }
        coordv[26][d] *= RealType(0.125);
    }
}

// Compute area vector by quad triangulation (Grandy algorithm)
template<typename RealType>
__device__ inline void quad_area_by_triangulation_generic(const RealType areacoords[4][3], RealType areaVec[3])
{
    areaVec[0] = 0.0;
    areaVec[1] = 0.0;
    areaVec[2] = 0.0;

    // Quad centroid
    RealType xmid[3] = {
        RealType(0.25) * (areacoords[0][0] + areacoords[1][0] + areacoords[2][0] + areacoords[3][0]),
        RealType(0.25) * (areacoords[0][1] + areacoords[1][1] + areacoords[2][1] + areacoords[3][1]),
        RealType(0.25) * (areacoords[0][2] + areacoords[1][2] + areacoords[2][2] + areacoords[3][2])
    };

    RealType r1[3] = {
        areacoords[0][0] - xmid[0],
        areacoords[0][1] - xmid[1],
        areacoords[0][2] - xmid[2]
    };

    // Sum cross products for 4 triangles
    for (int itri = 0; itri < 4; ++itri) {
        int t_index = (itri + 1) % 4;
        RealType r2[3] = {
            areacoords[t_index][0] - xmid[0],
            areacoords[t_index][1] - xmid[1],
            areacoords[t_index][2] - xmid[2]
        };

        // Cross product r1 x r2
        areaVec[0] += r1[1] * r2[2] - r2[1] * r1[2];
        areaVec[1] += r1[2] * r2[0] - r2[2] * r1[0];
        areaVec[2] += r1[0] * r2[1] - r2[0] * r1[1];

        r1[0] = r2[0];
        r1[1] = r2[1];
        r1[2] = r2[2];
    }

    areaVec[0] *= RealType(0.5);
    areaVec[1] *= RealType(0.5);
    areaVec[2] *= RealType(0.5);
}

// Pre-compute area vectors for all elements and all SCS
template<typename KeyType, typename RealType>
__global__ void precomputeAreaVectorsKernel(
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    const KeyType* __restrict__ d_conn4,
    const KeyType* __restrict__ d_conn5,
    const KeyType* __restrict__ d_conn6,
    const KeyType* __restrict__ d_conn7,
    size_t numElements,
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    RealType* __restrict__ d_areaVec_x,
    RealType* __restrict__ d_areaVec_y,
    RealType* __restrict__ d_areaVec_z)
{
    size_t elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Load element connectivity
    KeyType nodes[8];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];
    nodes[4] = d_conn4[elemIdx];
    nodes[5] = d_conn5[elemIdx];
    nodes[6] = d_conn6[elemIdx];
    nodes[7] = d_conn7[elemIdx];

    // Load coordinates
    RealType coords[8][3];
    for (int n = 0; n < 8; ++n) {
        KeyType node = nodes[n];
        coords[n][0] = d_x[node];
        coords[n][1] = d_y[node];
        coords[n][2] = d_z[node];
    }

    // Subdivide hex into 27 points
    RealType coordv[27][3];
    subdivide_hex_8_generic(coords, coordv);

    // Compute area vector for each SCS
    for (int scs = 0; scs < 12; ++scs) {
        // Get the 4 quad nodes for this SCS
        RealType scscoords[4][3];
        for (int inode = 0; inode < 4; ++inode) {
            int idx = d_hexEdgeFacetTable[scs][inode];
            for (int d = 0; d < 3; ++d) {
                scscoords[inode][d] = coordv[idx][d];
            }
        }

        // Compute area vector
        RealType areaVec[3];
        quad_area_by_triangulation_generic(scscoords, areaVec);

        // Store
        size_t offset = elemIdx * 12 + scs;
        d_areaVec_x[offset] = areaVec[0];
        d_areaVec_y[offset] = areaVec[1];
        d_areaVec_z[offset] = areaVec[2];
    }
}

// Host function to pre-compute area vectors
template<typename KeyType, typename RealType>
void precomputeAreaVectorsGpu(
    const KeyType* d_conn0, const KeyType* d_conn1,
    const KeyType* d_conn2, const KeyType* d_conn3,
    const KeyType* d_conn4, const KeyType* d_conn5,
    const KeyType* d_conn6, const KeyType* d_conn7,
    size_t numElements,
    const RealType* d_x, const RealType* d_y, const RealType* d_z,
    RealType* d_areaVec_x, RealType* d_areaVec_y, RealType* d_areaVec_z,
    cudaStream_t stream = 0)
{
    int blockSize = 256;
    int numBlocks = (numElements + blockSize - 1) / blockSize;
    precomputeAreaVectorsKernel<KeyType, RealType><<<numBlocks, blockSize, 0, stream>>>(
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_conn4, d_conn5, d_conn6, d_conn7,
        numElements,
        d_x, d_y, d_z,
        d_areaVec_x, d_areaVec_y, d_areaVec_z
    );
}

} // namespace fem
} // namespace mars
