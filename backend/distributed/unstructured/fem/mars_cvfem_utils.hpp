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

} // namespace fem
} // namespace mars
