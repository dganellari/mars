#pragma once

#include <cuda_runtime.h>

namespace mars {
namespace fem {

// =============================================================================
// AoS node data layout for CVFEM assembly.
//
// WHY AoS for nodes?
//   SoA wins when all threads access the SAME field for CONSECUTIVE elements:
//     d_x[e], d_x[e+1], ... → coalesced (1 sector per 4 threads for doubles)
//
//   AoS wins when each thread accesses ALL fields for ONE randomly-indexed node:
//     thread i reads phi, gamma, beta, x, y, z, gx, gy, gz at node nodes[n]
//     → with SoA: 9 arrays × 1 scattered sector each = 9 sectors per thread
//     → with AoS: 1 struct × 3 sectors = 3 sectors per thread  (3x less L2 traffic)
//
//   CVFEM node reads are the scatter/gather pattern → AoS wins.
//   Element-level reads (connectivity, mdot, areaVec) remain SoA (sequential).
//
// Layout: 9 doubles = 72 bytes.  Aligned to 32 bytes (= 1 L1 sector) so nodes
// never split across more sectors than necessary (ceil(72/32) = 3 sectors/node).
// =============================================================================
struct alignas(32) NodeData {
    double x, y, z;
    double phi, gamma, beta;
    double gx, gy, gz;     // grad_phi_x, grad_phi_y, grad_phi_z
};
// sizeof(NodeData) = 72 bytes = 3 L1 sectors.
// With AoS: 32 threads × random nodes → 96 sectors, 75% utilization.
// With SoA: 32 threads × 9 fields  → 288 sectors, 25% utilization.

// =============================================================================
// GPU kernel: pack 9 SoA arrays → NodeData AoS array (one thread per node).
// =============================================================================
template<typename RealType>
__global__ void packNodeDataKernel(
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    const RealType* __restrict__ d_phi,
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_beta,
    const RealType* __restrict__ d_gx,
    const RealType* __restrict__ d_gy,
    const RealType* __restrict__ d_gz,
    NodeData* __restrict__ d_nodeData,
    size_t numNodes)
{
    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    d_nodeData[i].x     = d_x[i];
    d_nodeData[i].y     = d_y[i];
    d_nodeData[i].z     = d_z[i];
    d_nodeData[i].phi   = d_phi[i];
    d_nodeData[i].gamma = d_gamma[i];
    d_nodeData[i].beta  = d_beta[i];
    d_nodeData[i].gx    = d_gx[i];
    d_nodeData[i].gy    = d_gy[i];
    d_nodeData[i].gz    = d_gz[i];
}

// Host wrapper — call after field data is ready on device.
// Returns allocated d_nodeData pointer (caller owns, must cudaFree).
template<typename RealType>
NodeData* packNodeData(
    const RealType* d_x,  const RealType* d_y,  const RealType* d_z,
    const RealType* d_phi, const RealType* d_gamma, const RealType* d_beta,
    const RealType* d_gx, const RealType* d_gy, const RealType* d_gz,
    size_t numNodes,
    cudaStream_t stream = 0)
{
    NodeData* d_nodeData = nullptr;
    cudaMalloc(&d_nodeData, numNodes * sizeof(NodeData));
    const int threads = 256;
    const int blocks  = static_cast<int>((numNodes + threads - 1) / threads);
    packNodeDataKernel<RealType><<<blocks, threads, 0, stream>>>(
        d_x, d_y, d_z, d_phi, d_gamma, d_beta, d_gx, d_gy, d_gz,
        d_nodeData, numNodes);
    return d_nodeData;
}

} // namespace fem
} // namespace mars
