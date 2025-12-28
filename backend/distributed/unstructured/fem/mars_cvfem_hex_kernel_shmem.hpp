#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

// =============================================================================
// Optimized CVFEM assembly kernel exploiting GH100 memory hierarchy
// =============================================================================
// Strategy: Thread-per-element for SCS integration, warp-cooperative assembly
//
// Phase 1: Each thread computes its element's local matrix (register-based)
// Phase 2: Warp-cooperative CSR column search using warp shuffle
// Phase 3: Coalesced global memory writes
// =============================================================================

// Warp-parallel column search: all 32 lanes search in parallel
// Returns the position of 'target' in colInd[start:end], or -1 if not found
__device__ __forceinline__ int warpParallelSearch(
    const int* __restrict__ colInd,
    int start, int end, int target)
{
    int laneId = threadIdx.x % 32;
    int len = end - start;

    // Each lane checks different positions
    int found_pos = -1;
    for (int base = 0; base < len; base += 32) {
        int idx = base + laneId;
        if (idx < len) {
            if (colInd[start + idx] == target) {
                found_pos = start + idx;
            }
        }
        // Reduce across warp to find if any lane found it
        unsigned mask = __ballot_sync(0xffffffff, found_pos >= 0);
        if (mask != 0) {
            // Someone found it - get the position from that lane
            int srcLane = __ffs(mask) - 1;
            found_pos = __shfl_sync(0xffffffff, found_pos, srcLane);
            return found_pos;
        }
    }
    return -1;
}

template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void cvfem_hex_assembly_kernel_shmem(
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
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_phi,
    const RealType* __restrict__ d_beta,
    const RealType* __restrict__ d_grad_phi_x,
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    const RealType* __restrict__ d_mdot,
    const RealType* __restrict__ d_areaVec_x,
    const RealType* __restrict__ d_areaVec_y,
    const RealType* __restrict__ d_areaVec_z,
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // ==========================================================================
    // Phase 1: Load data and compute local element matrix (all in registers)
    // ==========================================================================

    KeyType nodes[8];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];
    nodes[4] = d_conn4[elemIdx];
    nodes[5] = d_conn5[elemIdx];
    nodes[6] = d_conn6[elemIdx];
    nodes[7] = d_conn7[elemIdx];

    RealType coords[8][3];
    RealType phi[8], gamma[8], beta[8];
    RealType grad_phi[8][3];
    int dofs[8];
    uint8_t ownership[8];

    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        KeyType node = nodes[n];
        coords[n][0] = d_x[node];
        coords[n][1] = d_y[node];
        coords[n][2] = d_z[node];
        phi[n] = d_phi[node];
        gamma[n] = d_gamma[node];
        beta[n] = d_beta[node];
        grad_phi[n][0] = d_grad_phi_x[node];
        grad_phi[n][1] = d_grad_phi_y[node];
        grad_phi[n][2] = d_grad_phi_z[node];
        dofs[n] = d_node_to_dof[node];
        ownership[n] = d_ownership[node];
    }

    // Local element matrix and RHS (in registers)
    RealType lhs[64] = {0.0};
    RealType rhs[8] = {0.0};

    // SCS integration loop
    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        RealType mdot = d_mdot[elemIdx * 12 + ip];

        RealType areaVec[3];
        computeAreaVector(ip, coords, areaVec);

        RealType phi_ip = 0.0, gamma_ip = 0.0;
        RealType coords_ip[3] = {0.0, 0.0, 0.0};

        #pragma unroll
        for (int node = 0; node < 8; ++node) {
            RealType sf = hexShapeFcn[ip][node];
            phi_ip += sf * phi[node];
            gamma_ip += sf * gamma[node];
            coords_ip[0] += sf * coords[node][0];
            coords_ip[1] += sf * coords[node][1];
            coords_ip[2] += sf * coords[node][2];
        }

        RealType phi_upwind, beta_upwind;
        RealType dcorr = 0.0;

        if (mdot > 0.0) {
            phi_upwind = phi[nodeL];
            beta_upwind = beta[nodeL];
            dcorr = grad_phi[nodeL][0] * (coords_ip[0] - coords[nodeL][0]) +
                    grad_phi[nodeL][1] * (coords_ip[1] - coords[nodeL][1]) +
                    grad_phi[nodeL][2] * (coords_ip[2] - coords[nodeL][2]);
        } else {
            phi_upwind = phi[nodeR];
            beta_upwind = beta[nodeR];
            dcorr = grad_phi[nodeR][0] * (coords_ip[0] - coords[nodeR][0]) +
                    grad_phi[nodeR][1] * (coords_ip[1] - coords[nodeR][1]) +
                    grad_phi[nodeR][2] * (coords_ip[2] - coords[nodeR][2]);
        }
        dcorr *= beta_upwind;

        RealType adv_flux = mdot * (phi_upwind + dcorr);

        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        RealType lhsfac_L = 0.5 * (mdot + fabs(mdot));
        RealType lhsfac_R = 0.5 * (mdot - fabs(mdot));

        lhs[nodeL * 8 + nodeL] += lhsfac_L;
        lhs[nodeR * 8 + nodeL] -= lhsfac_L;
        lhs[nodeL * 8 + nodeR] += lhsfac_R;
        lhs[nodeR * 8 + nodeR] -= lhsfac_R;

        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int node = 0; node < 8; ++node) {
            RealType diff_coeff = -(gamma_ip * (dndx[node][0] * areaVec[0] +
                                                 dndx[node][1] * areaVec[1] +
                                                 dndx[node][2] * areaVec[2]));

            lhs[nodeL * 8 + node] += diff_coeff;
            lhs[nodeR * 8 + node] -= diff_coeff;

            rhs[nodeL] -= diff_coeff * phi[node];
            rhs[nodeR] += diff_coeff * phi[node];
        }
    }

    // ==========================================================================
    // Phase 2: Global assembly with single-pass vectorized column search
    // ==========================================================================

    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        int row_dof = dofs[i];
        if (ownership[i] == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        int row_start = matrix->rowPtr[row_dof];
        int row_end = matrix->rowPtr[row_dof + 1];
        int row_len = row_end - row_start;

        // Single-pass vectorized search: find all 8 column positions at once
        int positions[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

        for (int k = 0; k < row_len; ++k) {
            int col = matrix->colInd[row_start + k];
            #pragma unroll
            for (int j = 0; j < 8; ++j) {
                if (col == dofs[j]) {
                    positions[j] = row_start + k;
                }
            }
        }

        RealType diag_val = lhs[i * 8 + i];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            if (i == j) continue;
            if (dofs[j] < 0) continue;

            RealType value = lhs[i * 8 + j];

            if (positions[j] >= 0) {
                atomicAdd(&matrix->values[positions[j]], value);
            } else {
                diag_val += value;
            }
        }

        // Use found diagonal position
        if (positions[i] >= 0) {
            atomicAdd(&matrix->values[positions[i]], diag_val);
        } else {
            atomicAdd(matrix->diag(row_dof), diag_val);
        }
    }
}

} // namespace fem
} // namespace mars
