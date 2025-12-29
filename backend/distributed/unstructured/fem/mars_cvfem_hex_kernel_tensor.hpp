#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

// =============================================================================
// Full Matrix CVFEM assembly kernel
// =============================================================================
// Alternative approach: assemble full 8×8 local stiffness matrix then scatter
//
// Key differences from shmem kernel:
// - Assembles complete lhs[8][8] matrix locally (simpler logic)
// - Scatters all 64 entries to CSR with atomics
// - Uses 64 doubles for lhs[] vs 16 doubles (rhs[8] + diag[8]) in shmem
//
// Trade-offs:
// + Simpler assembly logic (no conditional atomic paths)
// + Better vectorization potential for matrix operations
// + Enables future tensor core optimizations
// - More register pressure (64 vs 16 accumulators)
// - More atomic operations (up to 64 vs ~32 per element)
//
// Best for: Architectures with abundant registers and fast atomics (GH200)
// =============================================================================

template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void cvfem_hex_assembly_kernel_tensor(
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
    // Load element data
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
    uint8_t own[8];

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
        own[n] = d_ownership[node];
    }

    // ==========================================================================
    // Full 8×8 local element matrix and RHS vector
    // ==========================================================================
    // Assemble full local matrix for better compute patterns
    // Trade-off: 64 doubles (48 more than shmem) but simpler assembly logic
    RealType lhs[64];
    RealType rhs[8];

    // Initialize to zero
    #pragma unroll
    for (int i = 0; i < 64; ++i) lhs[i] = 0.0;
    #pragma unroll
    for (int i = 0; i < 8; ++i) rhs[i] = 0.0;

    // ==========================================================================
    // SCS integration loop
    // ==========================================================================
    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        RealType mdot = d_mdot[elemIdx * 12 + ip];

        // Pre-computed area vectors
        RealType areaVec[3];
        areaVec[0] = d_areaVec_x[elemIdx * 12 + ip];
        areaVec[1] = d_areaVec_y[elemIdx * 12 + ip];
        areaVec[2] = d_areaVec_z[elemIdx * 12 + ip];

        // Interpolate to SCS
        RealType gamma_ip = 0.0;
        RealType coords_ip[3] = {0.0, 0.0, 0.0};

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType sf = hexShapeFcn[ip][n];
            gamma_ip += sf * gamma[n];
            coords_ip[0] += sf * coords[n][0];
            coords_ip[1] += sf * coords[n][1];
            coords_ip[2] += sf * coords[n][2];
        }

        // Advection
        RealType phi_upwind, beta_upwind, dcorr = 0.0;
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

        // Advection LHS contributions
        RealType lhsfac_L = 0.5 * (mdot + fabs(mdot));
        RealType lhsfac_R = 0.5 * (mdot - fabs(mdot));

        lhs[nodeL * 8 + nodeL] += lhsfac_L;
        lhs[nodeR * 8 + nodeL] -= lhsfac_L;
        lhs[nodeL * 8 + nodeR] += lhsfac_R;
        lhs[nodeR * 8 + nodeR] -= lhsfac_R;

        // Diffusion - compute shape derivatives on-the-fly
        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType diff_coeff = -(gamma_ip * (dndx[n][0] * areaVec[0] +
                                                 dndx[n][1] * areaVec[1] +
                                                 dndx[n][2] * areaVec[2]));

            // RHS contributions
            rhs[nodeL] -= diff_coeff * phi[n];
            rhs[nodeR] += diff_coeff * phi[n];

            // LHS contributions (full matrix)
            lhs[nodeL * 8 + n] += diff_coeff;
            lhs[nodeR * 8 + n] -= diff_coeff;
        }
    }

    // ==========================================================================
    // Potential tensor core optimization (future work)
    // ==========================================================================
    // WMMA API could accelerate certain matrix operations, but requires:
    // - Specific data layout (column-major fragments)
    // - Alignment constraints
    // - Warp-synchronous execution
    // Current scalar assembly is simpler and may be sufficient

    // ==========================================================================
    // Assemble into global CSR matrix
    // ==========================================================================
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        KeyType row_node = nodes[i];
        int row_dof = dofs[i];
        uint8_t ownership = own[i];

        if (ownership == 0 || row_dof < 0) continue;

        // Assemble RHS
        atomicAdd(&d_rhs[row_dof], rhs[i]);

        // Assemble LHS - find positions in CSR
        int row_start = matrix->rowPtr[row_dof];
        int row_end = matrix->rowPtr[row_dof + 1];
        int diag_pos = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            int col_dof = dofs[j];
            if (col_dof < 0) continue;

            RealType value = lhs[i * 8 + j];
            if (value == 0.0) continue; // Skip zeros

            if (i == j) {
                // Diagonal
                atomicAdd(&matrix->values[diag_pos], value);
            } else {
                // Off-diagonal - binary search
                int left = row_start;
                int right = row_end - 1;

                while (left <= right) {
                    int mid = (left + right) >> 1;
                    int col = matrix->colInd[mid];

                    if (col == col_dof) {
                        atomicAdd(&matrix->values[mid], value);
                        break;
                    } else if (col < col_dof) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                }
            }
        }
    }
}

} // namespace fem
} // namespace mars
