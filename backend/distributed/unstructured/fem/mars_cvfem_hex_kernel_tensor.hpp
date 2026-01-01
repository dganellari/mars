#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_common.hpp"

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

    // ==========================================================================
    // Optimized memory access: Software prefetching + loop unrolling
    // ==========================================================================
    // Strategy: Unroll loop by 4 to expose more ILP and help compiler optimize
    // memory accesses. The compiler can better schedule loads and hide latency.

    // Process nodes 0-3
    {
        KeyType n0 = nodes[0], n1 = nodes[1], n2 = nodes[2], n3 = nodes[3];

        // Coordinates - all loads issued together for better scheduling
        coords[0][0] = d_x[n0]; coords[0][1] = d_y[n0]; coords[0][2] = d_z[n0];
        coords[1][0] = d_x[n1]; coords[1][1] = d_y[n1]; coords[1][2] = d_z[n1];
        coords[2][0] = d_x[n2]; coords[2][1] = d_y[n2]; coords[2][2] = d_z[n2];
        coords[3][0] = d_x[n3]; coords[3][1] = d_y[n3]; coords[3][2] = d_z[n3];

        // Scalar fields
        phi[0] = d_phi[n0]; phi[1] = d_phi[n1]; phi[2] = d_phi[n2]; phi[3] = d_phi[n3];
        gamma[0] = d_gamma[n0]; gamma[1] = d_gamma[n1]; gamma[2] = d_gamma[n2]; gamma[3] = d_gamma[n3];
        beta[0] = d_beta[n0]; beta[1] = d_beta[n1]; beta[2] = d_beta[n2]; beta[3] = d_beta[n3];

        // Gradients
        grad_phi[0][0] = d_grad_phi_x[n0]; grad_phi[0][1] = d_grad_phi_y[n0]; grad_phi[0][2] = d_grad_phi_z[n0];
        grad_phi[1][0] = d_grad_phi_x[n1]; grad_phi[1][1] = d_grad_phi_y[n1]; grad_phi[1][2] = d_grad_phi_z[n1];
        grad_phi[2][0] = d_grad_phi_x[n2]; grad_phi[2][1] = d_grad_phi_y[n2]; grad_phi[2][2] = d_grad_phi_z[n2];
        grad_phi[3][0] = d_grad_phi_x[n3]; grad_phi[3][1] = d_grad_phi_y[n3]; grad_phi[3][2] = d_grad_phi_z[n3];

        // Metadata
        dofs[0] = d_node_to_dof[n0]; dofs[1] = d_node_to_dof[n1];
        dofs[2] = d_node_to_dof[n2]; dofs[3] = d_node_to_dof[n3];
        own[0] = d_ownership[n0]; own[1] = d_ownership[n1];
        own[2] = d_ownership[n2]; own[3] = d_ownership[n3];
    }

    // Process nodes 4-7
    {
        KeyType n4 = nodes[4], n5 = nodes[5], n6 = nodes[6], n7 = nodes[7];

        coords[4][0] = d_x[n4]; coords[4][1] = d_y[n4]; coords[4][2] = d_z[n4];
        coords[5][0] = d_x[n5]; coords[5][1] = d_y[n5]; coords[5][2] = d_z[n5];
        coords[6][0] = d_x[n6]; coords[6][1] = d_y[n6]; coords[6][2] = d_z[n6];
        coords[7][0] = d_x[n7]; coords[7][1] = d_y[n7]; coords[7][2] = d_z[n7];

        phi[4] = d_phi[n4]; phi[5] = d_phi[n5]; phi[6] = d_phi[n6]; phi[7] = d_phi[n7];
        gamma[4] = d_gamma[n4]; gamma[5] = d_gamma[n5]; gamma[6] = d_gamma[n6]; gamma[7] = d_gamma[n7];
        beta[4] = d_beta[n4]; beta[5] = d_beta[n5]; beta[6] = d_beta[n6]; beta[7] = d_beta[n7];

        grad_phi[4][0] = d_grad_phi_x[n4]; grad_phi[4][1] = d_grad_phi_y[n4]; grad_phi[4][2] = d_grad_phi_z[n4];
        grad_phi[5][0] = d_grad_phi_x[n5]; grad_phi[5][1] = d_grad_phi_y[n5]; grad_phi[5][2] = d_grad_phi_z[n5];
        grad_phi[6][0] = d_grad_phi_x[n6]; grad_phi[6][1] = d_grad_phi_y[n6]; grad_phi[6][2] = d_grad_phi_z[n6];
        grad_phi[7][0] = d_grad_phi_x[n7]; grad_phi[7][1] = d_grad_phi_y[n7]; grad_phi[7][2] = d_grad_phi_z[n7];

        dofs[4] = d_node_to_dof[n4]; dofs[5] = d_node_to_dof[n5];
        dofs[6] = d_node_to_dof[n6]; dofs[7] = d_node_to_dof[n7];
        own[4] = d_ownership[n4]; own[5] = d_ownership[n5];
        own[6] = d_ownership[n6]; own[7] = d_ownership[n7];
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
