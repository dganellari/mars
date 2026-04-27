#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_common.hpp"

namespace mars {
namespace fem {

// =============================================================================
// Per-integration-point (perip) CVFEM assembly kernel
// =============================================================================
// Key difference from tensor kernel:
//   tensor:  accumulate lhs[64] over 12 SCS → scatter 64 values at end
//   perip:   pre-lookup all 64 CSR positions once → scatter immediately per SCS
//
// Register trade-off:
//   Remove lhs[64] doubles  = -128 registers
//   Add    pos[64] ints     =  +64 registers
//   Net                     =  -64 registers  (~172 → ~108-140 estimated)
//
// Additionally: binary search is done ONCE per element (pre-SCS phase),
// not inside the 12-iteration hot loop — removes scan overhead from SCS loop.
//
// Two variants:
//   cvfem_hex_assembly_kernel_tensor_perip      — no launch_bounds (natural regs)
//   cvfem_hex_assembly_kernel_tensor_perip_lb2  — __launch_bounds__(256,2)
//     Forces 2 blocks/SM. Spilling burden: ~7 regs vs ~44 regs for tensor+lb2.
// =============================================================================

// =============================================================================
// Shared kernel body (inlined into both variants via __device__ __forceinline__)
// =============================================================================
template<typename KeyType, typename RealType>
__device__ __forceinline__ void cvfem_hex_perip_body(
    int elemIdx,
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    const KeyType* __restrict__ d_conn4,
    const KeyType* __restrict__ d_conn5,
    const KeyType* __restrict__ d_conn6,
    const KeyType* __restrict__ d_conn7,
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
    // =========================================================================
    // Load connectivity
    // =========================================================================
    KeyType nodes[8];
    nodes[0] = d_conn0[elemIdx]; nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx]; nodes[3] = d_conn3[elemIdx];
    nodes[4] = d_conn4[elemIdx]; nodes[5] = d_conn5[elemIdx];
    nodes[6] = d_conn6[elemIdx]; nodes[7] = d_conn7[elemIdx];

    RealType coords[8][3];
    RealType phi[8], gamma[8], beta[8];
    RealType grad_phi[8][3];
    int      dofs[8];
    uint8_t  own[8];

    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        KeyType ni = nodes[n];
        coords[n][0]   = d_x[ni];
        coords[n][1]   = d_y[ni];
        coords[n][2]   = d_z[ni];
        phi[n]         = d_phi[ni];
        gamma[n]       = d_gamma[ni];
        beta[n]        = d_beta[ni];
        grad_phi[n][0] = d_grad_phi_x[ni];
        grad_phi[n][1] = d_grad_phi_y[ni];
        grad_phi[n][2] = d_grad_phi_z[ni];
        dofs[n]        = d_node_to_dof[ni];
        own[n]         = d_ownership[ni];
    }

    // =========================================================================
    // Pre-lookup all 64 CSR positions (64 ints = 64 regs vs lhs[64] = 128 regs)
    // Done once per element, outside the SCS hot loop.
    // pos[i*8+j] = index in matrix->values for (dof[i], dof[j]), or -1.
    // =========================================================================
    int pos[64];

    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        int row_dof = dofs[i];

        if (own[i] == 0 || row_dof < 0) {
            #pragma unroll
            for (int j = 0; j < 8; ++j) pos[i * 8 + j] = -1;
            continue;
        }

        int row_start = matrix->rowPtr[row_dof];
        int row_end   = matrix->rowPtr[row_dof + 1];
        int diag_pos  = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            int col_dof = dofs[j];
            if (col_dof < 0) { pos[i * 8 + j] = -1; continue; }
            if (i == j)      { pos[i * 8 + j] = diag_pos; continue; }

            int left = row_start, right = row_end - 1, found = -1;
            while (left <= right) {
                int mid = (left + right) >> 1;
                int col = matrix->colInd[mid];
                if      (col == col_dof) { found = mid; break; }
                else if (col <  col_dof)   left  = mid + 1;
                else                       right = mid - 1;
            }
            pos[i * 8 + j] = found;
        }
    }

    // =========================================================================
    // SCS loop: scatter LHS immediately, accumulate RHS
    // =========================================================================
    RealType rhs[8] = {};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        RealType mdot     = d_mdot[elemIdx * 12 + ip];
        RealType areaVec0 = d_areaVec_x[elemIdx * 12 + ip];
        RealType areaVec1 = d_areaVec_y[elemIdx * 12 + ip];
        RealType areaVec2 = d_areaVec_z[elemIdx * 12 + ip];

        // Interpolate gamma and coords to SCS
        RealType gamma_ip = 0.0;
        RealType cx = 0.0, cy = 0.0, cz = 0.0;
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType sf = hexShapeFcn[ip][n];
            gamma_ip += sf * gamma[n];
            cx += sf * coords[n][0];
            cy += sf * coords[n][1];
            cz += sf * coords[n][2];
        }

        // Advection upwind
        RealType phi_up, beta_up, dcorr;
        if (mdot > 0.0) {
            phi_up  = phi[nodeL];
            beta_up = beta[nodeL];
            dcorr   = grad_phi[nodeL][0] * (cx - coords[nodeL][0]) +
                      grad_phi[nodeL][1] * (cy - coords[nodeL][1]) +
                      grad_phi[nodeL][2] * (cz - coords[nodeL][2]);
        } else {
            phi_up  = phi[nodeR];
            beta_up = beta[nodeR];
            dcorr   = grad_phi[nodeR][0] * (cx - coords[nodeR][0]) +
                      grad_phi[nodeR][1] * (cy - coords[nodeR][1]) +
                      grad_phi[nodeR][2] * (cz - coords[nodeR][2]);
        }
        dcorr *= beta_up;

        RealType adv_flux = mdot * (phi_up + dcorr);
        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        // Advection LHS — only two non-zero pairs per SCS (upwind side)
        RealType lhsfac_L = 0.5 * (mdot + fabs(mdot));
        RealType lhsfac_R = 0.5 * (mdot - fabs(mdot));

        // Row nodeL
        if (lhsfac_L != 0.0) {
            int p = pos[nodeL * 8 + nodeL];
            if (p >= 0) atomicAdd(&matrix->values[p],  lhsfac_L);
            p = pos[nodeR * 8 + nodeL];
            if (p >= 0) atomicAdd(&matrix->values[p], -lhsfac_L);
        }
        if (lhsfac_R != 0.0) {
            int p = pos[nodeL * 8 + nodeR];
            if (p >= 0) atomicAdd(&matrix->values[p],  lhsfac_R);
            p = pos[nodeR * 8 + nodeR];
            if (p >= 0) atomicAdd(&matrix->values[p], -lhsfac_R);
        }

        // Diffusion
        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType diff_coeff = -(gamma_ip * (dndx[n][0] * areaVec0 +
                                                dndx[n][1] * areaVec1 +
                                                dndx[n][2] * areaVec2));
            rhs[nodeL] -= diff_coeff * phi[n];
            rhs[nodeR] += diff_coeff * phi[n];

            int pL = pos[nodeL * 8 + n];
            int pR = pos[nodeR * 8 + n];
            if (pL >= 0) atomicAdd(&matrix->values[pL],  diff_coeff);
            if (pR >= 0) atomicAdd(&matrix->values[pR], -diff_coeff);
        }
    }

    // =========================================================================
    // Scatter RHS
    // =========================================================================
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        if (own[i] != 0 && dofs[i] >= 0)
            atomicAdd(&d_rhs[dofs[i]], rhs[i]);
    }
}

// =============================================================================
// Variant 1: no launch_bounds (natural register count, ~108-140 regs estimated)
// =============================================================================
template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void cvfem_hex_assembly_kernel_tensor_perip(
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
    if (elemIdx >= (int)numElements) return;
    cvfem_hex_perip_body<KeyType, RealType>(
        elemIdx,
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_conn4, d_conn5, d_conn6, d_conn7,
        d_x, d_y, d_z,
        d_gamma, d_phi, d_beta,
        d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
        d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
        d_node_to_dof, d_ownership,
        matrix, d_rhs);
}

// =============================================================================
// Variant 2: __launch_bounds__(256, 2) forces 2 blocks/SM.
// With lhs[64] removed, spilling burden is ~7 regs (vs ~44 regs for tensor+lb2).
// =============================================================================
template<typename KeyType, typename RealType, int BlockSize = 256>
__launch_bounds__(256, 2)
__global__ void cvfem_hex_assembly_kernel_tensor_perip_lb2(
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
    if (elemIdx >= (int)numElements) return;
    cvfem_hex_perip_body<KeyType, RealType>(
        elemIdx,
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_conn4, d_conn5, d_conn6, d_conn7,
        d_x, d_y, d_z,
        d_gamma, d_phi, d_beta,
        d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
        d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
        d_node_to_dof, d_ownership,
        matrix, d_rhs);
}

} // namespace fem
} // namespace mars
