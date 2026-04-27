#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_common.hpp"
#include "mars_cvfem_node_data.hpp"

namespace mars {
namespace fem {

// =============================================================================
// AoS-layout CVFEM assembly kernel
// =============================================================================
// Identical to the tensor kernel, but node data is read from a packed
// NodeData AoS array instead of 9 separate SoA arrays.
//
// Why AoS is better here:
//   Each thread reads ALL 9 fields (x,y,z,phi,gamma,beta,gx,gy,gz) for ONE
//   randomly-indexed node. SoA forces 9 scattered loads (9 sectors per thread
//   = 288 sectors/warp). AoS fetches all fields in one 72-byte contiguous
//   read (3 sectors per thread = 96 sectors/warp) → 3× less L2 traffic.
//
// Connectivity, mdot, and areaVec remain SoA (sequential element access).
// =============================================================================

template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void cvfem_hex_assembly_kernel_tensor_aos(
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    const KeyType* __restrict__ d_conn4,
    const KeyType* __restrict__ d_conn5,
    const KeyType* __restrict__ d_conn6,
    const KeyType* __restrict__ d_conn7,
    size_t numElements,
    const NodeData* __restrict__ d_nodeData,   // AoS: all 9 fields packed
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
    // Load element connectivity
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
    // AoS node data loads — each node read fetches 72 bytes (3 sectors)
    // vs 9 separate SoA loads (9 sectors per node).
    // ==========================================================================
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        const NodeData& nd = d_nodeData[nodes[n]];
        coords[n][0]   = nd.x;
        coords[n][1]   = nd.y;
        coords[n][2]   = nd.z;
        phi[n]         = nd.phi;
        gamma[n]       = nd.gamma;
        beta[n]        = nd.beta;
        grad_phi[n][0] = nd.gx;
        grad_phi[n][1] = nd.gy;
        grad_phi[n][2] = nd.gz;
        dofs[n]        = d_node_to_dof[nodes[n]];
        own[n]         = d_ownership[nodes[n]];
    }

    // ==========================================================================
    // Full 8×8 local element matrix and RHS vector
    // ==========================================================================
    RealType lhs[64];
    RealType rhs[8];

    #pragma unroll
    for (int i = 0; i < 64; ++i) lhs[i] = 0.0;
    #pragma unroll
    for (int i = 0; i < 8;  ++i) rhs[i] = 0.0;

    // ==========================================================================
    // SCS integration loop (12 sub-control-surface integration points)
    // ==========================================================================
    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        RealType mdot = d_mdot[elemIdx * 12 + ip];

        RealType areaVec[3];
        areaVec[0] = d_areaVec_x[elemIdx * 12 + ip];
        areaVec[1] = d_areaVec_y[elemIdx * 12 + ip];
        areaVec[2] = d_areaVec_z[elemIdx * 12 + ip];

        // Interpolate gamma and coords to SCS
        RealType gamma_ip = 0.0;
        RealType coords_ip[3] = {0.0, 0.0, 0.0};

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType sf = hexShapeFcn[ip][n];
            gamma_ip      += sf * gamma[n];
            coords_ip[0]  += sf * coords[n][0];
            coords_ip[1]  += sf * coords[n][1];
            coords_ip[2]  += sf * coords[n][2];
        }

        // Advection (upwind + deferred correction)
        RealType phi_upwind, beta_upwind, dcorr = 0.0;
        if (mdot > 0.0) {
            phi_upwind  = phi[nodeL];
            beta_upwind = beta[nodeL];
            dcorr = grad_phi[nodeL][0] * (coords_ip[0] - coords[nodeL][0]) +
                    grad_phi[nodeL][1] * (coords_ip[1] - coords[nodeL][1]) +
                    grad_phi[nodeL][2] * (coords_ip[2] - coords[nodeL][2]);
        } else {
            phi_upwind  = phi[nodeR];
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

        // Diffusion — compute shape derivatives on-the-fly
        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType diff_coeff = -(gamma_ip * (dndx[n][0] * areaVec[0] +
                                                dndx[n][1] * areaVec[1] +
                                                dndx[n][2] * areaVec[2]));
            rhs[nodeL] -= diff_coeff * phi[n];
            rhs[nodeR] += diff_coeff * phi[n];
            lhs[nodeL * 8 + n] += diff_coeff;
            lhs[nodeR * 8 + n] -= diff_coeff;
        }
    }

    // ==========================================================================
    // Scatter to global CSR matrix (atomics; same as tensor kernel)
    // ==========================================================================
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        int row_dof   = dofs[i];
        uint8_t ownership = own[i];

        if (ownership == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        int row_start = matrix->rowPtr[row_dof];
        int row_end   = matrix->rowPtr[row_dof + 1];
        int diag_pos  = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            int col_dof = dofs[j];
            if (col_dof < 0) continue;

            RealType value = lhs[i * 8 + j];
            if (value == 0.0) continue;

            if (i == j) {
                atomicAdd(&matrix->values[diag_pos], value);
            } else {
                int left  = row_start;
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
