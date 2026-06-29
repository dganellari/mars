#pragma once

#include "mars_cvfem_kernel.hpp"
#include "mars_cvfem_hex_kernel.hpp"   // CSRMatrix POD struct
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// Full CVFEM assembly kernel for linear tetrahedra, perip-style.
//
// Same math as cvfem_tet_assembly_kernel_full (mars_cvfem_tet_kernel_full.hpp),
// but with the same memory pattern that makes the hex perip kernel a winner:
//
//   * Pre-lookup all 16 CSR positions once per element (4x4 element matrix),
//     before the SCS advection loop. Binary search per (row, col) instead of
//     linear scan inside the addValue hot loop.
//   * Diagonal fast-path via matrix->diagPtr[row_dof] (no search at all).
//   * Scatter via direct atomicAdd(&values[pos[i*4+j]], ...). No repeated
//     row-range lookup per face.
//   * __restrict__ + __ldg on the connectivity/coords/fields reads.
//   * __launch_bounds__(256, 4): tet has 4 nodes vs hex's 8, so per-thread
//     state is ~half. 4 blocks/SM is the right occupancy target.
//
// Sparsity must contain all (i,j) pairs for the 4 element corners (16-NNZ
// pattern per element). buildElementGraphSparsity already produces this.

template<typename KeyType, typename RealType>
__device__ __forceinline__ void cvfem_tet_full_perip_body(
    int elemIdx,
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
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
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    constexpr int nnodes = 4;
    constexpr int nscs   = 6;
    constexpr int ndim   = 3;

    // ------------------------------------------------------------------
    // Connectivity + per-node data, all in one tight pass with __ldg.
    // ------------------------------------------------------------------
    KeyType nodes[nnodes];
    nodes[0] = __ldg(&d_conn0[elemIdx]);
    nodes[1] = __ldg(&d_conn1[elemIdx]);
    nodes[2] = __ldg(&d_conn2[elemIdx]);
    nodes[3] = __ldg(&d_conn3[elemIdx]);

    RealType coords[nnodes][ndim];
    RealType phi[nnodes], gamma[nnodes], beta[nnodes];
    RealType grad_phi[nnodes][ndim];
    int      dofs[nnodes];
    uint8_t  own[nnodes];

    #pragma unroll
    for (int n = 0; n < nnodes; ++n) {
        KeyType ni = nodes[n];
        coords[n][0]   = __ldg(&d_x[ni]);
        coords[n][1]   = __ldg(&d_y[ni]);
        coords[n][2]   = __ldg(&d_z[ni]);
        phi[n]         = __ldg(&d_phi[ni]);
        gamma[n]       = __ldg(&d_gamma[ni]);
        beta[n]        = __ldg(&d_beta[ni]);
        grad_phi[n][0] = __ldg(&d_grad_phi_x[ni]);
        grad_phi[n][1] = __ldg(&d_grad_phi_y[ni]);
        grad_phi[n][2] = __ldg(&d_grad_phi_z[ni]);
        dofs[n]        = __ldg(&d_node_to_dof[ni]);
        own[n]         = d_ownership[ni];
    }

    // ------------------------------------------------------------------
    // Pre-lookup all 16 CSR positions (one binary search per (i,j)).
    // pos[i*4+j] = index in matrix->values for (dofs[i], dofs[j]), or -1.
    //
    // -1 means "skip the scatter":
    //   - col_dof < 0          (column DOF eliminated, e.g. periodic slave)
    //   - row_dof < 0 or ghost (skip the whole row; we still scatter to
    //                          d_rhs only for owned + valid rows below)
    //
    // Diagonal hit comes from matrix->diagPtr[row_dof] directly (no search).
    // ------------------------------------------------------------------
    int pos[nnodes * nnodes];

    #pragma unroll
    for (int i = 0; i < nnodes; ++i) {
        int row_dof = dofs[i];

        if (own[i] == 0 || row_dof < 0) {
            #pragma unroll
            for (int j = 0; j < nnodes; ++j) pos[i * nnodes + j] = -1;
            continue;
        }

        const int row_start = matrix->rowPtr[row_dof];
        const int row_end   = matrix->rowPtr[row_dof + 1];
        const int diag_pos  = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int j = 0; j < nnodes; ++j) {
            int col_dof = dofs[j];
            if (col_dof < 0) { pos[i * nnodes + j] = -1; continue; }
            if (i == j)      { pos[i * nnodes + j] = diag_pos; continue; }

            // Binary search on the sorted row slice. The sparsity builder
            // sorts column indices per row at allocation, so this is safe.
            int left  = row_start;
            int right = row_end - 1;
            int found = -1;
            while (left <= right) {
                int mid    = (left + right) >> 1;
                int midCol = matrix->colInd[mid];
                if (midCol == col_dof) { found = mid; break; }
                if (midCol < col_dof) left = mid + 1;
                else                  right = mid - 1;
            }
            pos[i * nnodes + j] = found;
        }
    }

    RealType det;
    RealType dNdx[nnodes][ndim];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);

    RealType lhs[nnodes * nnodes] = {RealType(0)};
    RealType rhs[nnodes]          = {RealType(0)};

    const RealType volume = fabs(det) / RealType(6);
    const RealType gamma_elem = RealType(0.25) * (gamma[0] + gamma[1] + gamma[2] + gamma[3]);

    // Linear tet stiffness: K_ij = volume * gamma * (dNdx_i . dNdx_j).
    // Constant per-element since N is linear. 1-pt Gauss = exact.
    #pragma unroll
    for (int i = 0; i < nnodes; ++i) {
        #pragma unroll
        for (int j = 0; j < nnodes; ++j) {
            RealType k_ij = dNdx[i][0]*dNdx[j][0]
                          + dNdx[i][1]*dNdx[j][1]
                          + dNdx[i][2]*dNdx[j][2];
            lhs[i * nnodes + j] = volume * gamma_elem * k_ij;
        }
    }

    // rhs -= K * phi  (residual form, same convention as the hex full kernel).
    #pragma unroll
    for (int i = 0; i < nnodes; ++i) {
        RealType r = RealType(0);
        #pragma unroll
        for (int j = 0; j < nnodes; ++j) r += lhs[i * nnodes + j] * phi[j];
        rhs[i] -= r;
    }

    // ------------------------------------------------------------------
    // Advection on 6 SCS edges. mdot is per (element, ip).
    // Upwind value + deferred correction (dcorr) using grad(phi) at the
    // upwind node. Matrix contributions are the standard 2x2 upwind block.
    // ------------------------------------------------------------------
    #pragma unroll
    for (int ip = 0; ip < nscs; ++ip) {
        RealType mdot = __ldg(&d_mdot[elemIdx * nscs + ip]);
        if (mdot == RealType(0)) continue;

        int nodeL, nodeR;
        Tet4CVFEM::get_scs_nodes(ip, nodeL, nodeR);

        // scs_coords / shape_fcn take double& and double[]; keep these double
        // to match their signatures (the values are exact constants anyway),
        // then accumulate in RealType. Matches the _full / _graph kernels.
        double xi, eta, zeta;
        Tet4CVFEM::scs_coords(ip, xi, eta, zeta);
        double N[nnodes];
        Tet4CVFEM::shape_fcn(xi, eta, zeta, N);

        RealType coords_ip[ndim] = {RealType(0), RealType(0), RealType(0)};
        #pragma unroll
        for (int n = 0; n < nnodes; ++n) {
            #pragma unroll
            for (int d = 0; d < ndim; ++d) coords_ip[d] += RealType(N[n]) * coords[n][d];
        }

        int upw = (mdot > RealType(0)) ? nodeL : nodeR;
        RealType dcorr = RealType(0);
        #pragma unroll
        for (int d = 0; d < ndim; ++d) {
            dcorr += grad_phi[upw][d] * (coords_ip[d] - coords[upw][d]);
        }
        RealType adv_flux = mdot * (phi[upw] + beta[upw] * dcorr);

        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        RealType abs_mdot  = fabs(mdot);
        RealType lhsfac_L  = RealType(0.5) * (mdot + abs_mdot);
        RealType lhsfac_R  = RealType(0.5) * (mdot - abs_mdot);
        lhs[nodeL * nnodes + nodeL] += lhsfac_L;
        lhs[nodeR * nnodes + nodeL] -= lhsfac_L;
        lhs[nodeL * nnodes + nodeR] += lhsfac_R;
        lhs[nodeR * nnodes + nodeR] -= lhsfac_R;
    }

    // ------------------------------------------------------------------
    // Scatter via precomputed positions. RHS goes to owned + valid rows;
    // LHS goes everywhere pos[] is not -1.
    // ------------------------------------------------------------------
    #pragma unroll
    for (int i = 0; i < nnodes; ++i) {
        int row_dof = dofs[i];
        if (own[i] == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        #pragma unroll
        for (int j = 0; j < nnodes; ++j) {
            int p = pos[i * nnodes + j];
            if (p < 0) continue;
            atomicAdd(&matrix->values[p], lhs[i * nnodes + j]);
        }
    }
}

// Public kernel. __launch_bounds__(256, 4) targets 4 blocks/SM. Tet's
// per-thread state (4 corners) is half of hex (8 corners), so we can pack
// more blocks/SM than hex tensor_perip's (256, 2). Adjust if profiling
// shows register spilling.
template<typename KeyType, typename RealType>
__global__ __launch_bounds__(256, 4)
void cvfem_tet_assembly_kernel_full_perip(
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
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
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;
    cvfem_tet_full_perip_body<KeyType, RealType>(
        elemIdx,
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_x, d_y, d_z,
        d_gamma, d_phi, d_beta,
        d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
        d_mdot,
        d_node_to_dof, d_ownership,
        matrix, d_rhs);
}

} // namespace fem
} // namespace mars
