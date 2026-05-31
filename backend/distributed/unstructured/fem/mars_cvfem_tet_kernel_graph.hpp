#pragma once

#include "mars_cvfem_kernel.hpp"     // Tet4CVFEM reference data
#include "mars_cvfem_hex_kernel.hpp" // mars::fem::CSRMatrix lives here
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// CVFEM assembly kernel for linear tetrahedra with graph-based CSR sparsity
// and diagonal lumping for missing entries. Mirrors the hex graph kernel
// structure (mars_cvfem_hex_kernel_graph.hpp) so the dispatch in
// mars_cvfem_graph.cu can branch cleanly on element type.
//
// Tet4 specifics (from Tet4CVFEM in mars_cvfem_kernel.hpp):
//   - 6 sub-control surfaces (one per element edge)
//   - linear shape functions; dN/dxi constant per element
//   - SCS area vector: simplified as scaled edge vector
template<typename KeyType, typename RealType>
__global__ void cvfem_tet_assembly_kernel_graph(
    // Element connectivity (4 columns of local node IDs)
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    size_t numElements,
    // Node coordinates
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    // Field data (per node)
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_phi,
    const RealType* __restrict__ d_beta,
    const RealType* __restrict__ d_grad_phi_x,
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    // Element data: mdot per element per SCS (6 SCS for tets)
    const RealType* __restrict__ d_mdot,
    // DOF mapping
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    // Matrix assembly (CSR with graph-based sparsity)
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    constexpr int nnodes = 4;
    constexpr int nscs   = 6;
    constexpr int ndim   = 3;

    KeyType nodes[nnodes];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];

    RealType coords[nnodes][ndim];
    RealType phi[nnodes], gamma[nnodes], beta[nnodes];
    RealType grad_phi[nnodes][ndim];
    for (int n = 0; n < nnodes; ++n) {
        coords[n][0] = d_x[nodes[n]];
        coords[n][1] = d_y[nodes[n]];
        coords[n][2] = d_z[nodes[n]];
        phi[n]      = d_phi[nodes[n]];
        gamma[n]    = d_gamma[nodes[n]];
        beta[n]     = d_beta[nodes[n]];
        grad_phi[n][0] = d_grad_phi_x[nodes[n]];
        grad_phi[n][1] = d_grad_phi_y[nodes[n]];
        grad_phi[n][2] = d_grad_phi_z[nodes[n]];
    }

    RealType det;
    RealType dNdx[nnodes][ndim];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);

    RealType lhs[nnodes * nnodes] = {RealType(0)};
    RealType rhs[nnodes]          = {RealType(0)};

    RealType vol_scale = fabs(det) / RealType(6);   // tet volume

    // Diffusion: classical linear-tet stiffness K_ij = vol * gamma_elem *
    // (dNdx_i . dNdx_j). dNdx is element-constant for a linear tet, so this
    // is exact with a 1-point centroid rule (gamma at the centroid = mean of
    // the 4 nodal gammas). Same operator as the full kernel; assembled once
    // per element rather than per SCS. The graph CSR scatter below lumps any
    // (i,j) entry missing from the reduced sparsity onto the diagonal, so a
    // dense 4x4 stencil stays consistent (row-sum preserved).
    RealType gamma_elem = RealType(0.25) * (gamma[0] + gamma[1] + gamma[2] + gamma[3]);
    for (int i = 0; i < nnodes; ++i) {
        for (int j = 0; j < nnodes; ++j) {
            RealType k_ij = RealType(0);
            for (int d = 0; d < ndim; ++d) k_ij += dNdx[i][d] * dNdx[j][d];
            lhs[i * nnodes + j] += vol_scale * gamma_elem * k_ij;
        }
    }
    // Residual form: rhs -= K * phi (matches assembleFull; no-op when phi=0).
    for (int i = 0; i < nnodes; ++i) {
        RealType r = RealType(0);
        for (int j = 0; j < nnodes; ++j) r += lhs[i * nnodes + j] * phi[j];
        rhs[i] -= r;
    }

    for (int ip = 0; ip < nscs; ++ip) {
        int nodeL, nodeR;
        Tet4CVFEM::get_scs_nodes(ip, nodeL, nodeR);

        double xi, eta, zeta;
        Tet4CVFEM::scs_coords(ip, xi, eta, zeta);

        double N[nnodes];
        Tet4CVFEM::shape_fcn(xi, eta, zeta, N);

        // coords_ip: SCS integration point, used by the advection deferred
        // correction below. Diffusion is handled element-wise above, not here.
        RealType coords_ip[ndim] = {RealType(0), RealType(0), RealType(0)};
        for (int n = 0; n < nnodes; ++n) {
            for (int d = 0; d < ndim; ++d) coords_ip[d] += N[n] * coords[n][d];
        }

        RealType mdot = d_mdot[elemIdx * nscs + ip];

        // Advection (upwind + deferred-correction matching the hex kernel pattern).
        RealType phi_L = phi[nodeL];
        RealType phi_R = phi[nodeR];
        RealType phi_upwind, beta_upwind, dcorr = RealType(0);
        if (mdot > RealType(0)) {
            phi_upwind  = phi_L;
            beta_upwind = beta[nodeL];
            for (int d = 0; d < ndim; ++d) {
                RealType dx = coords_ip[d] - coords[nodeL][d];
                dcorr += grad_phi[nodeL][d] * dx;
            }
        } else {
            phi_upwind  = phi_R;
            beta_upwind = beta[nodeR];
            for (int d = 0; d < ndim; ++d) {
                RealType dx = coords_ip[d] - coords[nodeR][d];
                dcorr += grad_phi[nodeR][d] * dx;
            }
        }
        dcorr *= beta_upwind;
        RealType adv_flux = mdot * (phi_upwind + dcorr);

        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        RealType lhsfac_L = RealType(0.5) * (mdot + fabs(mdot));
        RealType lhsfac_R = RealType(0.5) * (mdot - fabs(mdot));
        lhs[nodeL * nnodes + nodeL] += lhsfac_L;
        lhs[nodeR * nnodes + nodeL] -= lhsfac_L;
        lhs[nodeL * nnodes + nodeR] += lhsfac_R;
        lhs[nodeR * nnodes + nodeR] -= lhsfac_R;
    }

    // Scatter into CSR with graph-sparsity + diagonal lumping for missing entries.
    for (int i = 0; i < nnodes; ++i) {
        KeyType row_node = nodes[i];
        int row_dof      = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];
        if (ownership == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        RealType diag_lump = RealType(0);
        int start = matrix->rowPtr[row_dof];
        int end   = matrix->rowPtr[row_dof + 1];

        for (int j = 0; j < nnodes; ++j) {
            KeyType col_node = nodes[j];
            int col_dof      = d_node_to_dof[col_node];
            if (col_dof < 0) continue;

            RealType value = lhs[i * nnodes + j];
            if (row_dof == col_dof) continue;

            bool found = false;
            for (int k = start; k < end; ++k) {
                if (matrix->colInd[k] == col_dof) {
                    atomicAdd(&matrix->values[k], value);
                    found = true;
                    break;
                }
            }
            if (!found) diag_lump += value;
        }

        RealType diag_value = lhs[i * nnodes + i] + diag_lump;
        for (int k = start; k < end; ++k) {
            if (matrix->colInd[k] == row_dof) {
                atomicAdd(&matrix->values[k], diag_value);
                break;
            }
        }
    }
}

} // namespace fem
} // namespace mars
