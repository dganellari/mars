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

    double coords[nnodes][ndim];
    double phi[nnodes], gamma[nnodes], beta[nnodes];
    double grad_phi[nnodes][ndim];
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

    // Jacobian: J[i][j] = sum_n coords[n][i] * dNdxi[n][j].
    // dNdxi is constant for linear tet:
    //   dNdxi[0] = (-1,-1,-1), dNdxi[1] = (1,0,0), dNdxi[2] = (0,1,0), dNdxi[3] = (0,0,1)
    // so J[i][j] = coords[1][i]-coords[0][i] for j=0, etc.
    double J[ndim][ndim];
    for (int i = 0; i < ndim; ++i) {
        J[i][0] = coords[1][i] - coords[0][i];
        J[i][1] = coords[2][i] - coords[0][i];
        J[i][2] = coords[3][i] - coords[0][i];
    }
    double det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
               - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
               + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    double inv_det = 1.0 / det;

    double Jinv[ndim][ndim];
    Jinv[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * inv_det;
    Jinv[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2]) * inv_det;
    Jinv[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * inv_det;
    Jinv[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2]) * inv_det;
    Jinv[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * inv_det;
    Jinv[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2]) * inv_det;
    Jinv[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * inv_det;
    Jinv[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1]) * inv_det;
    Jinv[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * inv_det;

    // Physical-space gradients: dN/dx = dN/dxi * Jinv.
    // For linear tet, dNdx is constant per element.
    double dNdx[nnodes][ndim];
    // dNdxi[0] = (-1,-1,-1)
    for (int d = 0; d < ndim; ++d) dNdx[0][d] = -Jinv[0][d] - Jinv[1][d] - Jinv[2][d];
    // dNdxi[1] = (1,0,0)
    for (int d = 0; d < ndim; ++d) dNdx[1][d] =  Jinv[0][d];
    // dNdxi[2] = (0,1,0)
    for (int d = 0; d < ndim; ++d) dNdx[2][d] =  Jinv[1][d];
    // dNdxi[3] = (0,0,1)
    for (int d = 0; d < ndim; ++d) dNdx[3][d] =  Jinv[2][d];

    double lhs[nnodes * nnodes] = {0.0};
    double rhs[nnodes]          = {0.0};

    double vol_scale = fabs(det) / 6.0;   // tet volume

    for (int ip = 0; ip < nscs; ++ip) {
        int nodeL, nodeR;
        Tet4CVFEM::get_scs_nodes(ip, nodeL, nodeR);

        double xi, eta, zeta;
        Tet4CVFEM::scs_coords(ip, xi, eta, zeta);

        double N[nnodes];
        Tet4CVFEM::shape_fcn(xi, eta, zeta, N);

        double gamma_ip = 0.0;
        double coords_ip[ndim] = {0.0, 0.0, 0.0};
        for (int n = 0; n < nnodes; ++n) {
            gamma_ip += N[n] * gamma[n];
            for (int d = 0; d < ndim; ++d) coords_ip[d] += N[n] * coords[n][d];
        }

        // Simplified SCS area vector: edge vector scaled by element volume / nscs.
        double areaVec[ndim];
        for (int d = 0; d < ndim; ++d) {
            areaVec[d] = (coords[nodeR][d] - coords[nodeL][d]) * vol_scale / nscs;
        }

        double mdot = d_mdot[elemIdx * nscs + ip];

        // Advection (upwind + deferred-correction matching the hex kernel pattern).
        double phi_L = phi[nodeL];
        double phi_R = phi[nodeR];
        double phi_upwind, beta_upwind, dcorr = 0.0;
        if (mdot > 0.0) {
            phi_upwind  = phi_L;
            beta_upwind = beta[nodeL];
            for (int d = 0; d < ndim; ++d) {
                double dx = coords_ip[d] - coords[nodeL][d];
                dcorr += grad_phi[nodeL][d] * dx;
            }
        } else {
            phi_upwind  = phi_R;
            beta_upwind = beta[nodeR];
            for (int d = 0; d < ndim; ++d) {
                double dx = coords_ip[d] - coords[nodeR][d];
                dcorr += grad_phi[nodeR][d] * dx;
            }
        }
        dcorr *= beta_upwind;
        double adv_flux = mdot * (phi_upwind + dcorr);

        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        double lhsfac_L = 0.5 * (mdot + fabs(mdot));
        double lhsfac_R = 0.5 * (mdot - fabs(mdot));
        lhs[nodeL * nnodes + nodeL] += lhsfac_L;
        lhs[nodeR * nnodes + nodeL] -= lhsfac_L;
        lhs[nodeL * nnodes + nodeR] += lhsfac_R;
        lhs[nodeR * nnodes + nodeR] -= lhsfac_R;

        // Diffusion: per-node contribution. dNdx is element-constant so
        // hoist outside the SCS loop conceptually, but kept inside for
        // structural parity with the hex kernel; cost is negligible (4 nodes).
        for (int n = 0; n < nnodes; ++n) {
            double diff_coeff = 0.0;
            for (int d = 0; d < ndim; ++d) {
                diff_coeff -= gamma_ip * dNdx[n][d] * areaVec[d];
            }
            lhs[nodeL * nnodes + n] += diff_coeff;
            lhs[nodeR * nnodes + n] -= diff_coeff;
            rhs[nodeL] -= diff_coeff * phi[n];
            rhs[nodeR] += diff_coeff * phi[n];
        }
    }

    // Scatter into CSR with graph-sparsity + diagonal lumping for missing entries.
    for (int i = 0; i < nnodes; ++i) {
        KeyType row_node = nodes[i];
        int row_dof      = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];
        if (ownership == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        double diag_lump = 0.0;
        int start = matrix->rowPtr[row_dof];
        int end   = matrix->rowPtr[row_dof + 1];

        for (int j = 0; j < nnodes; ++j) {
            KeyType col_node = nodes[j];
            int col_dof      = d_node_to_dof[col_node];
            if (col_dof < 0) continue;

            double value = lhs[i * nnodes + j];
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

        double diag_value = lhs[i * nnodes + i] + diag_lump;
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
