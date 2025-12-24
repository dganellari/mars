#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

// CVFEM assembly kernel with graph-based sparsity pattern and diagonal lumping
// This version matches STK's reduced sparsity approach where off-diagonal
// entries not in the graph are lumped to the diagonal
template<typename KeyType, typename RealType>
__global__ void cvfem_hex_assembly_kernel_graph(
    // Element connectivity (local node IDs)
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    const KeyType* __restrict__ d_conn4,
    const KeyType* __restrict__ d_conn5,
    const KeyType* __restrict__ d_conn6,
    const KeyType* __restrict__ d_conn7,
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
    // Element data
    const RealType* __restrict__ d_mdot,
    const RealType* __restrict__ d_areaVec_x,
    const RealType* __restrict__ d_areaVec_y,
    const RealType* __restrict__ d_areaVec_z,
    // DOF mapping
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    // Matrix assembly (CSR format with graph-based sparsity)
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Get element nodes
    KeyType nodes[8];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];
    nodes[4] = d_conn4[elemIdx];
    nodes[5] = d_conn5[elemIdx];
    nodes[6] = d_conn6[elemIdx];
    nodes[7] = d_conn7[elemIdx];

    // Gather nodal coordinates
    double coords[8][3];
    for (int n = 0; n < 8; ++n) {
        coords[n][0] = d_x[nodes[n]];
        coords[n][1] = d_y[nodes[n]];
        coords[n][2] = d_z[nodes[n]];
    }

    // Debug: print element 0 coordinates
    if (elemIdx == 0 && blockIdx.x == 0 && threadIdx.x == 0) {
        printf("\n=== MARS Element 0 Node Info ===\n");
        printf("Local node IDs: [%u, %u, %u, %u, %u, %u, %u, %u]\n",
               nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
        for (int ni = 0; ni < 8; ++ni) {
            printf("Node %d (local=%u): [%.10e, %.10e, %.10e]\n",
                   ni, nodes[ni], coords[ni][0], coords[ni][1], coords[ni][2]);
        }
        printf("\n");
    }

    // Gather nodal field values
    double phi[8], gamma[8], beta[8];
    double grad_phi[8][3];
    for (int n = 0; n < 8; ++n) {
        phi[n] = d_phi[nodes[n]];
        gamma[n] = d_gamma[nodes[n]];
        beta[n] = d_beta[nodes[n]];
        grad_phi[n][0] = d_grad_phi_x[nodes[n]];
        grad_phi[n][1] = d_grad_phi_y[nodes[n]];
        grad_phi[n][2] = d_grad_phi_z[nodes[n]];
    }

    // Local element matrix and RHS (8x8 matrix, 8-vector)
    double lhs[64] = {0.0};
    double rhs[8] = {0.0};

    // Loop over 12 sub-control surfaces (SCS)
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        double mdot = d_mdot[elemIdx * 12 + ip];

        // Compute area vector from geometry (like STK does)
        double areaVec[3];
        computeAreaVector(ip, coords, areaVec);

        // Debug: print element 0 area vectors
        if (elemIdx == 0 && blockIdx.x == 0 && threadIdx.x == 0) {
            printf("MARS Element 0, SCS %d: areaVec = [%.10e, %.10e, %.10e]\n",
                   ip, areaVec[0], areaVec[1], areaVec[2]);
        }

        // Interpolate fields to SCS integration point
        double phi_ip = 0.0, gamma_ip = 0.0;
        double grad_phi_ip[3] = {0.0, 0.0, 0.0};
        double coords_ip[3] = {0.0, 0.0, 0.0};
        for (int node = 0; node < 8; ++node) {
            double sf = hexShapeFcn[ip][node];
            phi_ip += sf * phi[node];
            gamma_ip += sf * gamma[node];
            for (int d = 0; d < 3; ++d) {
                grad_phi_ip[d] += sf * grad_phi[node][d];
                coords_ip[d] += sf * coords[node][d];
            }
        }

        // Advection flux
        double phi_L = phi[nodeL];
        double phi_R = phi[nodeR];
        double phi_upwind, beta_upwind;
        double dcorr = 0.0;

        if (mdot > 0.0) {
            phi_upwind = phi_L;
            beta_upwind = beta[nodeL];
            for (int d = 0; d < 3; ++d) {
                double dx = coords_ip[d] - coords[nodeL][d];
                dcorr += grad_phi[nodeL][d] * dx;
            }
        } else {
            phi_upwind = phi_R;
            beta_upwind = beta[nodeR];
            for (int d = 0; d < 3; ++d) {
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

        lhs[nodeL * 8 + nodeL] += lhsfac_L;
        lhs[nodeR * 8 + nodeL] -= lhsfac_L;
        lhs[nodeL * 8 + nodeR] += lhsfac_R;
        lhs[nodeR * 8 + nodeR] -= lhsfac_R;

        // Diffusion flux
        double dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        for (int node = 0; node < 8; ++node) {
            double diff_coeff = 0.0;
            for (int d = 0; d < 3; ++d) {
                diff_coeff -= gamma_ip * dndx[node][d] * areaVec[d];
            }

            lhs[nodeL * 8 + node] += diff_coeff;
            lhs[nodeR * 8 + node] -= diff_coeff;

            rhs[nodeL] -= diff_coeff * phi[node];
            rhs[nodeR] += diff_coeff * phi[node];
        }
    }

    // Assemble into global system with diagonal lumping for missing entries
    for (int i = 0; i < 8; ++i) {
        KeyType row_node = nodes[i];
        int row_dof = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];

        if (ownership == 0 || row_dof < 0) continue;

        // Assemble RHS
        atomicAdd(&d_rhs[row_dof], rhs[i]);

        // Diagonal accumulator for lumped entries
        double diag_lump = 0.0;
        int num_found = 0;
        int num_lumped = 0;

        // Assemble LHS with diagonal lumping
        for (int j = 0; j < 8; ++j) {
            KeyType col_node = nodes[j];
            int col_dof = d_node_to_dof[col_node];
            if (col_dof < 0) continue;

            double value = lhs[i * 8 + j];

            if (row_dof == col_dof) {
                // Diagonal entry - will add lumped contributions later
                continue;
            }

            // Try to find this entry in the sparsity pattern
            bool found = false;
            int start = matrix->rowPtr[row_dof];
            int end = matrix->rowPtr[row_dof + 1];

            for (int k = start; k < end; ++k) {
                if (matrix->colInd[k] == col_dof) {
                    atomicAdd(&matrix->values[k], value);
                    found = true;
                    num_found++;
                    break;
                }
            }

            // If entry not found in sparsity pattern, lump to diagonal (STK behavior)
            if (!found) {
                diag_lump += value;
                num_lumped++;
            }
        }

        // Debug: print stats for first element
        if (elemIdx == 0 && i == 0 && blockIdx.x == 0 && threadIdx.x == 0) {
            printf("Element 0, node 0: found=%d, lumped=%d, diag_lump=%.6e, diag_orig=%.6e\n",
                   num_found, num_lumped, diag_lump, lhs[i * 8 + i]);
            printf("  row_dof=%d, graph row has %d entries: [",
                   row_dof, matrix->rowPtr[row_dof + 1] - matrix->rowPtr[row_dof]);
            for (int k = matrix->rowPtr[row_dof]; k < matrix->rowPtr[row_dof + 1]; ++k) {
                printf("%d", matrix->colInd[k]);
                if (k < matrix->rowPtr[row_dof + 1] - 1) printf(", ");
            }
            printf("]\n");
            printf("  element nodes col_dofs: [");
            for (int jj = 0; jj < 8; ++jj) {
                printf("%d", d_node_to_dof[nodes[jj]]);
                if (jj < 7) printf(", ");
            }
            printf("]\n");
        }

        // Add diagonal entry + lumped contributions
        double diag_value = lhs[i * 8 + i] + diag_lump;

        // Find diagonal entry
        int start = matrix->rowPtr[row_dof];
        int end = matrix->rowPtr[row_dof + 1];
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
