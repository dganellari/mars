#pragma once

#include "mars_base.hpp"
#include "mars_sparse_matrix.hpp"
#include <cuda_runtime.h>

namespace mars {

// CVFEM reference element data for Tet4
struct Tet4CVFEM {
    static constexpr int num_scs = 6;  // 6 interior sub-control surfaces
    static constexpr int nodes_per_elem = 4;
    static constexpr int spatial_dim = 3;
    
    // SCS left-right node pairs
    __device__ __host__ static void get_scs_nodes(int scs_id, int& node_L, int& node_R) {
        const int lr_pairs[12] = {0,1, 1,2, 0,2, 0,3, 1,3, 2,3};
        node_L = lr_pairs[2*scs_id];
        node_R = lr_pairs[2*scs_id + 1];
    }
    
    // Linear tet shape functions at a point (xi, eta, zeta)
    __device__ __host__ static void shape_fcn(double xi, double eta, double zeta, double N[4]) {
        N[0] = 1.0 - xi - eta - zeta;
        N[1] = xi;
        N[2] = eta;
        N[3] = zeta;
    }
    
    // Shape function derivatives in parametric space (constant for linear tet)
    __device__ __host__ static void shape_deriv(double dNdxi[4][3]) {
        // dN/dxi
        dNdxi[0][0] = -1.0; dNdxi[0][1] = -1.0; dNdxi[0][2] = -1.0;
        dNdxi[1][0] =  1.0; dNdxi[1][1] =  0.0; dNdxi[1][2] =  0.0;
        dNdxi[2][0] =  0.0; dNdxi[2][1] =  1.0; dNdxi[2][2] =  0.0;
        dNdxi[3][0] =  0.0; dNdxi[3][1] =  0.0; dNdxi[3][2] =  1.0;
    }
    
    // SCS integration point coordinates (edge midpoints)
    __device__ __host__ static void scs_coords(int scs_id, double& xi, double& eta, double& zeta) {
        const double coords[6][3] = {
            {0.5, 0.5, 0.0},  // edge 0-1
            {0.0, 0.5, 0.5},  // edge 1-2
            {0.5, 0.0, 0.5},  // edge 0-2
            {0.5, 0.0, 0.0},  // edge 0-3
            {0.0, 0.5, 0.0},  // edge 1-3
            {0.0, 0.0, 0.5}   // edge 2-3
        };
        xi = coords[scs_id][0];
        eta = coords[scs_id][1];
        zeta = coords[scs_id][2];
    }
};

// CUDA kernel for CVFEM assembly
template<typename KeyType, typename RealType>
__global__ void cvfem_assembly_kernel(
    // Element connectivity (SFC keys)
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    // Local ID mapping
    const KeyType* __restrict__ d_sfc_to_local,
    size_t nodeCount,
    // Node coordinates
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    // Fields
    const RealType* __restrict__ d_gamma,      // diffusion coefficient [nodeCount]
    const RealType* __restrict__ d_phi,        // solution [nodeCount]
    const RealType* __restrict__ d_beta,       // upwind param [nodeCount]
    const RealType* __restrict__ d_grad_phi_x, // gradient [nodeCount]
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    const RealType* __restrict__ d_mdot,       // mass flow [elemCount * 6]
    // Sparse matrix assembly
    RealType* __restrict__ d_A_values,         // sparse matrix values
    const int* __restrict__ d_A_colidx,        // CSR column indices
    const int* __restrict__ d_A_rowptr,        // CSR row pointers
    const int* __restrict__ d_local_dofs,      // local DOF indices for nodes
    const int* __restrict__ d_node_ownership,  // node ownership (0=ghost, 1=owned, 2=shared)
    // RHS vector
    RealType* __restrict__ d_b,                // RHS vector
    size_t elementCount,
    bool include_advection)
{
    int elem_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (elem_id >= elementCount) return;
    
    using CVFEM = Tet4CVFEM;
    constexpr int nnodes = CVFEM::nodes_per_elem;
    constexpr int nscs = CVFEM::num_scs;
    constexpr int ndim = CVFEM::spatial_dim;
    
    // Get element nodes (local IDs via SFC lookup)
    KeyType sfc_nodes[nnodes] = {
        d_conn0[elem_id], d_conn1[elem_id], 
        d_conn2[elem_id], d_conn3[elem_id]
    };
    
    int elem_nodes[nnodes];
    RealType elem_coords[nnodes][ndim];
    RealType ws_gamma[nnodes], ws_phi[nnodes], ws_beta[nnodes];
    RealType ws_grad_phi[nnodes][ndim];
    
    // Gather data
    for (int n = 0; n < nnodes; ++n) {
        // Binary search for local ID
        KeyType target = sfc_nodes[n];
        int left = 0, right = nodeCount - 1;
        int local_id = -1;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (d_sfc_to_local[mid] == target) {
                local_id = mid;
                break;
            }
            if (d_sfc_to_local[mid] < target) left = mid + 1;
            else right = mid - 1;
        }
        
        elem_nodes[n] = local_id;
        elem_coords[n][0] = d_x[local_id];
        elem_coords[n][1] = d_y[local_id];
        elem_coords[n][2] = d_z[local_id];
        ws_gamma[n] = d_gamma[local_id];
        ws_phi[n] = d_phi[local_id];
        ws_beta[n] = d_beta[local_id];
        ws_grad_phi[n][0] = d_grad_phi_x[local_id];
        ws_grad_phi[n][1] = d_grad_phi_y[local_id];
        ws_grad_phi[n][2] = d_grad_phi_z[local_id];
    }
    
    // Compute Jacobian and inverse
    RealType dNdxi[nnodes][ndim];
    CVFEM::shape_deriv(dNdxi);
    
    // J = dx/dxi
    RealType J[ndim][ndim] = {0};
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
            for (int n = 0; n < nnodes; ++n) {
                J[i][j] += elem_coords[n][i] * dNdxi[n][j];
            }
        }
    }
    
    // Inverse Jacobian (3x3)
    RealType det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
                 - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
                 + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    RealType inv_det = 1.0 / det;
    
    RealType Jinv[ndim][ndim];
    Jinv[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * inv_det;
    Jinv[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2]) * inv_det;
    Jinv[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * inv_det;
    Jinv[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2]) * inv_det;
    Jinv[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * inv_det;
    Jinv[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2]) * inv_det;
    Jinv[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * inv_det;
    Jinv[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1]) * inv_det;
    Jinv[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * inv_det;
    
    // Physical space gradients: dN/dx = dN/dxi * dxi/dx
    RealType dNdx[nnodes][ndim];
    for (int n = 0; n < nnodes; ++n) {
        for (int i = 0; i < ndim; ++i) {
            dNdx[n][i] = 0.0;
            for (int j = 0; j < ndim; ++j) {
                dNdx[n][i] += dNdxi[n][j] * Jinv[j][i];
            }
        }
    }
    
    // Loop over SCS integration points
    for (int scs = 0; scs < nscs; ++scs) {
        int node_L, node_R;
        CVFEM::get_scs_nodes(scs, node_L, node_R);
        
        // Get SCS coordinates
        RealType xi, eta, zeta;
        CVFEM::scs_coords(scs, xi, eta, zeta);
        
        // Shape functions at SCS
        RealType N[nnodes];
        CVFEM::shape_fcn(xi, eta, zeta, N);
        
        // Interpolate to SCS
        RealType gamma_scs = 0.0;
        RealType coord_scs[ndim] = {0};
        for (int n = 0; n < nnodes; ++n) {
            gamma_scs += N[n] * ws_gamma[n];
            for (int d = 0; d < ndim; ++d) {
                coord_scs[d] += N[n] * elem_coords[n][d];
            }
        }
        
        // Area vector (simplified: edge vector scaled by volume)
        RealType area_vec[ndim];
        RealType vol_scale = fabs(det) / 6.0;  // tet volume
        for (int d = 0; d < ndim; ++d) {
            area_vec[d] = (elem_coords[node_R][d] - elem_coords[node_L][d]) * vol_scale / nscs;
        }
        
        // Get mdot
        RealType mdot = include_advection ? d_mdot[elem_id * nscs + scs] : 0.0;
        
        //=============
        // ADVECTION
        //=============
        if (include_advection && fabs(mdot) > 1e-14) {
            RealType phi_L = ws_phi[node_L];
            RealType phi_R = ws_phi[node_R];
            
            // Upwind with deferred correction
            RealType phi_upw, beta, dx[ndim];
            if (mdot > 0.0) {
                phi_upw = phi_L;
                beta = ws_beta[node_L];
                for (int d = 0; d < ndim; ++d) {
                    dx[d] = coord_scs[d] - elem_coords[node_L][d];
                }
            } else {
                phi_upw = phi_R;
                beta = ws_beta[node_R];
                for (int d = 0; d < ndim; ++d) {
                    dx[d] = coord_scs[d] - elem_coords[node_R][d];
                }
            }
            
            RealType deferred = 0.0;
            int upw_node = (mdot > 0.0) ? node_L : node_R;
            for (int d = 0; d < ndim; ++d) {
                deferred += beta * dx[d] * ws_grad_phi[upw_node][d];
            }
            
            RealType adv_flux = mdot * (phi_upw + deferred);
            
            // Get local DOF indices
            int local_L = d_local_dofs[elem_nodes[node_L]];
            int local_R = d_local_dofs[elem_nodes[node_R]];
            int ownership_L = d_node_ownership[elem_nodes[node_L]];
            int ownership_R = d_node_ownership[elem_nodes[node_R]];
            
            // RHS contribution
            if (ownership_L == 1 || ownership_L == 2) {  // owned or shared
                atomicAdd(&d_b[local_L], -adv_flux);
            }
            if (ownership_R == 1 || ownership_R == 2) {  // owned or shared
                atomicAdd(&d_b[local_R], adv_flux);
            }
            
            // LHS contribution (linearized upwind)
            RealType lhs_L = 0.5 * (mdot + fabs(mdot));
            RealType lhs_R = 0.5 * (mdot - fabs(mdot));
            
            // Matrix assembly for advection
            if (ownership_L == 1 || ownership_L == 2) {
                int rowStart = d_A_rowptr[local_L];
                int rowEnd = d_A_rowptr[local_L + 1];
                atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, local_L, lhs_L);
                atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, local_R, -lhs_L);
            }
            if (ownership_R == 1 || ownership_R == 2) {
                int rowStart = d_A_rowptr[local_R];
                int rowEnd = d_A_rowptr[local_R + 1];
                atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, local_L, -lhs_R);
                atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, local_R, lhs_R);
            }
        }
        
        //=============
        // DIFFUSION
        //=============
        for (int n = 0; n < nnodes; ++n) {
            // Diffusion: -gamma * grad(N) · area_vec
            RealType diff_contrib = 0.0;
            for (int d = 0; d < ndim; ++d) {
                diff_contrib += -gamma_scs * dNdx[n][d] * area_vec[d];
            }
            
            // Get local DOF indices
            int local_L = d_local_dofs[elem_nodes[node_L]];
            int local_R = d_local_dofs[elem_nodes[node_R]];
            int local_N = d_local_dofs[elem_nodes[n]];
            int ownership_L = d_node_ownership[elem_nodes[node_L]];
            int ownership_R = d_node_ownership[elem_nodes[node_R]];
            
            // RHS contribution: -K*u
            RealType phi_n = ws_phi[n];
            if (ownership_L == 1 || ownership_L == 2) {
                atomicAdd(&d_b[local_L], -diff_contrib * phi_n);
            }
            if (ownership_R == 1 || ownership_R == 2) {
                atomicAdd(&d_b[local_R], diff_contrib * phi_n);
            }
            
            // LHS contribution: add to matrix
            if (ownership_L == 1 || ownership_L == 2) {
                int rowStart = d_A_rowptr[local_L];
                int rowEnd = d_A_rowptr[local_L + 1];
                atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, local_N, diff_contrib);
            }
            if (ownership_R == 1 || ownership_R == 2) {
                int rowStart = d_A_rowptr[local_R];
                int rowEnd = d_A_rowptr[local_R + 1];
                atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, local_N, -diff_contrib);
            }
        }
    }
}

} // namespace mars
