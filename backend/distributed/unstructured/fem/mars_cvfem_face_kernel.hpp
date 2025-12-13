#pragma once

#include "mars_base.hpp"
#include "mars_sparse_matrix.hpp"
#include <cuda_runtime.h>

namespace mars {

// Face-based CVFEM assembly kernel (production implementation)
// Loops over faces instead of elements, computes fluxes across faces
template<typename KeyType, typename RealType>
__global__ void cvfem_face_assembly_kernel(
    // Face topology
    const KeyType* __restrict__ d_faceNodes,        // [numFaces * nodesPerFace] SFC keys
    const KeyType* __restrict__ d_faceToElemOffsets, // [numFaces + 1] CSR offsets
    const KeyType* __restrict__ d_faceToElemList,    // [numInteriorFaces*2 + numBoundaryFaces] element IDs
    const uint8_t* __restrict__ d_isBoundaryFace,   // [numFaces]
    const RealType* __restrict__ d_faceNormalX,     // [numFaces]
    const RealType* __restrict__ d_faceNormalY,     // [numFaces]
    const RealType* __restrict__ d_faceNormalZ,     // [numFaces]
    const RealType* __restrict__ d_faceArea,        // [numFaces]
    size_t numFaces,
    int nodesPerFace,
    // Node data
    const KeyType* __restrict__ d_sfc_to_local,
    size_t nodeCount,
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    // Field data
    const RealType* __restrict__ d_gamma,      // diffusion coefficient [nodeCount]
    const RealType* __restrict__ d_phi,        // solution [nodeCount]
    const RealType* __restrict__ d_beta,       // upwind param [nodeCount]
    const RealType* __restrict__ d_grad_phi_x, // gradient [nodeCount]
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    const RealType* __restrict__ d_mdot,       // mass flow [numFaces]
    // Sparse matrix assembly
    RealType* __restrict__ d_A_values,         // sparse matrix values
    const int* __restrict__ d_A_colidx,        // CSR column indices
    const int* __restrict__ d_A_rowptr,        // CSR row pointers
    const int* __restrict__ d_local_dofs,      // local DOF indices for nodes
    const int* __restrict__ d_node_ownership,  // node ownership (0=ghost, 1=owned, 2=shared)
    // RHS vector
    RealType* __restrict__ d_b,                // RHS vector
    bool include_advection)
{
    int face_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (face_id >= numFaces) return;
    
    // Get face nodes (3 for tet triangular face, 4 for hex quad face)
    // constexpr int nodesPerFace = 3;  // Triangular face for tetrahedra
    constexpr int ndim = 3;
    
    KeyType sfc_face_nodes[4];  // max 4
    int local_face_nodes[4];
    RealType face_node_coords[4][ndim];
    RealType phi_face[4], gamma_face[4], beta_face[4];
    RealType grad_phi_face[4][ndim];
    
    // Gather face node data
    for (int n = 0; n < nodesPerFace; ++n) {
        KeyType sfc_key = d_faceNodes[face_id * nodesPerFace + n];
        sfc_face_nodes[n] = sfc_key;
        
        // Binary search for local node ID
        int left = 0, right = nodeCount - 1;
        int local_id = -1;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (d_sfc_to_local[mid] == sfc_key) {
                local_id = mid;
                break;
            }
            if (d_sfc_to_local[mid] < sfc_key) left = mid + 1;
            else right = mid - 1;
        }
        
        local_face_nodes[n] = local_id;
        face_node_coords[n][0] = d_x[local_id];
        face_node_coords[n][1] = d_y[local_id];
        face_node_coords[n][2] = d_z[local_id];
        phi_face[n] = d_phi[local_id];
        gamma_face[n] = d_gamma[local_id];
        beta_face[n] = d_beta[local_id];
        grad_phi_face[n][0] = d_grad_phi_x[local_id];
        grad_phi_face[n][1] = d_grad_phi_y[local_id];
        grad_phi_face[n][2] = d_grad_phi_z[local_id];
    }
    
    // Get face geometry
    RealType normal[ndim] = {d_faceNormalX[face_id], d_faceNormalY[face_id], d_faceNormalZ[face_id]};
    RealType area = d_faceArea[face_id];
    
    // Face centroid
    RealType face_centroid[ndim] = {0, 0, 0};
    for (int n = 0; n < nodesPerFace; ++n) {
        for (int d = 0; d < ndim; ++d) {
            face_centroid[d] += face_node_coords[n][d] / nodesPerFace;
        }
    }
    
    // Interpolate fields to face centroid (simple average for now)
    RealType gamma_f = 0, phi_f = 0, beta_f = 0;
    RealType grad_phi_f[ndim] = {0, 0, 0};
    for (int n = 0; n < nodesPerFace; ++n) {
        gamma_f += gamma_face[n] / nodesPerFace;
        phi_f += phi_face[n] / nodesPerFace;
        beta_f += beta_face[n] / nodesPerFace;
        for (int d = 0; d < ndim; ++d) {
            grad_phi_f[d] += grad_phi_face[n][d] / nodesPerFace;
        }
    }
    
    // Get adjacent elements (left and right)
    int offset_start = d_faceToElemOffsets[face_id];
    int offset_end = d_faceToElemOffsets[face_id + 1];
    
    // Mass flow rate across face
    RealType mdot = d_mdot[face_id];
    
    //===============================================
    // DIFFUSION FLUX: -gamma * grad(phi) · normal
    //===============================================
    RealType diff_flux = 0;
    for (int d = 0; d < ndim; ++d) {
        diff_flux += -gamma_f * grad_phi_f[d] * normal[d] * area;
    }
    
    //==============================================
    // ADVECTION FLUX: mdot * phi_upwind
    //==============================================
    RealType adv_flux = 0;
    if (include_advection && fabs(mdot) > 1e-14) {
        // Upwind: use phi from upstream side
        RealType phi_upw = (mdot > 0) ? phi_f : phi_f;  // Simplified - need proper L/R element interpolation
        adv_flux = mdot * phi_upw;
    }
    
    RealType total_flux = diff_flux + adv_flux;
    
    //===============================================
    // ASSEMBLE INTO GLOBAL SYSTEM
    //===============================================
    
    // For each node on the face, contribute flux
    // In CVFEM, flux is distributed to face nodes based on shape functions
    for (int n = 0; n < nodesPerFace; ++n) {
        int local_node = local_face_nodes[n];
        int local_dof = d_local_dofs[local_node];
        int ownership = d_node_ownership[local_node];
        
        if (ownership != 1 && ownership != 2) continue;  // Skip ghost nodes
        
        // Shape function weight (equal distribution for simplicity)
        RealType weight = 1.0 / nodesPerFace;
        
        // RHS contribution
        atomicAdd(&d_b[local_dof], -total_flux * weight);
        
        // LHS contribution (diffusion matrix)
        int rowStart = d_A_rowptr[local_dof];
        int rowEnd = d_A_rowptr[local_dof + 1];
        
        // For each face node, add diffusion stencil
        for (int m = 0; m < nodesPerFace; ++m) {
            int local_node_m = local_face_nodes[m];
            int local_dof_m = d_local_dofs[local_node_m];
            
            // Diffusion contribution: grad(N_m) · normal
            // Simplified: use constant gradient approximation
            RealType diff_coeff = gamma_f * area * weight;
            
            fem::atomicAddSparseEntry(d_A_values, d_A_colidx, rowStart, rowEnd, 
                               local_dof_m, diff_coeff);
        }
    }
}

} // namespace mars
