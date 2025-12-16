#pragma once

#include "mars_base.hpp"
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// Simple CSR matrix wrapper for GPU assembly
template<typename RealType>
struct CSRMatrix {
    int* rowPtr;
    int* colInd;
    RealType* values;
    int numRows;
    int nnz;

    __device__ void addValue(int row, int col, RealType val) {
        int start = rowPtr[row];
        int end = rowPtr[row + 1];
        for (int i = start; i < end; ++i) {
            if (colInd[i] == col) {
                atomicAdd(&values[i], val);
                return;
            }
        }
    }
};

// Hex element SCS shape functions at face centers (6 faces for hex8)
__device__ __constant__ double hexSfConst[6][8] = {
    {0.25, 0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.25}, // ip 0: face ksi=-1
    {0.0, 0.25, 0.25, 0.0, 0.0, 0.25, 0.25, 0.0}, // ip 1: face ksi=1
    {0.25, 0.25, 0.0, 0.0, 0.25, 0.25, 0.0, 0.0}, // ip 2: face eta=-1
    {0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.25, 0.25}, // ip 3: face eta=1
    {0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0}, // ip 4: face zeta=-1
    {0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25}  // ip 5: face zeta=1
};

// Left/right sub-control volume pairs for each SCS integration point
__device__ __constant__ int hexLRSCV[24] = {
    0,1, 1,2, 2,3, 0,3,  // faces in ksi-eta plane
    4,5, 5,6, 6,7, 4,7,  // faces in ksi-eta plane (top)
    0,4, 1,5, 2,6, 3,7   // faces connecting bottom to top
};

// Hex element parametric derivatives at SCS face centers
__device__ __constant__ double hexDerivConst[6][8][3] = {
    // ip 0: ksi=-1, eta=0, zeta=0
    {{-0.125, -0.25, -0.125}, {-0.125, 0, 0}, {-0.125, 0, 0}, {-0.125, -0.25, -0.25},
     {-0.25, -0.5, 0.25}, {-0.25, 0, 0}, {-0.25, 0, 0}, {-0.25, -0.5, 0.25}},
    // ip 1: ksi=1, eta=0, zeta=0
    {{-0.125, 0, 0}, {-0.125, -0.25, -0.25}, {-0.125, -0.25, -0.25}, {-0.125, 0, 0},
     {-0.25, 0, 0}, {-0.25, -0.5, 0.25}, {-0.25, -0.5, 0.25}, {-0.25, 0, 0}},
    // ip 2: ksi=0, eta=-1, zeta=0
    {{-0.25, -0.125, -0.25}, {-0.25, -0.125, -0.25}, {0, -0.125, 0}, {0, -0.125, 0},
     {-0.5, -0.25, 0.25}, {-0.5, -0.25, 0.25}, {0, -0.25, 0}, {0, -0.25, 0}},
    // ip 3: ksi=0, eta=1, zeta=0
    {{0, -0.125, 0}, {0, -0.125, 0}, {-0.25, -0.125, -0.25}, {-0.25, -0.125, -0.25},
     {0, -0.25, 0}, {0, -0.25, 0}, {-0.5, -0.25, 0.25}, {-0.5, -0.25, 0.25}},
    // ip 4: ksi=0, eta=0, zeta=-1
    {{-0.25, -0.25, -0.125}, {-0.25, -0.25, -0.125}, {-0.25, -0.25, -0.125}, {-0.25, -0.25, -0.125},
     {0, 0, -0.125}, {0, 0, -0.125}, {0, 0, -0.125}, {0, 0, -0.125}},
    // ip 5: ksi=0, eta=0, zeta=1
    {{0, 0, -0.125}, {0, 0, -0.125}, {0, 0, -0.125}, {0, 0, -0.125},
     {-0.25, -0.25, -0.125}, {-0.25, -0.25, -0.125}, {-0.25, -0.25, -0.125}, {-0.25, -0.25, -0.125}}
};

// 3x3 matrix inversion
__device__ inline void invert3x3(double J[3][3], double invJ[3][3]) {
    double det = J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
               - J[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
               + J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    double invDet = 1.0 / det;
    invJ[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * invDet;
    invJ[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2]) * invDet;
    invJ[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * invDet;
    invJ[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2]) * invDet;
    invJ[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * invDet;
    invJ[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2]) * invDet;
    invJ[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * invDet;
    invJ[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1]) * invDet;
    invJ[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * invDet;
}

// Compute shape function derivatives in physical space (dN/dx, dN/dy, dN/dz)
__device__ inline void computeShapeDerivatives(
    int ip,
    const double coords[8][3],
    double dndx[8][3])
{
    // Compute Jacobian at integration point ip
    double J[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

    for (int node = 0; node < 8; ++node) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += hexDerivConst[ip][node][j] * coords[node][i];
            }
        }
    }

    // Invert Jacobian
    double invJ[3][3];
    invert3x3(J, invJ);

    // Transform parametric derivatives to physical derivatives
    for (int node = 0; node < 8; ++node) {
        for (int i = 0; i < 3; ++i) {
            dndx[node][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dndx[node][i] += invJ[i][j] * hexDerivConst[ip][node][j];
            }
        }
    }
}

// Subdivide hex8 into 27 points (8 corners + 12 edge mids + 6 face centers + 1 volume center)
__device__ inline void subdivide_hex_8(const double coords[8][3], double coordv[27][3])
{
    // Copy 8 corner nodes
    for (int n = 0; n < 8; ++n) {
        for (int d = 0; d < 3; ++d) {
            coordv[n][d] = coords[n][d];
        }
    }

    for (int d = 0; d < 3; ++d) {
        // 12 edge midpoints
        coordv[8][d] = 0.5 * (coords[0][d] + coords[1][d]);   // edge 1
        coordv[9][d] = 0.5 * (coords[1][d] + coords[2][d]);   // edge 2
        coordv[10][d] = 0.5 * (coords[2][d] + coords[3][d]);  // edge 3
        coordv[11][d] = 0.5 * (coords[3][d] + coords[0][d]);  // edge 4
        coordv[13][d] = 0.5 * (coords[4][d] + coords[5][d]);  // edge 5
        coordv[14][d] = 0.5 * (coords[5][d] + coords[6][d]);  // edge 6
        coordv[15][d] = 0.5 * (coords[6][d] + coords[7][d]);  // edge 7
        coordv[16][d] = 0.5 * (coords[7][d] + coords[4][d]);  // edge 8
        coordv[18][d] = 0.5 * (coords[1][d] + coords[5][d]);  // edge 9
        coordv[19][d] = 0.5 * (coords[0][d] + coords[4][d]);  // edge 10
        coordv[21][d] = 0.5 * (coords[3][d] + coords[7][d]);  // edge 11
        coordv[22][d] = 0.5 * (coords[2][d] + coords[6][d]);  // edge 12

        // 6 face centers
        coordv[12][d] = 0.25 * (coords[0][d] + coords[1][d] + coords[2][d] + coords[3][d]); // face 0
        coordv[17][d] = 0.25 * (coords[4][d] + coords[5][d] + coords[6][d] + coords[7][d]); // face 1
        coordv[20][d] = 0.25 * (coords[0][d] + coords[1][d] + coords[4][d] + coords[5][d]); // face 2
        coordv[23][d] = 0.25 * (coords[2][d] + coords[3][d] + coords[6][d] + coords[7][d]); // face 3
        coordv[24][d] = 0.25 * (coords[1][d] + coords[2][d] + coords[5][d] + coords[6][d]); // face 4
        coordv[25][d] = 0.25 * (coords[0][d] + coords[3][d] + coords[4][d] + coords[7][d]); // face 5

        // Volume centroid
        coordv[26][d] = 0.0;
        for (int n = 0; n < 8; ++n) {
            coordv[26][d] += coords[n][d];
        }
        coordv[26][d] *= 0.125;
    }
}

// Compute area vector by quad triangulation (Grandy algorithm)
__device__ inline void quad_area_by_triangulation(const double areacoords[4][3], double areaVec[3])
{
    areaVec[0] = 0.0;
    areaVec[1] = 0.0;
    areaVec[2] = 0.0;

    // Quad centroid
    double xmid[3] = {
        0.25 * (areacoords[0][0] + areacoords[1][0] + areacoords[2][0] + areacoords[3][0]),
        0.25 * (areacoords[0][1] + areacoords[1][1] + areacoords[2][1] + areacoords[3][1]),
        0.25 * (areacoords[0][2] + areacoords[1][2] + areacoords[2][2] + areacoords[3][2])
    };

    double r1[3] = {
        areacoords[0][0] - xmid[0],
        areacoords[0][1] - xmid[1],
        areacoords[0][2] - xmid[2]
    };

    // Sum cross products for 4 triangles
    for (int itri = 0; itri < 4; ++itri) {
        int t_index = (itri + 1) % 4;
        double r2[3] = {
            areacoords[t_index][0] - xmid[0],
            areacoords[t_index][1] - xmid[1],
            areacoords[t_index][2] - xmid[2]
        };

        // Cross product r1 x r2
        areaVec[0] += r1[1] * r2[2] - r2[1] * r1[2];
        areaVec[1] += r1[2] * r2[0] - r2[2] * r1[0];
        areaVec[2] += r1[0] * r2[1] - r2[0] * r1[1];

        r1[0] = r2[0];
        r1[1] = r2[1];
        r1[2] = r2[2];
    }

    areaVec[0] *= 0.5;
    areaVec[1] *= 0.5;
    areaVec[2] *= 0.5;
}

// SCS face quad node indices (into 27-point subdivision)
__device__ __constant__ int hexEdgeFacetTable[12][4] = {
    {20, 8, 12, 26},  {24, 9, 12, 26},  {10, 12, 26, 23}, {11, 25, 26, 12},
    {13, 20, 26, 17}, {17, 14, 24, 26}, {17, 15, 23, 26}, {16, 17, 26, 25},
    {19, 20, 26, 25}, {20, 18, 24, 26}, {22, 23, 26, 24}, {21, 25, 26, 23}
};

// Compute area vector for a sub-control surface using STK method
__device__ inline void computeAreaVector(
    int scs_id,
    const double coords[8][3],
    double areaVec[3])
{
    // Subdivide hex into 27 points
    double coordv[27][3];
    subdivide_hex_8(coords, coordv);

    // Get the 4 quad nodes for this SCS
    double scscoords[4][3];
    for (int inode = 0; inode < 4; ++inode) {
        int idx = hexEdgeFacetTable[scs_id][inode];
        for (int d = 0; d < 3; ++d) {
            scscoords[inode][d] = coordv[idx][d];
        }
    }

    // Compute area vector by triangulation
    quad_area_by_triangulation(scscoords, areaVec);
}

// CVFEM assembly kernel for hex elements
// Loops over elements, then over 12 sub-control surfaces (SCS) per element
template<typename KeyType, typename RealType>
__global__ void cvfem_hex_assembly_kernel(
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
    const RealType* __restrict__ d_gamma,      // diffusion coefficient
    const RealType* __restrict__ d_phi,        // solution field
    const RealType* __restrict__ d_beta,       // upwind limiting parameter
    const RealType* __restrict__ d_grad_phi_x, // gradient of phi
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    // Element data (per element, per face)
    const RealType* __restrict__ d_mdot,       // mass flow rate [numElements * 12]
    const RealType* __restrict__ d_areaVec_x,  // area vector [numElements * 12]
    const RealType* __restrict__ d_areaVec_y,
    const RealType* __restrict__ d_areaVec_z,
    // DOF mapping
    const int* __restrict__ d_node_to_dof,     // local node ID -> DOF ID
    const uint8_t* __restrict__ d_ownership,   // 0=ghost, 1=owned, 2=shared
    // Matrix assembly
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
    double lhs[64] = {0.0};  // 8*8
    double rhs[8] = {0.0};

    // Loop over 12 sub-control surfaces (SCS)
    for (int ip = 0; ip < 12; ++ip) {
        // Get left and right sub-control volumes
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        // Get mass flow rate (use input if available, else dummy)
        double mdot = d_mdot[elemIdx * 12 + ip];

        // Compute area vector from geometry using STK method
        double areaVec[3];
        computeAreaVector(ip, coords, areaVec);

        int face_id = ip / 2;  // 6 faces, 2 SCS per face (for shape functions)

        // Interpolate fields to SCS integration point
        double phi_ip = 0.0, gamma_ip = 0.0;
        double grad_phi_ip[3] = {0.0, 0.0, 0.0};
        double coords_ip[3] = {0.0, 0.0, 0.0};

        for (int node = 0; node < 8; ++node) {
            double sf = hexSfConst[face_id][node];
            phi_ip += sf * phi[node];
            gamma_ip += sf * gamma[node];
            for (int d = 0; d < 3; ++d) {
                grad_phi_ip[d] += sf * grad_phi[node][d];
                coords_ip[d] += sf * coords[node][d];
            }
        }

        //======================================
        // ADVECTION FLUX
        //======================================
        double phi_L = phi[nodeL];
        double phi_R = phi[nodeR];
        double phi_upwind, beta_upwind;
        double dcorr = 0.0;

        if (mdot > 0.0) {
            phi_upwind = phi_L;
            beta_upwind = beta[nodeL];
            // Higher-order correction
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

        // RHS contribution from advection
        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        // LHS contribution from advection (linearization)
        double lhsfac_L = 0.5 * (mdot + fabs(mdot));
        double lhsfac_R = 0.5 * (mdot - fabs(mdot));

        lhs[nodeL * 8 + nodeL] += lhsfac_L;
        lhs[nodeR * 8 + nodeL] -= lhsfac_L;
        lhs[nodeL * 8 + nodeR] += lhsfac_R;
        lhs[nodeR * 8 + nodeR] -= lhsfac_R;

        //======================================
        // DIFFUSION FLUX
        //======================================
        // Compute shape function derivatives at this SCS
        double dndx[8][3];
        computeShapeDerivatives(face_id, coords, dndx);

        // Diffusion flux contribution for each node
        for (int node = 0; node < 8; ++node) {
            double diff_coeff = 0.0;
            for (int d = 0; d < 3; ++d) {
                diff_coeff -= gamma_ip * dndx[node][d] * areaVec[d];
            }

            // Assemble into local matrix and RHS
            lhs[nodeL * 8 + node] += diff_coeff;
            lhs[nodeR * 8 + node] -= diff_coeff;

            rhs[nodeL] -= diff_coeff * phi[node];
            rhs[nodeR] += diff_coeff * phi[node];
        }
    }

    //======================================
    // ASSEMBLE INTO GLOBAL SYSTEM
    //======================================
    for (int i = 0; i < 8; ++i) {
        KeyType row_node = nodes[i];
        int row_dof = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];

        // Only assemble owned and shared nodes
        if (ownership == 0 || row_dof < 0) continue;  // Skip pure ghost nodes

        // Assemble RHS
        atomicAdd(&d_rhs[row_dof], rhs[i]);

        // Assemble LHS
        for (int j = 0; j < 8; ++j) {
            KeyType col_node = nodes[j];
            int col_dof = d_node_to_dof[col_node];

            // Skip ghost columns
            if (col_dof < 0) continue;

            // Find col_dof in row's sparsity pattern
            matrix->addValue(row_dof, col_dof, lhs[i * 8 + j]);
        }
    }
}

} // namespace fem
} // namespace mars
