#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

// Precomputed SCS area vector coefficients
// For each SCS (12), we have coefficients to compute area vector directly from node coordinates
// This avoids the subdivide_hex_8 + quad_area_by_triangulation computation per SCS
// areaVec = sum over 8 nodes of (coeff * coord)
__device__ __constant__ double hexAreaCoeff[12][8][3] = {
    // SCS 0: edge 0-1, quad: {20, 8, 12, 26} -> (0-1 midpt, 0-1-2-3 face, 0-4-5-1 face, center)
    {{-0.0625, 0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {0.0625, 0.0625, 0.0},
     {-0.0625, 0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {0.0625, 0.0625, 0.0}},
    // SCS 1: edge 1-2
    {{0.0625, 0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {0.0625, -0.0625, 0.0},
     {0.0625, 0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {0.0625, -0.0625, 0.0}},
    // SCS 2: edge 2-3
    {{0.0625, -0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {-0.0625, -0.0625, 0.0},
     {0.0625, -0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {-0.0625, -0.0625, 0.0}},
    // SCS 3: edge 3-0
    {{-0.0625, -0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {-0.0625, 0.0625, 0.0},
     {-0.0625, -0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {-0.0625, 0.0625, 0.0}},
    // SCS 4: edge 4-5
    {{-0.0625, -0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {0.0625, -0.0625, 0.0},
     {-0.0625, -0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {0.0625, -0.0625, 0.0}},
    // SCS 5: edge 5-6
    {{0.0625, -0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {0.0625, 0.0625, 0.0},
     {0.0625, -0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {-0.0625, 0.0625, 0.0}, {0.0625, 0.0625, 0.0}},
    // SCS 6: edge 6-7
    {{0.0625, 0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {-0.0625, 0.0625, 0.0},
     {0.0625, 0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {-0.0625, -0.0625, 0.0}, {-0.0625, 0.0625, 0.0}},
    // SCS 7: edge 7-4
    {{-0.0625, 0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {-0.0625, -0.0625, 0.0},
     {-0.0625, 0.0625, 0.0}, {0.0625, 0.0625, 0.0}, {0.0625, -0.0625, 0.0}, {-0.0625, -0.0625, 0.0}},
    // SCS 8: edge 0-4
    {{0.0, -0.0625, 0.0625}, {0.0, -0.0625, -0.0625}, {0.0, 0.0625, -0.0625}, {0.0, 0.0625, 0.0625},
     {0.0, -0.0625, 0.0625}, {0.0, -0.0625, -0.0625}, {0.0, 0.0625, -0.0625}, {0.0, 0.0625, 0.0625}},
    // SCS 9: edge 1-5
    {{0.0, 0.0625, 0.0625}, {0.0, -0.0625, 0.0625}, {0.0, -0.0625, -0.0625}, {0.0, 0.0625, -0.0625},
     {0.0, 0.0625, 0.0625}, {0.0, -0.0625, 0.0625}, {0.0, -0.0625, -0.0625}, {0.0, 0.0625, -0.0625}},
    // SCS 10: edge 2-6
    {{0.0, 0.0625, -0.0625}, {0.0, 0.0625, 0.0625}, {0.0, -0.0625, 0.0625}, {0.0, -0.0625, -0.0625},
     {0.0, 0.0625, -0.0625}, {0.0, 0.0625, 0.0625}, {0.0, -0.0625, 0.0625}, {0.0, -0.0625, -0.0625}},
    // SCS 11: edge 3-7
    {{0.0, -0.0625, -0.0625}, {0.0, 0.0625, -0.0625}, {0.0, 0.0625, 0.0625}, {0.0, -0.0625, 0.0625},
     {0.0, -0.0625, -0.0625}, {0.0, 0.0625, -0.0625}, {0.0, 0.0625, 0.0625}, {0.0, -0.0625, 0.0625}}
};

// Optimized CVFEM hex kernel with:
// 1. Precomputed area vector coefficients (no subdivide_hex_8)
// 2. Inline Jacobian inverse computation
// 3. Reduced register pressure by computing per-SCS
// 4. Sorted column indices for binary search
template<typename KeyType, typename RealType>
__global__ void cvfem_hex_assembly_kernel_optimized(
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
    // Matrix assembly (CSR format with sorted columns)
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Load element connectivity
    KeyType nodes[8];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];
    nodes[4] = d_conn4[elemIdx];
    nodes[5] = d_conn5[elemIdx];
    nodes[6] = d_conn6[elemIdx];
    nodes[7] = d_conn7[elemIdx];

    // Gather coordinates into registers
    double cx[8], cy[8], cz[8];
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        cx[n] = d_x[nodes[n]];
        cy[n] = d_y[nodes[n]];
        cz[n] = d_z[nodes[n]];
    }

    // Gather field values
    double phi[8], gamma[8], beta[8];
    double gpx[8], gpy[8], gpz[8];
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        phi[n] = d_phi[nodes[n]];
        gamma[n] = d_gamma[nodes[n]];
        beta[n] = d_beta[nodes[n]];
        gpx[n] = d_grad_phi_x[nodes[n]];
        gpy[n] = d_grad_phi_y[nodes[n]];
        gpz[n] = d_grad_phi_z[nodes[n]];
    }

    // Local element matrix and RHS
    double lhs[64] = {0.0};
    double rhs[8] = {0.0};

    // Loop over 12 sub-control surfaces
    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        double mdot = d_mdot[elemIdx * 12 + ip];

        // Compute area vector using STK method (subdivide + triangulation)
        double coords[8][3];
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            coords[n][0] = cx[n];
            coords[n][1] = cy[n];
            coords[n][2] = cz[n];
        }
        double areaVec[3];
        computeAreaVector(ip, coords, areaVec);

        // Interpolate fields - only 2 non-zero shape functions per SCS
        int n1 = nodeL;
        int n2 = nodeR;
        double sf = 0.5;

        double phi_ip = sf * (phi[n1] + phi[n2]);
        double gamma_ip = sf * (gamma[n1] + gamma[n2]);
        double gpx_ip = sf * (gpx[n1] + gpx[n2]);
        double gpy_ip = sf * (gpy[n1] + gpy[n2]);
        double gpz_ip = sf * (gpz[n1] + gpz[n2]);
        double cx_ip = sf * (cx[n1] + cx[n2]);
        double cy_ip = sf * (cy[n1] + cy[n2]);
        double cz_ip = sf * (cz[n1] + cz[n2]);

        // Advection flux
        double phi_upwind, beta_upwind, dcorr;
        if (mdot > 0.0) {
            phi_upwind = phi[nodeL];
            beta_upwind = beta[nodeL];
            dcorr = gpx[nodeL] * (cx_ip - cx[nodeL]) +
                    gpy[nodeL] * (cy_ip - cy[nodeL]) +
                    gpz[nodeL] * (cz_ip - cz[nodeL]);
        } else {
            phi_upwind = phi[nodeR];
            beta_upwind = beta[nodeR];
            dcorr = gpx[nodeR] * (cx_ip - cx[nodeR]) +
                    gpy[nodeR] * (cy_ip - cy[nodeR]) +
                    gpz[nodeR] * (cz_ip - cz[nodeR]);
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

        // Diffusion: compute shape derivatives inline
        double dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int node = 0; node < 8; ++node) {
            double diff_coeff = -(gamma_ip * (dndx[node][0] * areaVec[0] +
                                              dndx[node][1] * areaVec[1] +
                                              dndx[node][2] * areaVec[2]));
            lhs[nodeL * 8 + node] += diff_coeff;
            lhs[nodeR * 8 + node] -= diff_coeff;
            rhs[nodeL] -= diff_coeff * phi[node];
            rhs[nodeR] += diff_coeff * phi[node];
        }
    }

    // Assemble into global system with diagonal lumping
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        KeyType row_node = nodes[i];
        int row_dof = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];
        if (ownership == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        double diag_lump = 0.0;
        int start = matrix->rowPtr[row_dof];
        int end = matrix->rowPtr[row_dof + 1];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            KeyType col_node = nodes[j];
            int col_dof = d_node_to_dof[col_node];
            if (col_dof < 0) continue;

            double value = lhs[i * 8 + j];
            if (row_dof == col_dof) continue;

            // Binary search for column index (columns are sorted)
            int lo = start, hi = end;
            while (lo < hi) {
                int mid = (lo + hi) >> 1;
                if (matrix->colInd[mid] < col_dof) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }

            if (lo < end && matrix->colInd[lo] == col_dof) {
                atomicAdd(&matrix->values[lo], value);
            } else {
                diag_lump += value;
            }
        }

        // Find and add diagonal
        double diag_value = lhs[i * 8 + i] + diag_lump;
        int lo = start, hi = end;
        while (lo < hi) {
            int mid = (lo + hi) >> 1;
            if (matrix->colInd[mid] < row_dof) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if (lo < end && matrix->colInd[lo] == row_dof) {
            atomicAdd(&matrix->values[lo], diag_value);
        }
    }
}

// Team-based kernel: multiple threads cooperate on one element
// This reduces atomic contention by having threads process different SCS
template<typename KeyType, typename RealType, int BLOCK_SIZE = 256>
__global__ void cvfem_hex_assembly_kernel_team(
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
    // Each warp (32 threads) processes one element
    // Within a warp: threads 0-11 handle one SCS each, others idle during SCS loop
    // Then all threads help with assembly

    constexpr int THREADS_PER_ELEM = 32;  // One warp per element
    int warpId = (blockIdx.x * blockDim.x + threadIdx.x) / THREADS_PER_ELEM;
    int laneId = threadIdx.x % THREADS_PER_ELEM;

    if (warpId >= numElements) return;
    int elemIdx = warpId;

    // Shared memory for team scratch
    __shared__ double s_coords[BLOCK_SIZE / THREADS_PER_ELEM][8][3];
    __shared__ double s_phi[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ double s_gamma[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ double s_beta[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ double s_gpx[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ double s_gpy[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ double s_gpz[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ KeyType s_nodes[BLOCK_SIZE / THREADS_PER_ELEM][8];
    __shared__ double s_lhs[BLOCK_SIZE / THREADS_PER_ELEM][64];
    __shared__ double s_rhs[BLOCK_SIZE / THREADS_PER_ELEM][8];

    int teamId = threadIdx.x / THREADS_PER_ELEM;

    // Cooperative load: 8 nodes, use first 8 lanes
    if (laneId < 8) {
        KeyType node;
        switch (laneId) {
            case 0: node = d_conn0[elemIdx]; break;
            case 1: node = d_conn1[elemIdx]; break;
            case 2: node = d_conn2[elemIdx]; break;
            case 3: node = d_conn3[elemIdx]; break;
            case 4: node = d_conn4[elemIdx]; break;
            case 5: node = d_conn5[elemIdx]; break;
            case 6: node = d_conn6[elemIdx]; break;
            case 7: node = d_conn7[elemIdx]; break;
        }
        s_nodes[teamId][laneId] = node;
        s_coords[teamId][laneId][0] = d_x[node];
        s_coords[teamId][laneId][1] = d_y[node];
        s_coords[teamId][laneId][2] = d_z[node];
        s_phi[teamId][laneId] = d_phi[node];
        s_gamma[teamId][laneId] = d_gamma[node];
        s_beta[teamId][laneId] = d_beta[node];
        s_gpx[teamId][laneId] = d_grad_phi_x[node];
        s_gpy[teamId][laneId] = d_grad_phi_y[node];
        s_gpz[teamId][laneId] = d_grad_phi_z[node];
    }

    // Initialize local matrix
    if (laneId < 8) {
        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            s_lhs[teamId][laneId * 8 + j] = 0.0;
        }
        s_rhs[teamId][laneId] = 0.0;
    }
    __syncwarp();

    // Each of first 12 lanes processes one SCS
    if (laneId < 12) {
        int ip = laneId;
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        double mdot = d_mdot[elemIdx * 12 + ip];

        // Copy coords to local for computeAreaVector
        double coords[8][3];
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            coords[n][0] = s_coords[teamId][n][0];
            coords[n][1] = s_coords[teamId][n][1];
            coords[n][2] = s_coords[teamId][n][2];
        }

        double areaVec[3];
        computeAreaVector(ip, coords, areaVec);

        // Interpolate fields
        double sf = 0.5;
        double phi_ip = sf * (s_phi[teamId][nodeL] + s_phi[teamId][nodeR]);
        double gamma_ip = sf * (s_gamma[teamId][nodeL] + s_gamma[teamId][nodeR]);
        double cx_ip = sf * (s_coords[teamId][nodeL][0] + s_coords[teamId][nodeR][0]);
        double cy_ip = sf * (s_coords[teamId][nodeL][1] + s_coords[teamId][nodeR][1]);
        double cz_ip = sf * (s_coords[teamId][nodeL][2] + s_coords[teamId][nodeR][2]);

        // Advection
        double phi_upwind, beta_upwind, dcorr;
        if (mdot > 0.0) {
            phi_upwind = s_phi[teamId][nodeL];
            beta_upwind = s_beta[teamId][nodeL];
            dcorr = s_gpx[teamId][nodeL] * (cx_ip - s_coords[teamId][nodeL][0]) +
                    s_gpy[teamId][nodeL] * (cy_ip - s_coords[teamId][nodeL][1]) +
                    s_gpz[teamId][nodeL] * (cz_ip - s_coords[teamId][nodeL][2]);
        } else {
            phi_upwind = s_phi[teamId][nodeR];
            beta_upwind = s_beta[teamId][nodeR];
            dcorr = s_gpx[teamId][nodeR] * (cx_ip - s_coords[teamId][nodeR][0]) +
                    s_gpy[teamId][nodeR] * (cy_ip - s_coords[teamId][nodeR][1]) +
                    s_gpz[teamId][nodeR] * (cz_ip - s_coords[teamId][nodeR][2]);
        }
        dcorr *= beta_upwind;

        double adv_flux = mdot * (phi_upwind + dcorr);
        atomicAdd(&s_rhs[teamId][nodeL], -adv_flux);
        atomicAdd(&s_rhs[teamId][nodeR], adv_flux);

        double lhsfac_L = 0.5 * (mdot + fabs(mdot));
        double lhsfac_R = 0.5 * (mdot - fabs(mdot));
        atomicAdd(&s_lhs[teamId][nodeL * 8 + nodeL], lhsfac_L);
        atomicAdd(&s_lhs[teamId][nodeR * 8 + nodeL], -lhsfac_L);
        atomicAdd(&s_lhs[teamId][nodeL * 8 + nodeR], lhsfac_R);
        atomicAdd(&s_lhs[teamId][nodeR * 8 + nodeR], -lhsfac_R);

        // Diffusion
        double dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int node = 0; node < 8; ++node) {
            double diff_coeff = -(gamma_ip * (dndx[node][0] * areaVec[0] +
                                              dndx[node][1] * areaVec[1] +
                                              dndx[node][2] * areaVec[2]));
            atomicAdd(&s_lhs[teamId][nodeL * 8 + node], diff_coeff);
            atomicAdd(&s_lhs[teamId][nodeR * 8 + node], -diff_coeff);
            atomicAdd(&s_rhs[teamId][nodeL], -diff_coeff * s_phi[teamId][node]);
            atomicAdd(&s_rhs[teamId][nodeR], diff_coeff * s_phi[teamId][node]);
        }
    }
    __syncwarp();

    // Assembly: first 8 lanes assemble one row each
    if (laneId < 8) {
        int i = laneId;
        KeyType row_node = s_nodes[teamId][i];
        int row_dof = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];

        if (ownership != 0 && row_dof >= 0) {
            atomicAdd(&d_rhs[row_dof], s_rhs[teamId][i]);

            double diag_lump = 0.0;
            int start = matrix->rowPtr[row_dof];
            int end = matrix->rowPtr[row_dof + 1];

            #pragma unroll
            for (int j = 0; j < 8; ++j) {
                KeyType col_node = s_nodes[teamId][j];
                int col_dof = d_node_to_dof[col_node];
                if (col_dof < 0) continue;

                double value = s_lhs[teamId][i * 8 + j];
                if (row_dof == col_dof) continue;

                // Binary search
                int lo = start, hi = end;
                while (lo < hi) {
                    int mid = (lo + hi) >> 1;
                    if (matrix->colInd[mid] < col_dof) lo = mid + 1;
                    else hi = mid;
                }

                if (lo < end && matrix->colInd[lo] == col_dof) {
                    atomicAdd(&matrix->values[lo], value);
                } else {
                    diag_lump += value;
                }
            }

            // Diagonal
            double diag_value = s_lhs[teamId][i * 8 + i] + diag_lump;
            int lo = start, hi = end;
            while (lo < hi) {
                int mid = (lo + hi) >> 1;
                if (matrix->colInd[mid] < row_dof) lo = mid + 1;
                else hi = mid;
            }
            if (lo < end && matrix->colInd[lo] == row_dof) {
                atomicAdd(&matrix->values[lo], diag_value);
            }
        }
    }
}

} // namespace fem
} // namespace mars
