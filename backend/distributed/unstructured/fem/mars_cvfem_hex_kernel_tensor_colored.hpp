#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_common.hpp"

namespace mars {
namespace fem {

// =============================================================================
// Atomic-free CVFEM tensor kernel for graph-colored assembly.
//
// Identical to cvfem_hex_assembly_kernel_tensor except:
//   1. Element index comes from d_elemList[linearIdx] (coloring permutation).
//   2. All atomicAdd calls are replaced by direct += writes.
//
// Correctness: within a single color, no two elements share a node, so the
// same DOF is never written by two threads simultaneously.  Direct writes
// are therefore race-free without atomics.
//
// Usage: launch once per color with d_elemList = colorPerm + colorOffset[c],
//        numThisColor = colorOffset[c+1] - colorOffset[c].
//        Matrix/RHS must be zeroed before the first color's launch.
// =============================================================================
template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void cvfem_hex_assembly_kernel_tensor_colored(
    const int*     __restrict__ d_elemList,   // permutation: linearIdx → elemIdx
    size_t numThisColor,
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
    const int*     __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    const int linearIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (static_cast<size_t>(linearIdx) >= numThisColor) return;

    // Remap to actual element index via coloring permutation
    const int elemIdx = d_elemList[linearIdx];

    // ------------------------------------------------------------------
    // Load element connectivity
    // ------------------------------------------------------------------
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

    // Nodes 0-3
    {
        KeyType n0=nodes[0], n1=nodes[1], n2=nodes[2], n3=nodes[3];
        coords[0][0]=d_x[n0]; coords[0][1]=d_y[n0]; coords[0][2]=d_z[n0];
        coords[1][0]=d_x[n1]; coords[1][1]=d_y[n1]; coords[1][2]=d_z[n1];
        coords[2][0]=d_x[n2]; coords[2][1]=d_y[n2]; coords[2][2]=d_z[n2];
        coords[3][0]=d_x[n3]; coords[3][1]=d_y[n3]; coords[3][2]=d_z[n3];
        phi[0]=d_phi[n0]; phi[1]=d_phi[n1]; phi[2]=d_phi[n2]; phi[3]=d_phi[n3];
        gamma[0]=d_gamma[n0]; gamma[1]=d_gamma[n1]; gamma[2]=d_gamma[n2]; gamma[3]=d_gamma[n3];
        beta[0]=d_beta[n0]; beta[1]=d_beta[n1]; beta[2]=d_beta[n2]; beta[3]=d_beta[n3];
        grad_phi[0][0]=d_grad_phi_x[n0]; grad_phi[0][1]=d_grad_phi_y[n0]; grad_phi[0][2]=d_grad_phi_z[n0];
        grad_phi[1][0]=d_grad_phi_x[n1]; grad_phi[1][1]=d_grad_phi_y[n1]; grad_phi[1][2]=d_grad_phi_z[n1];
        grad_phi[2][0]=d_grad_phi_x[n2]; grad_phi[2][1]=d_grad_phi_y[n2]; grad_phi[2][2]=d_grad_phi_z[n2];
        grad_phi[3][0]=d_grad_phi_x[n3]; grad_phi[3][1]=d_grad_phi_y[n3]; grad_phi[3][2]=d_grad_phi_z[n3];
        dofs[0]=d_node_to_dof[n0]; dofs[1]=d_node_to_dof[n1];
        dofs[2]=d_node_to_dof[n2]; dofs[3]=d_node_to_dof[n3];
        own[0]=d_ownership[n0]; own[1]=d_ownership[n1];
        own[2]=d_ownership[n2]; own[3]=d_ownership[n3];
    }
    // Nodes 4-7
    {
        KeyType n4=nodes[4], n5=nodes[5], n6=nodes[6], n7=nodes[7];
        coords[4][0]=d_x[n4]; coords[4][1]=d_y[n4]; coords[4][2]=d_z[n4];
        coords[5][0]=d_x[n5]; coords[5][1]=d_y[n5]; coords[5][2]=d_z[n5];
        coords[6][0]=d_x[n6]; coords[6][1]=d_y[n6]; coords[6][2]=d_z[n6];
        coords[7][0]=d_x[n7]; coords[7][1]=d_y[n7]; coords[7][2]=d_z[n7];
        phi[4]=d_phi[n4]; phi[5]=d_phi[n5]; phi[6]=d_phi[n6]; phi[7]=d_phi[n7];
        gamma[4]=d_gamma[n4]; gamma[5]=d_gamma[n5]; gamma[6]=d_gamma[n6]; gamma[7]=d_gamma[n7];
        beta[4]=d_beta[n4]; beta[5]=d_beta[n5]; beta[6]=d_beta[n6]; beta[7]=d_beta[n7];
        grad_phi[4][0]=d_grad_phi_x[n4]; grad_phi[4][1]=d_grad_phi_y[n4]; grad_phi[4][2]=d_grad_phi_z[n4];
        grad_phi[5][0]=d_grad_phi_x[n5]; grad_phi[5][1]=d_grad_phi_y[n5]; grad_phi[5][2]=d_grad_phi_z[n5];
        grad_phi[6][0]=d_grad_phi_x[n6]; grad_phi[6][1]=d_grad_phi_y[n6]; grad_phi[6][2]=d_grad_phi_z[n6];
        grad_phi[7][0]=d_grad_phi_x[n7]; grad_phi[7][1]=d_grad_phi_y[n7]; grad_phi[7][2]=d_grad_phi_z[n7];
        dofs[4]=d_node_to_dof[n4]; dofs[5]=d_node_to_dof[n5];
        dofs[6]=d_node_to_dof[n6]; dofs[7]=d_node_to_dof[n7];
        own[4]=d_ownership[n4]; own[5]=d_ownership[n5];
        own[6]=d_ownership[n6]; own[7]=d_ownership[n7];
    }

    // ------------------------------------------------------------------
    // Assemble local 8×8 matrix and RHS
    // ------------------------------------------------------------------
    RealType lhs[64];
    RealType rhs[8];
    #pragma unroll
    for (int i = 0; i < 64; ++i) lhs[i] = 0.0;
    #pragma unroll
    for (int i = 0; i <  8; ++i) rhs[i] = 0.0;

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        const int nodeL = hexLRSCV[ip * 2];
        const int nodeR = hexLRSCV[ip * 2 + 1];

        const RealType mdot = d_mdot[elemIdx * 12 + ip];
        const RealType aVx  = d_areaVec_x[elemIdx * 12 + ip];
        const RealType aVy  = d_areaVec_y[elemIdx * 12 + ip];
        const RealType aVz  = d_areaVec_z[elemIdx * 12 + ip];

        // Interpolate gamma and coords to SCS
        RealType gamma_ip = 0.0;
        RealType cx = 0.0, cy = 0.0, cz = 0.0;
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            const RealType sf = hexShapeFcn[ip][n];
            gamma_ip += sf * gamma[n];
            cx += sf * coords[n][0];
            cy += sf * coords[n][1];
            cz += sf * coords[n][2];
        }

        // Advection
        RealType phi_up, beta_up, dcorr;
        if (mdot > 0.0) {
            phi_up  = phi[nodeL];  beta_up = beta[nodeL];
            dcorr = grad_phi[nodeL][0]*(cx-coords[nodeL][0])
                  + grad_phi[nodeL][1]*(cy-coords[nodeL][1])
                  + grad_phi[nodeL][2]*(cz-coords[nodeL][2]);
        } else {
            phi_up  = phi[nodeR];  beta_up = beta[nodeR];
            dcorr = grad_phi[nodeR][0]*(cx-coords[nodeR][0])
                  + grad_phi[nodeR][1]*(cy-coords[nodeR][1])
                  + grad_phi[nodeR][2]*(cz-coords[nodeR][2]);
        }
        dcorr *= beta_up;
        const RealType adv = mdot * (phi_up + dcorr);
        rhs[nodeL] -= adv;
        rhs[nodeR] += adv;

        const RealType fL = 0.5*(mdot + fabs(mdot));
        const RealType fR = 0.5*(mdot - fabs(mdot));
        lhs[nodeL*8+nodeL] += fL;
        lhs[nodeR*8+nodeL] -= fL;
        lhs[nodeL*8+nodeR] += fR;
        lhs[nodeR*8+nodeR] -= fR;

        // Diffusion
        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            const RealType dc = -(gamma_ip * (dndx[n][0]*aVx +
                                              dndx[n][1]*aVy +
                                              dndx[n][2]*aVz));
            rhs[nodeL] -= dc * phi[n];
            rhs[nodeR] += dc * phi[n];
            lhs[nodeL*8+n] += dc;
            lhs[nodeR*8+n] -= dc;
        }
    }

    // ------------------------------------------------------------------
    // CSR scatter — NO ATOMICS (same-color elements never share a node)
    // ------------------------------------------------------------------
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        if (own[i] == 0 || dofs[i] < 0) continue;

        const int row_dof = dofs[i];
        d_rhs[row_dof] += rhs[i];  // direct write, no race

        const int row_start = matrix->rowPtr[row_dof];
        const int row_end   = matrix->rowPtr[row_dof + 1];
        const int diag_pos  = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            const int col_dof = dofs[j];
            if (col_dof < 0) continue;
            const RealType value = lhs[i*8+j];
            if (value == 0.0) continue;

            if (i == j) {
                matrix->values[diag_pos] += value;  // direct write
            } else {
                // Binary search for off-diagonal position
                int left = row_start, right = row_end - 1;
                while (left <= right) {
                    const int mid = (left + right) >> 1;
                    const int col = matrix->colInd[mid];
                    if (col == col_dof) {
                        matrix->values[mid] += value;  // direct write
                        break;
                    } else if (col < col_dof) { left  = mid + 1; }
                    else                       { right = mid - 1; }
                }
            }
        }
    }
}

} // namespace fem
} // namespace mars
