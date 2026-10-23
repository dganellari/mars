#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

template <int Value> struct TensorEarlyIndex {
    static constexpr int value = Value;
};

__host__ __device__ constexpr int tensor_early_face(int stage) {
    constexpr int order[12] = {0, 3, 8, 1, 9, 2, 10, 11, 4, 7, 5, 6};
    return order[stage];
}
__host__ __device__ constexpr int tensor_early_last(int row) {
    constexpr int last[8] = {2, 4, 6, 7, 9, 10, 11, 11};
    return last[row];
}

constexpr bool tensor_early_schedule_valid() {
    bool seen[12] = {};
    for (int stage = 0; stage < 12; ++stage) {
        const int ip = tensor_early_face(stage);
        if (ip < 0 || ip >= 12 || seen[ip])
            return false;
        seen[ip] = true;
    }
    for (int row = 0; row < 8; ++row) {
        int last = -1;
        for (int stage = 0; stage < 12; ++stage) {
            const int ip = tensor_early_face(stage);
            if (hexLeftNode(ip) == row || hexRightNode(ip) == row)
                last = stage;
        }
        if (last != tensor_early_last(row))
            return false;
    }
    return true;
}
static_assert(tensor_early_schedule_valid());

// The launch bounds give the compiler the real block size: without them it must allow 1024-thread
// blocks, caps the kernel at 128 VGPRs and spills. With (256, 2) the assembly on MI300A takes 14.1 instead
// of 25.4 ms (hex192, same matrix). Launch with at most 256 threads.
template <typename KeyType, typename RealType, int BlockSize = 256>
__global__ void __launch_bounds__(BlockSize, 2) cvfem_hex_assembly_kernel_tensor_early(
    const KeyType *__restrict__ d_conn0, const KeyType *__restrict__ d_conn1,
    const KeyType *__restrict__ d_conn2, const KeyType *__restrict__ d_conn3,
    const KeyType *__restrict__ d_conn4, const KeyType *__restrict__ d_conn5,
    const KeyType *__restrict__ d_conn6, const KeyType *__restrict__ d_conn7, size_t numElements,
    const RealType *__restrict__ d_x, const RealType *__restrict__ d_y, const RealType *__restrict__ d_z,
    const RealType *__restrict__ d_gamma, const RealType *__restrict__ d_phi,
    const RealType *__restrict__ d_beta, const RealType *__restrict__ d_grad_phi_x,
    const RealType *__restrict__ d_grad_phi_y, const RealType *__restrict__ d_grad_phi_z,
    const RealType *__restrict__ d_mdot, const RealType *__restrict__ d_areaVec_x,
    const RealType *__restrict__ d_areaVec_y, const RealType *__restrict__ d_areaVec_z,
    const int *__restrict__ d_node_to_dof, const uint8_t *__restrict__ d_ownership,
    CSRMatrix<RealType> *matrix, RealType *__restrict__ d_rhs) {
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements)
        return;

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
    RealType phi[8], gamma[8];

    {
        KeyType n0 = nodes[0], n1 = nodes[1], n2 = nodes[2], n3 = nodes[3];

        coords[0][0] = d_x[n0];
        coords[0][1] = d_y[n0];
        coords[0][2] = d_z[n0];
        coords[1][0] = d_x[n1];
        coords[1][1] = d_y[n1];
        coords[1][2] = d_z[n1];
        coords[2][0] = d_x[n2];
        coords[2][1] = d_y[n2];
        coords[2][2] = d_z[n2];
        coords[3][0] = d_x[n3];
        coords[3][1] = d_y[n3];
        coords[3][2] = d_z[n3];

        phi[0] = d_phi[n0];
        phi[1] = d_phi[n1];
        phi[2] = d_phi[n2];
        phi[3] = d_phi[n3];
        gamma[0] = d_gamma[n0];
        gamma[1] = d_gamma[n1];
        gamma[2] = d_gamma[n2];
        gamma[3] = d_gamma[n3];
    }

    {
        KeyType n4 = nodes[4], n5 = nodes[5], n6 = nodes[6], n7 = nodes[7];

        coords[4][0] = d_x[n4];
        coords[4][1] = d_y[n4];
        coords[4][2] = d_z[n4];
        coords[5][0] = d_x[n5];
        coords[5][1] = d_y[n5];
        coords[5][2] = d_z[n5];
        coords[6][0] = d_x[n6];
        coords[6][1] = d_y[n6];
        coords[6][2] = d_z[n6];
        coords[7][0] = d_x[n7];
        coords[7][1] = d_y[n7];
        coords[7][2] = d_z[n7];

        phi[4] = d_phi[n4];
        phi[5] = d_phi[n5];
        phi[6] = d_phi[n6];
        phi[7] = d_phi[n7];
        gamma[4] = d_gamma[n4];
        gamma[5] = d_gamma[n5];
        gamma[6] = d_gamma[n6];
        gamma[7] = d_gamma[n7];
    }

    RealType lhs[64];
    RealType rhs[8];

#pragma unroll
    for (int i = 0; i < 64; ++i)
        lhs[i] = 0.0;
#pragma unroll
    for (int i = 0; i < 8; ++i)
        rhs[i] = 0.0;

    // Finish rows after their last incident face; each face is evaluated once.
    auto integrate = [&] __device__(auto stage_tag) __attribute__((always_inline)) {
        constexpr int stage = decltype(stage_tag)::value;
        constexpr int ip = tensor_early_face(stage);
        const int nodeL = hexLeftNode(ip);
        const int nodeR = hexRightNode(ip);

        RealType mdot = d_mdot[elemIdx * 12 + ip];

        RealType areaVec[3];
        areaVec[0] = d_areaVec_x[elemIdx * 12 + ip];
        areaVec[1] = d_areaVec_y[elemIdx * 12 + ip];
        areaVec[2] = d_areaVec_z[elemIdx * 12 + ip];

        RealType gamma_ip = 0.0;
        RealType coords_ip[3] = {0.0, 0.0, 0.0};

#pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType sf = hexShapeFcn[ip][n];
            gamma_ip += sf * gamma[n];
            coords_ip[0] += sf * coords[n][0];
            coords_ip[1] += sf * coords[n][1];
            coords_ip[2] += sf * coords[n][2];
        }

        RealType phi_upwind, beta_upwind, dcorr = 0.0;
        if (mdot > 0.0) {
            phi_upwind = phi[nodeL];
            beta_upwind = d_beta[nodes[nodeL]];
            dcorr = d_grad_phi_x[nodes[nodeL]] * (coords_ip[0] - coords[nodeL][0]) +
                    d_grad_phi_y[nodes[nodeL]] * (coords_ip[1] - coords[nodeL][1]) +
                    d_grad_phi_z[nodes[nodeL]] * (coords_ip[2] - coords[nodeL][2]);
        } else {
            phi_upwind = phi[nodeR];
            beta_upwind = d_beta[nodes[nodeR]];
            dcorr = d_grad_phi_x[nodes[nodeR]] * (coords_ip[0] - coords[nodeR][0]) +
                    d_grad_phi_y[nodes[nodeR]] * (coords_ip[1] - coords[nodeR][1]) +
                    d_grad_phi_z[nodes[nodeR]] * (coords_ip[2] - coords[nodeR][2]);
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

        RealType J[3][3] = {};
        for (int node = 0; node < 8; ++node)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    J[i][j] += hexDerivConst[ip][node][j] * coords[node][i];
        RealType invJ[3][3];
        invert3x3(J, invJ);
        // grad_x(N) dot area = grad_reference(N) dot (J^-1 area).
        RealType reference_area[3] = {};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                reference_area[j] += invJ[j][i] * areaVec[i];

#pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType diff_coeff = -(gamma_ip * (hexDerivConst[ip][n][0] * reference_area[0] +
                                                hexDerivConst[ip][n][1] * reference_area[1] +
                                                hexDerivConst[ip][n][2] * reference_area[2]));

            rhs[nodeL] -= diff_coeff * phi[n];
            rhs[nodeR] += diff_coeff * phi[n];

            lhs[nodeL * 8 + n] += diff_coeff;
            lhs[nodeR * 8 + n] -= diff_coeff;
        }

        const int numOwnedRows = (matrix->numOwnedRows >= 0) ? matrix->numOwnedRows : matrix->numRows;
        auto scatter_row = [&] __device__(auto row_tag) __attribute__((always_inline)) {
            constexpr int i = decltype(row_tag)::value;
            if constexpr (tensor_early_last(i) == stage) {
                int row_dof = d_node_to_dof[nodes[i]];
                uint8_t ownership = d_ownership[nodes[i]];

                if (ownership == 0 || row_dof < 0 || row_dof >= numOwnedRows)
                    return;

                atomicAdd(&d_rhs[row_dof], rhs[i]);

                int row_start = matrix->rowPtr[row_dof];
                int row_end = matrix->rowPtr[row_dof + 1];
                int diag_pos = matrix->diagPtr[row_dof];

#pragma unroll
                for (int j = 0; j < 8; ++j) {
                    int col_dof = d_node_to_dof[nodes[j]];
                    if (col_dof < 0)
                        continue;

                    RealType value = lhs[i * 8 + j];
                    if (value == 0.0)
                        continue;

                    // Distinct periodic corners can contribute to the same diagonal.
                    if (col_dof == row_dof) {
                        atomicAdd(&matrix->values[diag_pos], value);
                    } else {
                        int left = row_start;
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
        };
        scatter_row(TensorEarlyIndex<0>{});
        scatter_row(TensorEarlyIndex<1>{});
        scatter_row(TensorEarlyIndex<2>{});
        scatter_row(TensorEarlyIndex<3>{});
        scatter_row(TensorEarlyIndex<4>{});
        scatter_row(TensorEarlyIndex<5>{});
        scatter_row(TensorEarlyIndex<6>{});
        scatter_row(TensorEarlyIndex<7>{});
    };
    integrate(TensorEarlyIndex<0>{});
    integrate(TensorEarlyIndex<1>{});
    integrate(TensorEarlyIndex<2>{});
    integrate(TensorEarlyIndex<3>{});
    integrate(TensorEarlyIndex<4>{});
    integrate(TensorEarlyIndex<5>{});
    integrate(TensorEarlyIndex<6>{});
    integrate(TensorEarlyIndex<7>{});
    integrate(TensorEarlyIndex<8>{});
    integrate(TensorEarlyIndex<9>{});
    integrate(TensorEarlyIndex<10>{});
    integrate(TensorEarlyIndex<11>{});
}
} // namespace fem
} // namespace mars
