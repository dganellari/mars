#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

// The tensor kernel's element matrix computed by 16 threads per element, in three phases through shared
// memory: lanes 0-15 load the element's 8 nodes (two lanes per node), lanes 0-11 integrate one
// sub-control surface each, and lanes 0-15 add the faces into the 8 rows (two lanes per row, 4 columns
// each) and scatter them. Each face is evaluated once, and no thread holds a 64-entry element matrix.
//
// Face f from node L to node R gives the 8 coefficients c_fj and the flux q_f; row L gains c_f and
// -q_f, row R loses them. Diffusion and advection are summed per face here, before the rows, so the
// entries differ from the tensor kernel's in the last bits (floating-point association).

// The face's shifted integration point (s1, s2, s3) in the reference cube [-1/2, 1/2]^3, and its nodes.
// Lane-dependent, so computed with selects: a constexpr table indexed at run time goes to memory.
__host__ __device__ inline void tensorFacesPoint(int ip, double& s1, double& s2, double& s3)
{
    const int k = ip & 3;
    if (ip < 8)
    {
        s1 = (k == 1) ? 0.5 : (k == 3) ? -0.5 : 0.0;
        s2 = (k == 0) ? -0.5 : (k == 2) ? 0.5 : 0.0;
        s3 = (ip < 4) ? -0.5 : 0.5;
    }
    else
    {
        s1 = (k == 1 || k == 2) ? 0.5 : -0.5;
        s2 = (k >= 2) ? 0.5 : -0.5;
        s3 = 0.0;
    }
}

__host__ __device__ constexpr int tensorFacesLeft(int ip)
{
    return ip < 8 ? (ip & 4) + ((ip & 3) == 3 ? 0 : (ip & 3)) : ip - 8;
}
__host__ __device__ constexpr int tensorFacesRight(int ip)
{
    return ip < 8 ? (ip & 4) + ((ip & 3) == 3 ? 3 : (ip & 3) + 1) : ip - 4;
}

// node i's k-th face (k = 0, 1, 2): the two faces in its plane, then the vertical one
__host__ __device__ constexpr int tensorFacesOfNode(int i, int k)
{
    const int base = i & 4, m = i & 3;
    return k == 0 ? base + ((m + 3) & 3) : k == 1 ? base + m : 8 + m;
}

constexpr bool tensorFacesTablesValid()
{
    for (int ip = 0; ip < 12; ++ip)
        if (tensorFacesLeft(ip) != hexLeftNode(ip) || tensorFacesRight(ip) != hexRightNode(ip)) return false;
    for (int i = 0; i < 8; ++i)
        for (int k = 0; k < 3; ++k)
        {
            const int ip = tensorFacesOfNode(i, k);
            if (ip < 0 || ip >= 12 || (hexLeftNode(ip) != i && hexRightNode(ip) != i)) return false;
            for (int k2 = 0; k2 < k; ++k2)
                if (tensorFacesOfNode(i, k2) == ip) return false;
        }
    return true;
}
static_assert(tensorFacesTablesValid(), "face/node maps must match hexLeftNode/hexRightNode");

template<typename RealType>
struct TensorFacesElement
{
    // the node data is dead once the faces are integrated, so the face results reuse its storage
    union
    {
        struct
        {
            RealType x[3][8], phi[8], gamma[8], beta[8], grad[3][8];
        } node;
        struct
        {
            RealType c[8][12], q[12]; // c[j][f]: lanes of consecutive faces write consecutive words
        } face;
    };
    int dof[8];
    uint8_t own[8];
};

template<typename KeyType, typename RealType, int BlockSize = 256, bool LumpMissing = false>
__global__ void __launch_bounds__(BlockSize) cvfem_hex_assembly_kernel_tensor_faces(
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
    static_assert(BlockSize % 16 == 0, "16 threads per element");
    constexpr int ElemsPerBlock = BlockSize / 16;
    __shared__ TensorFacesElement<RealType> shared[ElemsPerBlock];

    const int lane        = threadIdx.x % 16;
    const size_t elemIdx  = size_t(blockIdx.x) * ElemsPerBlock + threadIdx.x / 16;
    const bool active     = elemIdx < numElements;
    TensorFacesElement<RealType>& s = shared[threadIdx.x / 16];

    // Phase 1: lane i loads coordinates and phi, gamma of node i; lane 8 + i its beta, gradient, DOF.
    // Inactive tail elements skip the work but take part in every barrier.
    if (active)
    {
        const int i      = lane & 7;
        const KeyType* c = d_conn0;
        c                = (i == 1) ? d_conn1 : c;
        c                = (i == 2) ? d_conn2 : c;
        c                = (i == 3) ? d_conn3 : c;
        c                = (i == 4) ? d_conn4 : c;
        c                = (i == 5) ? d_conn5 : c;
        c                = (i == 6) ? d_conn6 : c;
        c                = (i == 7) ? d_conn7 : c;
        const KeyType n  = c[elemIdx];
        if (lane < 8)
        {
            s.node.x[0][i]  = d_x[n];
            s.node.x[1][i]  = d_y[n];
            s.node.x[2][i]  = d_z[n];
            s.node.phi[i]   = d_phi[n];
            s.node.gamma[i] = d_gamma[n];
        }
        else
        {
            s.node.beta[i]    = d_beta[n];
            s.node.grad[0][i] = d_grad_phi_x[n];
            s.node.grad[1][i] = d_grad_phi_y[n];
            s.node.grad[2][i] = d_grad_phi_z[n];
            s.dof[i]          = d_node_to_dof[n];
            s.own[i]          = d_ownership[n];
        }
    }
    __syncthreads();

    // Phase 2: lane ip integrates face ip, in registers; the results go to shared memory after a barrier,
    // as they overwrite the node data that the other faces still read.
    const int ip       = lane;
    const bool produce = active && ip < 12;
    RealType coef[8], q = 0;
    if (produce)
    {
        const int nodeL = tensorFacesLeft(ip), nodeR = tensorFacesRight(ip);
        const size_t f  = elemIdx * 12 + ip;
        const RealType mdot = d_mdot[f];
        const RealType area[3] = {d_areaVec_x[f], d_areaVec_y[f], d_areaVec_z[f]};

        // the integration point lies on the edge L-R: only N_L = N_R = 1/2 are nonzero
        const RealType gamma_ip = 0.5 * s.node.gamma[nodeL] + 0.5 * s.node.gamma[nodeR];
        RealType x_ip[3];
#pragma unroll
        for (int a = 0; a < 3; ++a)
            x_ip[a] = 0.5 * s.node.x[a][nodeL] + 0.5 * s.node.x[a][nodeR];

        const int up = (mdot > 0.0) ? nodeL : nodeR;
        RealType dcorr = s.node.grad[0][up] * (x_ip[0] - s.node.x[0][up]) +
                         s.node.grad[1][up] * (x_ip[1] - s.node.x[1][up]) +
                         s.node.grad[2][up] * (x_ip[2] - s.node.x[2][up]);
        dcorr *= s.node.beta[up];
        q = mdot * (s.node.phi[up] + dcorr);

        // dN_n/ds of the trilinear hex at (s1, s2, s3), node n at the corner signs (a, b, c): products of
        // the factors 1/2 +- s, which are 0, 1/2 or 1, so exactly hexDerivConst
        RealType s1, s2, s3;
        tensorFacesPoint(ip, s1, s2, s3);
        const RealType f1[2] = {0.5 - s1, 0.5 + s1}, f2[2] = {0.5 - s2, 0.5 + s2}, f3[2] = {0.5 - s3, 0.5 + s3};
        RealType dN[8][3];
#pragma unroll
        for (int n = 0; n < 8; ++n)
        {
            // hex node order: (-,-,-) (+,-,-) (+,+,-) (-,+,-), then the same at s3 = +
            const int a = (n & 3) == 1 || (n & 3) == 2, b = (n & 3) >= 2, c = n >= 4;
            const RealType sa = a ? 1.0 : -1.0, sb = b ? 1.0 : -1.0, sc = c ? 1.0 : -1.0;
            dN[n][0] = sa * (f2[b] * f3[c]);
            dN[n][1] = sb * (f1[a] * f3[c]);
            dN[n][2] = sc * (f1[a] * f2[b]);
        }

        RealType J[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
#pragma unroll
        for (int n = 0; n < 8; ++n)
#pragma unroll
            for (int a = 0; a < 3; ++a)
#pragma unroll
                for (int b = 0; b < 3; ++b)
                    J[a][b] += dN[n][b] * s.node.x[a][n];
        RealType invJ[3][3];
        invert3x3(J, invJ);
        // grad_x(N) dot area = grad_reference(N) dot (J^-1 area)
        RealType ref_area[3] = {0, 0, 0};
#pragma unroll
        for (int b = 0; b < 3; ++b)
#pragma unroll
            for (int a = 0; a < 3; ++a)
                ref_area[b] += invJ[b][a] * area[a];

        const RealType lhsfac_L = 0.5 * (mdot + fabs(mdot));
        const RealType lhsfac_R = 0.5 * (mdot - fabs(mdot));
#pragma unroll
        for (int n = 0; n < 8; ++n)
        {
            const RealType diff_coeff =
                -(gamma_ip * (dN[n][0] * ref_area[0] + dN[n][1] * ref_area[1] + dN[n][2] * ref_area[2]));
            q += diff_coeff * s.node.phi[n];
            coef[n] = diff_coeff + ((n == nodeL) ? lhsfac_L : 0.0) + ((n == nodeR) ? lhsfac_R : 0.0);
        }
    }
    __syncthreads();
    if (produce)
    {
#pragma unroll
        for (int n = 0; n < 8; ++n)
            s.face.c[n][ip] = coef[n];
        s.face.q[ip] = q;
    }
    __syncthreads();
    if (!active) return;

    // Phase 3: lanes i and 8 + i gather row i, columns 0-3 and 4-7, from the 3 faces of node i
    const int i    = lane & 7;
    const int col0 = (lane >> 3) * 4;
    RealType lhs[4] = {0, 0, 0, 0}, rhs = 0;
#pragma unroll
    for (int k = 0; k < 3; ++k)
    {
        const int face    = tensorFacesOfNode(i, k);
        const RealType sg = (tensorFacesLeft(face) == i) ? RealType(1) : RealType(-1); // L gains, R loses
#pragma unroll
        for (int j = 0; j < 4; ++j)
            lhs[j] += sg * s.face.c[col0 + j][face];
        rhs -= sg * s.face.q[face];
    }

    const int numOwnedRows = (matrix->numOwnedRows >= 0) ? matrix->numOwnedRows : matrix->numRows;
    const int row_dof      = s.dof[i];
    if (s.own[i] == 0 || row_dof < 0 || row_dof >= numOwnedRows) return;
    if (col0 == 0) atomicAdd(&d_rhs[row_dof], rhs);
    const int row_start = matrix->rowPtr[row_dof];
    const int row_end   = matrix->rowPtr[row_dof + 1];
    const int diag_pos  = matrix->diagPtr[row_dof];
    RealType diagonal = 0.0;
#pragma unroll
    for (int j = 0; j < 4; ++j)
    {
        const int col_dof = s.dof[col0 + j];
        if (col_dof < 0) continue;
        const RealType value = lhs[j];
        if (value == 0.0) continue;
        // DOF equality: under periodic collapse two corners can share the diagonal
        if (col_dof == row_dof)
        {
            if constexpr (LumpMissing) diagonal += value;
            else atomicAdd(&matrix->values[diag_pos], value);
            continue;
        }
        int left = row_start, right = row_end - 1;
        while (left <= right)
        {
            const int mid = (left + right) >> 1;
            const int col = matrix->colInd[mid];
            if (col == col_dof)
            {
                atomicAdd(&matrix->values[mid], value);
                break;
            }
            if (col < col_dof) left = mid + 1;
            else right = mid - 1;
        }
        if constexpr (LumpMissing) {
            if (left > right) diagonal += value;
        }
    }
    if constexpr (LumpMissing) {
        // Preserve the full row sum when the graph omits a coupling.
        if (diagonal != 0.0) atomicAdd(&matrix->values[diag_pos], diagonal);
    }
}

} // namespace fem
} // namespace mars
