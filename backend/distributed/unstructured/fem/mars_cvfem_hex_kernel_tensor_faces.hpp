#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "../utils/mars_gpu_warp.cuh"

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

// Loads and atomics through a pointer known to be in device global memory. The CSR arrays reached through
// the matrix struct are generic pointers, so on AMD the compiler emits flat instructions, which also count
// against the shared-memory wait counter; as global instructions the aggregated kernel is 1% faster on MI300A.
template<class T>
__device__ __forceinline__ T tensorFacesLoadGlobal(const T* p)
{
#if defined(__HIP_DEVICE_COMPILE__)
    return *(const __attribute__((address_space(1))) T*)p;
#else
    return *p;
#endif
}
template<class T>
__device__ __forceinline__ void tensorFacesAtomicAddGlobal(T* p, T v)
{
#if defined(__HIP_DEVICE_COMPILE__)
    __hip_atomic_fetch_add((__attribute__((address_space(1))) T*)p, v, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_AGENT);
#else
    atomicAdd(p, v);
#endif
}

// the largest power of 2 not above n, and its log2: the hash table of the aggregated scatter
constexpr int tensorFacesTableSlots(size_t n)
{
    int slots = 1;
    while (size_t(slots) * 2 <= n) slots *= 2;
    return slots;
}
constexpr int tensorFacesLog2(int n)
{
    int bits = 0;
    while ((1 << bits) < n) ++bits;
    return bits;
}

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

// Phases 1-3 of one element on its 16 lanes: each lane gets columns (lane / 8) * 4 to + 3 of row lane % 8,
// and the row's right-hand side. The kernels that look up the CSR rows themselves need LoadDofs.
template<typename KeyType, typename RealType, bool LoadDofs>
__device__ __forceinline__ void tensorFacesRow(TensorFacesElement<RealType>& s,
    int lane,
    size_t elemIdx,
    bool active,
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
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    RealType (&lhs)[4],
    RealType& rhs)
{
    // Phase 1: lane i loads coordinates and phi, gamma of node i; lane 8 + i its beta, gradient and, with
    // LoadDofs, DOF and ownership. Inactive tail elements skip the work but take part in every barrier.
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
            if constexpr (LoadDofs)
            {
                s.dof[i] = d_node_to_dof[n];
                s.own[i] = d_ownership[n];
            }
        }
    }
    // An element's 16 lanes share memory only among themselves and always sit in one warp or wavefront,
    // so warp barriers suffice: 12.2 against 12.4 ms on MI300A, 7.3 against 7.8 ms on H100 (hex192).
    syncWarp();

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
    syncWarp();
    if (produce)
    {
#pragma unroll
        for (int n = 0; n < 8; ++n)
            s.face.c[n][ip] = coef[n];
        s.face.q[ip] = q;
    }
    syncWarp();

    // Phase 3: lanes i and 8 + i gather row i, columns 0-3 and 4-7, from the 3 faces of node i
    const int i    = lane & 7;
    const int col0 = (lane >> 3) * 4;
#pragma unroll
    for (int j = 0; j < 4; ++j)
        lhs[j] = 0;
    rhs = 0;
    if (active)
    {
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
    }
}

// Aggregate (graph mode only): the block adds up the contributions of its elements to shared rows in a
// shared-memory table first, then sends one global atomic per destination. On MI300A every global atomic
// is forwarded past the L2 to the fabric; with 16 elements per block this cuts them 234 to 156 million
// and the assembly 12.3 to 10.1 ms (hex192). On H100, where atomics stay in the L2, it doubles the time.
// This path also loads a short CSR row whole, with independent loads, and matches the columns in
// registers instead of a binary search of dependent loads per column: 10.1 to 8.6 ms on MI300A.
template<typename KeyType, typename RealType, int BlockSize = 256, bool LumpMissing = false, bool Aggregate = false>
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

    RealType lhs[4], rhs;
    tensorFacesRow<KeyType, RealType, true>(s, lane, elemIdx, active, d_conn0, d_conn1, d_conn2, d_conn3, d_conn4,
                                            d_conn5, d_conn6, d_conn7, d_x, d_y, d_z, d_gamma, d_phi, d_beta,
                                            d_grad_phi_x, d_grad_phi_y, d_grad_phi_z, d_mdot, d_areaVec_x,
                                            d_areaVec_y, d_areaVec_z, d_node_to_dof, d_ownership, lhs, rhs);
    const int i    = lane & 7;
    const int col0 = (lane >> 3) * 4;

    const int numOwnedRows = (matrix->numOwnedRows >= 0) ? matrix->numOwnedRows : matrix->numRows;
    const int row_dof      = active ? s.dof[i] : -1;
    const bool ownRow      = active && s.own[i] != 0 && row_dof >= 0 && row_dof < numOwnedRows;

    if constexpr (!Aggregate)
    {
        if (!ownRow) return;
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
    else
    {
        static_assert(LumpMissing, "aggregation bounds the destinations per row by the graph pattern");
        auto ld  = [](const int* p) { return tensorFacesLoadGlobal(p); };
        auto add = [&](int k, RealType v) { tensorFacesAtomicAddGlobal(k >= 0 ? &matrix->values[k] : &d_rhs[-1 - k], v); };
        // A destination is a matrix value index, or -1 - row for the right-hand side. This lane's go to
        // fixed slots (compile-time indices keep them in registers): right-hand side, diagonal, 4 columns.
        constexpr int EmptyKey = -2147483647 - 1;
        int key[6];
        RealType val[6];
#pragma unroll
        for (int t = 0; t < 6; ++t)
        {
            key[t] = EmptyKey;
            val[t] = 0;
        }
        if (ownRow)
        {
            if (col0 == 0)
            {
                key[0] = -1 - row_dof;
                val[0] = rhs;
            }
            const int row_start = ld(&matrix->rowPtr[row_dof]);
            const int row_end   = ld(&matrix->rowPtr[row_dof + 1]);
            RealType diagonal   = 0.0;
            constexpr int RowMax = 8;
            const bool wholeRow  = row_end - row_start <= RowMax;
            int rowCols[RowMax];
#pragma unroll
            for (int t = 0; t < RowMax; ++t)
                rowCols[t] = (wholeRow && row_start + t < row_end) ? ld(&matrix->colInd[row_start + t]) : -1;
#pragma unroll
            for (int j = 0; j < 4; ++j)
            {
                const int col_dof = s.dof[col0 + j];
                const RealType value = lhs[j];
                if (col_dof < 0 || value == 0.0) continue;
                if (col_dof == row_dof)
                {
                    diagonal += value;
                    continue;
                }
                if (wholeRow)
                {
                    int pos = -1;
#pragma unroll
                    for (int t = 0; t < RowMax; ++t)
                        if (rowCols[t] == col_dof) pos = row_start + t;
                    if (pos >= 0)
                    {
                        key[2 + j] = pos;
                        val[2 + j] = value;
                    }
                    else diagonal += value;
                    continue;
                }
                int left = row_start, right = row_end - 1;
                while (left <= right)
                {
                    const int mid = (left + right) >> 1;
                    const int col = ld(&matrix->colInd[mid]);
                    if (col == col_dof)
                    {
                        key[2 + j] = mid;
                        val[2 + j] = value;
                        break;
                    }
                    if (col < col_dof) left = mid + 1;
                    else right = mid - 1;
                }
                if (left > right) diagonal += value;
            }
            if (diagonal != 0.0)
            {
                key[1] = ld(&matrix->diagPtr[row_dof]);
                val[1] = diagonal;
            }
        }

        // Neighbouring elements share rows, so the block adds up its contributions in a shared hash table
        // and sends one global atomic per distinct destination. The table reuses the element storage, dead
        // now. A row has at most 5 destinations in the graph pattern (rhs, diagonal, 3 edge neighbours in
        // the element), so the table is at most 5/8 full; should periodic aliases ever fill it, the
        // contribution goes straight to global memory.
        constexpr int Slots    = tensorFacesTableSlots(sizeof(shared) / (sizeof(RealType) + sizeof(int)));
        constexpr int SlotBits = tensorFacesLog2(Slots);
        RealType* tableVal = reinterpret_cast<RealType*>(shared);
        int* tableKey      = reinterpret_cast<int*>(tableVal + Slots);
        __syncthreads();
        for (int t = threadIdx.x; t < Slots; t += BlockSize)
        {
            tableKey[t] = EmptyKey;
            tableVal[t] = 0;
        }
        __syncthreads();
#pragma unroll
        for (int t = 0; t < 6; ++t)
        {
            if (key[t] == EmptyKey) continue;
            unsigned h  = (unsigned(key[t]) * 2654435761u) >> (32 - SlotBits);
            bool stored = false;
            for (int probe = 0; probe < Slots && !stored; ++probe)
            {
                const int prev = atomicCAS(&tableKey[h], EmptyKey, key[t]);
                if (prev == EmptyKey || prev == key[t])
                {
                    atomicAdd(&tableVal[h], val[t]);
                    stored = true;
                }
                h = (h + 1) & (Slots - 1);
            }
            if (!stored) add(key[t], val[t]);
        }
        __syncthreads();
        for (int t = threadIdx.x; t < Slots; t += BlockSize)
        {
            const int k = tableKey[t];
            if (k == EmptyKey) continue;
            const RealType v = tableVal[t];
            if (v != 0.0) add(k, v);
        }
    }
}

} // namespace fem
} // namespace mars
