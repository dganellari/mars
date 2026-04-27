#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_common.hpp"

namespace mars {
namespace fem {

// =============================================================================
// Smem node-cache CVFEM assembly kernel
// =============================================================================
// Problem with tensor_perip_lb2: node reads hit L2 at ~30 cycles each.
// Each block of 256 Morton-ordered elements references ~300-400 unique nodes
// out of 256×8=2048 total node references (neighbors share nodes heavily).
//
// Strategy: deduplicate + cooperative load unique node data into smem once,
// then serve all 2048 node reads from smem at ~5 cycles instead of L2.
//
// Uses DYNAMIC shared memory to exceed the 48 KB static smem limit.
// cudaFuncSetAttribute(cudaFuncAttributeMaxDynamicSharedMemorySize) must be
// called before the first launch to unlock the larger allocation.
//
// Dynamic smem layout (~53 KB per block, 2 blocks/SM → 106 KB < 228 KB GH200):
//   hash_keys [HTSIZE ints]   =  4 KB  (node IDs, -1=empty)
//   hash_slots[HTSIZE ints]   =  4 KB  (smem slot per unique node)
//   slot_counter[1 int]       =  4 B   (atomic slot allocator)
//   <8-byte padding>
//   node_data [MAX_NODES×9 doubles] = 36 KB  (x,y,z,phi,gamma,beta,gx,gy,gz)
//   dof       [MAX_NODES ints]      =  2 KB
//   own       [MAX_NODES uint8_t]   =  0.5 KB
//   <4-byte padding>
//   row_start [MAX_NODES ints]      =  2 KB  (rowPtr[dof]   cached)
//   row_end   [MAX_NODES ints]      =  2 KB  (rowPtr[dof+1] cached)
//   diag_pos  [MAX_NODES ints]      =  2 KB  (diagPtr[dof]  cached)
//   Total                           = ~53 KB
// =============================================================================

// Compile-time smem size for a given HTSIZE and MAX_NODES
template<int HTSIZE, int MAX_NODES>
constexpr size_t smemCacheBytes() {
    return HTSIZE * sizeof(int)          // hash_keys
         + HTSIZE * sizeof(int)          // hash_slots
         + sizeof(int)                   // slot_counter
         + 4                             // padding: align double to 8 bytes
         + MAX_NODES * 9 * sizeof(double)// node_data
         + MAX_NODES * sizeof(int)       // dof
         + MAX_NODES * sizeof(uint8_t)   // own
         + (MAX_NODES % 4 != 0 ? (4 - MAX_NODES % 4) : 0) // pad to int-align
         + 3 * MAX_NODES * sizeof(int);  // row_start, row_end, diag_pos
}

__device__ __forceinline__ int smem_hash(unsigned int key, int mask) {
    key ^= key >> 16;
    key *= 0x45d9f3bu;
    key ^= key >> 16;
    return (int)(key & (unsigned int)mask);
}

// Helper: compute pointer layout from a raw dynamic smem buffer
// All pointers reference into the same dynamically-allocated smem region.
template<int HTSIZE, int MAX_NODES>
struct SmemPtrs {
    int*     hash_keys;
    int*     hash_slots;
    int*     slot_counter;
    double*  node_data;    // [MAX_NODES * 9]
    int*     dof;
    uint8_t* own;
    int*     row_start;
    int*     row_end;
    int*     diag_pos;

    __device__ explicit SmemPtrs(char* base) {
        hash_keys    = reinterpret_cast<int*>(base);
        hash_slots   = hash_keys    + HTSIZE;
        slot_counter = hash_slots   + HTSIZE;
        // align to 8 bytes for double
        char* p8     = reinterpret_cast<char*>(slot_counter + 1);
        p8           = reinterpret_cast<char*>((reinterpret_cast<uintptr_t>(p8) + 7) & ~7u);
        node_data    = reinterpret_cast<double*>(p8);
        dof          = reinterpret_cast<int*>(node_data + MAX_NODES * 9);
        own          = reinterpret_cast<uint8_t*>(dof + MAX_NODES);
        // align to 4 bytes for int
        char* p4     = reinterpret_cast<char*>(own + MAX_NODES);
        p4           = reinterpret_cast<char*>((reinterpret_cast<uintptr_t>(p4) + 3) & ~3u);
        row_start    = reinterpret_cast<int*>(p4);
        row_end      = row_start + MAX_NODES;
        diag_pos     = row_end   + MAX_NODES;
    }
};

template<typename KeyType, typename RealType,
         int BlockSize = 256,
         int HTSIZE    = 1024,   // power of 2, > 2× expected unique nodes per block
         int MAX_NODES = 512>
__launch_bounds__(256, 2)
__global__ void cvfem_hex_assembly_kernel_smem_cache(
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
    // Dynamic smem — size set at launch via cudaFuncSetAttribute
    extern __shared__ char s_smem[];
    SmemPtrs<HTSIZE, MAX_NODES> s(s_smem);

    const int tid     = threadIdx.x;
    const int elemIdx = blockIdx.x * BlockSize + tid;
    const bool valid  = (elemIdx < (int)numElements);
    const int HT_MASK = HTSIZE - 1;

    // =========================================================================
    // Phase 1: Initialize hash table (all threads)
    // =========================================================================
    for (int i = tid; i < HTSIZE; i += BlockSize)
        s.hash_keys[i] = -1;
    if (tid == 0) *s.slot_counter = 0;
    __syncthreads();

    // =========================================================================
    // Phase 2: Load connectivity + CAS-insert node IDs into hash table.
    // Hopper hardware smem atomics: ~1 cycle throughput.
    // If two threads insert the same nid: the second sees prev==nid → breaks.
    // =========================================================================
    KeyType my_nodes[8];
    if (valid) {
        my_nodes[0] = d_conn0[elemIdx]; my_nodes[1] = d_conn1[elemIdx];
        my_nodes[2] = d_conn2[elemIdx]; my_nodes[3] = d_conn3[elemIdx];
        my_nodes[4] = d_conn4[elemIdx]; my_nodes[5] = d_conn5[elemIdx];
        my_nodes[6] = d_conn6[elemIdx]; my_nodes[7] = d_conn7[elemIdx];

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            int nid = (int)my_nodes[n];
            int h   = smem_hash((unsigned)nid, HT_MASK);
            while (true) {
                int prev = atomicCAS(&s.hash_keys[h], -1, nid);
                if (prev == -1 || prev == nid) break;
                h = (h + 1) & HT_MASK;
            }
        }
    }
    __syncthreads();

    // =========================================================================
    // Phase 3: Assign sequential smem slots to occupied hash entries
    // =========================================================================
    for (int i = tid; i < HTSIZE; i += BlockSize)
        if (s.hash_keys[i] != -1)
            s.hash_slots[i] = atomicAdd(s.slot_counter, 1);
    __syncthreads();

    // =========================================================================
    // Phase 4: Cooperatively load unique node data + CSR metadata → smem.
    // ~400 unique nodes / 256 threads ≈ 1-2 nodes per thread.
    // Each unique node's data is loaded exactly once per block (not per element).
    // =========================================================================
    for (int i = tid; i < HTSIZE; i += BlockSize) {
        int nid = s.hash_keys[i];
        if (nid < 0) continue;
        int sl = s.hash_slots[i];
        if (sl >= MAX_NODES) continue;  // safety: very poor spatial locality

        double* base = s.node_data + sl * 9;
        base[0] = d_x[nid];           base[1] = d_y[nid];           base[2] = d_z[nid];
        base[3] = d_phi[nid];         base[4] = d_gamma[nid];       base[5] = d_beta[nid];
        base[6] = d_grad_phi_x[nid];  base[7] = d_grad_phi_y[nid];  base[8] = d_grad_phi_z[nid];

        int dof         = d_node_to_dof[nid];
        s.dof[sl]       = dof;
        s.own[sl]       = d_ownership[nid];

        // Cache rowPtr/diagPtr so pre-lookup phase hits smem not L2
        if (dof >= 0) {
            s.row_start[sl] = matrix->rowPtr[dof];
            s.row_end[sl]   = matrix->rowPtr[dof + 1];
            s.diag_pos[sl]  = matrix->diagPtr[dof];
        } else {
            s.row_start[sl] = 0;
            s.row_end[sl]   = 0;
            s.diag_pos[sl]  = -1;
        }
    }
    __syncthreads();

    if (!valid) return;

    // =========================================================================
    // Phase 5: Look up smem slot for each of my 8 nodes
    // =========================================================================
    int slot_map[8];
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        int nid = (int)my_nodes[n];
        int h   = smem_hash((unsigned)nid, HT_MASK);
        while (s.hash_keys[h] != nid) h = (h + 1) & HT_MASK;
        slot_map[n] = s.hash_slots[h];
    }

    // =========================================================================
    // Phase 6a: Load node data from smem → registers (~5 cycle latency)
    // =========================================================================
    RealType coords[8][3], phi[8], gamma[8], beta[8], grad_phi[8][3];
    int dofs[8]; uint8_t own[8];

    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        int sl             = slot_map[n];
        const double* base = s.node_data + sl * 9;
        coords[n][0]   = base[0]; coords[n][1] = base[1]; coords[n][2] = base[2];
        phi[n]         = base[3];
        gamma[n]       = base[4];
        beta[n]        = base[5];
        grad_phi[n][0] = base[6]; grad_phi[n][1] = base[7]; grad_phi[n][2] = base[8];
        dofs[n]        = s.dof[sl];
        own[n]         = s.own[sl];
    }

    // =========================================================================
    // Phase 6b: Pre-lookup all 64 CSR positions.
    // rowPtr/diagPtr reads now hit smem; colInd binary search still hits L2.
    // =========================================================================
    int pos[64];
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        if (own[i] == 0 || dofs[i] < 0) {
            #pragma unroll
            for (int j = 0; j < 8; ++j) pos[i * 8 + j] = -1;
            continue;
        }
        int sl_i      = slot_map[i];
        int row_start = s.row_start[sl_i];  // smem
        int row_end   = s.row_end[sl_i];    // smem
        int diag_pos  = s.diag_pos[sl_i];   // smem

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            int col_dof = dofs[j];
            if (col_dof < 0) { pos[i * 8 + j] = -1;        continue; }
            if (i == j)       { pos[i * 8 + j] = diag_pos;  continue; }

            int left = row_start, right = row_end - 1, found = -1;
            while (left <= right) {
                int mid = (left + right) >> 1;
                int col = matrix->colInd[mid];
                if      (col == col_dof) { found = mid; break; }
                else if (col <  col_dof)   left  = mid + 1;
                else                       right = mid - 1;
            }
            pos[i * 8 + j] = found;
        }
    }

    // =========================================================================
    // Phase 6c: SCS integration — scatter LHS per SCS, accumulate RHS
    // =========================================================================
    RealType rhs[8] = {};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip * 2];
        int nodeR = hexLRSCV[ip * 2 + 1];

        RealType mdot     = d_mdot[elemIdx * 12 + ip];
        RealType areaVec0 = d_areaVec_x[elemIdx * 12 + ip];
        RealType areaVec1 = d_areaVec_y[elemIdx * 12 + ip];
        RealType areaVec2 = d_areaVec_z[elemIdx * 12 + ip];

        RealType gamma_ip = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType sf = hexShapeFcn[ip][n];
            gamma_ip += sf * gamma[n];
            cx += sf * coords[n][0];
            cy += sf * coords[n][1];
            cz += sf * coords[n][2];
        }

        RealType phi_up, beta_up, dcorr;
        if (mdot > 0.0) {
            phi_up  = phi[nodeL]; beta_up = beta[nodeL];
            dcorr   = grad_phi[nodeL][0] * (cx - coords[nodeL][0]) +
                      grad_phi[nodeL][1] * (cy - coords[nodeL][1]) +
                      grad_phi[nodeL][2] * (cz - coords[nodeL][2]);
        } else {
            phi_up  = phi[nodeR]; beta_up = beta[nodeR];
            dcorr   = grad_phi[nodeR][0] * (cx - coords[nodeR][0]) +
                      grad_phi[nodeR][1] * (cy - coords[nodeR][1]) +
                      grad_phi[nodeR][2] * (cz - coords[nodeR][2]);
        }
        dcorr *= beta_up;

        RealType adv_flux = mdot * (phi_up + dcorr);
        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        RealType lhsfac_L = 0.5 * (mdot + fabs(mdot));
        RealType lhsfac_R = 0.5 * (mdot - fabs(mdot));
        if (lhsfac_L != 0.0) {
            int p = pos[nodeL * 8 + nodeL]; if (p >= 0) atomicAdd(&matrix->values[p],  lhsfac_L);
                p = pos[nodeR * 8 + nodeL]; if (p >= 0) atomicAdd(&matrix->values[p], -lhsfac_L);
        }
        if (lhsfac_R != 0.0) {
            int p = pos[nodeL * 8 + nodeR]; if (p >= 0) atomicAdd(&matrix->values[p],  lhsfac_R);
                p = pos[nodeR * 8 + nodeR]; if (p >= 0) atomicAdd(&matrix->values[p], -lhsfac_R);
        }

        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);

        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType diff_coeff = -(gamma_ip * (dndx[n][0] * areaVec0 +
                                                dndx[n][1] * areaVec1 +
                                                dndx[n][2] * areaVec2));
            rhs[nodeL] -= diff_coeff * phi[n];
            rhs[nodeR] += diff_coeff * phi[n];
            int pL = pos[nodeL * 8 + n]; if (pL >= 0) atomicAdd(&matrix->values[pL],  diff_coeff);
            int pR = pos[nodeR * 8 + n]; if (pR >= 0) atomicAdd(&matrix->values[pR], -diff_coeff);
        }
    }

    #pragma unroll
    for (int i = 0; i < 8; ++i)
        if (own[i] != 0 && dofs[i] >= 0)
            atomicAdd(&d_rhs[dofs[i]], rhs[i]);
}

} // namespace fem
} // namespace mars
