#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_common.hpp"
#include <mma.h>
#include <type_traits>

namespace mars {
namespace fem {

// =============================================================================
// Shared-memory struct for one warp (one element) in the WMMA kernel.
// Defined at namespace scope so sizeof() is available on the host for launch.
//
// FP64 WMMA uses 8x8x4 tiles (SM80+, standard WMMA API via <mma.h>).
// SM90a (GH200) additionally supports WGMMA m64n8k4 via PTX / cuda::ptx —
// upgrade to that path when CUTLASS cute is available.
//
// Matrix decomposition: lhs_diff = S^T x D
//   S[12x8]  scatter:    S[ip][n] = +1 if n==nodeL[ip], -1 if n==nodeR[ip]
//   D[12x8]  diff-coeff: D[ip][n] = -gamma_ip * (dndx[ip][n] . areaVec[ip])
//   S^T[8x12] x D[12x8] = lhs_diff[8x8]
//
// S^T is element-independent (same hexLRSCV for every hex element).
//
// Tiled: K=12 = 3 tiles of K=4 → 3 x (A[8x4] * B[4x8]) accumulated into C[8x8]
//   Tile 0: K = 0..3  (SCS ip = 0..3)
//   Tile 1: K = 4..7  (SCS ip = 4..7)
//   Tile 2: K = 8..11 (SCS ip = 8..11)
// =============================================================================
struct alignas(128) CvfemWmmaWarpData {
    // WMMA 8x8x4 tiles — 3 tiles, each stored col-major
    //   smem_A tile t: A[8x4] col-major (ldm=8) → 32 doubles
    //   smem_B tile t: B[4x8] col-major (ldm=4) → 32 doubles
    double smem_A[3 * 32];    //  768 bytes — S^T tiles (element-independent)
    double smem_B[3 * 32];    //  768 bytes — D   tiles (element-specific)
    double smem_C[8 * 8];     //  512 bytes — C[8x8] result, col-major (ldm=8)

    // Node data (loaded by lanes 0-7, read by all)
    double coords[8][3];      //  192 bytes
    double phi[8];            //   64 bytes
    double gamma_n[8];        //   64 bytes
    double beta[8];           //   64 bytes
    double grad_phi[8][3];    //  192 bytes

    // SCS data (loaded by lanes 0-11)
    double mdot[12];          //   96 bytes
    double areaVec[12][3];    //  288 bytes

    // DOF metadata (loaded by lanes 0-7)
    int     dofs[8];          //   32 bytes
    uint8_t own[8];           //    8 bytes
    uint8_t _pad[7];          //    7 bytes — keep rhs[] 8-byte aligned

    // Accumulated local element vectors / matrix
    // rhs[i] and lhs[i][j] are zeroed in Phase 2 then updated in Phase 3
    // via shared-memory atomicAdd (FP64 smem atomics: hardware on SM80+).
    double rhs[8];            //   64 bytes
    double lhs[8][8];         //  512 bytes — advection (from Phase 3) + diffusion (from WMMA)
};
// Total per warp: ~3.6 KB → 8 warps/block ≈ 29 KB (well within 228 KB on GH200)

// =============================================================================
// FP64 WMMA CVFEM Assembly  (requires SM80+ / CUDA 11.0+)
// Optimised for Hopper GH200 (SM90a): smem atomics are hardware-native,
// smem = 228 KB, FP64 tensor throughput ~67 TFLOPS.
// =============================================================================
// Thread organisation (BlockSize=256 = 8 warps = 8 elements/block):
//
//   Phase 2: lanes 0-7  load node data; all 32 zero smem_A/B + rhs + lhs
//            lanes 0-11 load SCS data
//   Phase 3: lanes 0-11 (one per SCS) — PARALLEL:
//              Jacobian → dndx → fill smem_A (S^T tile) + smem_B (D tile)
//              diffusion RHS + advection RHS  (smem atomicAdd)
//              advection LHS 4 entries        (smem atomicAdd)
//   Phase 4: all 32 lanes — WMMA × 3 tiles → smem_C = S^T * D = lhs_diff
//   Phase 5: lanes 0-7  add smem_C to w.lhs (lhs already has advection from Phase 3)
//   Phase 6: lanes 0-7  parallel CSR scatter (one lane per node row)
// =============================================================================

template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void __launch_bounds__(256, 4)
cvfem_hex_assembly_kernel_tensor_wmma(
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
    const int*     __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    static_assert(BlockSize % 32 == 0,
        "BlockSize must be a multiple of 32");
    static_assert(std::is_same<RealType, double>::value,
        "WMMA FP64 kernel requires RealType=double");

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800
    // FP64 8x8x4 WMMA available on SM80+ (A100, H100, GH200) via standard <mma.h>.
    // On SM90a (GH200) WGMMA m64n8k4 would give higher throughput — see notes.

    using namespace nvcuda;

    constexpr int WarpsPerBlock = BlockSize / 32;

    const int    warpId  = threadIdx.x / 32;
    const int    laneId  = threadIdx.x % 32;
    const size_t elemIdx = (size_t)blockIdx.x * WarpsPerBlock + warpId;

    if (elemIdx >= numElements) return;

    extern __shared__ CvfemWmmaWarpData warp_data[];
    CvfemWmmaWarpData& w = warp_data[warpId];

    // =========================================================================
    // Phase 1 — Load connectivity (all 32 lanes, broadcast from L1)
    // =========================================================================
    KeyType nodes[8];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];
    nodes[4] = d_conn4[elemIdx];
    nodes[5] = d_conn5[elemIdx];
    nodes[6] = d_conn6[elemIdx];
    nodes[7] = d_conn7[elemIdx];

    // =========================================================================
    // Phase 2 — Distribute loads + zero-init all accumulators
    //
    // smem_A[96] + smem_B[96] = 192 doubles → 6 doubles/lane (exact).
    // rhs[8] + lhs[64] = 72 doubles; lanes 0-7 zero rhs + lhs (9 doubles each).
    // smem_C is written by WMMA before it is read, so no init needed.
    // =========================================================================

    // Zero smem_A and smem_B: 192 doubles / 32 lanes = 6 doubles per lane
    {
        const int base = laneId * 6;
        w.smem_A[base]     = 0.0;  w.smem_A[base + 1] = 0.0;  w.smem_A[base + 2] = 0.0;
        w.smem_A[base + 3] = 0.0;  w.smem_A[base + 4] = 0.0;  w.smem_A[base + 5] = 0.0;
        w.smem_B[base]     = 0.0;  w.smem_B[base + 1] = 0.0;  w.smem_B[base + 2] = 0.0;
        w.smem_B[base + 3] = 0.0;  w.smem_B[base + 4] = 0.0;  w.smem_B[base + 5] = 0.0;
    }

    // Lanes 0-7: load node data + zero rhs and lhs row
    if (laneId < 8) {
        const int     lane = laneId;
        const KeyType ni   = nodes[lane];

        w.coords[lane][0]   = d_x[ni];
        w.coords[lane][1]   = d_y[ni];
        w.coords[lane][2]   = d_z[ni];
        w.phi[lane]         = d_phi[ni];
        w.gamma_n[lane]     = d_gamma[ni];
        w.beta[lane]        = d_beta[ni];
        w.grad_phi[lane][0] = d_grad_phi_x[ni];
        w.grad_phi[lane][1] = d_grad_phi_y[ni];
        w.grad_phi[lane][2] = d_grad_phi_z[ni];
        w.dofs[lane]        = d_node_to_dof[ni];
        w.own[lane]         = d_ownership[ni];

        w.rhs[lane] = 0.0;
        #pragma unroll
        for (int j = 0; j < 8; ++j) w.lhs[lane][j] = 0.0;
    }

    // Lanes 0-11: load SCS data
    if (laneId < 12) {
        const int ip      = laneId;
        w.mdot[ip]        = d_mdot[elemIdx * 12 + ip];
        w.areaVec[ip][0]  = d_areaVec_x[elemIdx * 12 + ip];
        w.areaVec[ip][1]  = d_areaVec_y[elemIdx * 12 + ip];
        w.areaVec[ip][2]  = d_areaVec_z[elemIdx * 12 + ip];
    }

    __syncwarp();

    // =========================================================================
    // Phase 3 — Parallel SCS compute: diffusion tiles + advection + RHS
    //           (lanes 0-11, one SCS per lane)
    //
    // All SCS work runs simultaneously:
    //   • smem_A / smem_B filled with no bank conflicts (each lane owns a unique
    //     k_local column within its tile)
    //   • Diffusion RHS and advection (RHS + LHS) accumulated via smem atomicAdd
    //     (FP64 shared-memory atomics are hardware instructions on SM80+;
    //      on GH200 SM90a they incur < 5 cycles per operation)
    //
    // hexShapeFcn[ip][*] is sparse: exactly two non-zero entries of 0.5, at
    // nodeL and nodeR.  Exploit this to avoid 8-element loops where possible.
    // =========================================================================
    if (laneId < 12) {
        const int ip      = laneId;
        const int t       = ip >> 2;       // tile index (0, 1, 2)
        const int k_local = ip & 3;        // position within tile (0..3)

        const int nodeL   = hexLRSCV[ip * 2];
        const int nodeR   = hexLRSCV[ip * 2 + 1];

        // gamma at SCS — hexShapeFcn is 0.5 at nodeL and nodeR, 0 elsewhere
        const double gamma_ip = 0.5 * (w.gamma_n[nodeL] + w.gamma_n[nodeR]);

        // Physical shape derivatives (Jacobian inversion for this SCS)
        double dndx[8][3];
        computeShapeDerivatives(ip, w.coords, dndx);

        // Diffusion coefficients: D[ip][n] = -gamma_ip * (dndx[n] . areaVec[ip])
        double diff_coeff[8];
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            diff_coeff[n] = -(gamma_ip * (
                dndx[n][0] * w.areaVec[ip][0] +
                dndx[n][1] * w.areaVec[ip][1] +
                dndx[n][2] * w.areaVec[ip][2]));
        }

        // Fill smem_A (S^T tile t): A[nodeL][k_local] = +1, A[nodeR][k_local] = -1
        // No conflicts: each lane writes its own k_local column.
        w.smem_A[t * 32 + k_local * 8 + nodeL] =  1.0;
        w.smem_A[t * 32 + k_local * 8 + nodeR] = -1.0;

        // Fill smem_B (D tile t): B[k_local][j] = diff_coeff[j]
        // No conflicts: each lane writes its own k_local row.
        #pragma unroll
        for (int j = 0; j < 8; ++j)
            w.smem_B[t * 32 + j * 4 + k_local] = diff_coeff[j];

        // ---- Diffusion RHS --------------------------------------------------
        // rhs[nodeL] -= sum_n D[ip][n]*phi[n],  rhs[nodeR] += same
        double diff_rhs = 0.0;
        #pragma unroll
        for (int n = 0; n < 8; ++n)
            diff_rhs += diff_coeff[n] * w.phi[n];
        atomicAdd(&w.rhs[nodeL], -diff_rhs);
        atomicAdd(&w.rhs[nodeR], +diff_rhs);

        // ---- Advection -------------------------------------------------------
        // SCS midpoint — sparse hexShapeFcn: 0.5*(nodeL + nodeR)
        const double cxL = w.coords[nodeL][0], cxR = w.coords[nodeR][0];
        const double cyL = w.coords[nodeL][1], cyR = w.coords[nodeR][1];
        const double czL = w.coords[nodeL][2], czR = w.coords[nodeR][2];
        const double cx  = 0.5 * (cxL + cxR);
        const double cy  = 0.5 * (cyL + cyR);
        const double cz  = 0.5 * (czL + czR);

        const double mdot = w.mdot[ip];

        double phi_upwind, beta_upwind, dcorr;
        if (mdot > 0.0) {
            phi_upwind  = w.phi[nodeL];
            beta_upwind = w.beta[nodeL];
            dcorr = w.grad_phi[nodeL][0] * (cx - cxL)
                  + w.grad_phi[nodeL][1] * (cy - cyL)
                  + w.grad_phi[nodeL][2] * (cz - czL);
        } else {
            phi_upwind  = w.phi[nodeR];
            beta_upwind = w.beta[nodeR];
            dcorr = w.grad_phi[nodeR][0] * (cx - cxR)
                  + w.grad_phi[nodeR][1] * (cy - cyR)
                  + w.grad_phi[nodeR][2] * (cz - czR);
        }
        dcorr *= beta_upwind;

        const double adv_flux = mdot * (phi_upwind + dcorr);
        atomicAdd(&w.rhs[nodeL], -adv_flux);
        atomicAdd(&w.rhs[nodeR], +adv_flux);

        // Advection LHS: 4 entries per SCS
        const double lhsfac_L = 0.5 * (mdot + fabs(mdot));
        const double lhsfac_R = 0.5 * (mdot - fabs(mdot));
        atomicAdd(&w.lhs[nodeL][nodeL], +lhsfac_L);
        atomicAdd(&w.lhs[nodeR][nodeL], -lhsfac_L);
        atomicAdd(&w.lhs[nodeL][nodeR], +lhsfac_R);
        atomicAdd(&w.lhs[nodeR][nodeR], -lhsfac_R);
    }

    __syncwarp();

    // =========================================================================
    // Phase 4 — WMMA: C[8x8] = sum_{t=0}^{2} A_t[8x4] * B_t[4x8]
    //   = S^T[8x12] * D[12x8] = lhs_diff[8x8]
    // All 32 lanes participate — WMMA requires the full warp.
    // =========================================================================
    wmma::fragment<wmma::accumulator, 8, 8, 4, double> c_frag;
    wmma::fill_fragment(c_frag, 0.0);

    #pragma unroll
    for (int t = 0; t < 3; ++t) {
        wmma::fragment<wmma::matrix_a, 8, 8, 4, double, wmma::col_major> a_frag;
        wmma::fragment<wmma::matrix_b, 8, 8, 4, double, wmma::col_major> b_frag;
        wmma::load_matrix_sync(a_frag, w.smem_A + t * 32, 8);  // A[8x4], ldm=8
        wmma::load_matrix_sync(b_frag, w.smem_B + t * 32, 4);  // B[4x8], ldm=4
        wmma::mma_sync(c_frag, a_frag, b_frag, c_frag);
    }
    // Store col-major: smem_C[j*8+i] = C[i][j] = lhs_diff[i][j]
    wmma::store_matrix_sync(w.smem_C, c_frag, 8, wmma::mem_col_major);

    __syncwarp();

    // =========================================================================
    // Phase 5 — Merge diffusion lhs (from smem_C) into w.lhs
    //           (lanes 0-7, one row each; w.lhs already has advection from Phase 3)
    // =========================================================================
    if (laneId < 8) {
        const int i = laneId;
        #pragma unroll
        for (int j = 0; j < 8; ++j)
            w.lhs[i][j] += w.smem_C[j * 8 + i];
    }

    __syncwarp();

    // =========================================================================
    // Phase 6 — Parallel CSR scatter (lanes 0-7, one node row each)
    // =========================================================================
    if (laneId < 8) {
        const int     i       = laneId;
        const int     row_dof = w.dofs[i];
        const uint8_t own_i   = w.own[i];

        if (own_i != 0 && row_dof >= 0) {
            atomicAdd(&d_rhs[row_dof], w.rhs[i]);

            const int row_start = matrix->rowPtr[row_dof];
            const int row_end   = matrix->rowPtr[row_dof + 1];
            const int diag_pos  = matrix->diagPtr[row_dof];

            #pragma unroll
            for (int j = 0; j < 8; ++j) {
                const int    col_dof = w.dofs[j];
                if (col_dof < 0) continue;

                const double value = w.lhs[i][j];
                if (value == 0.0) continue;

                if (i == j) {
                    atomicAdd(&matrix->values[diag_pos], value);
                } else {
                    int left = row_start, right = row_end - 1;
                    while (left <= right) {
                        const int mid = (left + right) >> 1;
                        const int col = matrix->colInd[mid];
                        if (col == col_dof) {
                            atomicAdd(&matrix->values[mid], value);
                            break;
                        } else if (col < col_dof) {
                            left  = mid + 1;
                        } else {
                            right = mid - 1;
                        }
                    }
                }
            }
        }
    }

#else  // __CUDA_ARCH__ < 800
    // FP64 WMMA (8x8x4) requires SM80+ (A100 / H100 / GH200).
    // Compile with -arch=sm_80 (or higher) to enable this kernel.
    (void)d_conn0; (void)d_conn1; (void)d_conn2; (void)d_conn3;
    (void)d_conn4; (void)d_conn5; (void)d_conn6; (void)d_conn7;
    (void)numElements;
    (void)d_x; (void)d_y; (void)d_z;
    (void)d_gamma; (void)d_phi; (void)d_beta;
    (void)d_grad_phi_x; (void)d_grad_phi_y; (void)d_grad_phi_z;
    (void)d_mdot; (void)d_areaVec_x; (void)d_areaVec_y; (void)d_areaVec_z;
    (void)d_node_to_dof; (void)d_ownership; (void)matrix; (void)d_rhs;
#endif
}

// =============================================================================
// Per-element shared data for the 2-elem/warp + WGMMA kernel.
// One instance per element within the warpgroup (8 elements/block).
// =============================================================================
struct alignas(128) CvfemWgmaElemData {
    double coords[8][3];   // 192 B
    double phi[8];         //  64 B
    double gamma_n[8];     //  64 B
    double beta[8];        //  64 B
    double grad_phi[8][3]; // 192 B
    double mdot[12];       //  96 B
    double areaVec[12][3]; // 288 B
    int    dofs[8];        //  32 B
    uint8_t own[8];        //   8 B
    uint8_t _pad[8];       //   8 B  (keep rhs[] double-aligned)
    double rhs[8];         //  64 B  (advection + diffusion RHS)
    double lhs_adv[8][8];  // 512 B  (advection LHS only; diffusion added from WGMMA C)
};
// sizeof(CvfemWgmaElemData) = 1,584 bytes

// =============================================================================
// Warpgroup-level smem layout (128 threads = 1 warpgroup = 8 elements).
//
// WGMMA m64n8k4 FP64 computes C[64×8] = A[64×4] × B[4×8] with K=12=3 tiles.
//   Formulation: lhs_diff = S^T × D → transposed: lhs_diff^T = D^T × S
//
//   A[64×4] per tile t (row-major, trans_a=1):
//     a[t][(8e+n)*4 + k] = dt[e][t*4+k][n]   ← D^T_t stacked for 8 elements
//
//   B[4×8] per tile t (col-major, trans_b=0):
//     s[t][i*4 + k]  = S[t*4+k][i]           ← scatter coeff (elem-independent)
//
//   C[64×8] output: c[(8e+n)*8 + i] = lhs_diff^T[n][i] = lhs_diff[i][n] for elem e
//
// WGMMA register layout (thread t in warpgroup, 0-127):
//   row      = t % 64
//   col_base = (t >= 64) ? 4 : 0
//   d{0..3} = C[row][col_base+{0..3}]
//
// All A/B/C smem regions are 128-byte aligned (verified by section sizes).
// =============================================================================
struct alignas(128) CvfemWgmaWarpgroupData {
    CvfemWgmaElemData elem[8];  // 8 × 1,584 = 12,672 B
    double dt[8][12][8];        //            =  6,144 B  — D^T[elem][ip][node]
    double a[3][64 * 4];        //            =  6,144 B  — A tiles (row-major [64×4])
    double s[3][4  * 8];        //            =    768 B  — B tiles (col-major [4×8])
    double c[64 * 8];           //            =  4,096 B  — WGMMA C output, lhs_diff^T
};
// Total: 29,824 B ≈ 29 KB → ≤ 7 blocks/SM from smem; register-limit typically 3-4.

// =============================================================================
// WGMMA + Option-A CVFEM Assembly Kernel (SM90+ / Hopper GH200)
//
// Thread organisation (BlockSize=128 = 1 warpgroup = 4 warps = 8 elements):
//   • 2 elements per warp: elemA = warpId*2, elemB = warpId*2+1
//   • Phase 2: lanes 0-7 load node data for elemA,
//              lanes 8-15 load node data for elemB,
//              lanes 0-11 load SCS for elemA, lanes 16-27 load SCS for elemB.
//   • Phase 3: lanes 0-11 (elemA) + 16-27 (elemB) — 75% lane utilisation —
//              compute Jacobians, fill dt[], accumulate advection via smem atomicAdd.
//   • Phase 4: all 128 threads cooperatively fill A tiles from dt[] and S tiles.
//   • Phase 5: WGMMA m64n8k4 × 3 tiles (wgmma.mma_async.sync.aligned).
//   • Phase 6: each thread stores d{0..3} to smem_C.
//   • Phase 7: threads 0-63 (1 per element-row) merge lhs_adv + lhs_diff → CSR.
// =============================================================================
template<typename KeyType, typename RealType, int BlockSize = 128>
__global__ void __launch_bounds__(128, 4)
cvfem_hex_assembly_kernel_wgmma(
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
    const int*     __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    static_assert(BlockSize == 128, "WGMMA kernel requires BlockSize=128 (1 warpgroup per block)");
    static_assert(std::is_same<RealType, double>::value, "WGMMA FP64 kernel requires RealType=double");

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800
    // Phases 1-3 are common to both WGMMA (SM90a) and WMMA (SM80+) fallback paths.
    // Phase 4 onward splits on __CUDA_ARCH_FEAT_SM90_ALL.

    // ---- thread / element indices ----
    const int    warpId     = threadIdx.x / 32;
    const int    laneId     = threadIdx.x % 32;
    const int    elemLocalA = warpId * 2;          // primary element (lanes 0-11 handle its SCS)
    const int    elemLocalB = warpId * 2 + 1;      // secondary element (lanes 16-27)
    const size_t elemGlobalA = (size_t)blockIdx.x * 8 + elemLocalA;
    const size_t elemGlobalB = (size_t)blockIdx.x * 8 + elemLocalB;
    const bool   validA = (elemGlobalA < numElements);
    const bool   validB = (elemGlobalB < numElements);

    extern __shared__ CvfemWgmaWarpgroupData wg_data[];
    CvfemWgmaWarpgroupData& g = wg_data[0];

    // =========================================================================
    // Phase 1 — Load connectivity (registers; all lanes redundantly load both)
    // =========================================================================
    KeyType nodesA[8], nodesB[8];
    if (validA) {
        nodesA[0]=d_conn0[elemGlobalA]; nodesA[1]=d_conn1[elemGlobalA];
        nodesA[2]=d_conn2[elemGlobalA]; nodesA[3]=d_conn3[elemGlobalA];
        nodesA[4]=d_conn4[elemGlobalA]; nodesA[5]=d_conn5[elemGlobalA];
        nodesA[6]=d_conn6[elemGlobalA]; nodesA[7]=d_conn7[elemGlobalA];
    }
    if (validB) {
        nodesB[0]=d_conn0[elemGlobalB]; nodesB[1]=d_conn1[elemGlobalB];
        nodesB[2]=d_conn2[elemGlobalB]; nodesB[3]=d_conn3[elemGlobalB];
        nodesB[4]=d_conn4[elemGlobalB]; nodesB[5]=d_conn5[elemGlobalB];
        nodesB[6]=d_conn6[elemGlobalB]; nodesB[7]=d_conn7[elemGlobalB];
    }

    // =========================================================================
    // Phase 2 — Zero accumulators + load node/SCS data into smem
    //
    // Threads 0-63 zero rhs[8] and lhs_adv[8][8] for their (element, row) pair.
    // Remaining load node/SCS data for both elements in each warp.
    // =========================================================================
    // Threads 0-63: zero rhs[row] and lhs_adv[row][0..7]
    if (threadIdx.x < 64) {
        const int e = threadIdx.x / 8;
        const int r = threadIdx.x % 8;
        g.elem[e].rhs[r] = 0.0;
        #pragma unroll
        for (int c = 0; c < 8; ++c) g.elem[e].lhs_adv[r][c] = 0.0;
    }

    // Lanes 0-7: node data for elemA
    if (laneId < 8 && validA) {
        const int n = laneId;
        const KeyType ni = nodesA[n];
        g.elem[elemLocalA].coords[n][0]   = d_x[ni];
        g.elem[elemLocalA].coords[n][1]   = d_y[ni];
        g.elem[elemLocalA].coords[n][2]   = d_z[ni];
        g.elem[elemLocalA].phi[n]         = d_phi[ni];
        g.elem[elemLocalA].gamma_n[n]     = d_gamma[ni];
        g.elem[elemLocalA].beta[n]        = d_beta[ni];
        g.elem[elemLocalA].grad_phi[n][0] = d_grad_phi_x[ni];
        g.elem[elemLocalA].grad_phi[n][1] = d_grad_phi_y[ni];
        g.elem[elemLocalA].grad_phi[n][2] = d_grad_phi_z[ni];
        g.elem[elemLocalA].dofs[n]        = d_node_to_dof[ni];
        g.elem[elemLocalA].own[n]         = d_ownership[ni];
    }

    // Lanes 8-15: node data for elemB (mapped: node = laneId-8)
    if (laneId >= 8 && laneId < 16 && validB) {
        const int n = laneId - 8;
        const KeyType ni = nodesB[n];
        g.elem[elemLocalB].coords[n][0]   = d_x[ni];
        g.elem[elemLocalB].coords[n][1]   = d_y[ni];
        g.elem[elemLocalB].coords[n][2]   = d_z[ni];
        g.elem[elemLocalB].phi[n]         = d_phi[ni];
        g.elem[elemLocalB].gamma_n[n]     = d_gamma[ni];
        g.elem[elemLocalB].beta[n]        = d_beta[ni];
        g.elem[elemLocalB].grad_phi[n][0] = d_grad_phi_x[ni];
        g.elem[elemLocalB].grad_phi[n][1] = d_grad_phi_y[ni];
        g.elem[elemLocalB].grad_phi[n][2] = d_grad_phi_z[ni];
        g.elem[elemLocalB].dofs[n]        = d_node_to_dof[ni];
        g.elem[elemLocalB].own[n]         = d_ownership[ni];
    }

    // Lanes 0-11: SCS data for elemA
    if (laneId < 12 && validA) {
        const int ip = laneId;
        g.elem[elemLocalA].mdot[ip]       = d_mdot[elemGlobalA * 12 + ip];
        g.elem[elemLocalA].areaVec[ip][0] = d_areaVec_x[elemGlobalA * 12 + ip];
        g.elem[elemLocalA].areaVec[ip][1] = d_areaVec_y[elemGlobalA * 12 + ip];
        g.elem[elemLocalA].areaVec[ip][2] = d_areaVec_z[elemGlobalA * 12 + ip];
    }

    // Lanes 16-27: SCS data for elemB
    if (laneId >= 16 && laneId < 28 && validB) {
        const int ip = laneId - 16;
        g.elem[elemLocalB].mdot[ip]       = d_mdot[elemGlobalB * 12 + ip];
        g.elem[elemLocalB].areaVec[ip][0] = d_areaVec_x[elemGlobalB * 12 + ip];
        g.elem[elemLocalB].areaVec[ip][1] = d_areaVec_y[elemGlobalB * 12 + ip];
        g.elem[elemLocalB].areaVec[ip][2] = d_areaVec_z[elemGlobalB * 12 + ip];
    }

    __syncthreads();

    // =========================================================================
    // Phase 3 — Parallel SCS: Jacobian → D^T fill + advection
    //   Lanes 0-11  → SCS 0-11 of elemA  (75% lane utilization per warp)
    //   Lanes 16-27 → SCS 0-11 of elemB
    //   Lanes 12-15, 28-31 idle.
    // =========================================================================
    const bool activeA = (laneId < 12) && validA;
    const bool activeB = (laneId >= 16 && laneId < 28) && validB;

    if (activeA || activeB) {
        const int ip      = activeA ? laneId : (laneId - 16);
        const int eLocal  = activeA ? elemLocalA : elemLocalB;
        CvfemWgmaElemData& e = g.elem[eLocal];

        const int nodeL = hexLRSCV[ip * 2];
        const int nodeR = hexLRSCV[ip * 2 + 1];

        const double gamma_ip = 0.5 * (e.gamma_n[nodeL] + e.gamma_n[nodeR]);

        double dndx[8][3];
        computeShapeDerivatives(ip, e.coords, dndx);

        // Diffusion coefficients: D[ip][n] = -gamma_ip * (dndx[n] · areaVec[ip])
        double diff_coeff[8];
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            diff_coeff[n] = -(gamma_ip * (
                dndx[n][0] * e.areaVec[ip][0] +
                dndx[n][1] * e.areaVec[ip][1] +
                dndx[n][2] * e.areaVec[ip][2]));
        }

        // Fill D^T: dt[eLocal][ip][n] = diff_coeff[n]  (each lane owns unique ip row)
        #pragma unroll
        for (int n = 0; n < 8; ++n)
            g.dt[eLocal][ip][n] = diff_coeff[n];

        // Diffusion RHS
        double diff_rhs = 0.0;
        #pragma unroll
        for (int n = 0; n < 8; ++n)
            diff_rhs += diff_coeff[n] * e.phi[n];
        atomicAdd(&e.rhs[nodeL], -diff_rhs);
        atomicAdd(&e.rhs[nodeR], +diff_rhs);

        // Advection (SCS midpoint from sparse hexShapeFcn = 0.5 at nodeL and nodeR)
        const double cxL = e.coords[nodeL][0], cxR = e.coords[nodeR][0];
        const double cyL = e.coords[nodeL][1], cyR = e.coords[nodeR][1];
        const double czL = e.coords[nodeL][2], czR = e.coords[nodeR][2];
        const double cx  = 0.5 * (cxL + cxR);
        const double cy  = 0.5 * (cyL + cyR);
        const double cz  = 0.5 * (czL + czR);

        const double mdot = e.mdot[ip];

        double phi_upwind, beta_upwind, dcorr;
        if (mdot > 0.0) {
            phi_upwind  = e.phi[nodeL];
            beta_upwind = e.beta[nodeL];
            dcorr = e.grad_phi[nodeL][0]*(cx-cxL) + e.grad_phi[nodeL][1]*(cy-cyL) + e.grad_phi[nodeL][2]*(cz-czL);
        } else {
            phi_upwind  = e.phi[nodeR];
            beta_upwind = e.beta[nodeR];
            dcorr = e.grad_phi[nodeR][0]*(cx-cxR) + e.grad_phi[nodeR][1]*(cy-cyR) + e.grad_phi[nodeR][2]*(cz-czR);
        }
        dcorr *= beta_upwind;

        const double adv_flux = mdot * (phi_upwind + dcorr);
        atomicAdd(&e.rhs[nodeL], -adv_flux);
        atomicAdd(&e.rhs[nodeR], +adv_flux);

        const double lhsfac_L = 0.5 * (mdot + fabs(mdot));
        const double lhsfac_R = 0.5 * (mdot - fabs(mdot));
        atomicAdd(&e.lhs_adv[nodeL][nodeL], +lhsfac_L);
        atomicAdd(&e.lhs_adv[nodeR][nodeL], -lhsfac_L);
        atomicAdd(&e.lhs_adv[nodeL][nodeR], +lhsfac_R);
        atomicAdd(&e.lhs_adv[nodeR][nodeR], -lhsfac_R);
    }

    __syncthreads();

#if defined(__CUDA_ARCH_FEAT_SM90_ALL)
    // =========================================================================
    // WGMMA path (SM90a / Hopper): async wgmma.mma_async m64n8k4 FP64
    // Rebuild with -DCMAKE_CUDA_ARCHITECTURES=90a to enable this path.
    // =========================================================================

    // Phase 4 — Fill WGMMA A tiles and S tiles cooperatively (all 128 threads)
    //
    // A[t][(8e+n)*4 + k] = dt[e][t*4+k][n]   row-major [64×4], stride=32B
    // S[t][i*4 + k]      = scatter(ip=t*4+k, node=i)  col-major [4×8], stride=32B
    // =========================================================================

    // Fill S tiles: 3 × 32 = 96 doubles.  Threads 0-95 each write 1 value.
    if (threadIdx.x < 96) {
        const int t       = threadIdx.x / 32;  // tile 0-2
        const int idx     = threadIdx.x % 32;  // flat index in B[4×8] col-major
        const int i       = idx / 4;            // column (node index) 0-7
        const int k       = idx % 4;            // row (K within tile) 0-3
        const int ip      = t * 4 + k;
        const int L = hexLRSCV[ip * 2];
        const int R = hexLRSCV[ip * 2 + 1];
        g.s[t][idx] = (i == L) ? 1.0 : (i == R) ? -1.0 : 0.0;
    }

    // Fill A tiles: 3 × 256 = 768 doubles.  Each thread fills 6 values (768/128=6).
    #pragma unroll
    for (int stride = 0; stride < 6; ++stride) {
        const int global_idx = threadIdx.x + stride * 128;
        // global_idx ∈ [0, 768): tile t, row (8e+n), col k_local
        const int t       = global_idx / 256;  // tile 0-2
        const int flat_A  = global_idx % 256;  // index within A[64×4]
        const int row     = flat_A / 4;         // 8e+n
        const int k       = flat_A % 4;         // k_local
        const int e       = row / 8;
        const int n       = row % 8;
        g.a[t][flat_A] = g.dt[e][t * 4 + k][n];
    }

    __syncthreads();

    // =========================================================================
    // Phase 5 — WGMMA m64n8k4 FP64 (SM90+)
    //   Three tiles: C[64×8] = A0[64×4]×B0[4×8] + A1×B1 + A2×B2
    //   A row-major (trans_a=1, stride=32B), B col-major (trans_b=0, stride=32B).
    //
    // Descriptor encoding (PTX ISA 8.5):
    //   bits [13:0]  = smem_addr >> 4
    //   bits [29:16] = leading_dim_stride_bytes >> 4
    //   (no swizzle, base_offset=0)
    // =========================================================================
    auto make_desc = [](const void* ptr, uint32_t stride_bytes) -> uint64_t {
        uint64_t desc = 0;
        const uint64_t addr = static_cast<uint64_t>(
                static_cast<uint32_t>(__cvta_generic_to_shared(ptr)));
        desc |= (addr >> 4) & 0x3FFFull;
        desc |= static_cast<uint64_t>((stride_bytes >> 4) & 0x3FFF) << 16;
        return desc;
    };

    // A tiles: row-major [64×4], row-stride = 4 doubles = 32 bytes
    const uint64_t desc_a0 = make_desc(g.a[0], 32u);
    const uint64_t desc_a1 = make_desc(g.a[1], 32u);
    const uint64_t desc_a2 = make_desc(g.a[2], 32u);
    // B tiles: col-major [4×8], col-stride = 4 doubles = 32 bytes
    const uint64_t desc_b0 = make_desc(g.s[0], 32u);
    const uint64_t desc_b1 = make_desc(g.s[1], 32u);
    const uint64_t desc_b2 = make_desc(g.s[2], 32u);

    double d0, d1, d2, d3;

    asm volatile("wgmma.fence.sync.aligned;\n" :::);

    // Tile 0 — scale_d=0 (initialise C, not accumulate)
    asm volatile(
        "wgmma.mma_async.sync.aligned.m64n8k4.f64.f64.f64 "
        "{%0,%1,%2,%3}, %4, %5, 0, 1, 1, 1, 0;\n"
        : "=d"(d0),"=d"(d1),"=d"(d2),"=d"(d3)
        : "l"(desc_a0),"l"(desc_b0));

    // Tile 1 — scale_d=1 (accumulate)
    asm volatile(
        "wgmma.mma_async.sync.aligned.m64n8k4.f64.f64.f64 "
        "{%0,%1,%2,%3}, %4, %5, 1, 1, 1, 1, 0;\n"
        : "+d"(d0),"+d"(d1),"+d"(d2),"+d"(d3)
        : "l"(desc_a1),"l"(desc_b1));

    // Tile 2 — scale_d=1 (accumulate)
    asm volatile(
        "wgmma.mma_async.sync.aligned.m64n8k4.f64.f64.f64 "
        "{%0,%1,%2,%3}, %4, %5, 1, 1, 1, 1, 0;\n"
        : "+d"(d0),"+d"(d1),"+d"(d2),"+d"(d3)
        : "l"(desc_a2),"l"(desc_b2));

    asm volatile("wgmma.commit_group.sync.aligned;\n" :::);
    asm volatile("wgmma.wait_group.sync.aligned 0;\n" :::);

    // =========================================================================
    // Phase 6 — Store C registers to smem_C
    //
    // Thread t (0-127) holds d{0..3} = C[t%64][col_base + 0..3]
    // where col_base = (t >= 64) ? 4 : 0.
    // c[(8e+n)*8 + i] = lhs_diff^T[n][i] = lhs_diff[i][n]
    // =========================================================================
    {
        const int row      = threadIdx.x % 64;          // 8e+n
        const int col_base = (threadIdx.x >= 64) ? 4 : 0;
        g.c[row * 8 + col_base + 0] = d0;
        g.c[row * 8 + col_base + 1] = d1;
        g.c[row * 8 + col_base + 2] = d2;
        g.c[row * 8 + col_base + 3] = d3;
    }

    __syncthreads();

    // =========================================================================
    // Phase 7 — Merge lhs + CSR scatter (threads 0-63, one per element-row)
    //   Thread 8e+r owns element e, node row r.
    //   lhs_diff[r][col_j] = g.c[(8e+col_j)*8 + r]  (transposed from C storage)
    // =========================================================================
    if (threadIdx.x < 64) {
        const int eLocal = threadIdx.x / 8;
        const int row_i  = threadIdx.x % 8;

        const size_t eGlobal = (size_t)blockIdx.x * 8 + eLocal;
        if (eGlobal >= numElements) return;

        const CvfemWgmaElemData& e = g.elem[eLocal];
        if (e.own[row_i] == 0 || e.dofs[row_i] < 0) return;

        const int row_dof   = e.dofs[row_i];
        atomicAdd(&d_rhs[row_dof], e.rhs[row_i]);

        const int row_start = matrix->rowPtr[row_dof];
        const int row_end   = matrix->rowPtr[row_dof + 1];
        const int diag_pos  = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int col_j = 0; col_j < 8; ++col_j) {
            const int col_dof = e.dofs[col_j];
            if (col_dof < 0) continue;

            // lhs_diff[row_i][col_j] stored transposed in smem_C
            const double value = e.lhs_adv[row_i][col_j]
                               + g.c[(8 * eLocal + col_j) * 8 + row_i];
            if (value == 0.0) continue;

            if (row_i == col_j) {
                atomicAdd(&matrix->values[diag_pos], value);
            } else {
                int left = row_start, right = row_end - 1;
                while (left <= right) {
                    const int mid = (left + right) >> 1;
                    const int col = matrix->colInd[mid];
                    if (col == col_dof) {
                        atomicAdd(&matrix->values[mid], value);
                        break;
                    } else if (col < col_dof) {
                        left  = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                }
            }
        }
    }

#else // !__CUDA_ARCH_FEAT_SM90_ALL  →  WMMA fallback (SM80-SM89 or sm_90 without sm_90a)
    // =========================================================================
    // WMMA fallback: same 2-elem/warp structure (75% Jacobian utilization),
    // but uses nvcuda::wmma 8×8×4 FP64 instead of WGMMA for the diffusion step.
    //
    // Reuses smem arrays with different interpretations:
    //   g.s[t][k*8+i]            = S^T[i][k]  WMMA A col-major ldm=8
    //   g.a[t][e*32 + col*4+row] = D[row][col] WMMA B col-major ldm=4 for elem e
    //   g.c[e*64 + col*8+row]    = lhs_diff[row][col] col-major ldm=8 for elem e
    // =========================================================================

    // Phase 4 — Fill WMMA tiles cooperatively
    using namespace nvcuda;

    // A tiles (S^T[8×4], element-independent, col-major ldm=8):
    // g.s[t][k*8+i] = S^T[i][k]
    if (threadIdx.x < 96) {
        const int t  = threadIdx.x / 32;
        const int k  = (threadIdx.x % 32) / 8;
        const int i  = (threadIdx.x % 32) % 8;
        const int ip = t * 4 + k;
        const int L  = hexLRSCV[ip * 2];
        const int R  = hexLRSCV[ip * 2 + 1];
        g.s[t][k * 8 + i] = (i == L) ? 1.0 : (i == R) ? -1.0 : 0.0;
    }

    // B tiles (D[4×8] per element, col-major ldm=4):
    // g.a[t][e*32 + col*4+row] = D[t*4+row][col] = g.dt[e][t*4+row][col]
    #pragma unroll
    for (int stride = 0; stride < 6; ++stride) {
        const int gidx = threadIdx.x + stride * 128;
        const int t    = gidx / 256;
        const int flat = gidx % 256;
        const int e    = flat / 32;
        const int col  = (flat % 32) / 4;
        const int row  = (flat % 32) % 4;
        g.a[t][flat] = g.dt[e][t * 4 + row][col];
    }

    __syncthreads();

    // Phase 5 — WMMA per warp: elemA then elemB (both uniform across warp → no divergence)
    if (validA) {
        wmma::fragment<wmma::accumulator, 8, 8, 4, double> c_frag;
        wmma::fill_fragment(c_frag, 0.0);
        #pragma unroll
        for (int t = 0; t < 3; ++t) {
            wmma::fragment<wmma::matrix_a, 8, 8, 4, double, wmma::col_major> a_frag;
            wmma::fragment<wmma::matrix_b, 8, 8, 4, double, wmma::col_major> b_frag;
            wmma::load_matrix_sync(a_frag, &g.s[t][0],              8);  // S^T ldm=8
            wmma::load_matrix_sync(b_frag, &g.a[t][elemLocalA * 32], 4); // D_A ldm=4
            wmma::mma_sync(c_frag, a_frag, b_frag, c_frag);
        }
        wmma::store_matrix_sync(&g.c[elemLocalA * 64], c_frag, 8, wmma::mem_col_major);
    }

    if (validB) {
        wmma::fragment<wmma::accumulator, 8, 8, 4, double> c_frag;
        wmma::fill_fragment(c_frag, 0.0);
        #pragma unroll
        for (int t = 0; t < 3; ++t) {
            wmma::fragment<wmma::matrix_a, 8, 8, 4, double, wmma::col_major> a_frag;
            wmma::fragment<wmma::matrix_b, 8, 8, 4, double, wmma::col_major> b_frag;
            wmma::load_matrix_sync(a_frag, &g.s[t][0],              8);
            wmma::load_matrix_sync(b_frag, &g.a[t][elemLocalB * 32], 4);
            wmma::mma_sync(c_frag, a_frag, b_frag, c_frag);
        }
        wmma::store_matrix_sync(&g.c[elemLocalB * 64], c_frag, 8, wmma::mem_col_major);
    }

    __syncthreads();

    // Phase 6 — CSR scatter (threads 0-63, one per element-row)
    // lhs_diff[row_i][col_j] = g.c[eLocal*64 + col_j*8 + row_i] (col-major ldm=8)
    if (threadIdx.x < 64) {
        const int eLocal = threadIdx.x / 8;
        const int row_i  = threadIdx.x % 8;

        const size_t eGlobal = (size_t)blockIdx.x * 8 + eLocal;
        if (eGlobal >= numElements) return;

        const CvfemWgmaElemData& e = g.elem[eLocal];
        if (e.own[row_i] == 0 || e.dofs[row_i] < 0) return;

        const int row_dof = e.dofs[row_i];
        atomicAdd(&d_rhs[row_dof], e.rhs[row_i]);

        const int row_start = matrix->rowPtr[row_dof];
        const int row_end   = matrix->rowPtr[row_dof + 1];
        const int diag_pos  = matrix->diagPtr[row_dof];

        #pragma unroll
        for (int col_j = 0; col_j < 8; ++col_j) {
            const int col_dof = e.dofs[col_j];
            if (col_dof < 0) continue;

            const double value = e.lhs_adv[row_i][col_j]
                               + g.c[eLocal * 64 + col_j * 8 + row_i];
            if (value == 0.0) continue;

            if (row_i == col_j) {
                atomicAdd(&matrix->values[diag_pos], value);
            } else {
                int left = row_start, right = row_end - 1;
                while (left <= right) {
                    const int mid = (left + right) >> 1;
                    const int col = matrix->colInd[mid];
                    if (col == col_dof) {
                        atomicAdd(&matrix->values[mid], value); break;
                    } else if (col < col_dof) { left  = mid + 1; }
                    else                      { right = mid - 1; }
                }
            }
        }
    }

#endif // __CUDA_ARCH_FEAT_SM90_ALL

#else  // < SM80: FP64 WMMA not available
    (void)d_conn0; (void)d_conn1; (void)d_conn2; (void)d_conn3;
    (void)d_conn4; (void)d_conn5; (void)d_conn6; (void)d_conn7;
    (void)numElements;
    (void)d_x; (void)d_y; (void)d_z;
    (void)d_gamma; (void)d_phi; (void)d_beta;
    (void)d_grad_phi_x; (void)d_grad_phi_y; (void)d_grad_phi_z;
    (void)d_mdot; (void)d_areaVec_x; (void)d_areaVec_y; (void)d_areaVec_z;
    (void)d_node_to_dof; (void)d_ownership; (void)matrix; (void)d_rhs;
#endif // __CUDA_ARCH__ >= 800
}

} // namespace fem
} // namespace mars
