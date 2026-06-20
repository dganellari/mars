#pragma once

// Matrix-free high-order Laplacian apply on FP64 tensor cores (DMMA).
//
// Drop-in numerical twin of the scalar ho_laplacian_apply_cg in
// mars_ho_laplacian.hpp: same elemDof gather/scatter, same per-point metric G
// (6 components Gxx,Gxy,Gxz,Gyy,Gyz,Gzz), same 1D differentiation matrix D,
// same contiguous node layout u[i*n*n + j*n + k] (i slowest, k fastest). The
// only allowed difference is floating-point reduction order, so the two paths
// agree to ~1e-15 relative on the same inputs.
//
// Why tensor cores: sum-factorization turns the element Laplacian into a small
// stack of dense 1D contractions C = D . U. At P=7 (n=8) each contraction is an
// 8x8x8 GEMM, which maps cleanly onto the SM80+ FP64 WMMA m8n8k4 fragment
// (K=8 = two k-tiles of 4). This first phase targets P=7 only -- the native
// 8-wide tile -- so there is no tile padding and the whole warp stays busy.
//
// Organisation: ONE WARP per element, all 32 lanes drive every WMMA op (the API
// requires the full warp). A block holds several warps. Per warp we stage in
// shared memory: D, D^T, the input tensor, the three derivative tensors, the
// output accumulator, and one reusable A/B/C tile triplet for the MMA.
//
// Launch contract (the caller must honor these; see ho_laplacian_dmma_launch):
//   - grid.x = ceil(numElements / WarpsPerBlock), block = BlockSize.
//   - d_y_dof zeroed before launch (scatter is additive via atomicAdd).
//   - per-warp smem is 22528 B, so a 4-warp block needs 88 KB DYNAMIC smem,
//     above the 48 KB static cap. The launcher must opt in once via
//     cudaFuncSetAttribute(..., cudaFuncAttributeMaxDynamicSharedMemorySize, ...)
//     and pass WarpsPerBlock*sizeof(HoLaplacianDmmaWarpData) as the 3rd <<<>>>
//     argument. ho_laplacian_dmma_launch below does exactly this.

#include <cuda_runtime.h>
#include <cstddef>
#include <mma.h>
#include <type_traits>

namespace mars {
namespace fem {

// Per-warp shared state for the DMMA Laplacian apply at P=7 (n=8, N3=512).
// One instance per warp in the block; sized so the host can pass the right
// dynamic-smem byte count to the launch.
//
// Buffer reuse over the element's lifetime:
//   gather      -> b_u
//   forward     -> b0 = D.u(dir0), b1 = D.u(dir1), b2 = D.u(dir2)   (read b_u)
//   metric      -> overwrite b0,b1,b2 in place with gr,gs,gt        (pointwise)
//   integrate   -> b_y += D^T.gr(dir0) + D^T.gs(dir1) + D^T.gt(dir2)
// Occupancy is smem-capped, so the working set is trimmed two ways:
//   - DT (the transpose of D) is NOT stored; the integrate-back reads D
//     transposed on the fly (D[s*8+a] instead of DT[a*8+s]). Saves 64 doubles.
//   - b_u (input) and b_y (output) never live at the same time -- the input is
//     dead once the three forward contractions are done, before the output is
//     zeroed/accumulated -- so they SHARE storage via a union. Saves 512 doubles.
// sizeof = 17920 B/warp (D + 3*512 + union 512 + A+B+C). Smaller blocks then fit
// more warps/SM (the 12.5%-occupancy cap was set by this footprint).
struct alignas(128) HoLaplacianDmmaWarpData {
    double D[64];     // 1D differentiation matrix, row-major D[a*8+s]
    double b0[512];   // dir0 derivative -> gr after metric
    double b1[512];   // dir1 derivative -> gs after metric
    double b2[512];   // dir2 derivative -> gt after metric
    union {
        double b_u[512];  // gathered input  u[i*64 + j*8 + k]  (forward pass)
        double b_y[512];  // output accumulator y[l]            (integrate-back)
    };

    // One reusable WMMA tile triplet. A is [8x4] col-major ldm=8, B is [4x8]
    // col-major ldm=4, C is [8x8] col-major ldm=8 (C[i][j] at smem_C[j*8+i]).
    double smem_A[32];
    double smem_B[32];
    double smem_C[64];
};

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800

// Run one batched 8x8x8 contraction tile: C[8x8] = A[8x8] . B[8x8] via two
// m8n8k4 MMAs (K = two tiles of 4). `fill_A(kt)` / `fill_B(kt)` are warp-wide
// lambdas that stage k-tile kt of the A and B operands into w.smem_A / w.smem_B
// (32 doubles each, written cooperatively across the 32 lanes). Result lands in
// w.smem_C, col-major ldm=8: C[i][j] = w.smem_C[j*8 + i].
template<typename FillA, typename FillB>
__device__ inline void dmma_tile_8x8x8(HoLaplacianDmmaWarpData& w,
                                        FillA fill_A, FillB fill_B)
{
    using namespace nvcuda;
    wmma::fragment<wmma::accumulator, 8, 8, 4, double> c_frag;
    wmma::fill_fragment(c_frag, 0.0);

    #pragma unroll
    for (int kt = 0; kt < 2; ++kt)
    {
        fill_A(kt);
        fill_B(kt);
        __syncwarp();

        wmma::fragment<wmma::matrix_a, 8, 8, 4, double, wmma::col_major> a_frag;
        wmma::fragment<wmma::matrix_b, 8, 8, 4, double, wmma::col_major> b_frag;
        wmma::load_matrix_sync(a_frag, w.smem_A, 8);  // A[8x4] col-major, ldm=8
        wmma::load_matrix_sync(b_frag, w.smem_B, 4);  // B[4x8] col-major, ldm=4
        wmma::mma_sync(c_frag, a_frag, b_frag, c_frag);
        __syncwarp();                                 // tile reused next kt
    }

    wmma::store_matrix_sync(w.smem_C, c_frag, 8, wmma::mem_col_major);
    __syncwarp();
}

// Scalar (CUDA-core) twin of dmma_tile_8x8x8: IDENTICAL staging (fill_A/fill_B
// into smem_A[8x4] / smem_B[4x8] per k-tile) and IDENTICAL smem_C col-major
// output, but the 8x8x8 product runs as scalar FMAs across the 32 lanes instead
// of the WMMA tensor-core op. This lets the benchmark isolate the tensor-core
// contribution: same warp/smem structure, only the contraction hardware differs.
template<typename FillA, typename FillB>
__device__ inline void scalar_tile_8x8x8(HoLaplacianDmmaWarpData& w,
                                         FillA fill_A, FillB fill_B)
{
    const int lane = threadIdx.x % 32;
    double acc[2] = {0.0, 0.0};   // this lane owns outputs o=lane and o=lane+32
    #pragma unroll
    for (int kt = 0; kt < 2; ++kt)
    {
        fill_A(kt);
        fill_B(kt);
        __syncwarp();
        #pragma unroll
        for (int t = 0; t < 2; ++t)
        {
            const int o = lane + t * 32;       // 0..63, col-major o = j*8 + i
            const int i = o & 7, j = o >> 3;
            double s = 0.0;
            #pragma unroll
            for (int k = 0; k < 4; ++k)
                s += w.smem_A[k * 8 + i] * w.smem_B[j * 4 + k];
            acc[t] += s;
        }
        __syncwarp();
    }
    w.smem_C[lane]      = acc[0];
    w.smem_C[lane + 32] = acc[1];
    __syncwarp();
}

// Dispatch the 8x8x8 contraction: TileMode 0 = WMMA tensor core, 1 = scalar.
template<int TileMode, typename FillA, typename FillB>
__device__ inline void tile_8x8x8(HoLaplacianDmmaWarpData& w, FillA fill_A, FillB fill_B)
{
    if constexpr (TileMode == 0) dmma_tile_8x8x8(w, fill_A, fill_B);
    else                         scalar_tile_8x8x8(w, fill_A, fill_B);
}

#endif // __CUDA_ARCH__ >= 800

// CG matrix-free Laplacian apply on FP64 tensor cores (DMMA). One WARP/element.
//   d_u_dof  : global input field
//   d_y_dof  : global output (must be zeroed before launch; scatter is additive)
//   d_elemDof: [numElements * N3] local-node -> global-DOF
//   d_D      : [n*n] 1D differentiation matrix, shared by all elements
//   d_G      : [numElements * N3 * 6] per-point symmetric metric
// BlockSize must be a multiple of 32 (one warp = one element).
// GMode selects how the metric G is stored/read:
//   0 (default): per-quad-point, d_G is [numElements * N3 * 6] (general/curved).
//   1 (affine):  per-ELEMENT constant geometric metric Ghat, d_G is
//                [numElements * 6], and the GLL quadrature weight w_i*w_j*w_k is
//                applied in-kernel from d_w[n]. For straight-sided (affine) hexes
//                Ghat = J^{-1}J^{-T}|J| is constant over the element, so this
//                reads 6 doubles/element instead of 6*N3 -- ~512x less metric
//                traffic at P=7, which is the dominant apply bandwidth. d_w is
//                unused (pass nullptr) when GMode==0.
template<typename RealType, int P, int BlockSize = 128, int GMode = 0, int TileMode = 0>
__global__ void __launch_bounds__(BlockSize)
ho_laplacian_apply_cg_dmma(const RealType* __restrict__ d_u_dof,
                           RealType* __restrict__ d_y_dof,
                           const int* __restrict__ d_elemDof,
                           const RealType* __restrict__ d_D,
                           const RealType* __restrict__ d_G,
                           const RealType* __restrict__ d_w,
                           size_t numElements)
{
    static_assert(std::is_same<RealType, double>::value,
        "DMMA Laplacian kernel requires RealType=double (FP64 tensor cores)");
    static_assert(P == 7,
        "DMMA Laplacian phase 1 targets P=7 (n=8): the native m8n8k4 tile size");
    static_assert(BlockSize % 32 == 0, "BlockSize must be a multiple of 32");

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800
    using namespace nvcuda;

    constexpr int n   = P + 1;   // 8
    constexpr int nn  = n * n;   // 64
    constexpr int N3  = nn * n;  // 512
    constexpr int WarpsPerBlock = BlockSize / 32;

    const int    warpId  = threadIdx.x / 32;
    const int    laneId  = threadIdx.x % 32;
    // e depends only on blockIdx/warpId (warp-uniform), so the early return
    // retires whole warps -- a partial warp never reaches a wmma collective op.
    const size_t e       = (size_t)blockIdx.x * WarpsPerBlock + warpId;
    if (e >= numElements) return;

    extern __shared__ HoLaplacianDmmaWarpData dmma_warp[];
    HoLaplacianDmmaWarpData& w = dmma_warp[warpId];

    const int*      edof  = d_elemDof + e * N3;
    const RealType* Gbase = d_G + e * (size_t)(GMode == 0 ? N3 * 6 : 6);

    // --- Stage D and gather u (all 32 lanes share the work). D^T is not stored;
    // the integrate-back reads D transposed on the fly. ---
    for (int i = laneId; i < nn; i += 32)
        w.D[i] = d_D[i];                     // D[a*8+s], a=i/8, s=i%8
    for (int l = laneId; l < N3; l += 32)
        w.b_u[l] = d_u_dof[edof[l]];
    __syncwarp();

    // =========================================================================
    // Forward contractions: b0 = D.u (dir0), b1 = D.u (dir1), b2 = D.u (dir2).
    //
    // Each direction is 8 independent 8x8x8 GEMMs (one per outer batch), and the
    // whole warp must drive each MMA, so we loop the 8 batches serially. The
    // ONLY thing that changes between directions is how the data operand is
    // gathered into the tile (its smem leading dimension / stride) and which
    // operand carries the matrix D. Output uses the col-major C layout
    // C[i][j] = w.smem_C[j*8 + i].
    // =========================================================================

    // dir0 (contract slowest index i): C[8x64] = D[8x8] . U[8x64],
    //   U[s][bc] = u[s*64 + bc], bc = j*8+k contiguous. Tile bc in 8 chunks of 8.
    //   A = D-tile:  A[a][k] = D[a*8 + (kt*4+k)]      -> smem_A[k*8 + a], ldm_a=8.
    //   B = U-slice: B[k][col] = u[(kt*4+k)*64 + bctile*8 + col] -> smem_B[col*4+k], ldm_b=4.
    //   store: out[a*64 + bctile*8 + col] = C[a][col] = smem_C[col*8 + a].
    for (int bctile = 0; bctile < n; ++bctile)
    {
        auto fillA = [&](int kt) {
            const int a = laneId % n, k = laneId / n;   // a in 0..7, k in 0..3
            w.smem_A[k * n + a] = w.D[a * n + (kt * 4 + k)];
        };
        auto fillB = [&](int kt) {
            const int k = laneId % 4, col = laneId / 4; // k in 0..3, col in 0..7
            w.smem_B[col * 4 + k] = w.b_u[(kt * 4 + k) * nn + bctile * n + col];
        };
        tile_8x8x8<TileMode>(w, fillA, fillB);
        for (int idx = laneId; idx < nn; idx += 32)
        {
            const int a = idx / n, col = idx % n;
            w.b0[a * nn + bctile * n + col] = w.smem_C[col * n + a];
        }
        __syncwarp();
    }

    // dir1 (contract middle index j): per slab a, C_a[8x8] = D[8x8] . Xa[8x8],
    //   Xa[s][c] = u[a*64 + s*8 + c]. 8 slabs (a), each one 8x8x8 GEMM.
    //   A = D-tile:   A[b][k] = D[b*8 + (kt*4+k)]          -> smem_A[k*8 + b], ldm_a=8.
    //   B = Xa-slice: B[k][c] = u[a*64 + (kt*4+k)*8 + c]   -> smem_B[c*4 + k], ldm_b=4.
    //   store: out[a*64 + b*8 + c] = C_a[b][c] = smem_C[c*8 + b].
    for (int a = 0; a < n; ++a)
    {
        auto fillA = [&](int kt) {
            const int b = laneId % n, k = laneId / n;
            w.smem_A[k * n + b] = w.D[b * n + (kt * 4 + k)];
        };
        auto fillB = [&](int kt) {
            const int k = laneId % 4, c = laneId / 4;
            w.smem_B[c * 4 + k] = w.b_u[a * nn + (kt * 4 + k) * n + c];
        };
        tile_8x8x8<TileMode>(w, fillA, fillB);
        for (int idx = laneId; idx < nn; idx += 32)
        {
            const int b = idx / n, c = idx % n;
            w.b1[a * nn + b * n + c] = w.smem_C[c * n + b];
        }
        __syncwarp();
    }

    // dir2 (contract fastest index k): C[64x8] = X[64x8] . Dt[8x8],
    //   X[ab][s] = u[ab*8 + s], ab = i*8+j contiguous. Tile ab in 8 chunks of 8.
    //   Dt[s][c] = D[c*8 + s] is the transposed READ of D (not D^T-of-operator):
    //   the scalar dir2 uses dval = D[c*n+s], so here the small matrix is the B
    //   operand and reads D row c, column s.
    //   A = X-slice: A[ar][k] = u[(abtile*8+ar)*8 + (kt*4+k)] -> smem_A[k*8 + ar], ldm_a=8.
    //   B = Dt-tile: B[k][c]  = D[c*8 + (kt*4+k)]             -> smem_B[c*4 + k], ldm_b=4.
    //   store: out[(abtile*8+i)*8 + j] = C[i][j] = smem_C[j*8 + i].
    for (int abtile = 0; abtile < n; ++abtile)
    {
        auto fillA = [&](int kt) {
            const int ar = laneId % n, k = laneId / n;
            w.smem_A[k * n + ar] = w.b_u[(abtile * n + ar) * n + (kt * 4 + k)];
        };
        auto fillB = [&](int kt) {
            const int k = laneId % 4, c = laneId / 4;
            w.smem_B[c * 4 + k] = w.D[c * n + (kt * 4 + k)];
        };
        tile_8x8x8<TileMode>(w, fillA, fillB);
        for (int idx = laneId; idx < nn; idx += 32)
        {
            const int i = idx / n, j = idx % n;
            const int ab = abtile * n + i;
            w.b2[ab * n + j] = w.smem_C[j * n + i];
        }
        __syncwarp();
    }

    // =========================================================================
    // Metric (pointwise per quad point, CUDA cores). 512 points / 32 lanes =
    // 16 points/lane. (gr,gs,gt) = G * (ur,us,ut), done in place: each point's
    // outputs depend only on that point's three inputs, so read-then-write is safe.
    // =========================================================================
    for (int p = laneId; p < N3; p += 32)
    {
        double Gxx, Gxy, Gxz, Gyy, Gyz, Gzz, wgt;
        if constexpr (GMode == 0)
        {
            // per-point: weight already folded into the stored G.
            Gxx = Gbase[p*6+0]; Gxy = Gbase[p*6+1]; Gxz = Gbase[p*6+2];
            Gyy = Gbase[p*6+3]; Gyz = Gbase[p*6+4]; Gzz = Gbase[p*6+5];
            wgt = 1.0;
        }
        else
        {
            // affine: element-constant Ghat, GLL weight applied here.
            Gxx = Gbase[0]; Gxy = Gbase[1]; Gxz = Gbase[2];
            Gyy = Gbase[3]; Gyz = Gbase[4]; Gzz = Gbase[5];
            const int i = p / nn, j = (p / n) % n, k = p % n;
            wgt = d_w[i] * d_w[j] * d_w[k];
        }
        const double ar = w.b0[p], as = w.b1[p], at = w.b2[p];
        w.b0[p] = wgt * (Gxx * ar + Gxy * as + Gxz * at);  // gr
        w.b1[p] = wgt * (Gxy * ar + Gyy * as + Gyz * at);  // gs
        w.b2[p] = wgt * (Gxz * ar + Gyz * as + Gzz * at);  // gt
    }
    __syncwarp();

    // Zero the accumulator before the three additive integrate-back contributions.
    for (int l = laneId; l < N3; l += 32) w.b_y[l] = 0.0;
    __syncwarp();

    // =========================================================================
    // Integrate back: y = D^T.gr (dir0) + D^T.gs (dir1) + D^T.gt (dir2).
    // Identical tile structure to the forward pass, with the matrix operand
    // replaced by D^T and the data operand being gr/gs/gt; results accumulate
    // into b_y. For dir2 the Dt tile reads D^T transposed: Dt[s][c] = DT[c*8 + s].
    // =========================================================================

    // dir0 with D^T: C[8x64] = D^T[8x8] . gr[8x64].
    for (int bctile = 0; bctile < n; ++bctile)
    {
        auto fillA = [&](int kt) {
            const int a = laneId % n, k = laneId / n;
            w.smem_A[k * n + a] = w.D[(kt * 4 + k) * n + a];   // D^T on the fly
        };
        auto fillB = [&](int kt) {
            const int k = laneId % 4, col = laneId / 4;
            w.smem_B[col * 4 + k] = w.b0[(kt * 4 + k) * nn + bctile * n + col];
        };
        tile_8x8x8<TileMode>(w, fillA, fillB);
        for (int idx = laneId; idx < nn; idx += 32)
        {
            const int a = idx / n, col = idx % n;
            w.b_y[a * nn + bctile * n + col] += w.smem_C[col * n + a];
        }
        __syncwarp();
    }

    // dir1 with D^T: per slab a, C_a[8x8] = D^T[8x8] . gs_a[8x8].
    for (int a = 0; a < n; ++a)
    {
        auto fillA = [&](int kt) {
            const int b = laneId % n, k = laneId / n;
            w.smem_A[k * n + b] = w.D[(kt * 4 + k) * n + b];   // D^T on the fly
        };
        auto fillB = [&](int kt) {
            const int k = laneId % 4, c = laneId / 4;
            w.smem_B[c * 4 + k] = w.b1[a * nn + (kt * 4 + k) * n + c];
        };
        tile_8x8x8<TileMode>(w, fillA, fillB);
        for (int idx = laneId; idx < nn; idx += 32)
        {
            const int b = idx / n, c = idx % n;
            w.b_y[a * nn + b * n + c] += w.smem_C[c * n + b];
        }
        __syncwarp();
    }

    // dir2 with D^T: C[64x8] = gt[64x8] . Dt[8x8], Dt[s][c] = DT[c*8 + s].
    for (int abtile = 0; abtile < n; ++abtile)
    {
        auto fillA = [&](int kt) {
            const int ar = laneId % n, k = laneId / n;
            w.smem_A[k * n + ar] = w.b2[(abtile * n + ar) * n + (kt * 4 + k)];
        };
        auto fillB = [&](int kt) {
            const int k = laneId % 4, c = laneId / 4;
            w.smem_B[c * 4 + k] = w.D[(kt * 4 + k) * n + c];   // D^T on the fly
        };
        tile_8x8x8<TileMode>(w, fillA, fillB);
        for (int idx = laneId; idx < nn; idx += 32)
        {
            const int i = idx / n, j = idx % n;
            const int ab = abtile * n + i;
            w.b_y[ab * n + j] += w.smem_C[j * n + i];
        }
        __syncwarp();
    }

    // --- Scatter (direct stiffness summation, additive into global DOFs) ---
    for (int l = laneId; l < N3; l += 32)
        atomicAdd(&d_y_dof[edof[l]], w.b_y[l]);

#else  // __CUDA_ARCH__ < 800
    // FP64 WMMA m8n8k4 requires SM80+ (A100 / H100 / GH200). Build with
    // -DCMAKE_CUDA_ARCHITECTURES=80 (or higher) to enable this kernel.
    (void)d_u_dof; (void)d_y_dof; (void)d_elemDof;
    (void)d_D; (void)d_G; (void)d_w; (void)numElements;
#endif
}

// Host launcher: opts into >48 KB dynamic smem (once, lazily) and launches with
// the correct grid / dynamic-smem byte count. d_y_dof MUST be zeroed by the
// caller first (scatter is additive). Returns the launch error so the caller can
// check it; does NOT synchronize.
// Shared launch body; GMode picks the metric mode, TileMode the contraction
// hardware (0 = WMMA tensor core, 1 = scalar twin with identical structure).
template<typename RealType, int P, int BlockSize, int GMode, int TileMode>
inline cudaError_t ho_laplacian_dmma_launch_impl(const RealType* d_u_dof,
                                                 RealType* d_y_dof,
                                                 const int* d_elemDof,
                                                 const RealType* d_D,
                                                 const RealType* d_G,
                                                 const RealType* d_w,
                                                 size_t numElements,
                                                 cudaStream_t stream)
{
    constexpr int  WarpsPerBlock = BlockSize / 32;
    constexpr int  smemBytes     = WarpsPerBlock * (int)sizeof(HoLaplacianDmmaWarpData);
    auto kernel = ho_laplacian_apply_cg_dmma<RealType, P, BlockSize, GMode, TileMode>;

    // The 88 KB/block we need is above the 48 KB static cap, so opt in once.
    // Idempotent and cheap; safe to call before every launch.
    cudaError_t attrErr = cudaFuncSetAttribute(
        kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, smemBytes);
    if (attrErr != cudaSuccess) return attrErr;

    const size_t grid = (numElements + WarpsPerBlock - 1) / WarpsPerBlock;
    kernel<<<(unsigned)grid, BlockSize, smemBytes, stream>>>(
        d_u_dof, d_y_dof, d_elemDof, d_D, d_G, d_w, numElements);
    return cudaGetLastError();
}

// Per-quad-point metric (general/curved). d_G is [numElements * N3 * 6].
template<typename RealType, int P, int BlockSize = 128>
inline cudaError_t ho_laplacian_dmma_launch(const RealType* d_u_dof,
                                            RealType* d_y_dof,
                                            const int* d_elemDof,
                                            const RealType* d_D,
                                            const RealType* d_G,
                                            size_t numElements,
                                            cudaStream_t stream = 0)
{
    return ho_laplacian_dmma_launch_impl<RealType, P, BlockSize, 0, 0>(
        d_u_dof, d_y_dof, d_elemDof, d_D, d_G, nullptr, numElements, stream);
}

// Affine metric (straight-sided hexes). d_Ghat is [numElements * 6] (constant
// geometric metric per element), d_w is the n GLL quadrature weights. ~512x less
// metric traffic than the per-point path; numerically identical for affine hexes.
template<typename RealType, int P, int BlockSize = 128>
inline cudaError_t ho_laplacian_dmma_launch_affine(const RealType* d_u_dof,
                                                   RealType* d_y_dof,
                                                   const int* d_elemDof,
                                                   const RealType* d_D,
                                                   const RealType* d_Ghat,
                                                   const RealType* d_w,
                                                   size_t numElements,
                                                   cudaStream_t stream = 0)
{
    return ho_laplacian_dmma_launch_impl<RealType, P, BlockSize, 1, 0>(
        d_u_dof, d_y_dof, d_elemDof, d_D, d_Ghat, d_w, numElements, stream);
}

// Matched-structure SCALAR baseline: identical warp-per-element + smem kernel as
// the affine DMMA path, but the 8x8x8 contraction runs on CUDA cores instead of
// tensor cores (TileMode=1). Differs from ho_laplacian_dmma_launch_affine in
// exactly one thing -- the tile op -- so (affine DMMA time)/(this time) isolates
// the tensor-core speedup, free of the memory-structure advantage that the
// thread-per-element ho_laplacian_apply_cg reference also carries.
template<typename RealType, int P, int BlockSize = 128>
inline cudaError_t ho_laplacian_smem_scalar_launch_affine(const RealType* d_u_dof,
                                                          RealType* d_y_dof,
                                                          const int* d_elemDof,
                                                          const RealType* d_D,
                                                          const RealType* d_Ghat,
                                                          const RealType* d_w,
                                                          size_t numElements,
                                                          cudaStream_t stream = 0)
{
    return ho_laplacian_dmma_launch_impl<RealType, P, BlockSize, 1, 1>(
        d_u_dof, d_y_dof, d_elemDof, d_D, d_Ghat, d_w, numElements, stream);
}

} // namespace fem
} // namespace mars
