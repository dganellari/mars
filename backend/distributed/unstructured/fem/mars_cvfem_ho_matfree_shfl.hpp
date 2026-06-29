#pragma once
// Register + warp-shuffle CVFEM PA apply  (MARS_HO_SHFL path, "the real engine").
//
// WHY: ncu shows the warp-per-element store-d_G apply is bound by the shared-memory INSTRUCTION
// pipe (LSU / Mem-Pipes-Busy ~99% at p=4 and p=7) while the FMA pipe is idle (SM ~30%). An earlier
// cuBLAS batched-GEMM probe (git history) confirmed the *formulation* recasts bit-exact but that
// relocating the contractions to GLOBAL memory only trades LSU pressure for DRAM -- it did NOT help.
// The fix that actually empties the LSU is the Nek/hipBone design: keep the per-element face data
// RESIDENT IN REGISTERS and exchange it across the element's threads with WARP SHUFFLES (__shfl) --
// the shuffle unit is NOT the LSU pipe. Measured p=4: LSU 98.7% -> 50.2%, 6.0 -> 10.9 GDOF/s (1.8x);
// p=1 +20%, p=3 +22%, all A.1 bit-identical.
//
// SCOPE: orders whose face fits one warp (NN <= 32 -> p<=4) AND whose warp-aligned lane group
// tpe = nextPow2(NN) doesn't waste too many lanes (so NOT p=2: NN=9 rounds to tpe=16 = 44% dead
// lanes, which costs more than the saved face round-trips -- measured -3% vs the baseline's
// perfectly-packed tpe=9, so p=2 falls back). p=1/3/4 route here; p=2 and p>=5 fall back to the
// reference warp kernel, so the full p=1..8 sweep still runs. my_u stays in shared (the step-1
// normal reads); only the FACE round-trips (2 tangential-D + 2 W contractions per (dir,l) -- the
// bulk of the shared traffic) move onto shuffles.

#include "mars_cvfem_ho_matfree.hpp"
#include <cstdlib>

namespace mars {
namespace fem {

constexpr int ho_next_pow2(int x) { int p = 1; while (p < x) p <<= 1; return p; }

// One warp-aligned lane group per element; lane == the (s,r) face slot for lane < NN. The step-1
// outputs interp(bi)/deriv(di) stay in registers; the tangential-D and both W contractions read
// other lanes' values via __shfl_sync(width=tpe) instead of a shared face buffer. Inactive lanes
// (lane >= NN) compute clamped garbage that is never scattered, but MUST take every __shfl_sync
// (group-collective) -- so the shuffles are unconditional and only the my_y scatter is guarded.
template<typename RealType, int P, int BlockSize, int ElemsPerBlock, typename GStore = double, bool Recompute = false>
__global__ void __launch_bounds__(BlockSize)
ho_cvfem_apply_kernel_shfl(const RealType* __restrict__ d_u,
                           RealType* __restrict__ d_y,
                           const int* __restrict__ d_elemDof,
                           const GStore* __restrict__ d_G,
                           size_t numElements,
                           const int* __restrict__ d_elemList = nullptr,
                           size_t count = 0,
                           const RealType* __restrict__ d_corners = nullptr)
{
    constexpr int N = P + 1, NN = N * N, N3 = NN * N, E = ElemsPerBlock;
    constexpr int tpe = BlockSize / E;
    static_assert(tpe <= 32 && (tpe & (tpe - 1)) == 0, "shfl apply: tpe must be a power of 2 <= 32");
    static_assert(NN <= tpe, "shfl apply: face must fit the element's lane group");
    static_assert(std::is_same<RealType, double>::value, "shfl apply: double only");

    // smem: [u_sh | y_sh] only -- no face buffers (registers + shuffles replace them).
    extern __shared__ char ho_smem_raw[];
    RealType* u_sh = reinterpret_cast<RealType*>(ho_smem_raw);
    RealType* y_sh = u_sh + (size_t)E * N3;

    const int t = threadIdx.x;
    const int localElem = t / tpe;
    const int lane      = t % tpe;            // == the (s,r) slot index for lane < NN
    const size_t slot_e = (size_t)blockIdx.x * E + localElem;
    const size_t e = (d_elemList != nullptr)
                       ? (slot_e < count ? (size_t)d_elemList[slot_e] : numElements)
                       : slot_e;
    RealType* my_u = u_sh + localElem * N3;
    RealType* my_y = y_sh + localElem * N3;
    RealType* my_cn = nullptr;
    if constexpr (Recompute) my_cn = y_sh + (size_t)E * N3 + localElem * 24;   // corners slab (E*24 doubles)
    const bool valid = (e < numElements);
    const int* edof  = valid ? (d_elemDof + e * N3) : nullptr;

    for (int l = lane; l < N3; l += tpe) {
        my_u[l] = (valid && edof[l] >= 0) ? d_u[edof[l]] : RealType(0);
        my_y[l] = RealType(0);
    }
    if constexpr (Recompute)        // 8 corners/element into smem, shared by the element's lanes
        for (int c = lane; c < 24; c += tpe) my_cn[c] = valid ? d_corners[e * 24 + c] : RealType(0);
    __syncwarp();

    // store-d_G base ptr -- under !Recompute compiles to the original (d_G + offset), no extra check;
    // stays null + unused under Recompute (the d_G arithmetic lives only in the !Recompute branch).
    const GStore* Gbase = nullptr;
    if constexpr (!Recompute) Gbase = d_G + (valid ? e : 0) * (size_t)(3 * P * NN * 3);
    // Clamp inactive lanes to slot (0,0): keeps their my_u / Gbase reads in-bounds (results are
    // garbage but never scattered) and lets them join every shuffle.
    const int s = (lane < NN) ? lane / N : 0;
    const int r = (lane < NN) ? lane % N : 0;
    const unsigned mask = 0xFFFFFFFFu;        // whole warp converged; width=tpe subdivides per element

    for (int dir = 0; dir < 3; ++dir) {
        for (int l = 0; l < P; ++l) {
            // step 1: normal interp + deriv -> registers (no shared face write)
            RealType bi = 0, di = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) {
                RealType uq = my_u[ho_idx<N, NN>(dir, q, s, r)];
                bi += c_Btil[l * N + q] * uq;
                di += c_Dtil[l * N + q] * uq;
            }
            // step 2: tangential D of interp via shuffles of bi (within the element group); flux -> reg
            RealType dt2 = 0, dt1 = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) dt2 += c_D[r * N + q] * __shfl_sync(mask, bi, s * N + q, tpe);
            #pragma unroll
            for (int q = 0; q < N; ++q) dt1 += c_D[s * N + q] * __shfl_sync(mask, bi, q * N + r, tpe);
            RealType fb;
            if constexpr (Recompute) {              // regenerate the metric inline from the element corners
                double gg[3];
                ho_cvfem_metric_point<P>(reinterpret_cast<const double(*)[3]>(my_cn), dir, l, s, r, gg);
                fb = RealType(gg[2]) * di + RealType(gg[0]) * dt2 + RealType(gg[1]) * dt1;
            } else {                                // store-d_G: byte-identical to the pre-MF-shuffle expression
                const GStore* g = Gbase + (size_t)(((dir * P + l) * N + s) * N + r) * 3;
                fb = RealType(g[2]) * di + RealType(g[0]) * dt2 + RealType(g[1]) * dt1;
            }
            // step 3: W along r via shuffle of flux -> register
            RealType v = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) v += c_W[r * N + q] * __shfl_sync(mask, fb, s * N + q, tpe);
            // step 4: W along s via shuffle of v -> scatter to the two SCS-bounding nodes
            RealType intf = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) intf += c_W[s * N + q] * __shfl_sync(mask, v, q * N + r, tpe);
            if (lane < NN) {
                my_y[ho_idx<N, NN>(dir, l,     s, r)] -= intf;
                my_y[ho_idx<N, NN>(dir, l + 1, s, r)] += intf;
            }
            __syncwarp();   // serialize the my_y accumulation across dir (column ownership is per-dir)
        }
    }

    if (e < numElements)
        for (int l = lane; l < N3; l += tpe) {
            int dof = edof[l];
            if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);
        }
}

// "Nek-style" register-line variant of the single-warp shuffle (MARS_HO_NEK): the normal column
// my_u[ho_idx(dir,q,s,r)] (q=0..N-1) is IDENTICAL for every l in the inner sweep, so it is HOISTED
// out of the l-loop into registers (regCol). Step 1 then reads registers instead of shared, cutting
// the per-(dir,l) shared my_u reads by a factor of P. Everything else mirrors the shuffle kernel
// exactly -- the math is unchanged, so A.1 stays bit-identical; this only measures whether moving the
// last shared reads (step 1) onto registers, Nek-style, buys anything on top of the shuffle.
template<typename RealType, int P, int BlockSize, int ElemsPerBlock, typename GStore = double>
__global__ void __launch_bounds__(BlockSize)
ho_cvfem_apply_kernel_shfl_nek(const RealType* __restrict__ d_u,
                               RealType* __restrict__ d_y,
                               const int* __restrict__ d_elemDof,
                               const GStore* __restrict__ d_G,
                               size_t numElements,
                               const int* __restrict__ d_elemList = nullptr,
                               size_t count = 0)
{
    constexpr int N = P + 1, NN = N * N, N3 = NN * N, E = ElemsPerBlock;
    constexpr int tpe = BlockSize / E;
    static_assert(tpe <= 32 && (tpe & (tpe - 1)) == 0, "nek apply: tpe must be a power of 2 <= 32");
    static_assert(NN <= tpe, "nek apply: face must fit the element's lane group");
    static_assert(std::is_same<RealType, double>::value, "nek apply: double only");

    extern __shared__ char ho_smem_raw[];
    RealType* u_sh = reinterpret_cast<RealType*>(ho_smem_raw);
    RealType* y_sh = u_sh + (size_t)E * N3;

    const int t = threadIdx.x;
    const int localElem = t / tpe;
    const int lane      = t % tpe;
    const size_t slot_e = (size_t)blockIdx.x * E + localElem;
    const size_t e = (d_elemList != nullptr)
                       ? (slot_e < count ? (size_t)d_elemList[slot_e] : numElements)
                       : slot_e;
    RealType* my_u = u_sh + localElem * N3;
    RealType* my_y = y_sh + localElem * N3;
    const bool valid = (e < numElements);
    const int* edof  = valid ? (d_elemDof + e * N3) : nullptr;

    for (int l = lane; l < N3; l += tpe) {
        my_u[l] = (valid && edof[l] >= 0) ? d_u[edof[l]] : RealType(0);
        my_y[l] = RealType(0);
    }
    __syncwarp();

    const GStore* Gbase = d_G + (valid ? e : 0) * (size_t)(3 * P * NN * 3);
    const int s = (lane < NN) ? lane / N : 0;
    const int r = (lane < NN) ? lane % N : 0;
    const unsigned mask = 0xFFFFFFFFu;

    for (int dir = 0; dir < 3; ++dir) {
        // HOIST: the normal column is l-independent -> read it once per dir into registers.
        RealType regCol[N];
        #pragma unroll
        for (int q = 0; q < N; ++q) regCol[q] = my_u[ho_idx<N, NN>(dir, q, s, r)];
        for (int l = 0; l < P; ++l) {
            RealType bi = 0, di = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) {
                bi += c_Btil[l * N + q] * regCol[q];
                di += c_Dtil[l * N + q] * regCol[q];
            }
            RealType dt2 = 0, dt1 = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) dt2 += c_D[r * N + q] * __shfl_sync(mask, bi, s * N + q, tpe);
            #pragma unroll
            for (int q = 0; q < N; ++q) dt1 += c_D[s * N + q] * __shfl_sync(mask, bi, q * N + r, tpe);
            const GStore* g = Gbase + (size_t)(((dir * P + l) * N + s) * N + r) * 3;
            const RealType fb = RealType(g[2]) * di + RealType(g[0]) * dt2 + RealType(g[1]) * dt1;
            RealType v = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) v += c_W[r * N + q] * __shfl_sync(mask, fb, s * N + q, tpe);
            RealType intf = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) intf += c_W[s * N + q] * __shfl_sync(mask, v, q * N + r, tpe);
            if (lane < NN) {
                my_y[ho_idx<N, NN>(dir, l,     s, r)] -= intf;
                my_y[ho_idx<N, NN>(dir, l + 1, s, r)] += intf;
            }
            __syncwarp();
        }
    }

    if (e < numElements)
        for (int l = lane; l < N3; l += tpe) {
            int dof = edof[l];
            if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);
        }
}

// Multi-warp HYBRID for the orders whose face spans >1 warp (NN > 32 -> p>=5). Following Nek's
// lesson: the cross-warp direction goes through SHARED (spans warps), the intra-warp direction stays
// on shuffles. Rows are PADDED to rowWidth = nextPow2(N) so each row is a power-of-2 lane block that
// stays inside one 32-lane warp (p=7's N=8 needs no padding; p=5/6 pad N=6/7 to 8-wide rows; p=8 pads
// N=9 to 16). The tangential-r and W-r contractions read within a row (intra-warp) -> __shfl; the
// tangential-s and W-s contractions read down a column (crosses warps) -> a thin shared buffer faceA
// (reused for bi then v). Cuts ~half the face shared reads vs the baseline. my_u stays shared. The
// padding wastes lanes (p=6 23%, p=5 44%, p=8 49%) -- whether it pays is measured per order.
template<typename RealType, int P, int BlockSize, int ElemsPerBlock, typename GStore = double, bool Hoist = false, int MinBlocks = 0, bool Recompute = false>
__global__ void __launch_bounds__(BlockSize, MinBlocks)
ho_cvfem_apply_kernel_shfl_mw(const RealType* __restrict__ d_u,
                              RealType* __restrict__ d_y,
                              const int* __restrict__ d_elemDof,
                              const GStore* __restrict__ d_G,
                              size_t numElements,
                              const int* __restrict__ d_elemList = nullptr,
                              size_t count = 0,
                              const RealType* __restrict__ d_corners = nullptr)
{
    constexpr int N = P + 1, NN = N * N, N3 = NN * N, E = ElemsPerBlock;
    constexpr int rowWidth = ho_next_pow2(N);          // padded row width: power of 2, divides 32
    constexpr int tpe = BlockSize / E;
    constexpr int rpw = 32 / rowWidth;                 // rows per warp
    static_assert(tpe == ((N * rowWidth + 31) / 32) * 32, "mw: tpe = N padded rows rounded to whole warps");
    static_assert(32 % rowWidth == 0, "mw: padded row width divides 32");
    static_assert(std::is_same<RealType, double>::value, "mw: double only");

    extern __shared__ char ho_smem_raw[];
    RealType* u_sh  = reinterpret_cast<RealType*>(ho_smem_raw);
    RealType* y_sh  = u_sh + (size_t)E * N3;
    RealType* fa_sh = y_sh + (size_t)E * N3;  // faceA column staging: tpe per element (padded lanes)

    const int t = threadIdx.x;
    const int localElem = t / tpe;
    const int lane      = t % tpe;
    const size_t slot_e = (size_t)blockIdx.x * E + localElem;
    const size_t e = (d_elemList != nullptr)
                       ? (slot_e < count ? (size_t)d_elemList[slot_e] : numElements)
                       : slot_e;
    RealType* my_u  = u_sh  + localElem * N3;
    RealType* my_y  = y_sh  + localElem * N3;
    RealType* faceA = fa_sh + localElem * tpe;
    RealType* my_cn = nullptr;
    if constexpr (Recompute) my_cn = fa_sh + (size_t)E * tpe + localElem * 24;   // corners slab (E*24)
    const bool valid = (e < numElements);
    const int* edof  = valid ? (d_elemDof + e * N3) : nullptr;

    for (int l = lane; l < N3; l += tpe) {
        my_u[l] = (valid && edof[l] >= 0) ? d_u[edof[l]] : RealType(0);
        my_y[l] = RealType(0);
    }
    if constexpr (Recompute)
        for (int c = lane; c < 24; c += tpe) my_cn[c] = valid ? d_corners[e * 24 + c] : RealType(0);
    __syncthreads();

    const GStore* Gbase = nullptr;     // store-d_G: original code under !Recompute (byte-identical)
    if constexpr (!Recompute) Gbase = d_G + (valid ? e : 0) * (size_t)(3 * P * NN * 3);
    // Padded lane layout: lane = sp*rowWidth + rp; a slot is active only where sp<N and rp<N.
    // Padding lanes compute clamped garbage and never scatter, but join every shuffle (collective).
    const int sp = lane / rowWidth, rp = lane % rowWidth;
    const bool active = (sp < N && rp < N);
    const int s = active ? sp : 0, r = active ? rp : 0;
    const int rowLane = (sp % rpw) * rowWidth;         // base within-warp lane of this lane's row
    const unsigned mask = 0xFFFFFFFFu;

    for (int dir = 0; dir < 3; ++dir) {
        // Hoist (MARS_HO_NEK): the normal column is l-independent -> read it once per dir into regs.
        RealType regCol[Hoist ? N : 1];
        if constexpr (Hoist) {
            #pragma unroll
            for (int q = 0; q < N; ++q) regCol[q] = my_u[ho_idx<N, NN>(dir, q, s, r)];
        }
        for (int l = 0; l < P; ++l) {
            RealType bi = 0, di = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) {
                RealType uq;
                if constexpr (Hoist) uq = regCol[q];
                else                 uq = my_u[ho_idx<N, NN>(dir, q, s, r)];
                bi += c_Btil[l * N + q] * uq;
                di += c_Dtil[l * N + q] * uq;
            }
            faceA[lane] = bi;
            __syncthreads();                  // bi visible for the cross-warp column read
            RealType dt2 = 0, dt1 = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) dt2 += c_D[r * N + q] * __shfl_sync(mask, bi, rowLane + q); // row
            #pragma unroll
            for (int q = 0; q < N; ++q) dt1 += c_D[s * N + q] * faceA[q * rowWidth + r];                   // col
            RealType fb;
            if constexpr (Recompute) {              // regenerate the metric inline from the element corners
                double gg[3];
                ho_cvfem_metric_point<P>(reinterpret_cast<const double(*)[3]>(my_cn), dir, l, s, r, gg);
                fb = RealType(gg[2]) * di + RealType(gg[0]) * dt2 + RealType(gg[1]) * dt1;
            } else {                                // store-d_G: byte-identical to the pre-MF-shuffle expression
                const GStore* g = Gbase + (size_t)(((dir * P + l) * N + s) * N + r) * 3;
                fb = RealType(g[2]) * di + RealType(g[0]) * dt2 + RealType(g[1]) * dt1;
            }
            RealType v = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) v += c_W[r * N + q] * __shfl_sync(mask, fb, rowLane + q);    // row
            __syncthreads();                  // finish all dt1 reads of bi before overwriting faceA
            faceA[lane] = v;
            __syncthreads();                  // v visible for the column read
            RealType intf = 0;
            #pragma unroll
            for (int q = 0; q < N; ++q) intf += c_W[s * N + q] * faceA[q * rowWidth + r];                   // col
            if (active) {
                my_y[ho_idx<N, NN>(dir, l,     s, r)] -= intf;
                my_y[ho_idx<N, NN>(dir, l + 1, s, r)] += intf;
            }
            __syncthreads();                  // faceA reused by next l
        }
    }

    if (e < numElements)
        for (int l = lane; l < N3; l += tpe) {
            int dof = edof[l];
            if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);
        }
}

// Launcher: NN<=32 (p=1,2,3,4) -> single-warp shuffle (tpe = nextPow2(NN)); NN>32 (p=5,6,7,8) -> the
// padded multi-warp hybrid above. Every order p=1..8 now runs on a shuffle kernel (p=2 is kept despite
// its 44% dead lanes for a uniform sweep). Drop-in for ho_cvfem_apply_launch; selected by MARS_HO_SHFL.
template<typename RealType, int P, typename GStore = double, bool Recompute = false,
         int BlockSize = HoCvfemLaunchDefault<P>::Block,
         int ElemsPerBlock = HoCvfemLaunchDefault<P>::Elems>
inline cudaError_t ho_cvfem_apply_launch_shfl(const RealType* d_u, RealType* d_y,
                                              const int* d_elemDof, const GStore* d_G, size_t numElements,
                                              cudaStream_t stream = 0,
                                              const int* d_elemList = nullptr, size_t count = 0,
                                              const RealType* d_corners = nullptr)
{
    constexpr int N = P + 1, NN = N * N, N3 = NN * N;
    if constexpr (NN <= 32 &&
                  std::is_same<RealType, double>::value && std::is_same<GStore, double>::value) {
        constexpr int tpe   = ho_next_pow2(NN);   // 4,16,16,32 for p=1,2,3,4 (p=2 = 44% dead lanes)
        constexpr int Block = 256;
        constexpr int Elems = Block / tpe;
        const size_t smemBytes = (size_t)(Elems * 2 * N3 + (Recompute ? Elems * 24 : 0)) * sizeof(RealType);
        auto kern = ho_cvfem_apply_kernel_shfl<RealType, P, Block, Elems, GStore, Recompute>;
        cudaFuncSetAttribute(kern, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smemBytes);
        const size_t nBlocks = (d_elemList != nullptr ? (count + Elems - 1) / Elems
                                                      : (numElements + Elems - 1) / Elems);
        kern<<<(unsigned)nBlocks, Block, smemBytes, stream>>>(
            d_u, d_y, d_elemDof, d_G, numElements, d_elemList, count, d_corners);
        return cudaGetLastError();
    } else if constexpr (NN > 32 &&
                         std::is_same<RealType, double>::value && std::is_same<GStore, double>::value) {
        constexpr int rowWidth  = ho_next_pow2(N);
        constexpr int elemLanes = N * rowWidth;
        constexpr int tpe       = ((elemLanes + 31) / 32) * 32;        // whole warps: p=5/6/7 -> 64, p=8 -> 160
        constexpr int Elems     = (256 / tpe > 0) ? (256 / tpe) : 1;   // p=5/6/7 -> 4, p=8 -> 1
        constexpr int Block     = Elems * tpe;
        const size_t smemBytes = (size_t)(Elems * (2 * N3 + tpe) + (Recompute ? Elems * 24 : 0)) * sizeof(RealType);
        const size_t nBlocks = (d_elemList != nullptr ? (count + Elems - 1) / Elems
                                                      : (numElements + Elems - 1) / Elems);
        // Per-order occupancy cap: min-blocks=3 caps registers -> 3 blocks/SM (37.5%), MEASURED
        // +14% on p=5/6 and +1% on p=8. p=7 (122 regs) spills under the cap and loses 17%, so it
        // keeps the default (min-blocks=0, 2 blocks/SM). All shuffle, one kernel -- just a constant.
        // (Recompute adds the inline metric's registers; mb is unchanged here -- re-tune on-device if it spills.)
        constexpr int mb = (P == 7) ? 0 : 3;
        auto kern = ho_cvfem_apply_kernel_shfl_mw<RealType, P, Block, Elems, GStore, false, mb, Recompute>;
        cudaFuncSetAttribute(kern, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smemBytes);
        kern<<<(unsigned)nBlocks, Block, smemBytes, stream>>>(
            d_u, d_y, d_elemDof, d_G, numElements, d_elemList, count, d_corners);
        return cudaGetLastError();
    } else {
        return ho_cvfem_apply_launch<RealType, P, GStore>(
            d_u, d_y, d_elemDof, d_G, numElements, stream, d_elemList, count);
    }
}

// Nek launcher: register-line (hoisted normal column) for the FULL sweep -- NN<=32 (p<=4) uses the
// hoisted single-warp kernel, NN>32 (p>=5) uses the padded multi-warp kernel with Hoist=true. Selected
// by MARS_HO_NEK; the shuffle path (Hoist=false default) is byte-identical, so MARS_HO_SHFL is a clean A/B.
template<typename RealType, int P, typename GStore = double,
         int BlockSize = HoCvfemLaunchDefault<P>::Block,
         int ElemsPerBlock = HoCvfemLaunchDefault<P>::Elems>
inline cudaError_t ho_cvfem_apply_launch_nek(const RealType* d_u, RealType* d_y,
                                             const int* d_elemDof, const GStore* d_G, size_t numElements,
                                             cudaStream_t stream = 0,
                                             const int* d_elemList = nullptr, size_t count = 0)
{
    constexpr int N = P + 1, NN = N * N, N3 = NN * N;
    if constexpr (NN <= 32 &&
                  std::is_same<RealType, double>::value && std::is_same<GStore, double>::value) {
        constexpr int tpe   = ho_next_pow2(NN);
        constexpr int Block = 256;
        constexpr int Elems = Block / tpe;
        const size_t smemBytes = (size_t)(Elems * 2 * N3) * sizeof(RealType);
        auto kern = ho_cvfem_apply_kernel_shfl_nek<RealType, P, Block, Elems, GStore>;
        cudaFuncSetAttribute(kern, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smemBytes);
        const size_t nBlocks = (d_elemList != nullptr ? (count + Elems - 1) / Elems
                                                      : (numElements + Elems - 1) / Elems);
        kern<<<(unsigned)nBlocks, Block, smemBytes, stream>>>(
            d_u, d_y, d_elemDof, d_G, numElements, d_elemList, count);
        return cudaGetLastError();
    } else if constexpr (NN > 32 &&
                         std::is_same<RealType, double>::value && std::is_same<GStore, double>::value) {
        // multi-warp register-line: the padded mw kernel with Hoist=true (the full-sweep nek path).
        constexpr int rowWidth  = ho_next_pow2(N);
        constexpr int elemLanes = N * rowWidth;
        constexpr int tpe       = ((elemLanes + 31) / 32) * 32;
        constexpr int Elems     = (256 / tpe > 0) ? (256 / tpe) : 1;
        constexpr int Block     = Elems * tpe;
        const size_t smemBytes = (size_t)(Elems * (2 * N3 + tpe)) * sizeof(RealType);
        auto kern = ho_cvfem_apply_kernel_shfl_mw<RealType, P, Block, Elems, GStore, true>;
        cudaFuncSetAttribute(kern, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smemBytes);
        const size_t nBlocks = (d_elemList != nullptr ? (count + Elems - 1) / Elems
                                                      : (numElements + Elems - 1) / Elems);
        kern<<<(unsigned)nBlocks, Block, smemBytes, stream>>>(
            d_u, d_y, d_elemDof, d_G, numElements, d_elemList, count, nullptr);   // d_corners unused (store-d_G nek)
        return cudaGetLastError();
    } else {
        return ho_cvfem_apply_launch_shfl<RealType, P, GStore>(
            d_u, d_y, d_elemDof, d_G, numElements, stream, d_elemList, count);
    }
}

} // namespace fem
} // namespace mars
