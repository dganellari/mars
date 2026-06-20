#pragma once

// Optimized GPU matrix-free apply for the high-order CVFEM diffusion operator
// (Knaus Alg 2). This is the device twin of the host reference
// applyHoCvfemElement (mars_cvfem_ho_apply.hpp): same elemDof gather/scatter,
// same per-element metric G layout, same per-direction sum-factorized sweep over
// subcontrol-surface (SCS) faces. The only allowed numerical difference is FP
// reduction order, so the two agree to ~1e-13 relative.
//
// WHY this design (the operator class is shared-memory-bandwidth + occupancy
// bound, NOT compute bound -- ncu prior from the sibling DMMA Galerkin kernel:
// 90KB smem/block capped it to 2 blocks/SM = 12.5% occupancy). The levers:
//
//   1. SMALL reference operators (Btil, Dtil, D, W) live in __constant__ memory,
//      NOT shared -- zero per-block smem for them. They are tiny ((p+1)^2-ish)
//      and broadcast-cached. This is the single biggest occupancy win.
//   2. Alg 2 never needs a full n^3 derivative TENSOR in smem. Each SCS face is
//      only nn=(p+1)^2 work. We map one CUDA thread per tangential (s,r) face
//      slot; that thread owns its entire normal column, so the normal-direction
//      contractions (Btil/Dtil) and the -/+ scatter to nodes l,l+1 need NO
//      cross-thread exchange and NO __syncthreads. Only the tangential D and W
//      contractions exchange data, via two small nn smem face buffers.
//   3. Multi-element tiling at low p (ElemsPerBlock): a p=1 element is 8 nodes /
//      4 face slots -- one block/element starves the SM. Pack E elements/block.
//   4. Per-block smem = E*(2*n3 for u_sh+y_sh) + E*(2 face buffers, each padded
//      nn -> N*(N+1) to break shared-memory bank conflicts -- see line 299). At p=7
//      this is ~8.7KB vs the DMMA 90KB -> >=4x occupancy headroom.
//   5. GMode { PerPoint, Affine }: PerPoint reads the full host metric layout
//      (curved/sheared hexes); Affine reads a per-element constant triple plus a
//      constant-mem quadrature weight table -> far less metric DRAM traffic.
//
// Index/layout contract (ported VERBATIM from applyHoCvfemElement):
//   idx(dir,nrm,t1,t2): dir0 -> nrm*nn + t1*n + t2
//                       dir1 -> t1*nn + nrm*n + t2
//                       dir2 -> t1*nn + t2*n + nrm
//   G layout: g = d_G[((dir*p + l)*n + s)*n + r], g[2]=normal, g[1]=tang1(s),
//             g[0]=tang2(r). metric = detJ * J^{-1} J^{-T} (NOT (JJ^T)^{-1}).
//
// Launch contract:
//   - d_y zeroed before launch (scatter is additive via atomicAdd).
//   - d_elemDof is [numElements * N3] int, l = i*nn + j*n + k (i slowest), built
//     host-side by HODofHandler and uploaded to device by the caller.
//   - reference operators uploaded once via ho_cvfem_upload_operators().
//
// LINKAGE: the __constant__ banks (c_Btil etc., c_hexCornerRef) are defined in
// this header. __constant__ cannot be `inline` for ODR, so this header must be
// included by EXACTLY ONE translation unit. Including it from a second .cu will
// produce duplicate-symbol link errors. If multi-TU use is ever needed, move the
// __constant__ definitions + ho_cvfem_upload_operators into one .cu and expose
// only declarations.

#include <cuda_runtime.h>
#include <cassert>
#include <cstddef>
#include <type_traits>

namespace mars {
namespace fem {

// --------------------------------------------------------------------------
// Reference operators in __constant__ memory. WHY constant (not smem): they are
// read-only, tiny, and identical for every element -> the constant cache
// broadcasts a single value to the whole warp in one transaction, and they cost
// ZERO shared memory, which is the binding occupancy resource. One flat bank per
// operator covers all P up to HO_CVFEM_MAX_P; a per-P offset is not needed
// because a single P is instantiated per kernel and we index [0, len(P)).
//
// Sizes at the max order: Btil,Dtil are p*(p+1); D,W are (p+1)^2.
// At P=8 (n=9): Btil/Dtil = 72 each, D/W = 81 each -> 306 doubles = 2.4KB,
// well under the 64KB constant bank.
// --------------------------------------------------------------------------
#ifndef HO_CVFEM_MAX_P
#define HO_CVFEM_MAX_P 8
#endif
namespace detail {
constexpr int kMaxN  = HO_CVFEM_MAX_P + 1;
constexpr int kBDlen = HO_CVFEM_MAX_P * kMaxN;   // Btil / Dtil
constexpr int kDWlen = kMaxN * kMaxN;            // D / W
}

__constant__ double c_Btil[detail::kBDlen];
__constant__ double c_Dtil[detail::kBDlen];
__constant__ double c_D[detail::kDWlen];
__constant__ double c_W[detail::kDWlen];
// Quadrature/solution node coordinates needed ONLY by the metric-precompute
// kernel (the apply itself never touches geometry). xi = P Gauss points (SCS
// normals), zeta = (P+1) GLL nodes (tangential slots).
__constant__ double c_xi[detail::kMaxN];
__constant__ double c_zeta[detail::kMaxN];

// Upload the host HoCvfemOperators to constant memory. Call ONCE per P before
// launching. Btil/Dtil have length P*(P+1); D/W have length (P+1)^2; xi has P;
// zeta has P+1. Returns the first CUDA error encountered.
inline cudaError_t ho_cvfem_upload_operators(int P,
                                             const double* h_Btil,
                                             const double* h_Dtil,
                                             const double* h_D,
                                             const double* h_W,
                                             const double* h_xi,
                                             const double* h_zeta)
{
    const int n = P + 1;
    cudaError_t err;
    err = cudaMemcpyToSymbol(c_Btil, h_Btil, sizeof(double) * (size_t)P * n); if (err) return err;
    err = cudaMemcpyToSymbol(c_Dtil, h_Dtil, sizeof(double) * (size_t)P * n); if (err) return err;
    err = cudaMemcpyToSymbol(c_D,    h_D,    sizeof(double) * (size_t)n * n); if (err) return err;
    err = cudaMemcpyToSymbol(c_W,    h_W,    sizeof(double) * (size_t)n * n); if (err) return err;
    err = cudaMemcpyToSymbol(c_xi,   h_xi,   sizeof(double) * (size_t)P);     if (err) return err;
    err = cudaMemcpyToSymbol(c_zeta, h_zeta, sizeof(double) * (size_t)n);     if (err) return err;
    return cudaSuccess;
}

// GMode selector mirrors the DMMA kernel.
enum HoCvfemGMode { HO_CVFEM_PERPOINT = 0, HO_CVFEM_AFFINE = 1 };

// NOTE: the kernel signatures and bodies below stay visible in BOTH nvcc passes
// (host + device). We deliberately do NOT wrap them in `#if __CUDA_ARCH__`: the
// launchers further down issue `kernel<<<...>>>` and are compiled in the HOST
// pass, where the kernel templates must exist to be addressable. This code uses
// only atomicAdd/__syncthreads (valid device functions, never instantiated for
// host), so no host-side trap body is needed -- mirrors mars_ho_laplacian_dmma.

// idx(dir,nrm,t1,t2) -- VERBATIM port of the host lambda in applyHoCvfemElement.
template<int N, int NN>
__device__ __forceinline__ int ho_idx(int dir, int nrm, int t1, int t2)
{
    if (dir == 0) return nrm * NN + t1 * N + t2;
    if (dir == 1) return t1 * NN + nrm * N + t2;
    return t1 * NN + t2 * N + nrm;
}

// Reference-cube corner signs, matching host hexCornerRef() exactly.
__constant__ int c_hexCornerRef[8][3] =
    {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};

// --------------------------------------------------------------------------
// Metric-precompute kernel (PerPoint). This is a PRECOMPUTE-ONCE buffer: d_G is
// built once and reused across every CG apply (the apply never touches geometry),
// so its register footprint / per-point corner reloads are amortized and not on
// the apply critical path. One thread per (e, dir, l, s, r) point.
// Ports computeElementMetric VERBATIM: metric = detJ * J^{-1} J^{-T}, with
// G[((dir*p+l)*n+s)*n+r] = (g[0]=tang2(r-axis), g[1]=tang1(s-axis), g[2]=normal).
// Output d_G is [numElements * 3*p*n*n] contiguous vec3 -- the EXACT host layout
// the apply reads. d_corners is [numElements][8][3] (8 corners x xyz per element).
// --------------------------------------------------------------------------
template<typename RealType, int P>
__global__ void ho_cvfem_metric_perpoint_kernel(const RealType* __restrict__ d_corners,
                                                RealType* __restrict__ d_G,
                                                size_t numElements)
{
    constexpr int N  = P + 1;
    constexpr int NN = N * N;
    constexpr int pointsPerElem = 3 * P * NN;   // (dir,l,s,r) over 3*p*n*n

    const size_t gid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= numElements * (size_t)pointsPerElem) return;

    const size_t e   = gid / pointsPerElem;
    int rem          = (int)(gid % pointsPerElem);
    const int dir    = rem / (P * NN);  rem %= (P * NN);
    const int l      = rem / NN;        rem %= NN;
    const int s      = rem / N;
    const int r      = rem % N;

    // Load 8 corners for this element (xyz contiguous per corner).
    double corners[8][3];
    const RealType* cp = d_corners + e * 24;
    #pragma unroll
    for (int c = 0; c < 8; ++c) {
        corners[c][0] = cp[c * 3 + 0];
        corners[c][1] = cp[c * 3 + 1];
        corners[c][2] = cp[c * 3 + 2];
    }

    // tangential axes per direction (host t1axis/t2axis).
    const int t1axis[3] = {1, 0, 0}, t2axis[3] = {2, 2, 1};
    double rf[3];
    rf[dir]          = c_xi[l];
    rf[t1axis[dir]]  = c_zeta[s];
    rf[t2axis[dir]]  = c_zeta[r];

    // Trilinear Jacobian J[a][b] = dx_a/dr_b -- VERBATIM port of the host
    // trilinearJacobian (inlined; the metric is computed once per point).
    double J[3][3];
    #pragma unroll
    for (int a = 0; a < 3; ++a)
        #pragma unroll
        for (int b = 0; b < 3; ++b) J[a][b] = 0.0;
    #pragma unroll
    for (int c = 0; c < 8; ++c) {
        const int S0 = c_hexCornerRef[c][0], S1 = c_hexCornerRef[c][1], S2 = c_hexCornerRef[c][2];
        double sh[3];
        sh[0] = 0.125 * S0 * (1 + S1 * rf[1]) * (1 + S2 * rf[2]);
        sh[1] = 0.125 * (1 + S0 * rf[0]) * S1 * (1 + S2 * rf[2]);
        sh[2] = 0.125 * (1 + S0 * rf[0]) * (1 + S1 * rf[1]) * S2;
        #pragma unroll
        for (int a = 0; a < 3; ++a)
            #pragma unroll
            for (int b = 0; b < 3; ++b) J[a][b] += sh[b] * corners[c][a];
    }

    double det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
               - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
               + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
    double Ji[3][3];
    Ji[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[2][1])/det; Ji[0][1]=(J[0][2]*J[2][1]-J[0][1]*J[2][2])/det; Ji[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/det;
    Ji[1][0]=(J[1][2]*J[2][0]-J[1][0]*J[2][2])/det; Ji[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[2][0])/det; Ji[1][2]=(J[0][2]*J[1][0]-J[0][0]*J[1][2])/det;
    Ji[2][0]=(J[1][0]*J[2][1]-J[1][1]*J[2][0])/det; Ji[2][1]=(J[0][1]*J[2][0]-J[0][0]*J[2][1])/det; Ji[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[1][0])/det;

    // Gvec = detJ * (J^{-1} J^{-T})[:,dir]
    double gvec[3];
    #pragma unroll
    for (int a = 0; a < 3; ++a) {
        double v = 0;
        #pragma unroll
        for (int k = 0; k < 3; ++k) v += Ji[a][k] * Ji[dir][k];
        gvec[a] = det * v;
    }

    RealType* g = d_G + (e * (size_t)pointsPerElem + (size_t)(((dir * P + l) * N + s) * N + r)) * 3;
    g[0] = gvec[t2axis[dir]];   // tang2 (r-axis)
    g[1] = gvec[t1axis[dir]];   // tang1 (s-axis)
    g[2] = gvec[dir];           // normal
}

// --------------------------------------------------------------------------
// The apply kernel. One CUDA thread per (s,r) tangential face slot; E elements
// packed per block (E = ElemsPerBlock). Thread t in [0, blockDim.x):
//   localElem = t / NN,  slot = t % NN,  s = slot / N,  r = slot % N.
// At high p NN can exceed the per-element thread count we want; we grid-stride
// the slot over the element's threads so a block can use fewer than NN threads
// per element (capped block) -- see threadsPerElem in the launcher.
//
// Per element, shared memory holds:
//   u_sh[E*N3] : gathered input, read-only across all 3 directions (re-strided
//                by ho_idx, NOT triple-allocated).
//   y_sh[E*N3] : accumulator, zeroed once. CANNOT alias u_sh: u is still read in
//                dir1/dir2 while y accumulates from dir0.
//   face[E*2*N*(N+1)]: two ping-pong tangential exchange buffers (interp, then
//                 flux, then tmp), reused every (dir,l). Row pitch is padded to
//                 N+1 (not N) to break shared-memory bank conflicts on the
//                 column-stride tangential reads -- pure layout, bit-identical.
//
// Each (dir,l) iteration:
//   step 1 (no sync): thread (s,r) reads its normal column u_sh[idx(dir,q,s,r)]
//          and forms interp_sr (Btil) and deriv_sr (Dtil). deriv stays in a
//          register (only this thread needs it). interp -> faceA[slot].
//   __syncthreads (interp visible to tangential D).
//   step 2 (no sync after, until flux staged): dt2 = sum_q D[r,q]*faceA[s,q],
//          dt1 = sum_q D[s,q]*faceA[q,r]; flux_sr = g2*deriv + g0*dt2 + g1*dt1.
//          flux -> faceB[slot].
//   __syncthreads.
//   step 3: tmp_sr = sum_q W[r,q]*faceB[s,q]; tmp -> faceA[slot] (reuse).
//   __syncthreads.
//   step 4: intf_sr = sum_q W[s,q]*faceA[q,r]; scatter into y_sh at nodes l,l+1.
//   __syncthreads (faceA reused by next l's interp).
//
// The scatter to y_sh is hazard-free WITHOUT atomics: thread (s,r) owns the
// full normal column for this dir, so only it writes y_sh[idx(dir,l,s,r)] and
// y_sh[idx(dir,l+1,s,r)]. Across the serial l-loop, node l+1 is touched by the
// "+=" of iteration l and the "-=" of iteration l+1, BOTH by the same thread in
// the same dir -> no cross-thread hazard, no sync needed for the y_sh write
// itself. (The trailing __syncthreads only guards the faceA reuse.)
// --------------------------------------------------------------------------
// d_elemList (optional): when non-null, the block/slot index addresses a SUBSET
// of elements -- element e = d_elemList[slot_e] for slot_e < count, else the
// invalid sentinel numElements (so tail threads still hit every block-collective
// __syncthreads). When null the index IS the element id (slot_e), bit-identical
// to the original full-range path. Used by the distributed overlap matvec to
// apply interior elements before the halo wait and boundary elements after.
template<typename RealType, int P, int BlockSize, int ElemsPerBlock, int GMode>
__global__ void __launch_bounds__(BlockSize)
ho_cvfem_apply_kernel(const RealType* __restrict__ d_u,
                      RealType* __restrict__ d_y,
                      const int* __restrict__ d_elemDof,
                      const RealType* __restrict__ d_G,
                      size_t numElements,
                      const int* __restrict__ d_elemList = nullptr,
                      size_t count = 0)
{
    // The host reference applyHoCvfemElement accumulates in double. Only the
    // double instantiation is gated to <1e-12. A float build would do the
    // sum-factorization in float (lower accuracy than the spec) and is
    // unvalidated -- block it until a float-tolerance gate exists.
    static_assert(std::is_same<RealType, double>::value,
        "HO-CVFEM matrix-free apply is validated for RealType=double only");
    constexpr int N   = P + 1;
    constexpr int NN  = N * N;
    constexpr int N3  = NN * N;
    constexpr int E   = ElemsPerBlock;
    constexpr int threadsPerElem = BlockSize / E;   // threads working one element
    // Face-buffer row pitch. N+1 padding was tried to break the column-stride bank
    // conflicts (ncu: ~4-way on ~55% of shared accesses) but REVERTED: on H100 it
    // raised registers 40->44, dropping theoretical occupancy 75%->62.5%, which
    // exactly canceled the conflict win -- the kernel is L1-throughput bound, not
    // occupancy bound, so neither moved throughput. NP=N is the cleaner baseline.
    constexpr int NP  = N;            // no padding (see note)
    constexpr int NNP = N * NP;       // face-buffer length (per buffer) = N*N

    // Shared layout: [u_sh | y_sh | face]. DYNAMIC smem (extern) so the launcher
    // can opt into the Hopper 228KB carveout via cudaFuncSetAttribute. Under the
    // 48KB static default the binding resource at p=2/3 was smem (E*2*N3 doubles),
    // capping blocks/SM to 2-3 = 25-37% occupancy; the carveout lifts those orders
    // back to thread-bound (~100%). The byte count passed at launch is
    // (E*2*N3 + E*2*NNP)*sizeof(RealType), NNP = N*(N+1) (padded face buffers).
    extern __shared__ char ho_smem_raw[];
    RealType* u_sh = reinterpret_cast<RealType*>(ho_smem_raw);
    RealType* y_sh = u_sh + (size_t)E * N3;
    RealType* face = y_sh + (size_t)E * N3;   // 2 ping-pong nn buffers per element

    const int t          = threadIdx.x;
    const int localElem  = t / threadsPerElem;
    const int laneInElem = t % threadsPerElem;   // 0..threadsPerElem-1
    // slot_e is the dense per-block element index (the original `e` expression).
    // With a list, it indexes d_elemList; without, it IS the element id. The full
    // -range path takes e = slot_e -- the EXACT original expression, so the only
    // added cost is one predicated select on a compile-time-null pointer per
    // element, changing no FP reduction, smem layout, or thread->slot mapping.
    const size_t slot_e  = (size_t)blockIdx.x * E + localElem;
    const size_t e       = (d_elemList != nullptr)
                             ? (slot_e < count ? (size_t)d_elemList[slot_e] : numElements)
                             : slot_e;

    // Per-element slabs.
    RealType* my_u    = u_sh + localElem * N3;
    RealType* my_y    = y_sh + localElem * N3;
    RealType* faceA   = face + localElem * 2 * NNP;
    RealType* faceB   = faceA + NNP;
    // A partial last block can have some localElems with e >= numElements. Those
    // threads MUST still execute every __syncthreads below (block-collective), so
    // we do NOT branch the loop on validity -- invalid threads just compute into
    // their own private smem slab (never read by a valid element) and skip the
    // global gather/scatter. `valid` guards only the global-memory touches.
    const bool valid = (e < numElements);
    const int* edof  = valid ? (d_elemDof + e * N3) : nullptr;

    // --- Gather u and zero y (grid-strided over the element's threads). dof<0
    //     means a constrained/absent node -> contributes 0, exactly like the
    //     host gather and the p=1 matfree kernel. Invalid elems zero u_sh so the
    //     contractions stay finite (NaN-free) even though their output is unused. ---
    // INVARIANT: edof is nullptr for invalid elems (e>=numElements). The `valid &&`
    // short-circuit MUST stay first so edof[l] is never dereferenced when invalid.
    for (int l = laneInElem; l < N3; l += threadsPerElem) {
        my_u[l] = (valid && edof[l] >= 0) ? d_u[edof[l]] : RealType(0);
        my_y[l] = RealType(0);
    }
    __syncthreads();

    // Metric base for this element. PerPoint stores a vec3 per (dir,l,s,r) point,
    // so the per-element stride is 3*P*NN*3 doubles -- this MUST match the metric
    // kernel's write index d_G[(e*pointsPerElem + idx)*3] with pointsPerElem =
    // 3*P*NN (line below). The within-element per-point *3 is added at the read in
    // step 2. Affine stores 3 doubles/element. Invalid elems point at element 0's
    // metric -- harmless, output is discarded.
    const RealType* Gbase = d_G + (valid ? e : 0) * (size_t)(GMode == HO_CVFEM_PERPOINT ? 3 * P * NN * 3 : 3);

    // Keep dir and l ROLLED. Full-unrolling them replicates the whole 4-step
    // contraction body 3*P times (21 copies at p=7), exploding registers and
    // I-cache on the very kernel whose goal is to RAISE occupancy by shrinking
    // the per-thread footprint. Only the short length-N inner q-loops are
    // unrolled (genuinely beneficial, bounded). Verify regs/thread with
    // -Xptxas -v on the cluster.
    for (int dir = 0; dir < 3; ++dir)
    {
            for (int l = 0; l < P; ++l)
            {
                // step 1: normal-direction interp(Btil) + deriv(Dtil). Each slot
                // owns its (s,r); we grid-stride slots so threadsPerElem may be
                // < NN at high p. interp goes to the shared face buffer (needed
                // by tangential D across threads); deriv is needed only by this
                // thread at the same (s,r), so it stays in a register cache
                // indexed by the grid-stride iteration. The cache size is the
                // worst-case slots handled by one thread = ceil(NN/threadsPerElem).
                constexpr int kDerivIters = (NN + threadsPerElem - 1) / threadsPerElem;
                RealType deriv_cache[kDerivIters];
                int it = 0;
                for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {
                    const int s = slot / N, r = slot % N;
                    RealType bi = 0, di = 0;
                    #pragma unroll
                    for (int q = 0; q < N; ++q) {
                        RealType uq = my_u[ho_idx<N, NN>(dir, q, s, r)];
                        bi += c_Btil[l * N + q] * uq;
                        di += c_Dtil[l * N + q] * uq;
                    }
                    faceA[s * NP + r] = bi;      // interp -> face buffer A (padded pitch)
                    deriv_cache[it]   = di;      // deriv stays register-resident
                }
                __syncthreads();                 // interp visible for tangential D

                // step 2: tangential D-derivatives of interp + flux assembly.
                // For Affine the cross-term coefficients g0,g1 are compile-time
                // zero, so the two tangential D-contractions (each 2*N faceA smem
                // reads/slot -- the largest step-2 smem-traffic term) are skipped
                // entirely: flux is just g2*deriv. This is the whole point of the
                // Affine specialization on a smem-bandwidth-bound kernel.
                it = 0;
                for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {
                    if constexpr (GMode == HO_CVFEM_PERPOINT) {
                        const int s = slot / N, r = slot % N;
                        RealType dt2 = 0, dt1 = 0;
                        #pragma unroll
                        for (int q = 0; q < N; ++q) dt2 += c_D[r * N + q] * faceA[s * NP + q]; // r-axis (t2)
                        #pragma unroll
                        for (int q = 0; q < N; ++q) dt1 += c_D[s * N + q] * faceA[q * NP + r]; // s-axis (t1)
                        const RealType* g = Gbase + (size_t)(((dir * P + l) * N + s) * N + r) * 3;
                        faceB[s * NP + r] = g[2] * deriv_cache[it] + g[0] * dt2 + g[1] * dt1;
                    } else {
                        // Affine: element-constant diagonal metric. For an
                        // orthogonal/axis-aligned hex the cross terms vanish and
                        // Gbase[dir] is the per-element constant normal coefficient
                        // -- the SAME role as `coeff` in applyHoCvfemElementCube
                        // (W already carries the face quadrature, so NO extra
                        // weight is folded here; folding one would over-count).
                        const int s = slot / N, r = slot % N;
                        faceB[s * NP + r] = Gbase[dir] * deriv_cache[it];
                    }
                }
                __syncthreads();                 // flux visible for W-integration

                // step 3: W-integration along r-axis (t2). tmp -> faceA (reuse).
                for (int slot = laneInElem; slot < NN; slot += threadsPerElem) {
                    const int s = slot / N, r = slot % N;
                    RealType v = 0;
                    #pragma unroll
                    for (int q = 0; q < N; ++q) v += c_W[r * N + q] * faceB[s * NP + q];
                    faceA[s * NP + r] = v;
                }
                __syncthreads();

                // step 4: W-integration along s-axis (t1) + distribute -/+ to the
                // two SCS-bounding nodes l, l+1. No atomics: this thread owns the
                // whole normal column for (dir,s,r).
                for (int slot = laneInElem; slot < NN; slot += threadsPerElem) {
                    const int s = slot / N, r = slot % N;
                    RealType intf = 0;
                    #pragma unroll
                    for (int q = 0; q < N; ++q) intf += c_W[s * N + q] * faceA[q * NP + r];
                    my_y[ho_idx<N, NN>(dir, l,     s, r)] -= intf;
                    my_y[ho_idx<N, NN>(dir, l + 1, s, r)] += intf;
                }
                __syncthreads();                 // faceA reused by next l's interp
            }
        }

    // --- Scatter to global (additive). dof<0 skipped, mirroring the gather and
    //     the p=1 matfree atomic-scatter to owned rows. The if(e<numElements)
    //     guard is required: edof is nullptr for invalid elems. ---
    if (e < numElements) {
        for (int l = laneInElem; l < N3; l += threadsPerElem) {
            int dof = edof[l];
            if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);
        }
    }
}

// --------------------------------------------------------------------------
// Default ElemsPerBlock / BlockSize per P (occupancy plan). The binding rule:
// threadsPerElem = Block/E should be the SMALLEST warp multiple >= NN, so no
// warp is more than partially idle and no whole warp does zero work. Decoupling
// threadsPerElem from NN was the source of the dead-lane waste (p>=4 had
// threadsPerElem=64 with NN=25 -> 39 idle lanes/elem, a half-empty block on an
// occupancy-bound kernel). Block is then E*threadsPerElem.
//   p=1: NN=4,  tpe=4  (8 elems/warp, full)              E=64  Block=256
//   p=2: NN=9,  tpe=9  (tpe==NN, no idle)                E=28  Block=252
//   p=3: NN=16, tpe=16 (tpe==NN)                         E=16  Block=256
//   p=4: NN=25, tpe=32 (1 warp, 7 idle)                  E=8   Block=256
//   p=5: NN=36, tpe=64 (2 warps, 28 idle)                E=2   Block=128
//   p=6: NN=49, tpe=64 (2 warps, 15 idle)                E=2   Block=128
//   p=7: NN=64, tpe=64 (2 warps, 0 idle)                 E=2   Block=128
//   p=8: NN=81, tpe=96 (3 warps, 15 idle)                E=1   Block=96
// kDerivIters = ceil(NN/tpe) is 1 for all of these. These are reasoned, not
// autotuned; the user should sweep E/Block on-device with -Xptxas -v + ncu.
// --------------------------------------------------------------------------
template<int P> struct HoCvfemLaunchDefault;
template<> struct HoCvfemLaunchDefault<1> { static constexpr int Block = 256; static constexpr int Elems = 64; };   // tpe=4
template<> struct HoCvfemLaunchDefault<2> { static constexpr int Block = 252; static constexpr int Elems = 28; };   // tpe=9
template<> struct HoCvfemLaunchDefault<3> { static constexpr int Block = 256; static constexpr int Elems = 16; };   // tpe=16
template<> struct HoCvfemLaunchDefault<4> { static constexpr int Block = 256; static constexpr int Elems = 8;  };   // tpe=32 (NN=25)
template<> struct HoCvfemLaunchDefault<5> { static constexpr int Block = 128; static constexpr int Elems = 2;  };   // tpe=64 (NN=36)
template<> struct HoCvfemLaunchDefault<6> { static constexpr int Block = 128; static constexpr int Elems = 2;  };   // tpe=64 (NN=49)
template<> struct HoCvfemLaunchDefault<7> { static constexpr int Block = 128; static constexpr int Elems = 2;  };   // tpe=64 (NN=64)
template<> struct HoCvfemLaunchDefault<8> { static constexpr int Block = 96;  static constexpr int Elems = 1;  };   // tpe=96 (NN=81)

// Metric-precompute launcher (PerPoint). d_corners is [numElements*8*3], d_G is
// [numElements * 3*P*(P+1)^2 * 3]. Operators (xi/zeta) MUST be uploaded first.
template<typename RealType, int P, int BlockSize = 256>
inline cudaError_t ho_cvfem_metric_perpoint_launch(const RealType* d_corners,
                                                   RealType* d_G,
                                                   size_t numElements,
                                                   cudaStream_t stream = 0)
{
    constexpr int N = P + 1;
    const size_t total = numElements * (size_t)(3 * P * N * N);
    const size_t grid  = (total + BlockSize - 1) / BlockSize;
    ho_cvfem_metric_perpoint_kernel<RealType, P>
        <<<(unsigned)grid, BlockSize, 0, stream>>>(d_corners, d_G, numElements);
    return cudaGetLastError();
}

// Generic launcher. grid = ceil(numElements / E). Uses DYNAMIC shared memory and
// opts into the Hopper 228KB carveout once (idempotent, cheap) so the mid orders
// p=2/3 are thread-bound rather than smem-capped. d_y MUST be pre-zeroed.
// Operators MUST be uploaded once via ho_cvfem_upload_operators before launch.
// d_elemList/count (optional): when d_elemList!=nullptr the launch applies only
// the `count` elements it names (grid = ceil(count/E)); numElements stays the DOF
// validity sentinel. When null the full [0,numElements) range runs, bit-identical
// to before (same grid, same kernel index expression).
template<typename RealType, int P, int BlockSize, int ElemsPerBlock, int GMode>
inline cudaError_t ho_cvfem_apply_launch_impl(const RealType* d_u,
                                              RealType* d_y,
                                              const int* d_elemDof,
                                              const RealType* d_G,
                                              size_t numElements,
                                              cudaStream_t stream,
                                              const int* d_elemList = nullptr,
                                              size_t count = 0)
{
    static_assert(BlockSize % ElemsPerBlock == 0,
        "BlockSize must be an integer multiple of ElemsPerBlock (threadsPerElem = BlockSize/E)");
    // A subset list must never name the invalid sentinel value numElements: a tail
    // thread in a partial block gets e=numElements and relies on valid=(e<numElements)
    // being false to skip its global gather/scatter. If a list entry equaled
    // numElements it would alias that sentinel and be silently dropped. Interior/
    // boundary lists only hold ids in [0,numElements) by construction, so this is a
    // contract guard for future callers, not a live bug.
    assert(d_elemList == nullptr || count <= numElements);
    constexpr int N   = P + 1;
    constexpr int NN  = N * N;
    constexpr int N3  = NN * N;
    // Must match the in-kernel layout: [u_sh | y_sh | face]. Face buffers use pitch
    // NP=N (the N+1 padding was reverted -- see kernel note); must track NNP=N*N.
    constexpr int NNP = N * N;
    constexpr int smemBytes =
        (int)(((size_t)ElemsPerBlock * 2 * N3 + (size_t)ElemsPerBlock * 2 * NNP) * sizeof(RealType));

    auto kernel = ho_cvfem_apply_kernel<RealType, P, BlockSize, ElemsPerBlock, GMode>;
    // Opt into >48KB dynamic smem; harmless when smemBytes <= 48KB. Idempotent.
    cudaError_t attrErr = cudaFuncSetAttribute(
        kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, smemBytes);
    if (attrErr != cudaSuccess) return attrErr;

    const size_t launchElems = (d_elemList != nullptr) ? count : numElements;
    const size_t grid = (launchElems + ElemsPerBlock - 1) / ElemsPerBlock;
    if (grid == 0) return cudaSuccess;   // empty subset: nothing to launch
    kernel<<<(unsigned)grid, BlockSize, smemBytes, stream>>>(
        d_u, d_y, d_elemDof, d_G, numElements, d_elemList, count);
    return cudaGetLastError();
}

// PerPoint (general/curved hexes). d_G = [numElements * 3*P*(P+1)^2] vec3.
template<typename RealType, int P,
         int BlockSize = HoCvfemLaunchDefault<P>::Block,
         int ElemsPerBlock = HoCvfemLaunchDefault<P>::Elems>
inline cudaError_t ho_cvfem_apply_launch(const RealType* d_u,
                                         RealType* d_y,
                                         const int* d_elemDof,
                                         const RealType* d_G,
                                         size_t numElements,
                                         cudaStream_t stream = 0,
                                         const int* d_elemList = nullptr,
                                         size_t count = 0)
{
    return ho_cvfem_apply_launch_impl<RealType, P, BlockSize, ElemsPerBlock, HO_CVFEM_PERPOINT>(
        d_u, d_y, d_elemDof, d_G, numElements, stream, d_elemList, count);
}

// Affine (axis-aligned/orthogonal hexes). d_Ghat = [numElements * 3] per-element
// constant normal coefficients (one per direction). Cross terms are zero.
//
// WARNING: UNVALIDATED. No gate exercises this path and there is no device
// helper that fills d_Ghat -- the caller must supply, for each element and
// direction dir, d_Ghat[e*3+dir] = the per-element constant gvec[dir] from
// computeElementMetric (the [2]=normal component, which is constant across all
// points only for an orthogonal/straight axis-aligned hex). For any sheared or
// curved hex the cross terms g0,g1 are nonzero and this path is SILENTLY WRONG;
// use ho_cvfem_apply_launch (PerPoint) instead. Add an affine gate to the test
// (compare vs PerPoint on an axis-aligned cube to <1e-12) before relying on it.
template<typename RealType, int P,
         int BlockSize = HoCvfemLaunchDefault<P>::Block,
         int ElemsPerBlock = HoCvfemLaunchDefault<P>::Elems>
inline cudaError_t ho_cvfem_apply_launch_affine(const RealType* d_u,
                                                RealType* d_y,
                                                const int* d_elemDof,
                                                const RealType* d_Ghat,
                                                size_t numElements,
                                                cudaStream_t stream = 0)
{
    return ho_cvfem_apply_launch_impl<RealType, P, BlockSize, ElemsPerBlock, HO_CVFEM_AFFINE>(
        d_u, d_y, d_elemDof, d_Ghat, numElements, stream);
}

} // namespace fem
} // namespace mars
