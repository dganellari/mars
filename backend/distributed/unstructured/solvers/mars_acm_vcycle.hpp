#ifndef MARS_ACM_VCYCLE_HPP
#define MARS_ACM_VCYCLE_HPP

// ACM Stage 4a: GPU multilevel V-cycle (Darwish-Saad-Hamdan 2008), single rank.
//
// Hierarchy = host-once greedy aggregation per level (the aggregation MAP is setup, not the hot
// cycle; the GPU DIRECTIONAL aggregation is Stage 5) + the Stage-3 GPU coarse operator
// (sum-of-fine == P^T A P) + order-zero injection / sum restriction. Smoother = damped point-Jacobi
// (Stage 1 showed block-Jacobi is unnecessary on the PSPG operator -- its pressure diagonal is O(0.1-1),
// not ~0). Coarsest = many Jacobi sweeps (the paper's serial V-cycle; the parallel master-gather
// direct solve is Stage 6). The whole apply path runs on the GPU (no D2H inside the cycle).
//
// This header also carries a HOST V-cycle replica (acmHostVcycle) built from the same hierarchy, used
// only as the Stage-4a validation oracle: one GPU cycle must equal one host cycle to round-off. We do
// NOT validate via standalone MG iteration -- stationary MG can diverge on an indefinite saddle
// operator for legitimate reasons; Krylov acceleration (FlexGMRES, Stage 4b) is the real harness.

#include "backend/distributed/unstructured/solvers/mars_acm_coarsen.hpp"
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/system/cuda/execution_policy.h>   // thrust::cuda::par.on(stream) for stream-ordered fills
#include <cusolverSp.h>
#include <cusparse.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

namespace mars {

template<class V> inline auto acmRaw(V& v) { return thrust::raw_pointer_cast(v.data()); }

template<typename RealType>
__global__ void acmSpmvKernel(const int* rowOff, const int* colInd, const RealType* vals,
                              const RealType* x, RealType* y, int ND)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x; if (r >= ND) return;
    RealType s = 0;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) s += vals[k] * x[colInd[k]];
    y[r] = s;
}

template<typename RealType>
__global__ void acmDinvKernel(const int* rowOff, const int* colInd, const RealType* vals,
                              RealType* dinv, int ND)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x; if (r >= ND) return;
    RealType d = 0;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) if (colInd[k] == r) { d = vals[k]; break; }
    dinv[r] = (fabs((double)d) > 1e-30) ? RealType(1) / d : RealType(0);
}

// one damped point-Jacobi sweep: xout = xin + omega * dinv * (b - A xin)   (proper Jacobi, ping-pong)
template<typename RealType>
__global__ void acmJacobiKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                const RealType* dinv, const RealType* b,
                                const RealType* xin, RealType* xout, RealType omega, int ND)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x; if (r >= ND) return;
    RealType s = 0;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) s += vals[k] * xin[colInd[k]];
    xout[r] = xin[r] + omega * dinv[r] * (b[r] - s);
}

// ---- COUPLED (block) smoother: invert the per-node 4x4 [u,v,w,p] block so each relaxation captures
// the local pressure-velocity coupling (point-Jacobi is segregated -- it never sees the G/D coupling).
// 4x4 Gauss-Jordan inverse with partial pivoting + tiny-pivot regularization (host+device).
template<typename T>
__host__ __device__ inline void acmInvert4x4(const T* A, T* Ai)
{
    T m[4][8];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) m[i][j] = A[4 * i + j];
        for (int j = 0; j < 4; ++j) m[i][4 + j] = (i == j) ? T(1) : T(0);
    }
    for (int col = 0; col < 4; ++col) {
        int piv = col; T best = fabs((double)m[col][col]);
        for (int r = col + 1; r < 4; ++r) { double v = fabs((double)m[r][col]); if (v > best) { best = v; piv = r; } }
        if (piv != col) for (int j = 0; j < 8; ++j) { T t = m[col][j]; m[col][j] = m[piv][j]; m[piv][j] = t; }
        T d = m[col][col];
        if (fabs((double)d) < 1e-30) d = (d >= T(0)) ? T(1e-30) : T(-1e-30);   // regularize singular block
        for (int j = 0; j < 8; ++j) m[col][j] /= d;
        for (int r = 0; r < 4; ++r) if (r != col) { T f = m[r][col]; for (int j = 0; j < 8; ++j) m[r][j] -= f * m[col][j]; }
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) Ai[4 * i + j] = m[i][4 + j];
}

// per node: extract the 4x4 diagonal [u,v,w,p] block from the CSR and store its inverse (16/node)
template<typename RealType>
__global__ void acmBlockInvKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                  RealType* binv, int nNodes)
{
    int node = blockIdx.x * blockDim.x + threadIdx.x; if (node >= nNodes) return;
    RealType B[16];
    for (int t = 0; t < 16; ++t) B[t] = RealType(0);
    for (int a = 0; a < 4; ++a) {
        int r = 4 * node + a;
        for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) {
            int c = colInd[k];
            if (c >= 4 * node && c < 4 * node + 4) B[4 * a + (c - 4 * node)] = vals[k];
        }
    }
    acmInvert4x4(B, binv + 16 * node);
}

// one damped BLOCK-Jacobi sweep: per node x_blk += omega * Binv * (b - A xin)_blk  (couples u,v,w,p)
template<typename RealType>
__global__ void acmBlockJacobiKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                     const RealType* binv, const RealType* b,
                                     const RealType* xin, RealType* xout, RealType omega, int nNodes)
{
    int node = blockIdx.x * blockDim.x + threadIdx.x; if (node >= nNodes) return;
    RealType res[4];
    for (int a = 0; a < 4; ++a) {
        int r = 4 * node + a; RealType s = 0;
        for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) s += vals[k] * xin[colInd[k]];
        res[a] = b[r] - s;
    }
    const RealType* Bi = binv + 16 * node;
    for (int a = 0; a < 4; ++a) {
        RealType dx = 0;
        for (int c = 0; c < 4; ++c) dx += Bi[4 * a + c] * res[c];
        xout[4 * node + a] = xin[4 * node + a] + omega * dx;
    }
}

// ---- block-ILU(0) smoother (darwish2008a's smoother) ----
// 4x4 block matmul C = A*B (row-major). Jacobi is weak on convection-dominated operators (smoothing
// factor -> 1 as Pe grows, host-validated in scripts/acm_ilu_smoother_study.py); ILU(0) sweeps error
// along the characteristics so it stays a strong smoother at high Pe. The triangular solves are made
// parallel by LEVEL-SCHEDULING on the node block-graph (a level = rows whose lower-deps are all earlier).
template<typename T>
__host__ __device__ inline void acmBlkMatmul4(const T* A, const T* B, T* C)
{
    for (int a = 0; a < 4; ++a)
        for (int c = 0; c < 4; ++c) { T s = 0; for (int b = 0; b < 4; ++b) s += A[4*a+b]*B[4*b+c]; C[4*a+c] = s; }
}

// forward block L-solve for one level: y_i = r_i - sum_{k<i stored} L_ik * y_k   (unit lower-block-diag)
template<typename RealType>
__global__ void acmIluLsolveKernel(const int* brow, const int* bcol, const int* bdiagPtr, const RealType* bilu,
                                   const int* nodesByLevel, int lstart, int lend, const RealType* r, RealType* y)
{
    int t = lstart + blockIdx.x * blockDim.x + threadIdx.x; if (t >= lend) return;
    int i = nodesByLevel[t];
    RealType yi[4]; for (int a = 0; a < 4; ++a) yi[a] = r[4*i+a];
    for (int kk = brow[i]; kk < bdiagPtr[i]; ++kk) {              // strictly-lower block-cols k<i
        int k = bcol[kk]; const RealType* L = bilu + 16*kk;
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) yi[a] -= L[4*a+b] * y[4*k+b];
    }
    for (int a = 0; a < 4; ++a) y[4*i+a] = yi[a];
}

// back block U-solve for one level: dx_i = inv(U_ii) * (y_i - sum_{j>i stored} U_ij * dx_j)
template<typename RealType>
__global__ void acmIluUsolveKernel(const int* brow, const int* bcol, const int* bdiagPtr, const RealType* bilu,
                                   const int* nodesByLevel, int lstart, int lend, const RealType* y, RealType* dx)
{
    int t = lstart + blockIdx.x * blockDim.x + threadIdx.x; if (t >= lend) return;
    int i = nodesByLevel[t];
    RealType ti[4]; for (int a = 0; a < 4; ++a) ti[a] = y[4*i+a];
    for (int jj = bdiagPtr[i] + 1; jj < brow[i+1]; ++jj) {        // strictly-upper block-cols j>i
        int j = bcol[jj]; const RealType* U = bilu + 16*jj;
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) ti[a] -= U[4*a+b] * dx[4*j+b];
    }
    const RealType* Di = bilu + 16*bdiagPtr[i];                   // stored inv(U_ii)
    for (int a = 0; a < 4; ++a) { RealType s = 0; for (int b = 0; b < 4; ++b) s += Di[4*a+b]*ti[b]; dx[4*i+a] = s; }
}

// ---- multicolor block-DILU solve kernels ----
// Split lower/upper by COLOR (not bdiagPtr/index): color order != index order, so a neighbor k<i can have
// color[k]>color[i] and belong to U. So we iterate the FULL row and branch on color[]. Off-diagonals are
// raw A (bblk); inv(D_i) is NODE-indexed (dinvD). One color per launch -> all nodes in it are independent.
template<typename RealType>
__global__ void acmDiluLsolveKernel(const int* brow, const int* bcol, const int* color,
        const RealType* bblk, const RealType* dinvD, const int* nodesByColor,
        int cstart, int cend, const RealType* r, RealType* y, RealType* w)
{
    int t = cstart + blockIdx.x * blockDim.x + threadIdx.x; if (t >= cend) return;
    int i = nodesByColor[t], ci = color[i];
    RealType wi[4]; for (int a = 0; a < 4; ++a) wi[a] = r[4*i+a];
    for (int kk = brow[i]; kk < brow[i+1]; ++kk) {               // strictly-lower-COLOR off-diags (full row scan)
        int k = bcol[kk]; if (k == i || color[k] >= ci) continue;
        const RealType* Aik = bblk + 16*kk;
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) wi[a] -= Aik[4*a+b] * y[4*k+b];
    }
    for (int a = 0; a < 4; ++a) w[4*i+a] = wi[a];                // w = pre-diagonal residual = z for the back solve
    const RealType* iD = dinvD + 16*i;
    for (int a = 0; a < 4; ++a) { RealType s = 0; for (int b = 0; b < 4; ++b) s += iD[4*a+b]*wi[b]; y[4*i+a] = s; }
}

template<typename RealType>
__global__ void acmDiluUsolveKernel(const int* brow, const int* bcol, const int* color,
        const RealType* bblk, const RealType* dinvD, const int* nodesByColor,
        int cstart, int cend, const RealType* z, RealType* x)
{
    int t = cstart + blockIdx.x * blockDim.x + threadIdx.x; if (t >= cend) return;
    int i = nodesByColor[t], ci = color[i];
    RealType ti[4]; for (int a = 0; a < 4; ++a) ti[a] = z[4*i+a];  // z = w STORED by the forward sweep (NOT recomputed as D_i*y_i; equal by construction since y_i=inv(D_i)w_i)
    for (int jj = brow[i]; jj < brow[i+1]; ++jj) {               // strictly-upper-COLOR off-diags (full row scan)
        int j = bcol[jj]; if (j == i || color[j] <= ci) continue;
        const RealType* Aij = bblk + 16*jj;
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) ti[a] -= Aij[4*a+b] * x[4*j+b];
    }
    const RealType* iD = dinvD + 16*i;
    for (int a = 0; a < 4; ++a) { RealType s = 0; for (int b = 0; b < 4; ++b) s += iD[4*a+b]*ti[b]; x[4*i+a] = s; }
}

template<typename RealType>
__global__ void acmAxpyKernel(RealType* x, const RealType* y, RealType a, int ND)
{ int r = blockIdx.x * blockDim.x + threadIdx.x; if (r < ND) x[r] += a * y[r]; }

// fused residual rtmp = b - A x in one pass (replaces a separate SpMV + subtract): one launch instead of
// two, one global round-trip instead of two. Same arithmetic (sum then subtract) -> bit-identical.
template<typename RealType>
__global__ void acmResidualKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                  const RealType* x, const RealType* b, RealType* rtmp, int ND)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x; if (r >= ND) return;
    RealType s = 0;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) s += vals[k] * x[colInd[k]];
    rtmp[r] = b[r] - s;
}

template<typename RealType>
__global__ void acmAddKernel(RealType* x, const RealType* y, int ND)
{ int r = blockIdx.x * blockDim.x + threadIdx.x; if (r < ND) x[r] += y[r]; }

// scatter a sparse CSR (row-major) into a dense COLUMN-MAJOR n x n matrix (pre-zeroed) for cuSOLVER-Dn
// getrf/getrs. Used to factor the small coarsest operator ONCE per Picard (vs the per-V-cycle cusolverSp QR).
template<typename RealType>
__global__ void acmDenseFromCsrKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                      RealType* dense, int n)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x; if (r >= n) return;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) dense[(size_t)colInd[k] * n + r] = vals[k];   // col-major [col*n+row]
}

template<typename RealType>
struct AcmLevel {
    int nNodes = 0, ND = 0, nCoarse = 0;
    thrust::device_vector<int> rowOff, colInd;
    thrust::device_vector<RealType> vals, dinv, binv;          // dinv: point-Jacobi; binv: 4x4 block-Jacobi
    thrust::device_vector<int> agg;                            // nNodes -> parent (empty on coarsest)
    thrust::device_vector<RealType> bvec, xvec, rtmp, xtmp, ytmp;  // cycle work vectors (ytmp = ILU L-solve scratch)
    // block-ILU(0) smoother (built only when MARS_ACM_SMOOTHER=ilu): node block-CSR + factored 4x4 blocks
    // (L below diag, inv(U_ii) on diag, U above) + level schedule on the lower block-graph.
    thrust::device_vector<int> brow, bcol, bdiagPtr, nodesByLevel;
    thrust::device_vector<RealType> bilu;
    std::vector<int> h_levStart;                               // host: per-level slices into nodesByLevel (launch sizing)
    int nLev = 0;
    // multicolor block-DILU smoother (built only when MARS_ACM_SMOOTHER=dilu): off-diagonals = raw A blocks
    // (bblk), diagonal = factored inv(D_i) (dilu_dinv, NODE-indexed); multicolor schedule (numColors << nLev
    // -> ~few launches/sweep vs ILU's thousands). Reuses brow/bcol/bdiagPtr.
    thrust::device_vector<int> color, nodesByColor;            // color: node->color id; nodesByColor: nodes grouped by color
    thrust::device_vector<RealType> bblk, dilu_dinv, ztmp;     // bblk: raw A 4x4 blocks; dilu_dinv: inv(D_i); ztmp: forward pre-diag residual
    std::vector<int> h_colStart;                               // host: per-color slices into nodesByColor
    int numColors = 0;
};

// host-once DIRECTIONAL aggregation (Mavriplis-style, paper Eqs 19-20) on the block graph of a coupled
// (4*node+comp) CSR. Strength S_IJ = ||A_IJ block||_F symmetrized -- the paper's algebraic "mutual
// coefficients" route: on a stretched mesh the strong couplings run across the thin cell direction, so
// aggregates elongate along the anisotropy. Two-sided strong test with factor beta (paper uses 0.5),
// strongest-first greedy seed-growth, then a cleanup pass that ABSORBS leftover nodes into the strongest
// neighboring aggregate (this is what keeps the coarsening ratio healthy -- avoids the singleton spray
// that made the isotropic greedy stall at ~1.3x on unstructured meshes). beta=0 -> isotropic.
template<typename RealType>
inline void acmDirectionalAggregate(const std::vector<int>& rowOff, const std::vector<int>& colInd,
                                    const std::vector<RealType>& vals, int nNodes, int kmax,
                                    RealType beta, std::vector<int>& agg, int& nCoarse)
{
    const int ND = 4 * nNodes;
    // block-Frobenius^2 per ordered block-pair (I,J), I!=J
    std::unordered_map<long long, RealType> sq;
    sq.reserve(static_cast<size_t>(colInd.size()));
    for (int r = 0; r < ND; ++r) {
        long long I = r >> 2;
        for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) {
            long long J = colInd[k] >> 2;
            if (I != J) sq[I * static_cast<long long>(nNodes) + J] += vals[k] * vals[k];
        }
    }
    // symmetric strength S_IJ = sqrt(sq_IJ + sq_JI), each undirected edge emitted once into both lists
    std::vector<std::vector<std::pair<int, RealType>>> nb(nNodes);
    for (const auto& kv : sq) {
        int I = static_cast<int>(kv.first / nNodes), J = static_cast<int>(kv.first % nNodes);
        if (I < J) {
            auto it = sq.find(J * static_cast<long long>(nNodes) + I);
            RealType s = std::sqrt(kv.second + (it != sq.end() ? it->second : RealType(0)));
            nb[I].push_back({J, s}); nb[J].push_back({I, s});
        }
    }
    std::vector<RealType> wmax(nNodes, 0);
    for (int i = 0; i < nNodes; ++i) {
        for (auto& e : nb[i]) wmax[i] = std::max(wmax[i], e.second);
        std::sort(nb[i].begin(), nb[i].end(), [](auto& a, auto& b) { return a.second > b.second; });  // strongest first
    }
    auto strong = [&](int i, int j, RealType s) { return s > beta * wmax[i] && s > beta * wmax[j]; };

    agg.assign(nNodes, -1); nCoarse = 0;
    // Pass 1: seed only where a strong unaggregated neighbor exists (no singleton-seeding), fuse strongest-first
    for (int seed = 0; seed < nNodes; ++seed) {
        if (agg[seed] != -1) continue;
        bool canSeed = false;
        for (auto& e : nb[seed]) if (agg[e.first] == -1 && strong(seed, e.first, e.second)) { canSeed = true; break; }
        if (!canSeed) continue;
        agg[seed] = nCoarse; int cnt = 1;
        for (auto& e : nb[seed]) {
            if (cnt >= kmax) break;
            if (agg[e.first] == -1 && strong(seed, e.first, e.second)) { agg[e.first] = nCoarse; ++cnt; }
        }
        ++nCoarse;
    }
    // Pass 2: absorb every leftover node into its strongest neighboring aggregate (true isolates -> own)
    for (int i = 0; i < nNodes; ++i) {
        if (agg[i] != -1) continue;
        int best = -1;
        for (auto& e : nb[i]) if (agg[e.first] != -1) { best = agg[e.first]; break; }   // nb sorted by strength
        agg[i] = (best >= 0) ? best : nCoarse++;
    }
}

// ILU smoother selected by env (A/B in one binary, no signature threading): MARS_ACM_SMOOTHER=ilu.
inline bool acmUseIlu()
{
    static int v = -1;
    if (v < 0) { const char* e = std::getenv("MARS_ACM_SMOOTHER"); v = (e && std::strcmp(e, "ilu") == 0) ? 1 : 0; }
    return v == 1;
}

// multicolor block-DILU smoother (= AMGX's MULTICOLOR_DILU): MARS_ACM_SMOOTHER=dilu. Multicolor-parallel
// (numColors launches/sweep) where ILU's level schedule is launch-bound (nLev ~ thousands on SFC order).
// Host study (acm_ilu_smoother_study.py): multicolor ILU == DILU, ~3x weaker than natural-order but beats Jacobi.
inline bool acmUseDilu()
{
    static int v = -1;
    if (v < 0) { const char* e = std::getenv("MARS_ACM_SMOOTHER"); v = (e && std::strcmp(e, "dilu") == 0) ? 1 : 0; }
    return v == 1;
}

// Build the block-ILU(0) factorization + level schedule for one level (host, once per hierarchy level;
// reuses the host CSR already copied for aggregation). Sparsity is fixed across Picard iters, so the
// factor cost amortizes. Block unit = 4x4 nodal [u,v,w,p] block (pivots on the whole block via
// acmInvert4x4 -- scalar ILU would pivot on the bare indefinite a_pp).
template<typename RealType>
inline void acmBuildLevelIlu(AcmLevel<RealType>& L, const std::vector<int>& rowOff,
                             const std::vector<int>& colInd, const std::vector<RealType>& vals)
{
    const int nN = L.nNodes;
    // 1. node block-CSR: UNION of all 4 comp-rows' block-cols + guaranteed diagonal, then STRUCTURALLY
    // SYMMETRIZED (add every transpose edge). Symmetry makes ONE lower-graph level schedule valid for
    // both the L and U solves -- robust even if a convective/upwind stencil makes A structurally
    // nonsymmetric (the schedule reuse would otherwise read stale upper unknowns in the back solve).
    // On the current ring-symmetric operator the symmetrization is a no-op (zero extra blocks).
    std::vector<std::set<int>> cset(nN);
    for (int I = 0; I < nN; ++I) {
        cset[I].insert(I);                                       // diagonal block always present
        for (int a = 0; a < 4; ++a) { int r = 4 * I + a; for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) cset[I].insert(colInd[k] >> 2); }
    }
    for (int I = 0; I < nN; ++I) { std::vector<int> row(cset[I].begin(), cset[I].end()); for (int J : row) cset[J].insert(I); }
    std::vector<int> brow(nN + 1, 0), bcol, bdiagPtr(nN, -1);
    for (int I = 0; I < nN; ++I) brow[I + 1] = brow[I] + (int)cset[I].size();
    bcol.resize(brow[nN]);
    for (int I = 0; I < nN; ++I) { int p = brow[I]; for (int J : cset[I]) { bcol[p] = J; if (J == I) bdiagPtr[I] = p; ++p; } }
    auto slot = [&](int i, int j) -> int {
        int lo = brow[i], hi = brow[i + 1];
        while (lo < hi) { int m = (lo + hi) / 2; if (bcol[m] < j) lo = m + 1; else if (bcol[m] > j) hi = m; else return m; }
        return -1;
    };
    // 2. extract the 4x4 blocks into bilu (16/block, brow/bcol order)
    std::vector<RealType> bilu(16 * brow[nN], RealType(0));
    for (int I = 0; I < nN; ++I)
        for (int a = 0; a < 4; ++a) {
            int r = 4 * I + a;
            for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) {
                int J = colInd[k] >> 2, b = colInd[k] & 3, s = slot(I, J);
                if (s >= 0) bilu[16 * s + 4 * a + b] = vals[k];
            }
        }
    // 3. block-IKJ ILU(0) in place: L_ik = A_ik*inv(U_kk) below, inv(U_ii) on the diagonal, U above (no fill)
    RealType Lik[16], tmp[16], Aii[16];
    for (int i = 0; i < nN; ++i) {
        for (int kk = brow[i]; kk < bdiagPtr[i]; ++kk) {                  // k < i
            int k = bcol[kk];
            acmBlkMatmul4(&bilu[16 * kk], &bilu[16 * bdiagPtr[k]], Lik);  // L_ik = A_ik * inv(U_kk)
            for (int t = 0; t < 16; ++t) bilu[16 * kk + t] = Lik[t];
            for (int jj = bdiagPtr[k] + 1; jj < brow[k + 1]; ++jj) {      // upper cols j>k of row k
                int j = bcol[jj], sij = slot(i, j);
                if (sij >= 0) { acmBlkMatmul4(Lik, &bilu[16 * jj], tmp); for (int t = 0; t < 16; ++t) bilu[16 * sij + t] -= tmp[t]; }
            }
        }
        acmInvert4x4(&bilu[16 * bdiagPtr[i]], Aii);                       // store inv(U_ii)
        for (int t = 0; t < 16; ++t) bilu[16 * bdiagPtr[i] + t] = Aii[t];
    }
    // 4. level schedule on the lower block-graph: level[i] = 1 + max_{k<i stored} level[k]
    std::vector<int> level(nN, 0); int maxLev = 0;
    for (int i = 0; i < nN; ++i) {
        int lv = 0;
        for (int kk = brow[i]; kk < bdiagPtr[i]; ++kk) lv = std::max(lv, level[bcol[kk]] + 1);
        level[i] = lv; maxLev = std::max(maxLev, lv);
    }
    int nLev = maxLev + 1;
    // diagnostic: nLev = the serial depth of the level-scheduled triangular solve (2*nLev kernel launches
    // per sweep). Large nLev (SFC ordering) -> launch-bound, the case for multicolor reordering. Prints
    // per build, so it also shows whether the hierarchy is being reused across Picard iters.
    std::fprintf(stderr, "[acm-ilu] build nN=%d nLev=%d (avg level width %.1f)\n", nN, nLev, double(nN) / nLev);
    std::vector<int> levStart(nLev + 1, 0);
    for (int i = 0; i < nN; ++i) levStart[level[i] + 1]++;
    for (int l = 0; l < nLev; ++l) levStart[l + 1] += levStart[l];
    std::vector<int> nodesByLevel(nN), pos(levStart.begin(), levStart.end());
    for (int i = 0; i < nN; ++i) nodesByLevel[pos[level[i]]++] = i;
    // 5. upload
    L.brow.assign(brow.begin(), brow.end());
    L.bcol.assign(bcol.begin(), bcol.end());
    L.bdiagPtr.assign(bdiagPtr.begin(), bdiagPtr.end());
    L.bilu.assign(bilu.begin(), bilu.end());
    L.nodesByLevel.assign(nodesByLevel.begin(), nodesByLevel.end());
    L.h_levStart = levStart; L.nLev = nLev;
    L.ytmp.assign(L.ND, RealType(0));
}

// device binary search for block (i,j) in the symmetrized block-CSR row i -> slot, or -1.
__device__ __forceinline__ int acmBlockSlot(const int* brow, const int* bcol, int i, int j)
{
    int lo = brow[i], hi = brow[i + 1];
    while (lo < hi) { int m = (lo + hi) >> 1; int c = bcol[m]; if (c < j) lo = m + 1; else if (c > j) hi = m; else return m; }
    return -1;
}

// GPU value extract: scatter the level's scalar CSR values into the 4x4 block array bblk (one thread/node).
// Replaces the host slot()-loop over all scalar nnz (the dominant chunk of the 2.7s setup). bblk pre-zeroed.
template<typename RealType>
__global__ void acmDiluExtractKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                     const int* brow, const int* bcol, RealType* bblk, int nN)
{
    int I = blockIdx.x * blockDim.x + threadIdx.x; if (I >= nN) return;
    for (int a = 0; a < 4; ++a) {
        int r = 4 * I + a;
        for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) {
            int J = colInd[k] >> 2, b = colInd[k] & 3, s = acmBlockSlot(brow, bcol, I, J);
            if (s >= 0) bblk[16 * s + 4 * a + b] = vals[k];
        }
    }
}

// GPU block-DILU diagonal factor for ONE color: D_i = A_ii - sum_{color(k)<color(i)} A_ik inv(D_k) A_ki,
// store inv(D_i). Launched per color in ascending order on one stream -> every lower-color inv(D_k) is
// already written (the exact dependency the host serial color loop relied on). bblk read-only (raw A).
template<typename RealType>
__global__ void acmDiluFactorColorKernel(const int* brow, const int* bcol, const int* bdiagPtr,
        const int* color, const int* nodesByColor, const RealType* bblk, RealType* dinvD, int cstart, int cend)
{
    int t = cstart + blockIdx.x * blockDim.x + threadIdx.x; if (t >= cend) return;
    int i = nodesByColor[t], ci = color[i];
    RealType Di[16]; for (int q = 0; q < 16; ++q) Di[q] = bblk[16 * bdiagPtr[i] + q];   // A_ii
    for (int kk = brow[i]; kk < brow[i + 1]; ++kk) {
        int k = bcol[kk]; if (k == i || color[k] >= ci) continue;                       // strictly lower-color
        int kT = acmBlockSlot(brow, bcol, k, i); if (kT < 0) continue;                  // A_ki (sym -> present)
        RealType T1[16], T2[16];
        acmBlkMatmul4(&bblk[16 * kk], &dinvD[16 * k], T1);                              // A_ik inv(D_k)
        acmBlkMatmul4(T1, &bblk[16 * kT], T2);                                          // (A_ik inv(D_k)) A_ki
        for (int q = 0; q < 16; ++q) Di[q] -= T2[q];
    }
    acmInvert4x4(Di, &dinvD[16 * i]);                                                   // inv(D_i)
}

// Build the multicolor block-DILU smoother for one level. STRUCTURE (symmetrized node block-CSR + greedy
// coloring) is built host-side (value-independent given the sparsity); the VALUES (bblk extract + the
// color-order block-diagonal factor D_i = A_ii - sum_{color(k)<color(i)} A_ik inv(D_k) A_ki) run on the
// GPU above -- the host value loops were the dominant part of the 2.7s/Picard setup.
template<typename RealType>
inline void acmBuildLevelDilu(AcmLevel<RealType>& L, const std::vector<int>& rowOff,
                              const std::vector<int>& colInd, const std::vector<RealType>& vals)
{
    const int nN = L.nNodes;
    // STRUCTURE (host, value-independent given the sparsity): symmetrized node block-CSR (union of the 4
    // comp-rows + guaranteed diagonal + every transpose edge), then a greedy first-fit coloring.
    std::vector<std::set<int>> cset(nN);
    for (int I = 0; I < nN; ++I) {
        cset[I].insert(I);
        for (int a = 0; a < 4; ++a) { int r = 4 * I + a; for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) cset[I].insert(colInd[k] >> 2); }
    }
    for (int I = 0; I < nN; ++I) { std::vector<int> row(cset[I].begin(), cset[I].end()); for (int J : row) cset[J].insert(I); }
    std::vector<int> brow(nN + 1, 0), bcol, bdiagPtr(nN, -1);
    for (int I = 0; I < nN; ++I) brow[I + 1] = brow[I] + (int)cset[I].size();
    bcol.resize(brow[nN]);
    for (int I = 0; I < nN; ++I) { int p = brow[I]; for (int J : cset[I]) { bcol[p] = J; if (J == I) bdiagPtr[I] = p; ++p; } }
    std::vector<int> color(nN, -1), forbidden(nN, -1);
    int numColors = 0;
    for (int i = 0; i < nN; ++i) {
        for (int kk = brow[i]; kk < brow[i + 1]; ++kk) { int j = bcol[kk]; if (j != i && color[j] >= 0) forbidden[color[j]] = i; }
        int c = 0; while (c < numColors && forbidden[c] == i) ++c;
        if (c == numColors) ++numColors;
        color[i] = c;
    }
    std::vector<int> colStart(numColors + 1, 0);
    for (int i = 0; i < nN; ++i) colStart[color[i] + 1]++;
    for (int c = 0; c < numColors; ++c) colStart[c + 1] += colStart[c];
    std::vector<int> nodesByColor(nN), pos(colStart.begin(), colStart.end());
    for (int i = 0; i < nN; ++i) nodesByColor[pos[color[i]]++] = i;
    std::fprintf(stderr, "[acm-dilu] build nN=%d numColors=%d (avg color width %.1f)\n", nN, numColors, double(nN) / numColors);
    // upload the structure
    L.brow.assign(brow.begin(), brow.end());
    L.bcol.assign(bcol.begin(), bcol.end());
    L.bdiagPtr.assign(bdiagPtr.begin(), bdiagPtr.end());
    L.color.assign(color.begin(), color.end());
    L.nodesByColor.assign(nodesByColor.begin(), nodesByColor.end());
    L.h_colStart = colStart; L.numColors = numColors;
    L.ytmp.assign(L.ND, RealType(0)); L.ztmp.assign(L.ND, RealType(0));
    // VALUES on the GPU (the host loops were the 2.7s/Picard bottleneck): extract raw A blocks into bblk
    // from the level's device CSR, then color-parallel block-diagonal factor. The ascending-color factor
    // launches MUST stay on the default stream -- the cross-color inv(D_k)-before-use dependency (and
    // extract-before-factor) is enforced purely by default-stream serialization; an async stream without
    // explicit deps would let a color read unwritten lower-color inv(D_k) and silently corrupt it. bblk pre-zeroed.
    L.bblk.assign(16 * brow[nN], RealType(0));
    L.dilu_dinv.resize(16 * nN);
    const int blk = 256, gAll = (nN + blk - 1) / blk;
    acmDiluExtractKernel<RealType><<<gAll, blk>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                                                  acmRaw(L.brow), acmRaw(L.bcol), acmRaw(L.bblk), nN);
    for (int c = 0; c < numColors; ++c) {
        int cs = colStart[c], ce = colStart[c + 1], g = (ce - cs + blk - 1) / blk;
        if (g > 0) acmDiluFactorColorKernel<RealType><<<g, blk>>>(acmRaw(L.brow), acmRaw(L.bcol), acmRaw(L.bdiagPtr),
            acmRaw(L.color), acmRaw(L.nodesByColor), acmRaw(L.bblk), acmRaw(L.dilu_dinv), cs, ce);
    }
}

// one or more block-ILU(0) sweeps: x += omega * (LU)^-1 (b - A x), triangular solves level-scheduled.
// All launches on stream `s` (default 0) so the cycle can be captured into a CUDA graph (Step 1).
template<typename RealType>
inline void acmIluSmoothGpu(AcmLevel<RealType>& L, int sweeps, RealType omega, cudaStream_t s = 0)
{
    const int blk = 256, grd = (L.ND + blk - 1) / blk;
    for (int sw = 0; sw < sweeps; ++sw) {
        acmResidualKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                                                  acmRaw(L.xvec), acmRaw(L.bvec), acmRaw(L.rtmp), L.ND);   // rtmp = b - A x (fused)
        for (int lv = 0; lv < L.nLev; ++lv) {                            // forward L-solve -> ytmp (ascending)
            int ls = L.h_levStart[lv], le = L.h_levStart[lv + 1], g = (le - ls + blk - 1) / blk;
            if (g > 0) acmIluLsolveKernel<RealType><<<g, blk, 0, s>>>(acmRaw(L.brow), acmRaw(L.bcol), acmRaw(L.bdiagPtr),
                acmRaw(L.bilu), acmRaw(L.nodesByLevel), ls, le, acmRaw(L.rtmp), acmRaw(L.ytmp));
        }
        for (int lv = L.nLev - 1; lv >= 0; --lv) {                       // back U-solve -> xtmp = dx (descending)
            int ls = L.h_levStart[lv], le = L.h_levStart[lv + 1], g = (le - ls + blk - 1) / blk;
            if (g > 0) acmIluUsolveKernel<RealType><<<g, blk, 0, s>>>(acmRaw(L.brow), acmRaw(L.bcol), acmRaw(L.bdiagPtr),
                acmRaw(L.bilu), acmRaw(L.nodesByLevel), ls, le, acmRaw(L.ytmp), acmRaw(L.xtmp));
        }
        acmAxpyKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.ND);            // x += omega*dx
    }
}

// one or more multicolor block-DILU sweeps: x += omega * (LU)^-1 (b - A x), triangular solves color-parallel.
// M = (D+L) D^-1 (D+U). Forward (D+L)y=r: w_i = r_i - sum_{lower-color} A_ik y_k, then y_i = inv(D_i) w_i.
// Since z_i = D_i y_i = w_i, the back solve (D+U)x=z reads w (ztmp) directly -- exact ONLY because the
// forward step multiplies by EXACTLY inv(D_i) (a future forward-diagonal damping would break the z:=w shortcut).
template<typename RealType>
inline void acmDiluSmoothGpu(AcmLevel<RealType>& L, int sweeps, RealType omega, cudaStream_t s = 0)
{
    const int blk = 256, grd = (L.ND + blk - 1) / blk;
    for (int sw = 0; sw < sweeps; ++sw) {
        acmResidualKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                                                  acmRaw(L.xvec), acmRaw(L.bvec), acmRaw(L.rtmp), L.ND);   // rtmp = b - A x (fused)
        for (int c = 0; c < L.numColors; ++c) {                            // forward L-solve -> ytmp (y), pre-diag w -> ztmp
            int cs = L.h_colStart[c], ce = L.h_colStart[c + 1], g = (ce - cs + blk - 1) / blk;
            if (g > 0) acmDiluLsolveKernel<RealType><<<g, blk, 0, s>>>(acmRaw(L.brow), acmRaw(L.bcol), acmRaw(L.color),
                acmRaw(L.bblk), acmRaw(L.dilu_dinv), acmRaw(L.nodesByColor), cs, ce, acmRaw(L.rtmp), acmRaw(L.ytmp), acmRaw(L.ztmp));
        }
        for (int c = L.numColors - 1; c >= 0; --c) {                       // back U-solve -> xtmp (dx), RHS = ztmp (w)
            int cs = L.h_colStart[c], ce = L.h_colStart[c + 1], g = (ce - cs + blk - 1) / blk;
            if (g > 0) acmDiluUsolveKernel<RealType><<<g, blk, 0, s>>>(acmRaw(L.brow), acmRaw(L.bcol), acmRaw(L.color),
                acmRaw(L.bblk), acmRaw(L.dilu_dinv), acmRaw(L.nodesByColor), cs, ce, acmRaw(L.ztmp), acmRaw(L.xtmp));
        }
        acmAxpyKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.ND);            // x += omega*dx
    }
}

template<typename RealType>
void acmBuildHierarchy(const thrust::device_vector<int>& rowOff0,
                       const thrust::device_vector<int>& colInd0,
                       const thrust::device_vector<RealType>& vals0, int nNodes0,
                       int kmax, int maxCoarseND, std::vector<AcmLevel<RealType>>& levels,
                       RealType beta = RealType(0.5))
{
    levels.clear();
    {
        AcmLevel<RealType> L;
        L.nNodes = nNodes0; L.ND = 4 * nNodes0;
        L.rowOff = rowOff0; L.colInd = colInd0; L.vals = vals0;
        levels.push_back(std::move(L));
    }
    for (int idx = 0; ; ++idx) {
        const int blk = 256, grd = (levels[idx].ND + blk - 1) / blk;
        levels[idx].dinv.resize(levels[idx].ND);
        acmDinvKernel<RealType><<<grd, blk>>>(acmRaw(levels[idx].rowOff), acmRaw(levels[idx].colInd),
                                              acmRaw(levels[idx].vals), acmRaw(levels[idx].dinv), levels[idx].ND);
        levels[idx].binv.resize(16 * levels[idx].nNodes);     // per-node 4x4 inverse (coupled block smoother)
        { const int gN = (levels[idx].nNodes + blk - 1) / blk;
          acmBlockInvKernel<RealType><<<gN, blk>>>(acmRaw(levels[idx].rowOff), acmRaw(levels[idx].colInd),
              acmRaw(levels[idx].vals), acmRaw(levels[idx].binv), levels[idx].nNodes); }
        levels[idx].bvec.assign(levels[idx].ND, 0); levels[idx].xvec.assign(levels[idx].ND, 0);
        levels[idx].rtmp.assign(levels[idx].ND, 0); levels[idx].xtmp.assign(levels[idx].ND, 0);
        if (levels[idx].ND <= maxCoarseND) break;

        std::vector<int> hro(levels[idx].rowOff.size()), hci(levels[idx].colInd.size());
        std::vector<RealType> hva(levels[idx].vals.size());
        thrust::copy(levels[idx].rowOff.begin(), levels[idx].rowOff.end(), hro.begin());
        thrust::copy(levels[idx].colInd.begin(), levels[idx].colInd.end(), hci.begin());
        thrust::copy(levels[idx].vals.begin(), levels[idx].vals.end(), hva.begin());
        std::vector<int> agg; int nCoarse;
        acmDirectionalAggregate<RealType>(hro, hci, hva, levels[idx].nNodes, kmax, beta, agg, nCoarse);
        if (nCoarse >= levels[idx].nNodes) break;            // no coarsening progress -> coarsest (no ILU, nLev=0)
        if (acmUseIlu()) acmBuildLevelIlu<RealType>(levels[idx], hro, hci, hva);        // smoothed level only -> build ILU
        else if (acmUseDilu()) acmBuildLevelDilu<RealType>(levels[idx], hro, hci, hva); // or multicolor block-DILU

        levels[idx].agg.assign(agg.begin(), agg.end());
        levels[idx].nCoarse = nCoarse;
        AcmLevel<RealType> C; C.nNodes = nCoarse; C.ND = 4 * nCoarse;
        buildCoarseOperator<RealType>(levels[idx].rowOff, levels[idx].colInd, levels[idx].vals,
                                      levels[idx].agg, levels[idx].nNodes, nCoarse,
                                      C.rowOff, C.colInd, C.vals);
        levels.push_back(std::move(C));
    }
}

// useBlock=true -> coupled per-node 4x4 block-Jacobi (captures the local [u,v,w,p] coupling);
// false -> segregated point-Jacobi (kept for the Stage-4a host-match gate).
template<typename RealType>
inline void acmSmoothGpu(AcmLevel<RealType>& L, int sweeps, RealType omega, bool useBlock, cudaStream_t s = 0)
{
    if (acmUseIlu() && L.nLev > 0) { acmIluSmoothGpu(L, sweeps, omega, s); return; }        // block-ILU(0) (coarsest nLev=0 -> Jacobi)
    if (acmUseDilu() && L.numColors > 0) { acmDiluSmoothGpu(L, sweeps, omega, s); return; } // multicolor block-DILU (coarsest -> Jacobi)
    const int blk = 256, grd = (L.ND + blk - 1) / blk, gN = (L.nNodes + blk - 1) / blk;
    for (int sw = 0; sw < sweeps; ++sw) {
        if (useBlock)
            acmBlockJacobiKernel<RealType><<<gN, blk, 0, s>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                acmRaw(L.binv), acmRaw(L.bvec), acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.nNodes);
        else
            acmJacobiKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                acmRaw(L.dinv), acmRaw(L.bvec), acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.ND);
        L.xvec.swap(L.xtmp);                                 // O(1) host pointer swap -- NOT in the captured DILU path (returns above)
    }
}

// coarsest direct solve (Algorithm A): A_c x = b via cusolverSp QR. Replaces the non-contractive Jacobi
// sweeps that poison a deep hierarchy. The coarse op is indefinite -> QR (not LU). cs==nullptr keeps the
// Jacobi-sweep coarsest (so the Stage-4a host-vs-GPU V-cycle match still holds bit-for-bit).
template<typename RealType>
inline void acmCoarseSolve(AcmLevel<RealType>& L, int coarse, RealType omega, bool useBlock,
                           cusolverSpHandle_t cs, cusparseMatDescr_t descr, cudaStream_t s = 0)
{
    if (!(cs && descr)) { acmSmoothGpu(L, coarse, omega, useBlock, s); return; }
    cusolverSpSetStream(cs, s);                            // QR on the cycle stream -> orders with the captured graphs
    int singularity = -1;                                  // <0 = full rank; >=0 = first deficient pivot
    const int nnzc = static_cast<int>(L.vals.size());
    if constexpr (std::is_same_v<RealType, double>)
        cusolverSpDcsrlsvqr(cs, L.ND, nnzc, descr, acmRaw(L.vals), acmRaw(L.rowOff), acmRaw(L.colInd),
                            acmRaw(L.bvec), 1e-12, 1, acmRaw(L.xvec), &singularity);   // reorder=1 (match ref, stabilize)
    else
        cusolverSpScsrlsvqr(cs, L.ND, nnzc, descr, acmRaw(L.vals), acmRaw(L.rowOff), acmRaw(L.colInd),
                            acmRaw(L.bvec), 1e-6f, 1, acmRaw(L.xvec), &singularity);
    cudaStreamSynchronize(s);                              // QR is non-capturable + needs the singularity readback
    // a rank-deficient coarse saddle block (near-null pressure mode) makes the QR least-squares solve
    // unreliable -> don't inject garbage into the cycle; reset and fall back to smoothing.
    if (singularity >= 0) {
        thrust::fill(thrust::cuda::par.on(s), L.xvec.begin(), L.xvec.end(), RealType(0));
        acmSmoothGpu(L, coarse, omega, useBlock, s);
    }
}

// V-cycle DOWN half: pre-smooth + residual + restrict, level Lidx down to just above the coarsest.
// Pure capturable kernel sequence (no host sync, DILU smoother has no buffer swap) -> CUDA-graph safe.
template<typename RealType>
void acmVcycleDownGpu(std::vector<AcmLevel<RealType>>& levels, int Lidx, int pre, RealType omega,
                      bool useBlock, cudaStream_t s)
{
    const int blk = 256, nL = (int)levels.size();
    for (int l = Lidx; l < nL - 1; ++l) {
        AcmLevel<RealType>& L = levels[l];
        AcmLevel<RealType>& C = levels[l + 1];
        const int grd = (L.ND + blk - 1) / blk, gN = (L.nNodes + blk - 1) / blk;
        acmSmoothGpu(L, pre, omega, useBlock, s);
        acmResidualKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                                                  acmRaw(L.xvec), acmRaw(L.bvec), acmRaw(L.rtmp), L.ND);   // rtmp = b - A x (fused)
        cudaMemsetAsync(acmRaw(C.bvec), 0, sizeof(RealType) * C.ND, s);   // capture-safe zero (all-zero bytes == 0.0)
        cudaMemsetAsync(acmRaw(C.xvec), 0, sizeof(RealType) * C.ND, s);
        acmRestrictAddKernel<RealType><<<gN, blk, 0, s>>>(acmRaw(L.rtmp), acmRaw(C.bvec), acmRaw(L.agg), L.nNodes);
    }
}

// V-cycle UP half: inject + add + post-smooth, just above the coarsest back up to level Lidx. Capturable.
template<typename RealType>
void acmVcycleUpGpu(std::vector<AcmLevel<RealType>>& levels, int Lidx, int post, RealType omega,
                    bool useBlock, cudaStream_t s)
{
    const int blk = 256, nL = (int)levels.size();
    for (int l = nL - 2; l >= Lidx; --l) {
        AcmLevel<RealType>& L = levels[l];
        AcmLevel<RealType>& C = levels[l + 1];
        const int grd = (L.ND + blk - 1) / blk, gN = (L.nNodes + blk - 1) / blk;
        acmInjectKernel<RealType><<<gN, blk, 0, s>>>(acmRaw(C.xvec), acmRaw(L.xtmp), acmRaw(L.agg), L.nNodes);
        acmAddKernel<RealType><<<grd, blk, 0, s>>>(acmRaw(L.xvec), acmRaw(L.xtmp), L.ND);
        acmSmoothGpu(L, post, omega, useBlock, s);
    }
}

// one V-cycle: down -> coarse QR -> up (equivalent to the old recursion, same ops/order). The down/up
// halves are separate so the preconditioner can CUDA-graph-capture each around the non-capturable QR.
template<typename RealType>
void acmVcycleGpu(std::vector<AcmLevel<RealType>>& levels, int Lidx,
                  int pre, int post, int coarse, RealType omega, bool useBlock = false,
                  cusolverSpHandle_t cs = nullptr, cusparseMatDescr_t descr = nullptr, cudaStream_t s = 0)
{
    const int nL = (int)levels.size();
    acmVcycleDownGpu(levels, Lidx, pre, omega, useBlock, s);
    acmCoarseSolve(levels[nL - 1], coarse, omega, useBlock, cs, descr, s);   // eager (non-capturable)
    acmVcycleUpGpu(levels, Lidx, post, omega, useBlock, s);
}

// ---- host replica (validation oracle): same hierarchy, same arithmetic, on the CPU ----
template<typename RealType>
struct AcmLevelHost {
    int nNodes, ND, nCoarse;
    std::vector<int> rowOff, colInd, agg;
    std::vector<RealType> vals, dinv;
};

template<typename RealType>
void acmHierarchyToHost(const std::vector<AcmLevel<RealType>>& g, std::vector<AcmLevelHost<RealType>>& h)
{
    h.resize(g.size());
    for (size_t l = 0; l < g.size(); ++l) {
        h[l].nNodes = g[l].nNodes; h[l].ND = g[l].ND; h[l].nCoarse = g[l].nCoarse;
        h[l].rowOff.resize(g[l].rowOff.size()); h[l].colInd.resize(g[l].colInd.size());
        h[l].vals.resize(g[l].vals.size());     h[l].dinv.resize(g[l].dinv.size());
        h[l].agg.resize(g[l].agg.size());
        thrust::copy(g[l].rowOff.begin(), g[l].rowOff.end(), h[l].rowOff.begin());
        thrust::copy(g[l].colInd.begin(), g[l].colInd.end(), h[l].colInd.begin());
        thrust::copy(g[l].vals.begin(),   g[l].vals.end(),   h[l].vals.begin());
        thrust::copy(g[l].dinv.begin(),   g[l].dinv.end(),   h[l].dinv.begin());
        thrust::copy(g[l].agg.begin(),    g[l].agg.end(),    h[l].agg.begin());
    }
}

template<typename RealType>
void acmHostVcycle(const std::vector<AcmLevelHost<RealType>>& Ls, int L,
                   std::vector<RealType>& x, const std::vector<RealType>& b,
                   int pre, int post, int coarse, RealType omega)
{
    const AcmLevelHost<RealType>& lv = Ls[L];
    auto smooth = [&](std::vector<RealType>& xx, const std::vector<RealType>& bb, int sw) {
        std::vector<RealType> t(lv.ND);
        for (int s = 0; s < sw; ++s) {
            for (int r = 0; r < lv.ND; ++r) {
                RealType sm = 0;
                for (int k = lv.rowOff[r]; k < lv.rowOff[r + 1]; ++k) sm += lv.vals[k] * xx[lv.colInd[k]];
                t[r] = xx[r] + omega * lv.dinv[r] * (bb[r] - sm);
            }
            xx.swap(t);
        }
    };
    if (L == (int)Ls.size() - 1) { smooth(x, b, coarse); return; }
    smooth(x, b, pre);
    std::vector<RealType> r(lv.ND);
    for (int i = 0; i < lv.ND; ++i) {
        RealType sm = 0;
        for (int k = lv.rowOff[i]; k < lv.rowOff[i + 1]; ++k) sm += lv.vals[k] * x[lv.colInd[k]];
        r[i] = b[i] - sm;
    }
    const AcmLevelHost<RealType>& cv = Ls[L + 1];
    std::vector<RealType> bc(cv.ND, 0), xc(cv.ND, 0);
    for (int i = 0; i < lv.nNodes; ++i) for (int c = 0; c < 4; ++c) bc[4 * lv.agg[i] + c] += r[4 * i + c];
    acmHostVcycle(Ls, L + 1, xc, bc, pre, post, coarse, omega);
    for (int i = 0; i < lv.nNodes; ++i) for (int c = 0; c < 4; ++c) x[4 * i + c] += xc[4 * lv.agg[i] + c];
    smooth(x, b, post);
}

}  // namespace mars

#endif  // MARS_ACM_VCYCLE_HPP
