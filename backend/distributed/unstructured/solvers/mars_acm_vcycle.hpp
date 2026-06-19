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
#include <cusolverSp.h>
#include <cusparse.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#include <cmath>

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

template<typename RealType>
__global__ void acmSubKernel(const RealType* a, const RealType* b, RealType* y, int ND)
{ int r = blockIdx.x * blockDim.x + threadIdx.x; if (r < ND) y[r] = a[r] - b[r]; }

template<typename RealType>
__global__ void acmAddKernel(RealType* x, const RealType* y, int ND)
{ int r = blockIdx.x * blockDim.x + threadIdx.x; if (r < ND) x[r] += y[r]; }

template<typename RealType>
struct AcmLevel {
    int nNodes = 0, ND = 0, nCoarse = 0;
    thrust::device_vector<int> rowOff, colInd;
    thrust::device_vector<RealType> vals, dinv, binv;          // dinv: point-Jacobi; binv: 4x4 block-Jacobi
    thrust::device_vector<int> agg;                            // nNodes -> parent (empty on coarsest)
    thrust::device_vector<RealType> bvec, xvec, rtmp, xtmp;    // per-level cycle work vectors
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
        if (nCoarse >= levels[idx].nNodes) break;            // no coarsening progress

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
inline void acmSmoothGpu(AcmLevel<RealType>& L, int sweeps, RealType omega, bool useBlock)
{
    const int blk = 256, grd = (L.ND + blk - 1) / blk, gN = (L.nNodes + blk - 1) / blk;
    for (int s = 0; s < sweeps; ++s) {
        if (useBlock)
            acmBlockJacobiKernel<RealType><<<gN, blk>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                acmRaw(L.binv), acmRaw(L.bvec), acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.nNodes);
        else
            acmJacobiKernel<RealType><<<grd, blk>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                acmRaw(L.dinv), acmRaw(L.bvec), acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.ND);
        L.xvec.swap(L.xtmp);                                 // O(1) pointer swap; xvec holds the new iterate
    }
}

// coarsest direct solve (Algorithm A): A_c x = b via cusolverSp QR. Replaces the non-contractive Jacobi
// sweeps that poison a deep hierarchy. The coarse op is indefinite -> QR (not LU). cs==nullptr keeps the
// Jacobi-sweep coarsest (so the Stage-4a host-vs-GPU V-cycle match still holds bit-for-bit).
template<typename RealType>
inline void acmCoarseSolve(AcmLevel<RealType>& L, int coarse, RealType omega, bool useBlock,
                           cusolverSpHandle_t cs, cusparseMatDescr_t descr)
{
    if (!(cs && descr)) { acmSmoothGpu(L, coarse, omega, useBlock); return; }
    int singularity = -1;                                  // <0 = full rank; >=0 = first deficient pivot
    const int nnzc = static_cast<int>(L.vals.size());
    if constexpr (std::is_same_v<RealType, double>)
        cusolverSpDcsrlsvqr(cs, L.ND, nnzc, descr, acmRaw(L.vals), acmRaw(L.rowOff), acmRaw(L.colInd),
                            acmRaw(L.bvec), 1e-12, 1, acmRaw(L.xvec), &singularity);   // reorder=1 (match ref, stabilize)
    else
        cusolverSpScsrlsvqr(cs, L.ND, nnzc, descr, acmRaw(L.vals), acmRaw(L.rowOff), acmRaw(L.colInd),
                            acmRaw(L.bvec), 1e-6f, 1, acmRaw(L.xvec), &singularity);
    cudaDeviceSynchronize();
    // a rank-deficient coarse saddle block (near-null pressure mode) makes the QR least-squares solve
    // unreliable -> don't inject garbage into the cycle; reset and fall back to smoothing.
    if (singularity >= 0) {
        thrust::fill(L.xvec.begin(), L.xvec.end(), RealType(0));
        acmSmoothGpu(L, coarse, omega, useBlock);
    }
}

// one V-cycle on levels[Lidx], operating on its bvec (in) / xvec (in-out). useBlock -> coupled block-Jacobi.
template<typename RealType>
void acmVcycleGpu(std::vector<AcmLevel<RealType>>& levels, int Lidx,
                  int pre, int post, int coarse, RealType omega, bool useBlock = false,
                  cusolverSpHandle_t cs = nullptr, cusparseMatDescr_t descr = nullptr)
{
    AcmLevel<RealType>& L = levels[Lidx];
    const int blk = 256, grd = (L.ND + blk - 1) / blk;
    if (Lidx == (int)levels.size() - 1) { acmCoarseSolve(L, coarse, omega, useBlock, cs, descr); return; }

    acmSmoothGpu(L, pre, omega, useBlock);
    acmSpmvKernel<RealType><<<grd, blk>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                                          acmRaw(L.xvec), acmRaw(L.rtmp), L.ND);
    acmSubKernel<RealType><<<grd, blk>>>(acmRaw(L.bvec), acmRaw(L.rtmp), acmRaw(L.rtmp), L.ND);

    AcmLevel<RealType>& C = levels[Lidx + 1];
    thrust::fill(C.bvec.begin(), C.bvec.end(), RealType(0));
    thrust::fill(C.xvec.begin(), C.xvec.end(), RealType(0));
    const int gN = (L.nNodes + blk - 1) / blk;
    acmRestrictAddKernel<RealType><<<gN, blk>>>(acmRaw(L.rtmp), acmRaw(C.bvec), acmRaw(L.agg), L.nNodes);

    acmVcycleGpu(levels, Lidx + 1, pre, post, coarse, omega, useBlock, cs, descr);

    acmInjectKernel<RealType><<<gN, blk>>>(acmRaw(C.xvec), acmRaw(L.xtmp), acmRaw(L.agg), L.nNodes);
    acmAddKernel<RealType><<<grd, blk>>>(acmRaw(L.xvec), acmRaw(L.xtmp), L.ND);
    acmSmoothGpu(L, post, omega, useBlock);
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
