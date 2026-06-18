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
#include <vector>
#include <set>
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
    thrust::device_vector<RealType> vals, dinv;
    thrust::device_vector<int> agg;                            // nNodes -> parent (empty on coarsest)
    thrust::device_vector<RealType> bvec, xvec, rtmp, xtmp;    // per-level cycle work vectors
};

// host-once greedy aggregation on the block graph of a coupled (4*node+comp) CSR
inline void acmHostAggregate(const std::vector<int>& rowOff, const std::vector<int>& colInd,
                             int nNodes, int kmax, std::vector<int>& agg, int& nCoarse)
{
    const int ND = 4 * nNodes;
    std::vector<std::set<int>> nb(nNodes);
    for (int r = 0; r < ND; ++r) {
        int I = r >> 2;
        for (int k = rowOff[r]; k < rowOff[r + 1]; ++k) nb[I].insert(colInd[k] >> 2);
    }
    agg.assign(nNodes, -1); nCoarse = 0;
    for (int i = 0; i < nNodes; ++i) {
        if (agg[i] != -1) continue;
        agg[i] = nCoarse; int cnt = 1;
        for (int j : nb[i]) { if (cnt >= kmax) break; if (j != i && agg[j] == -1) { agg[j] = nCoarse; ++cnt; } }
        ++nCoarse;
    }
}

template<typename RealType>
void acmBuildHierarchy(const thrust::device_vector<int>& rowOff0,
                       const thrust::device_vector<int>& colInd0,
                       const thrust::device_vector<RealType>& vals0, int nNodes0,
                       int kmax, int maxCoarseND, std::vector<AcmLevel<RealType>>& levels)
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
        levels[idx].bvec.assign(levels[idx].ND, 0); levels[idx].xvec.assign(levels[idx].ND, 0);
        levels[idx].rtmp.assign(levels[idx].ND, 0); levels[idx].xtmp.assign(levels[idx].ND, 0);
        if (levels[idx].ND <= maxCoarseND) break;

        std::vector<int> hro(levels[idx].rowOff.size()), hci(levels[idx].colInd.size());
        thrust::copy(levels[idx].rowOff.begin(), levels[idx].rowOff.end(), hro.begin());
        thrust::copy(levels[idx].colInd.begin(), levels[idx].colInd.end(), hci.begin());
        std::vector<int> agg; int nCoarse;
        acmHostAggregate(hro, hci, levels[idx].nNodes, kmax, agg, nCoarse);
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

template<typename RealType>
inline void acmSmoothGpu(AcmLevel<RealType>& L, int sweeps, RealType omega)
{
    const int blk = 256, grd = (L.ND + blk - 1) / blk;
    for (int s = 0; s < sweeps; ++s) {
        acmJacobiKernel<RealType><<<grd, blk>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
            acmRaw(L.dinv), acmRaw(L.bvec), acmRaw(L.xvec), acmRaw(L.xtmp), omega, L.ND);
        L.xvec.swap(L.xtmp);                                 // O(1) pointer swap; xvec holds the new iterate
    }
}

// one V-cycle on levels[Lidx], operating on its bvec (in) / xvec (in-out)
template<typename RealType>
void acmVcycleGpu(std::vector<AcmLevel<RealType>>& levels, int Lidx,
                  int pre, int post, int coarse, RealType omega)
{
    AcmLevel<RealType>& L = levels[Lidx];
    const int blk = 256, grd = (L.ND + blk - 1) / blk;
    if (Lidx == (int)levels.size() - 1) { acmSmoothGpu(L, coarse, omega); return; }

    acmSmoothGpu(L, pre, omega);
    acmSpmvKernel<RealType><<<grd, blk>>>(acmRaw(L.rowOff), acmRaw(L.colInd), acmRaw(L.vals),
                                          acmRaw(L.xvec), acmRaw(L.rtmp), L.ND);
    acmSubKernel<RealType><<<grd, blk>>>(acmRaw(L.bvec), acmRaw(L.rtmp), acmRaw(L.rtmp), L.ND);

    AcmLevel<RealType>& C = levels[Lidx + 1];
    thrust::fill(C.bvec.begin(), C.bvec.end(), RealType(0));
    thrust::fill(C.xvec.begin(), C.xvec.end(), RealType(0));
    const int gN = (L.nNodes + blk - 1) / blk;
    acmRestrictAddKernel<RealType><<<gN, blk>>>(acmRaw(L.rtmp), acmRaw(C.bvec), acmRaw(L.agg), L.nNodes);

    acmVcycleGpu(levels, Lidx + 1, pre, post, coarse, omega);

    acmInjectKernel<RealType><<<gN, blk>>>(acmRaw(C.xvec), acmRaw(L.xtmp), acmRaw(L.agg), L.nNodes);
    acmAddKernel<RealType><<<grd, blk>>>(acmRaw(L.xvec), acmRaw(L.xtmp), L.ND);
    acmSmoothGpu(L, post, omega);
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
