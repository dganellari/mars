#ifndef MARS_ACM_COARSEN_HPP
#define MARS_ACM_COARSEN_HPP

// ACM Stage 3: GPU coarse-operator + grid-transfer kernels for the agglomeration
// additive-correction multigrid (Darwish-Saad-Hamdan 2008).
//
// The coarse operator is the sum-of-fine Galerkin operator (paper Eqs 27-28): for a 0/1 block
// aggregate-membership P, A_coarse == P^T A P. The host replica (scripts/acm_stage0_galerkin_gate.py)
// proved this to machine precision. We build it the GPU way: map every fine nonzero (r,c) to its
// coarse (cr,cc), sort by a linear key, reduce_by_key to sum duplicates, split back to CSR. Internal
// faces (both endpoints in the same aggregate) fold into the coarse diagonal for free, because their
// (cr,cc) lands on a diagonal block. This is the same sort/reduce_by_key idiom the unstructured
// backend already uses for CSR construction.
//
// DOF layout is interleaved 4*node+comp (the coupled [u,v,w,p] block). Aggregation is NODE-granular
// (all 4 components of a node share an aggregate), so the coarse dof of fine dof r is
// 4*agg[r/4] + (r%4). The aggregation map agg[] is an INPUT here (host-once for now; the GPU
// directional aggregation is Stage 5).

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/counting_iterator.h>

namespace mars {

// per fine row r: coarse row cr = 4*agg[r/4] + r%4; emit key = NDc*cr + cc for each nz in the row.
template<typename RealType>
__global__ void acmEmitCoarseKeyKernel(const int* frowOff, const int* fcolInd, const int* agg,
                                        int NDf, long long NDc, long long* key)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= NDf) return;
    long long cr = 4LL * agg[r >> 2] + (r & 3);
    for (int k = frowOff[r]; k < frowOff[r + 1]; ++k) {
        int c = fcolInd[k];
        long long cc = 4LL * agg[c >> 2] + (c & 3);
        key[k] = NDc * cr + cc;
    }
}

// templated only for weak linkage (header included in multiple TUs); RealType is unused
template<typename RealType>
__global__ void acmSplitKeyKernel(const long long* ukey, long long NDc, int cnnz,
                                   int* crow, int* ccol)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= cnnz) return;
    crow[i] = static_cast<int>(ukey[i] / NDc);
    ccol[i] = static_cast<int>(ukey[i] % NDc);
}

// order-zero (injection) prolongation: xf[4*node+c] = xc[4*agg[node]+c]
template<typename RealType>
__global__ void acmInjectKernel(const RealType* xc, RealType* xf, const int* agg, int nFine)
{
    int node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= nFine) return;
    int I = agg[node];
    for (int c = 0; c < 4; ++c) xf[4 * node + c] = xc[4 * I + c];
}

// sum restriction (P^T): xc[4*agg[node]+c] += xf[4*node+c]  (xc must be pre-zeroed)
template<typename RealType>
__global__ void acmRestrictAddKernel(const RealType* xf, RealType* xc, const int* agg, int nFine)
{
    int node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= nFine) return;
    int I = agg[node];
    for (int c = 0; c < 4; ++c) atomicAdd(&xc[4 * I + c], xf[4 * node + c]);
}

// Build the coarse CSR (NDc=4*nCoarse) from the fine CSR (NDf=4*nFine) + aggregation map.
// Output CSR is column-sorted within each row on device (the linear key sorts cc ascending
// within a fixed cr), so no host sortColumns round-trip is needed.
template<typename RealType>
void buildCoarseOperator(const thrust::device_vector<int>& d_frowOff,
                         const thrust::device_vector<int>& d_fcolInd,
                         const thrust::device_vector<RealType>& d_fvals,
                         const thrust::device_vector<int>& d_agg,   // size nFine
                         int nFine, int nCoarse,
                         thrust::device_vector<int>& d_crowOff,
                         thrust::device_vector<int>& d_ccolInd,
                         thrust::device_vector<RealType>& d_cvals)
{
    const int NDf = 4 * nFine;
    const long long NDc = 4LL * nCoarse;
    const int nnzF = static_cast<int>(d_fvals.size());

    // 1) one key per fine nonzero = NDc*coarse_row + coarse_col
    thrust::device_vector<long long> d_key(nnzF);
    thrust::device_vector<RealType>  d_val(d_fvals);            // values move with their key
    {
        int blk = 256, grd = (NDf + blk - 1) / blk;
        acmEmitCoarseKeyKernel<RealType><<<grd, blk>>>(
            thrust::raw_pointer_cast(d_frowOff.data()),
            thrust::raw_pointer_cast(d_fcolInd.data()),
            thrust::raw_pointer_cast(d_agg.data()), NDf, NDc,
            thrust::raw_pointer_cast(d_key.data()));
        cudaDeviceSynchronize();
    }

    // 2) sort by key, 3) sum duplicates (the segmented reduce = Eqs 27-28)
    thrust::sort_by_key(d_key.begin(), d_key.end(), d_val.begin());
    thrust::device_vector<long long> d_ukey(nnzF);
    thrust::device_vector<RealType>  d_uval(nnzF);
    auto end = thrust::reduce_by_key(d_key.begin(), d_key.end(), d_val.begin(),
                                     d_ukey.begin(), d_uval.begin());
    int cnnz = static_cast<int>(end.first - d_ukey.begin());
    d_ukey.resize(cnnz); d_uval.resize(cnnz);

    // 4) split unique keys -> (coarse_row, coarse_col); values are the coarse vals
    d_ccolInd.resize(cnnz);
    d_cvals = d_uval;
    thrust::device_vector<int> d_crow(cnnz);
    {
        int blk = 256, grd = (cnnz + blk - 1) / blk;
        acmSplitKeyKernel<RealType><<<grd, blk>>>(thrust::raw_pointer_cast(d_ukey.data()), NDc, cnnz,
                                        thrust::raw_pointer_cast(d_crow.data()),
                                        thrust::raw_pointer_cast(d_ccolInd.data()));
        cudaDeviceSynchronize();
    }

    // rowOff via lower_bound of each row id in the sorted coarse-row array (keys are sorted, so
    // d_crow is non-decreasing): rowOff[r] = #entries with crow < r.
    d_crowOff.resize(NDc + 1);
    thrust::counting_iterator<int> rows(0);
    thrust::lower_bound(d_crow.begin(), d_crow.end(), rows, rows + (NDc + 1), d_crowOff.begin());
}

}  // namespace mars

#endif  // MARS_ACM_COARSEN_HPP
