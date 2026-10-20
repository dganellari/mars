#pragma once

#include "../utils/mars_gpu_runtime.hpp"
#include "../utils/mars_gpu_mma.cuh"
#include <cstddef>

// AMD path for P=7: the same D^T G D operator as the staged DMMA kernel.
// Two elements share a wavefront. Keeping only u/y and three gradients uses exactly 32 KiB
// per block; MFMA operands and results go directly between registers and those arrays.
namespace mars::fem
{
struct HoMfmaDirectData
{
    double u[512];
    double g[3][512];
};
static_assert(2 * sizeof(HoMfmaDirectData) == 32768);

// Logical output coordinates for each of the three contractions, k fastest.
template<int Direction>
__host__ __device__ constexpr int mfmaTensorIndex(int batch, int row, int col)
{
    static_assert(Direction >= 0 && Direction < 3);
    if constexpr (Direction == 0) return row * 64 + batch * 8 + col;
    else return batch * 64 + row * 8 + col;
}

#if defined(MARS_HAVE_FP64_MFMA)
template<int Direction, bool Transpose>
__device__ __forceinline__ void mfmaDirectContract(const double* D, const double* in0, const double* in1,
                                                   double* out0, double* out1, int batch)
{
    static_assert(Direction >= 0 && Direction < 3);
    const int lane = threadIdx.x % 64;
    const int q = lane / 16, block = (lane % 16) / 4, t = lane % 4;
    const int row = 4 * (block / 2) + t, col = 4 * (block % 2) + t;
    double acc0 = 0, acc1 = 0;
    #pragma unroll
    for (int kt = 0; kt < 2; ++kt)
    {
        const int s = 4 * kt + q;
        double a0, a1, b0, b1;
        if constexpr (Direction < 2)
        {
            a0 = a1 = D[Transpose ? s * 8 + row : row * 8 + s];
            b0 = in0[mfmaTensorIndex<Direction>(batch, s, col)];
            b1 = in1[mfmaTensorIndex<Direction>(batch, s, col)];
        }
        else
        {
            a0 = in0[mfmaTensorIndex<Direction>(batch, row, s)];
            a1 = in1[mfmaTensorIndex<Direction>(batch, row, s)];
            b0 = b1 = D[Transpose ? s * 8 + col : col * 8 + s];
        }
        acc0 = __builtin_amdgcn_mfma_f64_4x4x4f64(a0, b0, acc0, 0, 0, 0);
        acc1 = __builtin_amdgcn_mfma_f64_4x4x4f64(a1, b1, acc1, 0, 0, 0);
    }
    const int i = mfmaTensorIndex<Direction>(batch, 4 * (block / 2) + q, col);
    if constexpr (Transpose)
    {
        out0[i] += acc0;
        out1[i] += acc1;
    }
    else
    {
        out0[i] = acc0;
        out1[i] = acc1;
    }
    syncWarp();
}
#endif

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
template<int GMode>
__global__ void __launch_bounds__(64)
ho_laplacian_apply_mfma_direct(const double* __restrict__ u, double* __restrict__ y,
                               const int* __restrict__ elem_dof, const double* __restrict__ D,
                               const double* __restrict__ G, const double* __restrict__ weights,
                               size_t num_elements)
{
#if defined(MARS_HAVE_FP64_MFMA)
    __shared__ HoMfmaDirectData pair[2];
    const int group = threadIdx.x / 32, lane = threadIdx.x % 32;
    const size_t first = size_t(blockIdx.x) * 2;
    if (first >= num_elements) return;
    const size_t element = first + group;
    const bool active = element < num_elements;
    // Both half-waves must issue each MFMA; the unused partner repeats the last element.
    const size_t data_element = active ? element : num_elements - 1;
    const int* dof = elem_dof + data_element * 512;
    auto& w = pair[group];
    for (int p = lane; p < 512; p += 32) w.u[p] = u[dof[p]];
    syncWarp();
    for (int batch = 0; batch < 8; ++batch)
    {
        mfmaDirectContract<0, false>(D, pair[0].u, pair[1].u, pair[0].g[0], pair[1].g[0], batch);
        mfmaDirectContract<1, false>(D, pair[0].u, pair[1].u, pair[0].g[1], pair[1].g[1], batch);
        mfmaDirectContract<2, false>(D, pair[0].u, pair[1].u, pair[0].g[2], pair[1].g[2], batch);
    }
    const double* metric = G + data_element * (GMode == 0 ? 512 * 6 : 6);
    for (int p = lane; p < 512; p += 32)
    {
        const double* gp = metric + (GMode == 0 ? p * 6 : 0);
        double weight = 1;
        if constexpr (GMode == 1) weight = weights[p / 64] * weights[(p / 8) % 8] * weights[p % 8];
        const double ar = w.g[0][p], as = w.g[1][p], at = w.g[2][p];
        w.g[0][p] = weight * (gp[0] * ar + gp[1] * as + gp[2] * at);
        w.g[1][p] = weight * (gp[1] * ar + gp[3] * as + gp[4] * at);
        w.g[2][p] = weight * (gp[2] * ar + gp[4] * as + gp[5] * at);
        w.u[p] = 0;
    }
    syncWarp();
    for (int batch = 0; batch < 8; ++batch)
        mfmaDirectContract<0, true>(D, pair[0].g[0], pair[1].g[0], pair[0].u, pair[1].u, batch);
    for (int batch = 0; batch < 8; ++batch)
        mfmaDirectContract<1, true>(D, pair[0].g[1], pair[1].g[1], pair[0].u, pair[1].u, batch);
    for (int batch = 0; batch < 8; ++batch)
        mfmaDirectContract<2, true>(D, pair[0].g[2], pair[1].g[2], pair[0].u, pair[1].u, batch);
    if (active)
        for (int p = lane; p < 512; p += 32) atomicAdd(y + dof[p], w.u[p]);
#else
    __builtin_trap(); // This path requires gfx90a or gfx942 FP64 MFMA.
#endif
}

// Asynchronous additive apply; the caller owns output initialization and stream ordering.
template<int GMode = 0>
inline cudaError_t ho_laplacian_mfma_direct_launch(const double* u, double* y, const int* elem_dof,
                                                  const double* D, const double* G, const double* weights,
                                                  size_t num_elements, cudaStream_t stream = 0)
{
    static_assert(GMode == 0 || GMode == 1);
    if (num_elements == 0) return cudaSuccess;
    ho_laplacian_apply_mfma_direct<GMode><<<unsigned(num_elements / 2 + num_elements % 2), 64, 0, stream>>>(
        u, y, elem_dof, D, G, weights, num_elements);
    return cudaGetLastError();
}
#endif
} // namespace mars::fem
