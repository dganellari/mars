#pragma once

#include "../utils/mars_gpu_runtime.hpp"
#include "../utils/mars_gpu_mma.cuh"
#include <cstddef>

// AMD path for P=7: the same D^T G D operator as the staged DMMA kernel. Two elements share a wavefront.
namespace mars::fem
{
// Plane schedule. Keeping u and the three gradients of both elements in LDS takes 32 KiB per block: only
// two wavefronts then fit on a CU, and on MI300A they waited 80-87% of the time. Here each lane owns one
// (j, k) column of both elements and keeps u and y along i in registers, which is also where the in-plane
// MFMA results land. The i contraction and its transpose become lane-local FMAs, and LDS holds only the
// current i-plane of u and of the two in-plane fluxes: 3 KiB per block.
struct HoMfmaPlaneData
{
    double u[64];
    double w1[64];
    double w2[64];
};

#if defined(MARS_HAVE_FP64_MFMA)
// Contraction over j (Direction 1) or k (Direction 2) of one 8x8 plane [j][k] with D or D^T, for both
// elements. Results are added to acc0/acc1 at this lane's plane position (mfmaPlanePosition).
template<int Direction, bool Transpose>
__device__ __forceinline__ void mfmaPlaneContract(const double* D, const double* in0, const double* in1,
                                                  double& acc0, double& acc1)
{
    static_assert(Direction == 1 || Direction == 2);
    const int lane = threadIdx.x % 64;
    const int q = lane / 16, block = (lane % 16) / 4, t = lane % 4;
    const int row = 4 * (block / 2) + t, col = 4 * (block % 2) + t;
    #pragma unroll
    for (int kt = 0; kt < 2; ++kt)
    {
        const int s = 4 * kt + q;
        double a0, a1, b0, b1;
        if constexpr (Direction == 1)
        {
            a0 = a1 = D[Transpose ? s * 8 + row : row * 8 + s];
            b0 = in0[s * 8 + col];
            b1 = in1[s * 8 + col];
        }
        else
        {
            a0 = in0[row * 8 + s];
            a1 = in1[row * 8 + s];
            b0 = b1 = D[Transpose ? s * 8 + col : col * 8 + s];
        }
        acc0 = __builtin_amdgcn_mfma_f64_4x4x4f64(a0, b0, acc0, 0, 0, 0);
        acc1 = __builtin_amdgcn_mfma_f64_4x4x4f64(a1, b1, acc1, 0, 0, 0);
    }
}

__device__ __forceinline__ int mfmaPlanePosition()
{
    const int lane = threadIdx.x % 64;
    const int q = lane / 16, block = (lane % 16) / 4, t = lane % 4;
    return (4 * (block / 2) + q) * 8 + 4 * (block % 2) + t;
}
#endif

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
template<int GMode>
__global__ void __launch_bounds__(64)
ho_laplacian_apply_mfma_plane(const double* __restrict__ u, double* __restrict__ y,
                              const int* __restrict__ elem_dof, const double* __restrict__ D,
                              const double* __restrict__ G, const double* __restrict__ weights,
                              size_t num_elements)
{
#if defined(MARS_HAVE_FP64_MFMA)
    __shared__ HoMfmaPlaneData pair[2];
    const size_t first = size_t(blockIdx.x) * 2;
    if (first >= num_elements) return;
    // Both elements must issue every MFMA; a missing second element repeats the first and does not scatter.
    const bool second = first + 1 < num_elements;
    const int* dof0 = elem_dof + first * 512;
    const int* dof1 = elem_dof + (second ? first + 1 : first) * 512;
    const double* metric0 = G + first * (GMode == 0 ? 512 * 6 : 6);
    const double* metric1 = G + (second ? first + 1 : first) * (GMode == 0 ? 512 * 6 : 6);
    const int jk = mfmaPlanePosition();

    double u0[8], u1[8], y0[8], y1[8];
    #pragma unroll
    for (int i = 0; i < 8; ++i)
    {
        u0[i] = u[dof0[i * 64 + jk]];
        u1[i] = u[dof1[i * 64 + jk]];
        y0[i] = y1[i] = 0;
    }
    double plane_weight = 1;
    if constexpr (GMode == 1) plane_weight = weights[jk / 8] * weights[jk % 8];

    #pragma unroll
    for (int i = 0; i < 8; ++i)
    {
        pair[0].u[jk] = u0[i];
        pair[1].u[jk] = u1[i];
        syncWarp();
        double r0 = 0, r1 = 0, s0 = 0, s1 = 0, t0 = 0, t1 = 0;
        mfmaPlaneContract<1, false>(D, pair[0].u, pair[1].u, s0, s1);
        mfmaPlaneContract<2, false>(D, pair[0].u, pair[1].u, t0, t1);
        #pragma unroll
        for (int m = 0; m < 8; ++m)
        {
            r0 += D[i * 8 + m] * u0[m];
            r1 += D[i * 8 + m] * u1[m];
        }

        const double* g0 = metric0 + (GMode == 0 ? (i * 64 + jk) * 6 : 0);
        const double* g1 = metric1 + (GMode == 0 ? (i * 64 + jk) * 6 : 0);
        const double weight = GMode == 1 ? weights[i] * plane_weight : 1;
        const double w00 = weight * (g0[0] * r0 + g0[1] * s0 + g0[2] * t0);
        const double w01 = weight * (g0[1] * r0 + g0[3] * s0 + g0[4] * t0);
        const double w02 = weight * (g0[2] * r0 + g0[4] * s0 + g0[5] * t0);
        const double w10 = weight * (g1[0] * r1 + g1[1] * s1 + g1[2] * t1);
        const double w11 = weight * (g1[1] * r1 + g1[3] * s1 + g1[4] * t1);
        const double w12 = weight * (g1[2] * r1 + g1[4] * s1 + g1[5] * t1);
        #pragma unroll
        for (int m = 0; m < 8; ++m)
        {
            y0[m] += D[i * 8 + m] * w00;
            y1[m] += D[i * 8 + m] * w10;
        }

        pair[0].w1[jk] = w01;
        pair[1].w1[jk] = w11;
        pair[0].w2[jk] = w02;
        pair[1].w2[jk] = w12;
        syncWarp();
        mfmaPlaneContract<1, true>(D, pair[0].w1, pair[1].w1, y0[i], y1[i]);
        mfmaPlaneContract<2, true>(D, pair[0].w2, pair[1].w2, y0[i], y1[i]);
    }
    #pragma unroll
    for (int i = 0; i < 8; ++i)
    {
        atomicAdd(y + dof0[i * 64 + jk], y0[i]);
        if (second) atomicAdd(y + dof1[i * 64 + jk], y1[i]);
    }
#else
    __builtin_trap(); // This path requires gfx90a or gfx942 FP64 MFMA.
#endif
}

// Asynchronous additive apply; the caller owns output initialization and stream ordering.
template<int GMode = 0>
inline cudaError_t ho_laplacian_mfma_plane_launch(const double* u, double* y, const int* elem_dof,
                                                 const double* D, const double* G, const double* weights,
                                                 size_t num_elements, cudaStream_t stream = 0)
{
    static_assert(GMode == 0 || GMode == 1);
    if (num_elements == 0) return cudaSuccess;
    ho_laplacian_apply_mfma_plane<GMode><<<unsigned(num_elements / 2 + num_elements % 2), 64, 0, stream>>>(
        u, y, elem_dof, D, G, weights, num_elements);
    return cudaGetLastError();
}
#endif
} // namespace mars::fem
