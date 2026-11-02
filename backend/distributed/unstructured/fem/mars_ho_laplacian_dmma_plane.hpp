#pragma once

#include "../utils/mars_gpu_runtime.hpp"
#include <cstddef>

// CUDA path for P=7: the plane schedule of the AMD MFMA kernel (mars_ho_laplacian_mfma_direct.hpp) on
// FP64 mma.sync m8n8k4, one warp per element. The staged kernel keeps the element's tensors in shared
// memory (17.5 KiB per warp), which caps an H100 SM at 12 warps; it then issues in ~10% of its cycles.
// Here lane L owns the plane positions (j, k) = (L/4, 2(L%4)) and (L/4, 2(L%4)+1), which are exactly its
// m8n8k4 C-fragment entries, and keeps u and y along i for both in registers. The i contraction, its
// transpose and the metric are lane-local; shared memory holds the current i-plane of u, w1 and w2.
#if defined(MARS_DMMA_PLANE_DEVICE) || (!defined(__HIP_PLATFORM_AMD__) && !defined(__HIPCC__))
// FP64 mma.sync needs SM80+. Device passes for older targets compile the kernel with a trapping body;
// the public launchers select the staged kernel there.
#if defined(MARS_DMMA_PLANE_DEVICE) || !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 800
#define MARS_DMMA_PLANE_BODY 1
#endif

namespace mars::fem
{
struct HoDmmaPlaneData
{
    double u[64];
    double w1[64];
    double w2[64];
};

// C += A B with A 8x4 and B 4x8. Lane L holds A[L/4][L%4], B[L%4][L/4], C[L/4][2(L%4)] and C[L/4][2(L%4)+1].
__device__ __forceinline__ void dmmaM8n8k4(double a, double b, double& c0, double& c1)
{
#if defined(MARS_DMMA_PLANE_DEVICE)
    emulate_dmma(a, b, c0, c1);
#elif defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800
    asm volatile("mma.sync.aligned.m8n8k4.row.col.f64.f64.f64.f64 {%0, %1}, {%2}, {%3}, {%0, %1};"
                 : "+d"(c0), "+d"(c1)
                 : "d"(a), "d"(b));
#endif
}

// Contraction over j (Direction 1) or k (Direction 2) of one 8x8 plane [j][k] with D or D^T. Results are
// added to c0/c1 at this lane's plane positions.
template<int Direction, bool Transpose>
__device__ __forceinline__ void dmmaPlaneContract(const double* D, const double* in, double& c0, double& c1)
{
    static_assert(Direction == 1 || Direction == 2);
    const int lane = threadIdx.x % 32, r = lane / 4, q = lane % 4;
    #pragma unroll
    for (int kt = 0; kt < 2; ++kt)
    {
        const int s = 4 * kt + q;
        const double d = D[Transpose ? s * 8 + r : r * 8 + s];
        if constexpr (Direction == 1) dmmaM8n8k4(d, in[s * 8 + r], c0, c1);
        else dmmaM8n8k4(in[r * 8 + s], d, c0, c1);
    }
}

template<int GMode, int WarpsPerBlock>
__global__ void __launch_bounds__(WarpsPerBlock * 32)
ho_laplacian_apply_dmma_plane(const double* __restrict__ u, double* __restrict__ y,
                              const int* __restrict__ elem_dof, const double* __restrict__ D,
                              const double* __restrict__ G, const double* __restrict__ weights,
                              size_t num_elements)
{
#if defined(MARS_DMMA_PLANE_BODY)
    __shared__ HoDmmaPlaneData planes[WarpsPerBlock];
    const int warp = threadIdx.x / 32, lane = threadIdx.x % 32;
    const size_t e = size_t(blockIdx.x) * WarpsPerBlock + warp;
    // Warp-uniform, and the kernel only synchronizes within a warp.
    if (e >= num_elements) return;
    HoDmmaPlaneData& p = planes[warp];
    const int jk0 = (lane / 4) * 8 + 2 * (lane % 4), jk1 = jk0 + 1;
    const int* dof = elem_dof + e * 512;
    const double* metric = G + e * (GMode == 0 ? 512 * 6 : 6);

    double u0[8], u1[8], y0[8], y1[8];
    #pragma unroll
    for (int i = 0; i < 8; ++i)
    {
        u0[i] = u[dof[i * 64 + jk0]];
        u1[i] = u[dof[i * 64 + jk1]];
        y0[i] = y1[i] = 0;
    }
    double plane_weight0 = 1, plane_weight1 = 1;
    if constexpr (GMode == 1)
    {
        plane_weight0 = weights[jk0 / 8] * weights[jk0 % 8];
        plane_weight1 = weights[jk1 / 8] * weights[jk1 % 8];
    }

    #pragma unroll
    for (int i = 0; i < 8; ++i)
    {
        p.u[jk0] = u0[i];
        p.u[jk1] = u1[i];
        __syncwarp();
        double r0 = 0, r1 = 0, s0 = 0, s1 = 0, t0 = 0, t1 = 0;
        dmmaPlaneContract<1, false>(D, p.u, s0, s1);
        dmmaPlaneContract<2, false>(D, p.u, t0, t1);
        #pragma unroll
        for (int m = 0; m < 8; ++m)
        {
            r0 += D[i * 8 + m] * u0[m];
            r1 += D[i * 8 + m] * u1[m];
        }

        const double* g0 = metric + (GMode == 0 ? (i * 64 + jk0) * 6 : 0);
        const double* g1 = metric + (GMode == 0 ? (i * 64 + jk1) * 6 : 0);
        const double weight0 = GMode == 1 ? weights[i] * plane_weight0 : 1;
        const double weight1 = GMode == 1 ? weights[i] * plane_weight1 : 1;
        const double w00 = weight0 * (g0[0] * r0 + g0[1] * s0 + g0[2] * t0);
        const double w01 = weight0 * (g0[1] * r0 + g0[3] * s0 + g0[4] * t0);
        const double w02 = weight0 * (g0[2] * r0 + g0[4] * s0 + g0[5] * t0);
        const double w10 = weight1 * (g1[0] * r1 + g1[1] * s1 + g1[2] * t1);
        const double w11 = weight1 * (g1[1] * r1 + g1[3] * s1 + g1[4] * t1);
        const double w12 = weight1 * (g1[2] * r1 + g1[4] * s1 + g1[5] * t1);
        #pragma unroll
        for (int m = 0; m < 8; ++m)
        {
            y0[m] += D[i * 8 + m] * w00;
            y1[m] += D[i * 8 + m] * w10;
        }

        p.w1[jk0] = w01;
        p.w1[jk1] = w11;
        p.w2[jk0] = w02;
        p.w2[jk1] = w12;
        __syncwarp();
        dmmaPlaneContract<1, true>(D, p.w1, y0[i], y1[i]);
        dmmaPlaneContract<2, true>(D, p.w2, y0[i], y1[i]);
    }
    #pragma unroll
    for (int i = 0; i < 8; ++i)
    {
        atomicAdd(y + dof[i * 64 + jk0], y0[i]);
        atomicAdd(y + dof[i * 64 + jk1], y1[i]);
    }
#else
    __trap();
#endif
}

#if !defined(MARS_DMMA_PLANE_DEVICE)
// Asynchronous additive apply on SM80+ (FP64 mma.sync); the caller owns output initialization, stream
// ordering and the compute-capability check.
template<int GMode = 0, int WarpsPerBlock = 4>
inline cudaError_t ho_laplacian_dmma_plane_launch(const double* u, double* y, const int* elem_dof,
                                                 const double* D, const double* G, const double* weights,
                                                 size_t num_elements, cudaStream_t stream = 0)
{
    static_assert(GMode == 0 || GMode == 1);
    if (num_elements == 0) return cudaSuccess;
    const size_t grid = (num_elements + WarpsPerBlock - 1) / WarpsPerBlock;
    ho_laplacian_apply_dmma_plane<GMode, WarpsPerBlock><<<unsigned(grid), WarpsPerBlock * 32, 0, stream>>>(
        u, y, elem_dof, D, G, weights, num_elements);
    return cudaGetLastError();
}
#endif
} // namespace mars::fem
#endif
