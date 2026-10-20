#pragma once

// FP64 8x8 tile products for the tensor-core kernels on AMD. CUDA's FP64 WMMA tile is m8n8k4, one
// element per 32-lane warp. AMD has no 8x8 FP64 MFMA; v_mfma_f64_4x4x4f64 (gfx90a, gfx942) runs four
// 4x4x4 blocks per 64-lane wavefront, which is one element's 8x8 output (2x2 blocks) for one k-step of
// 4. A wavefront holds two elements, one per 32-lane group, so every k-step issues one instruction per
// element and every lane reads that element's tile.
//
// Layout measured on MI300A, lane l = 16k + 4b + t: A of block b holds [t][k], B holds [k][t], and the
// result in lane 16r + 4b + c is [r][c]. Block b covers rows 4*(b/2).. and columns 4*(b%2).. here.
// Tiles use the WMMA layouts: A[8x4] col-major (ldm 8), B[4x8] col-major (ldm 4), C[8x8] col-major
// (ldm 8). Both 32-lane groups must call these together.
#include "mars_gpu_warp.cuh"

#if defined(__HIP_DEVICE_COMPILE__) && (defined(__gfx90a__) || defined(__gfx942__))
#define MARS_HAVE_FP64_MFMA 1

namespace mars
{
// acc[e] += A_e . B_e for the two elements e of the wavefront
__device__ __forceinline__ void mfmaPairStep(double acc[2], const double* A0, const double* B0, const double* A1,
                                             const double* B1)
{
    const int l = threadIdx.x % 64;
    const int q = l / 16, b = (l % 16) / 4, t = l % 4;
    const int R = b >> 1, C = b & 1;
    acc[0] = __builtin_amdgcn_mfma_f64_4x4x4f64(A0[q * 8 + R * 4 + t], B0[(C * 4 + t) * 4 + q], acc[0], 0, 0, 0);
    acc[1] = __builtin_amdgcn_mfma_f64_4x4x4f64(A1[q * 8 + R * 4 + t], B1[(C * 4 + t) * 4 + q], acc[1], 0, 0, 0);
}

__device__ __forceinline__ void mfmaPairStore(const double acc[2], double* C0, double* C1)
{
    const int l = threadIdx.x % 64;
    const int q = l / 16, b = (l % 16) / 4, t = l % 4;
    const int R = b >> 1, C = b & 1;
    C0[(C * 4 + t) * 8 + R * 4 + q] = acc[0];
    C1[(C * 4 + t) * 8 + R * 4 + q] = acc[1];
}
} // namespace mars

#endif
