#pragma once

// Matrix-free high-order Laplacian operator apply for tensor-product hex
// spectral elements (GLL collocation). The element operator y_e = K_e u_e is
// evaluated via sum-factorization: the only real cost is the 1D
// differentiation-matrix contractions (D applied along each tensor direction),
// which is where FP64 tensor cores (DMMA m8n8k4) pay off at p>=3.
//
// This header has the SCALAR (CUDA-core) reference path. It is the correctness
// baseline and the speedup denominator for the DMMA path that slots into the
// same contraction structure (mars_ho_laplacian_dmma.hpp, next step).
//
// Conventions:
//   - order P, n = P+1 nodes per direction, N3 = n^3 nodes per element
//   - GLL collocation: quadrature points == nodes, so interpolation is the
//     identity and the mass matrix is diagonal. Only the derivative
//     contractions remain. This is the standard spectral-element structure.
//   - Per-element geometric factor G is the symmetric 3x3 metric tensor
//     w * J^{-1} J^{-T} |J| at each quad point, precomputed on the host/GPU
//     and passed in as 6 unique components per quad point (Gxx,Gxy,Gxz,Gyy,Gyz,Gzz).
//     For an affine (straight-sided) element it is constant; we still store
//     per-point to keep the path general for curved elements later.

#include <cuda_runtime.h>
#include <cstddef>

namespace mars {
namespace fem {

// Apply the 1D differentiation matrix D (n x n) along direction `dir` of a
// rank-3 tensor u[n][n][n] flattened as u[i*n*n + j*n + k] (i=slowest).
//   dir 0: contract index i ->  out[a][j][k] = sum_i D[a][i] u[i][j][k]
//   dir 1: contract index j
//   dir 2: contract index k (fastest)
// Scalar reference; one thread evaluates the whole element (clear, not fast).
template<typename RealType, int n>
__device__ inline void contract1d(const RealType* __restrict__ D,
                                   const RealType* __restrict__ in,
                                   RealType* __restrict__ out,
                                   int dir)
{
    constexpr int nn = n * n;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c)
            {
                RealType acc = RealType(0);
                for (int s = 0; s < n; ++s)
                {
                    RealType dval;
                    RealType uval;
                    if (dir == 0) { dval = D[a * n + s]; uval = in[s * nn + b * n + c]; }
                    else if (dir == 1) { dval = D[b * n + s]; uval = in[a * nn + s * n + c]; }
                    else { dval = D[c * n + s]; uval = in[a * nn + b * n + s]; }
                    acc += dval * uval;
                }
                // (a,b,c) enumerate the output point directly; `dir` only
                // selects which input index is contracted and which D row.
                out[a * nn + b * n + c] = acc;
            }
}

// Scalar matrix-free Laplacian apply. One thread per element.
//   d_u    : input nodal field, [numElements * N3]
//   d_y    : output (K_e u), [numElements * N3]
//   d_D    : 1D differentiation matrix, [n*n], shared by all elements
//   d_G    : per-element-per-point metric, [numElements * N3 * 6]
template<typename RealType, int P>
__global__ void ho_laplacian_apply_scalar(const RealType* __restrict__ d_u,
                                           RealType* __restrict__ d_y,
                                           const RealType* __restrict__ d_D,
                                           const RealType* __restrict__ d_G,
                                           size_t numElements)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;

    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    // Load D into registers/local (small: n*n).
    RealType D[n * n];
    for (int i = 0; i < n * n; ++i) D[i] = d_D[i];

    const RealType* u = d_u + e * N3;
    const RealType* G = d_G + e * N3 * 6;
    RealType*       y = d_y + e * N3;

    // Derivatives at quad points (collocation): ur = D u along each dir.
    RealType ur[N3], us[N3], ut[N3];
    contract1d<RealType, n>(D, u, ur, 0);
    contract1d<RealType, n>(D, u, us, 1);
    contract1d<RealType, n>(D, u, ut, 2);

    // Apply symmetric metric G at each point: (gr,gs,gt) = G * (ur,us,ut).
    RealType gr[N3], gs[N3], gt[N3];
    for (int p = 0; p < N3; ++p)
    {
        RealType Gxx = G[p * 6 + 0], Gxy = G[p * 6 + 1], Gxz = G[p * 6 + 2];
        RealType Gyy = G[p * 6 + 3], Gyz = G[p * 6 + 4], Gzz = G[p * 6 + 5];
        RealType a = ur[p], b = us[p], c = ut[p];
        gr[p] = Gxx * a + Gxy * b + Gxz * c;
        gs[p] = Gxy * a + Gyy * b + Gyz * c;
        gt[p] = Gxz * a + Gyz * b + Gzz * c;
    }

    // Integrate back: y = D^T gr (dir0) + D^T gs (dir1) + D^T gt (dir2).
    // D^T contraction: reuse contract1d with the transpose by swapping index.
    RealType yr[N3], ys[N3], yt[N3];
    // Build D^T once.
    RealType DT[n * n];
    for (int a = 0; a < n; ++a)
        for (int s = 0; s < n; ++s)
            DT[a * n + s] = D[s * n + a];

    contract1d<RealType, n>(DT, gr, yr, 0);
    contract1d<RealType, n>(DT, gs, ys, 1);
    contract1d<RealType, n>(DT, gt, yt, 2);

    for (int p = 0; p < N3; ++p) y[p] = yr[p] + ys[p] + yt[p];
}

// Continuous-Galerkin matrix-free Laplacian apply: y = A u over a GLOBAL DOF
// vector. Same per-element sum-factorization as the scalar kernel, wrapped in
//   gather:  u_local[l] = u_dof[elemDof[e,l]]
//   apply:   y_local = K_e u_local  (sum-factorization)
//   scatter: atomicAdd(y_dof[elemDof[e,l]], y_local[l])   (direct stiffness sum)
// CG continuity is whatever the elemDof map encodes: shared edge/face/corner
// DOFs that map to the same global id get their element contributions summed.
// d_elemDof: [numElements * N3] local-node -> global-DOF; d_y_dof must be
// zeroed before launch (scatter is additive).
template<typename RealType, int P>
__global__ void ho_laplacian_apply_cg(const RealType* __restrict__ d_u_dof,
                                      RealType* __restrict__ d_y_dof,
                                      const int* __restrict__ d_elemDof,
                                      const RealType* __restrict__ d_D,
                                      const RealType* __restrict__ d_G,
                                      size_t numElements)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;

    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    RealType D[n * n];
    for (int i = 0; i < n * n; ++i) D[i] = d_D[i];

    const int*      edof = d_elemDof + e * N3;
    const RealType* G    = d_G + e * N3 * 6;

    RealType u[N3];
    for (int l = 0; l < N3; ++l) u[l] = d_u_dof[edof[l]];

    RealType ur[N3], us[N3], ut[N3];
    contract1d<RealType, n>(D, u, ur, 0);
    contract1d<RealType, n>(D, u, us, 1);
    contract1d<RealType, n>(D, u, ut, 2);

    RealType gr[N3], gs[N3], gt[N3];
    for (int p = 0; p < N3; ++p)
    {
        RealType Gxx = G[p * 6 + 0], Gxy = G[p * 6 + 1], Gxz = G[p * 6 + 2];
        RealType Gyy = G[p * 6 + 3], Gyz = G[p * 6 + 4], Gzz = G[p * 6 + 5];
        RealType a = ur[p], b = us[p], c = ut[p];
        gr[p] = Gxx * a + Gxy * b + Gxz * c;
        gs[p] = Gxy * a + Gyy * b + Gyz * c;
        gt[p] = Gxz * a + Gyz * b + Gzz * c;
    }

    RealType DT[n * n];
    for (int a = 0; a < n; ++a)
        for (int s = 0; s < n; ++s)
            DT[a * n + s] = D[s * n + a];

    RealType yr[N3], ys[N3], yt[N3];
    contract1d<RealType, n>(DT, gr, yr, 0);
    contract1d<RealType, n>(DT, gs, ys, 1);
    contract1d<RealType, n>(DT, gt, yt, 2);

    for (int l = 0; l < N3; ++l)
        atomicAdd(&d_y_dof[edof[l]], yr[l] + ys[l] + yt[l]);
}

} // namespace fem
} // namespace mars
