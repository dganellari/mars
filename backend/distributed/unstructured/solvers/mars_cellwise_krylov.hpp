#pragma once
// Krylov solves on element-local (cell-wise) vectors: each element stores its own
// copy of every node it touches, and no assembled vector is ever formed
// (Wichrowski, "Coalesced Matrix-Free Finite Elements in Cell-Wise Storage",
// internal-notes/FlexibleCG.pdf).
//
// The operator maps a continuous element-local field to unassembled element
// residuals; it needs no communication. Elements talk only inside the
// preconditioner, through direct stiffness summation (DSS). The CVFEM operator is
// nonsymmetric, so instead of the paper's flexible CG this runs BiCGStab on the
// left-preconditioned operator P A. Every vector is then continuous, and the inner
// product that weights each element-local copy by 1/(number of copies) equals the
// assembled one, so the iterates are those of the assembled solve. Reference:
// marsir-mlir/test/cellwise_krylov_ref.py.
//
// First version: one structured block of hexahedra at p = 7. The DSS is the paper's
// dimensionally split cascade: three axis passes of one-to-one face sums, which
// complete the edge and vertex sums as a side effect. No atomics, no coloring, no
// index map.

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace mars {
namespace cellwise {

constexpr int kN = 8, kNN = 64, kN3 = 512;   // nodes per edge, face, element (p = 7)

#define MARS_CELLWISE_CK(call)                                                   \
    do {                                                                         \
        cudaError_t e_ = (call);                                                 \
        if (e_ != cudaSuccess) {                                                 \
            fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_),  \
                    __FILE__, __LINE__);                                         \
            std::abort();                                                        \
        }                                                                        \
    } while (0)

// nx * ny * nz elements; element (ex, ey, ez) is ((ex * ny) + ey) * nz + ez, and its
// local node (a, b, c) is (a * 8 + b) * 8 + c with local axis a along x.
struct Block {
    int nx, ny, nz;
    long long elements() const { return (long long)nx * ny * nz; }
    long long values() const { return elements() * kN3; }
};

// One thread per face node of each element pair adjacent along AXIS: sum the two
// copies, write the sum to both. Each copy belongs to exactly one pair in a pass,
// so plain stores are race-free; the passes run in sequence so that each sees the
// partial sums of the previous one.
template <int AXIS>
__global__ void dss_axis_kernel(double* __restrict__ v, int nx, int ny, int nz)
{
    int dims[3] = {nx, ny, nz};
    dims[AXIS] -= 1;
    const long long pairs = (long long)dims[0] * dims[1] * dims[2];
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= pairs * kNN) return;
    const int face = (int)(t % kNN);
    const long long pair = t / kNN;
    const int ez = (int)(pair % dims[2]);
    const long long q = pair / dims[2];
    const int ey = (int)(q % dims[1]);
    const int ex = (int)(q / dims[1]);
    const long long e = ((long long)ex * ny + ey) * nz + ez;
    const long long f = e + (AXIS == 0 ? (long long)ny * nz : AXIS == 1 ? nz : 1);
    const int u = face / kN, w = face % kN;
    int lo, hi;   // the shared plane: last in e, first in f
    if (AXIS == 0) { lo = (kN - 1) * kNN + u * kN + w; hi = u * kN + w; }
    else if (AXIS == 1) { lo = u * kNN + (kN - 1) * kN + w; hi = u * kNN + w; }
    else { lo = u * kNN + w * kN + (kN - 1); hi = u * kNN + w * kN; }
    const double s = v[e * kN3 + lo] + v[f * kN3 + hi];
    v[e * kN3 + lo] = s;
    v[f * kN3 + hi] = s;
}

inline void dss(double* d_v, const Block& b, cudaStream_t stream = 0)
{
    constexpr int threads = 256;
    auto blocks = [&](long long pairs) {
        return (unsigned)((pairs * kNN + threads - 1) / threads);
    };
    if (b.nx > 1)
        dss_axis_kernel<0><<<blocks((long long)(b.nx - 1) * b.ny * b.nz), threads, 0, stream>>>(
            d_v, b.nx, b.ny, b.nz);
    if (b.ny > 1)
        dss_axis_kernel<1><<<blocks((long long)b.nx * (b.ny - 1) * b.nz), threads, 0, stream>>>(
            d_v, b.nx, b.ny, b.nz);
    if (b.nz > 1)
        dss_axis_kernel<2><<<blocks((long long)b.nx * b.ny * (b.nz - 1)), threads, 0, stream>>>(
            d_v, b.nx, b.ny, b.nz);
    MARS_CELLWISE_CK(cudaGetLastError());
}

__device__ inline bool on_boundary(long long t, int nx, int ny, int nz)
{
    const long long e = t / kN3;
    const int l = (int)(t % kN3);
    const int a = l / kNN, b = (l / kN) % kN, c = l % kN;
    const int ez = (int)(e % nz);
    const long long q = e / nz;
    const int ey = (int)(q % ny);
    const int ex = (int)(q / ny);
    return (ex == 0 && a == 0) || (ex == nx - 1 && a == kN - 1) ||
           (ey == 0 && b == 0) || (ey == ny - 1 && b == kN - 1) ||
           (ez == 0 && c == 0) || (ez == nz - 1 && c == kN - 1);
}

// After the DSS: homogeneous Dirichlet copies to zero, the rest scaled by the
// assembled diagonal (Jacobi).
__global__ void jacobi_mask_kernel(double* __restrict__ z, const double* __restrict__ diag,
                                   long long n, int nx, int ny, int nz)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    z[t] = on_boundary(t, nx, ny, nz) ? 0.0 : z[t] / diag[t];
}

// z = P r: assemble the unassembled residual, then the Jacobi-and-Dirichlet step.
// The result is continuous, as the equivalence with the assembled solve requires.
inline void precondition(const double* d_r, double* d_z, const double* d_diag, const Block& b,
                         cudaStream_t stream = 0)
{
    MARS_CELLWISE_CK(cudaMemcpyAsync(d_z, d_r, b.values() * sizeof(double),
                                     cudaMemcpyDeviceToDevice, stream));
    dss(d_z, b, stream);
    constexpr int threads = 256;
    jacobi_mask_kernel<<<(unsigned)((b.values() + threads - 1) / threads), threads, 0, stream>>>(
        d_z, d_diag, b.values(), b.nx, b.ny, b.nz);
    MARS_CELLWISE_CK(cudaGetLastError());
}

// Deterministic two-pass reduction of sum a * b / copies: fixed per-block partial
// sums, then one block adds them in a fixed order, so repeated solves give
// bit-identical scalars.
constexpr int kDotThreads = 256, kDotBlocks = 1024;

__global__ void wdot_partial_kernel(const double* __restrict__ a, const double* __restrict__ b,
                                    const double* __restrict__ copies, long long n,
                                    double* __restrict__ partial)
{
    __shared__ double sh[kDotThreads];
    double s = 0.0;
    for (long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x; t < n;
         t += (long long)gridDim.x * blockDim.x)
        s += a[t] * b[t] / copies[t];
    sh[threadIdx.x] = s;
    __syncthreads();
    for (int w = kDotThreads / 2; w > 0; w >>= 1) {
        if ((int)threadIdx.x < w) sh[threadIdx.x] += sh[threadIdx.x + w];
        __syncthreads();
    }
    if (threadIdx.x == 0) partial[blockIdx.x] = sh[0];
}

__global__ void sum_partials_kernel(const double* __restrict__ partial, double* __restrict__ out)
{
    __shared__ double sh[kDotThreads];
    double s = 0.0;
    for (int i = threadIdx.x; i < kDotBlocks; i += kDotThreads) s += partial[i];
    sh[threadIdx.x] = s;
    __syncthreads();
    for (int w = kDotThreads / 2; w > 0; w >>= 1) {
        if ((int)threadIdx.x < w) sh[threadIdx.x] += sh[threadIdx.x + w];
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = sh[0];
}

// Krylov scalars live in device memory: the dot products write them, one-thread
// kernels combine them, and the vector updates read them, so no step waits for a
// round trip to the host. Slots of Workspace::scalars:
enum Scalar { kRho, kRhoOld, kAlpha, kOmega, kBeta, kRhV, kTS, kTT, kRR, kScalars };

__global__ void begin_iteration_kernel(double* __restrict__ sc)
{
    sc[kBeta] = (sc[kRho] / sc[kRhoOld]) * (sc[kAlpha] / sc[kOmega]);
    sc[kRhoOld] = sc[kRho];
}

__global__ void set_alpha_kernel(double* __restrict__ sc) { sc[kAlpha] = sc[kRho] / sc[kRhV]; }

__global__ void set_omega_kernel(double* __restrict__ sc) { sc[kOmega] = sc[kTS] / sc[kTT]; }

__global__ void init_scalars_kernel(double* __restrict__ sc)
{
    sc[kRhoOld] = 1.0;
    sc[kAlpha] = 1.0;
    sc[kOmega] = 1.0;
}

// BiCGStab vector updates, fused so each vector is read once per update.
__global__ void update_p_kernel(double* __restrict__ p, const double* __restrict__ r,
                                const double* __restrict__ v, const double* __restrict__ sc,
                                long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) p[t] = r[t] + sc[kBeta] * (p[t] - sc[kOmega] * v[t]);
}

__global__ void update_s_kernel(double* __restrict__ s, const double* __restrict__ r,
                                const double* __restrict__ v, const double* __restrict__ sc,
                                long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) s[t] = r[t] - sc[kAlpha] * v[t];
}

__global__ void update_xr_kernel(double* __restrict__ x, double* __restrict__ r,
                                 const double* __restrict__ p, const double* __restrict__ s,
                                 const double* __restrict__ tt, const double* __restrict__ sc,
                                 long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) {
        x[t] += sc[kAlpha] * p[t] + sc[kOmega] * s[t];
        r[t] = s[t] - sc[kOmega] * tt[t];
    }
}

__global__ void fill_kernel(double* __restrict__ v, double value, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) v[t] = value;
}

struct Workspace {
    double *r = nullptr, *rh = nullptr, *p = nullptr, *v = nullptr, *s = nullptr, *t = nullptr,
           *q = nullptr, *partial = nullptr, *scalars = nullptr;
    explicit Workspace(long long n)
    {
        for (double** a : {&r, &rh, &p, &v, &s, &t, &q})
            MARS_CELLWISE_CK(cudaMalloc(a, n * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&partial, kDotBlocks * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&scalars, kScalars * sizeof(double)));
    }
    ~Workspace()
    {
        for (double* a : {r, rh, p, v, s, t, q, partial, scalars}) cudaFree(a);
    }
    Workspace(const Workspace&) = delete;
    Workspace& operator=(const Workspace&) = delete;
};

// sum a * b / copies into the device slot `out`.
inline void wdot(const double* d_a, const double* d_b, const double* d_copies, long long n,
                 Workspace& ws, double* out, cudaStream_t stream = 0)
{
    wdot_partial_kernel<<<kDotBlocks, kDotThreads, 0, stream>>>(d_a, d_b, d_copies, n, ws.partial);
    sum_partials_kernel<<<1, kDotThreads, 0, stream>>>(ws.partial, out);
    MARS_CELLWISE_CK(cudaGetLastError());
}

// The residual norm is the one value the host needs: it decides when to stop.
inline double read_norm(const Workspace& ws, cudaStream_t stream)
{
    double h = 0.0;
    MARS_CELLWISE_CK(cudaMemcpyAsync(&h, ws.scalars + kRR, sizeof(double),
                                     cudaMemcpyDeviceToHost, stream));
    MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
    return std::sqrt(h);
}

// Number of element-local copies of every node: the DSS of ones.
inline void copies(double* d_copies, const Block& b, cudaStream_t stream = 0)
{
    constexpr int threads = 256;
    fill_kernel<<<(unsigned)((b.values() + threads - 1) / threads), threads, 0, stream>>>(
        d_copies, 1.0, b.values());
    dss(d_copies, b, stream);
}

struct SolveResult {
    int iterations = 0;
    double residual0 = 0.0, residual = 0.0;
    std::vector<double> history;   // ||P r||_w after each iteration
};

// Left-preconditioned BiCGStab on P A x = P b, x starting at zero. `apply(u, y)`
// computes y = A u on element-local vectors (u continuous, y unassembled).
template <typename Apply>
SolveResult bicgstab(Apply&& apply, const Block& b, const double* d_rhs, double* d_x,
                     const double* d_diag, const double* d_copies, double tol, int max_iterations,
                     cudaStream_t stream = 0)
{
    const long long n = b.values();
    constexpr int threads = 256;
    const unsigned grid = (unsigned)((n + threads - 1) / threads);
    Workspace ws(n);
    double* sc = ws.scalars;
    auto dot = [&](const double* x, const double* y, Scalar slot) {
        wdot(x, y, d_copies, n, ws, sc + slot, stream);
    };

    fill_kernel<<<grid, threads, 0, stream>>>(d_x, 0.0, n);
    precondition(d_rhs, ws.r, d_diag, b, stream);   // r = P (b - A 0)
    MARS_CELLWISE_CK(cudaMemcpyAsync(ws.rh, ws.r, n * sizeof(double), cudaMemcpyDeviceToDevice, stream));
    fill_kernel<<<grid, threads, 0, stream>>>(ws.p, 0.0, n);
    fill_kernel<<<grid, threads, 0, stream>>>(ws.v, 0.0, n);
    init_scalars_kernel<<<1, 1, 0, stream>>>(sc);
    MARS_CELLWISE_CK(cudaGetLastError());

    SolveResult res;
    dot(ws.r, ws.r, kRR);
    res.residual0 = read_norm(ws, stream);
    res.residual = res.residual0;
    for (int k = 0; k < max_iterations; ++k) {
        dot(ws.rh, ws.r, kRho);
        begin_iteration_kernel<<<1, 1, 0, stream>>>(sc);
        update_p_kernel<<<grid, threads, 0, stream>>>(ws.p, ws.r, ws.v, sc, n);
        apply(ws.p, ws.q);
        precondition(ws.q, ws.v, d_diag, b, stream);
        dot(ws.rh, ws.v, kRhV);
        set_alpha_kernel<<<1, 1, 0, stream>>>(sc);
        update_s_kernel<<<grid, threads, 0, stream>>>(ws.s, ws.r, ws.v, sc, n);
        apply(ws.s, ws.q);
        precondition(ws.q, ws.t, d_diag, b, stream);
        dot(ws.t, ws.s, kTS);
        dot(ws.t, ws.t, kTT);
        set_omega_kernel<<<1, 1, 0, stream>>>(sc);
        update_xr_kernel<<<grid, threads, 0, stream>>>(d_x, ws.r, ws.p, ws.s, ws.t, sc, n);
        MARS_CELLWISE_CK(cudaGetLastError());
        dot(ws.r, ws.r, kRR);
        res.residual = read_norm(ws, stream);
        res.history.push_back(res.residual);
        res.iterations = k + 1;
        if (res.residual <= tol * res.residual0) break;
    }
    return res;
}

}  // namespace cellwise
}  // namespace mars
