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
// First version: one structured block of hexahedra at p = 7. The DSS gathers every
// copy of a node from the neighbouring elements and adds them in the order of the
// paper's dimensionally split cascade (x pairs, then y pairs of those, then z), so
// all copies get the same bits and the result equals the cascade's exactly. Each
// element writes only its own values: no atomics, no coloring, no index map. Unlike
// the cascade, the gather is one pass, so it fuses with the Jacobi step and the dot
// products after it: a preconditioner call reads the residual once and writes the
// result once. Copy counts and the Dirichlet boundary come from the element
// position, not from stored arrays.

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace mars {
namespace cellwise {

constexpr int kN = 8, kNN = 64, kN3 = 512;   // nodes per edge, face, element (p = 7)
constexpr int kThreads = 256;                // threads per block of every kernel here

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
    __host__ __device__ long long elements() const { return (long long)nx * ny * nz; }
    __host__ __device__ long long values() const { return elements() * kN3; }
};

// The paper's cascade, in place: one thread per face node of each element pair
// adjacent along AXIS sums the two copies and writes the sum to both. Kept to time
// against the gather below.
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

inline void dss_cascade(double* d_v, const Block& b, cudaStream_t stream = 0)
{
    auto blocks = [&](long long pairs) {
        return (unsigned)((pairs * kNN + kThreads - 1) / kThreads);
    };
    if (b.nx > 1)
        dss_axis_kernel<0><<<blocks((long long)(b.nx - 1) * b.ny * b.nz), kThreads, 0, stream>>>(
            d_v, b.nx, b.ny, b.nz);
    if (b.ny > 1)
        dss_axis_kernel<1><<<blocks((long long)b.nx * (b.ny - 1) * b.nz), kThreads, 0, stream>>>(
            d_v, b.nx, b.ny, b.nz);
    if (b.nz > 1)
        dss_axis_kernel<2><<<blocks((long long)b.nx * b.ny * (b.nz - 1)), kThreads, 0, stream>>>(
            d_v, b.nx, b.ny, b.nz);
    MARS_CELLWISE_CK(cudaGetLastError());
}

struct Node {
    long long e, t;   // element, and value index e * kN3 + local index
    int a, b, c;      // local coordinates
    int ex, ey, ez;   // element coordinates
};

// Along one axis, local coordinate a of element index ei among n: -1 when the node
// is shared with element ei - 1, +1 with ei + 1, 0 when only this element has it.
__device__ inline int share_side(int a, int ei, int n)
{
    if (a == 0) return ei > 0 ? -1 : 0;
    if (a == kN - 1) return ei < n - 1 ? 1 : 0;
    return 0;
}

__device__ inline bool on_outer_boundary(const Node& nd, const Block& blk)
{
    auto outer = [](int a, int ei, int n) {
        return (a == 0 && ei == 0) || (a == kN - 1 && ei == n - 1);
    };
    return outer(nd.a, nd.ex, blk.nx) || outer(nd.b, nd.ey, blk.ny) || outer(nd.c, nd.ez, blk.nz);
}

// 1 / (number of copies): a power of two, so scaling by it is exact.
__device__ inline double weight(const Node& nd, const Block& blk)
{
    double w = 1.0;
    if (share_side(nd.a, nd.ex, blk.nx)) w *= 0.5;
    if (share_side(nd.b, nd.ey, blk.ny)) w *= 0.5;
    if (share_side(nd.c, nd.ez, blk.nz)) w *= 0.5;
    return w;
}

// Sum of every copy of the node, in the cascade's order: x pairs first, then y pairs
// of those sums, then z pairs. On a shared axis the pair is (low element, high
// element); the low one holds the node on its last plane, the high one on plane 0.
__device__ inline double gather_sum(const double* __restrict__ q, const Node& nd,
                                    const Block& blk)
{
    const int sx = share_side(nd.a, nd.ex, blk.nx);
    const int sy = share_side(nd.b, nd.ey, blk.ny);
    const int sz = share_side(nd.c, nd.ez, blk.nz);
    const long long stride_x = (long long)blk.ny * blk.nz, stride_y = blk.nz;
    const long long low = nd.e - (sx < 0 ? stride_x : 0) - (sy < 0 ? stride_y : 0) - (sz < 0 ? 1 : 0);
    const int a0 = sx ? kN - 1 : nd.a, b0 = sy ? kN - 1 : nd.b, c0 = sz ? kN - 1 : nd.c;
    double zsum = 0.0;
#pragma unroll
    for (int k = 0; k < 2; ++k) {
        if (k == 1 && !sz) break;
        double ysum = 0.0;
#pragma unroll
        for (int j = 0; j < 2; ++j) {
            if (j == 1 && !sy) break;
            const long long base = low + j * stride_y + k;
            const int bc = (j ? 0 : b0) * kN + (k ? 0 : c0);
            double xsum = q[base * kN3 + a0 * kNN + bc];
            if (sx) xsum += q[(base + stride_x) * kN3 + bc];
            ysum = j == 0 ? xsum : ysum + xsum;
        }
        zsum = k == 0 ? ysum : zsum + ysum;
    }
    return zsum;
}

// Calls f(node) for each node of this thread. Each thread block owns a contiguous
// range of elements, so the element coordinates step along without a division per
// element, and a gather's z neighbours (adjacent indices) are read by the same
// block close in time.
template <typename F>
__device__ inline void for_each_node(const Block& blk, F&& f)
{
    const long long E = blk.elements();
    const long long chunk = (E + gridDim.x - 1) / gridDim.x;
    const long long begin = (long long)blockIdx.x * chunk;
    const long long end = begin + chunk < E ? begin + chunk : E;
    if (begin >= end) return;
    Node nd;
    nd.ez = (int)(begin % blk.nz);
    nd.ey = (int)((begin / blk.nz) % blk.ny);
    nd.ex = (int)(begin / ((long long)blk.ny * blk.nz));
    for (nd.e = begin; nd.e < end; ++nd.e) {
        for (int l = threadIdx.x; l < kN3; l += kThreads) {
            nd.t = nd.e * kN3 + l;
            nd.a = l / kNN;
            nd.b = (l / kN) % kN;
            nd.c = l % kN;
            f(nd);
        }
        if (++nd.ez == blk.nz) {
            nd.ez = 0;
            if (++nd.ey == blk.ny) {
                nd.ey = 0;
                ++nd.ex;
            }
        }
    }
}

// Sums NV per-thread values over the thread block in a fixed order, so repeated
// solves give bit-identical scalars; thread 0 gets the totals.
template <int NV>
__device__ inline void block_sum(double (&v)[NV])
{
    __shared__ double sh[NV][kThreads / 32];
    const int lane = threadIdx.x % 32, warp = threadIdx.x / 32;
#pragma unroll
    for (int i = 0; i < NV; ++i) {
        for (int o = 16; o > 0; o >>= 1) v[i] += __shfl_down_sync(0xffffffffu, v[i], o);
        if (lane == 0) sh[i][warp] = v[i];
    }
    __syncthreads();
    if (threadIdx.x != 0) return;
#pragma unroll
    for (int i = 0; i < NV; ++i) {
        v[i] = 0.0;
        for (int w = 0; w < kThreads / 32; ++w) v[i] += sh[i][w];
    }
}

template <int NV>
__device__ inline void store_partial(double (&v)[NV], double* __restrict__ partial)
{
    block_sum<NV>(v);
    if (threadIdx.x == 0)
#pragma unroll
        for (int i = 0; i < NV; ++i) partial[blockIdx.x * NV + i] = v[i];
}

__global__ void __launch_bounds__(kThreads)
dss_kernel(const double* __restrict__ in, double* __restrict__ out, Block blk)
{
    for_each_node(blk, [&](const Node& nd) { out[nd.t] = gather_sum(in, nd, blk); });
}

// z = P q: gather-DSS, then zero on the outer boundary (homogeneous Dirichlet) and
// divide by the assembled diagonal elsewhere. With AZ / ZZ it also forms the
// weighted dot products <a, z> / <z, z> of the result while z is in registers.
template <bool AZ, bool ZZ>
__global__ void __launch_bounds__(kThreads)
precondition_kernel(const double* __restrict__ q, double* __restrict__ z,
                    const double* __restrict__ diag, const double* __restrict__ a, Block blk,
                    double* __restrict__ partial)
{
    constexpr int NV = AZ + ZZ;
    [[maybe_unused]] double acc[NV > 0 ? NV : 1] = {};
    for_each_node(blk, [&](const Node& nd) {
        const double v = on_outer_boundary(nd, blk) ? 0.0 : gather_sum(q, nd, blk) / diag[nd.t];
        z[nd.t] = v;
        if constexpr (NV > 0) {
            const double w = weight(nd, blk);
            if constexpr (AZ) acc[0] += a[nd.t] * v * w;
            if constexpr (ZZ) acc[NV - 1] += v * v * w;
        }
    });
    if constexpr (NV > 0) store_partial<NV>(acc, partial);
}

__global__ void __launch_bounds__(kThreads)
wdot_kernel(const double* __restrict__ a, const double* __restrict__ b, Block blk,
            double* __restrict__ partial)
{
    double acc[1] = {0.0};
    for_each_node(blk, [&](const Node& nd) { acc[0] += a[nd.t] * b[nd.t] * weight(nd, blk); });
    store_partial<1>(acc, partial);
}

// Krylov scalars live in device memory: the reductions write them and the vector
// updates read them, so no step waits for a round trip to the host.
enum Scalar { kRho, kAlpha, kOmega, kBeta, kRR, kScalars };

// What finish_kernel does with the totals of a reduction: store one, start the
// solve, or form alpha, omega, or beta (with the new rho) for the next update.
enum Step { kStepStore, kStepStart, kStepAlpha, kStepOmega, kStepRho };

// x += alpha p + omega s and r = s - omega t, plus the weighted <r, r> (the stopping
// test) and <rh, r> (the next rho).
__global__ void __launch_bounds__(kThreads)
update_xr_kernel(double* __restrict__ x, double* __restrict__ r, const double* __restrict__ p,
                 const double* __restrict__ s, const double* __restrict__ tt,
                 const double* __restrict__ rh, const double* __restrict__ sc, Block blk,
                 double* __restrict__ partial)
{
    const double alpha = sc[kAlpha], omega = sc[kOmega];
    double acc[2] = {0.0, 0.0};
    for_each_node(blk, [&](const Node& nd) {
        const long long i = nd.t;
        x[i] += alpha * p[i] + omega * s[i];
        const double ri = s[i] - omega * tt[i];
        r[i] = ri;
        const double w = weight(nd, blk);
        acc[0] += ri * ri * w;
        acc[1] += rh[i] * ri * w;
    });
    store_partial<2>(acc, partial);
}

template <Step STEP, int NV>
__global__ void __launch_bounds__(kThreads)
finish_kernel(const double* __restrict__ partial, int blocks, double* __restrict__ sc, int slot)
{
    double v[NV] = {};
    for (int i = threadIdx.x; i < blocks; i += kThreads)
#pragma unroll
        for (int j = 0; j < NV; ++j) v[j] += partial[i * NV + j];
    block_sum<NV>(v);
    if (threadIdx.x != 0) return;
    if constexpr (STEP == kStepStore) {
        sc[slot] = v[0];
    } else if constexpr (STEP == kStepStart) {   // rh = r0, so rho = <r0, r0>
        sc[kRR] = v[0];
        sc[kRho] = v[0];
        sc[kAlpha] = 1.0;
        sc[kOmega] = 1.0;
        sc[kBeta] = 0.0;
    } else if constexpr (STEP == kStepAlpha) {
        sc[kAlpha] = sc[kRho] / v[0];
    } else if constexpr (STEP == kStepOmega) {
        sc[kOmega] = v[0] / v[1];
    } else {
        sc[kRR] = v[0];
        sc[kBeta] = (v[1] / sc[kRho]) * (sc[kAlpha] / sc[kOmega]);
        sc[kRho] = v[1];
    }
}

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

__global__ void fill_kernel(double* __restrict__ v, double value, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) v[t] = value;
}

// One wave of thread blocks: as many as the GPU keeps resident for this kernel.
// A second, partial wave would leave most of the GPU idle while it runs.
template <typename Kernel>
int resident_grid(Kernel kernel)
{
    int dev = 0, sms = 0, per_sm = 0;
    MARS_CELLWISE_CK(cudaGetDevice(&dev));
    MARS_CELLWISE_CK(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev));
    MARS_CELLWISE_CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, kernel, kThreads, 0));
    return sms * per_sm;
}

// Scratch for the weighted dot products: per-block partial sums, the Krylov scalars,
// and the grid of each reducing kernel.
struct Reduction {
    double* partial = nullptr;
    double* scalars = nullptr;
    int grid_dot, grid_start, grid_alpha, grid_omega, grid_xr;
    Reduction()
        : grid_dot(resident_grid(wdot_kernel)),
          grid_start(resident_grid(precondition_kernel<false, true>)),
          grid_alpha(resident_grid(precondition_kernel<true, false>)),
          grid_omega(resident_grid(precondition_kernel<true, true>)),
          grid_xr(resident_grid(update_xr_kernel))
    {
        const int most = std::max({grid_dot, grid_start, grid_alpha, grid_omega, grid_xr});
        MARS_CELLWISE_CK(cudaMalloc(&partial, 2 * most * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&scalars, kScalars * sizeof(double)));
    }
    ~Reduction()
    {
        cudaFree(partial);
        cudaFree(scalars);
    }
    Reduction(const Reduction&) = delete;
    Reduction& operator=(const Reduction&) = delete;
};

// out = DSS(in): every copy of a node gets the sum of all its copies.
inline void dss(const double* d_in, double* d_out, const Block& b, cudaStream_t stream = 0)
{
    static const int grid = resident_grid(dss_kernel);
    dss_kernel<<<grid, kThreads, 0, stream>>>(d_in, d_out, b);
    MARS_CELLWISE_CK(cudaGetLastError());
}

// z = P q. The result is continuous, as the equivalence with the assembled solve requires.
inline void precondition(const double* d_q, double* d_z, const double* d_diag, const Block& b,
                         cudaStream_t stream = 0)
{
    static const int grid = resident_grid(precondition_kernel<false, false>);
    precondition_kernel<false, false><<<grid, kThreads, 0, stream>>>(d_q, d_z, d_diag, nullptr, b,
                                                                     nullptr);
    MARS_CELLWISE_CK(cudaGetLastError());
}

// Weighted (= assembled) dot product of two continuous fields into scalar `slot`.
inline void wdot(const double* d_a, const double* d_b, const Block& b, Reduction& red, int slot,
                 cudaStream_t stream = 0)
{
    wdot_kernel<<<red.grid_dot, kThreads, 0, stream>>>(d_a, d_b, b, red.partial);
    finish_kernel<kStepStore, 1><<<1, kThreads, 0, stream>>>(red.partial, red.grid_dot,
                                                            red.scalars, slot);
    MARS_CELLWISE_CK(cudaGetLastError());
}

// The square root of scalar kRR: the one value the host needs, to decide when to stop.
inline double read_norm(const Reduction& red, cudaStream_t stream)
{
    double h = 0.0;
    MARS_CELLWISE_CK(cudaMemcpyAsync(&h, red.scalars + kRR, sizeof(double),
                                     cudaMemcpyDeviceToHost, stream));
    MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
    return std::sqrt(h);
}

struct Workspace {
    double *r = nullptr, *rh = nullptr, *p = nullptr, *v = nullptr, *s = nullptr, *t = nullptr,
           *q = nullptr;
    Reduction red;
    explicit Workspace(long long n)
    {
        for (double** a : {&r, &rh, &p, &v, &s, &t, &q})
            MARS_CELLWISE_CK(cudaMalloc(a, n * sizeof(double)));
    }
    ~Workspace()
    {
        for (double* a : {r, rh, p, v, s, t, q}) cudaFree(a);
    }
    Workspace(const Workspace&) = delete;
    Workspace& operator=(const Workspace&) = delete;
};

struct SolveResult {
    int iterations = 0;
    double residual0 = 0.0, residual = 0.0;
    std::vector<double> history;   // ||P r||_w after each iteration
};

// Left-preconditioned BiCGStab on P A x = P b, x starting at zero. `apply(u, y)`
// computes y = A u on element-local vectors (u continuous, y unassembled). Per
// iteration: two operator calls, two fused preconditioner-and-dot passes, and three
// vector updates, one of which also forms the two dot products the next iteration
// needs.
template <typename Apply>
SolveResult bicgstab(Apply&& apply, const Block& b, const double* d_rhs, double* d_x,
                     const double* d_diag, double tol, int max_iterations, cudaStream_t stream = 0)
{
    const long long n = b.values();
    const unsigned grid = (unsigned)((n + kThreads - 1) / kThreads);
    Workspace ws(n);
    Reduction& red = ws.red;
    double* sc = red.scalars;

    fill_kernel<<<grid, kThreads, 0, stream>>>(d_x, 0.0, n);
    fill_kernel<<<grid, kThreads, 0, stream>>>(ws.p, 0.0, n);
    fill_kernel<<<grid, kThreads, 0, stream>>>(ws.v, 0.0, n);
    precondition_kernel<false, true><<<red.grid_start, kThreads, 0, stream>>>(
        d_rhs, ws.r, d_diag, nullptr, b, red.partial);   // r = P (b - A 0)
    finish_kernel<kStepStart, 1><<<1, kThreads, 0, stream>>>(red.partial, red.grid_start, sc, 0);
    MARS_CELLWISE_CK(cudaMemcpyAsync(ws.rh, ws.r, n * sizeof(double), cudaMemcpyDeviceToDevice, stream));
    MARS_CELLWISE_CK(cudaGetLastError());

    SolveResult res;
    res.residual0 = read_norm(red, stream);
    res.residual = res.residual0;
    for (int k = 0; k < max_iterations; ++k) {
        update_p_kernel<<<grid, kThreads, 0, stream>>>(ws.p, ws.r, ws.v, sc, n);
        apply(ws.p, ws.q);
        precondition_kernel<true, false><<<red.grid_alpha, kThreads, 0, stream>>>(
            ws.q, ws.v, d_diag, ws.rh, b, red.partial);
        finish_kernel<kStepAlpha, 1><<<1, kThreads, 0, stream>>>(red.partial, red.grid_alpha, sc, 0);
        update_s_kernel<<<grid, kThreads, 0, stream>>>(ws.s, ws.r, ws.v, sc, n);
        apply(ws.s, ws.q);
        precondition_kernel<true, true><<<red.grid_omega, kThreads, 0, stream>>>(
            ws.q, ws.t, d_diag, ws.s, b, red.partial);
        finish_kernel<kStepOmega, 2><<<1, kThreads, 0, stream>>>(red.partial, red.grid_omega, sc, 0);
        update_xr_kernel<<<red.grid_xr, kThreads, 0, stream>>>(d_x, ws.r, ws.p, ws.s, ws.t, ws.rh,
                                                               sc, b, red.partial);
        finish_kernel<kStepRho, 2><<<1, kThreads, 0, stream>>>(red.partial, red.grid_xr, sc, 0);
        MARS_CELLWISE_CK(cudaGetLastError());
        res.residual = read_norm(red, stream);
        res.history.push_back(res.residual);
        res.iterations = k + 1;
        if (res.residual <= tol * res.residual0) break;
    }
    return res;
}

}  // namespace cellwise
}  // namespace mars
