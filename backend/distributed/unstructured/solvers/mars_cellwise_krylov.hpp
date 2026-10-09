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
// The DSS is a one-pass gather (mars_cellwise_layout.hpp) fused with the Jacobi step
// and the dot products after it: a preconditioner call reads the residual once and
// writes the result once. On GH200 the fused preconditioner took 1.6 ms against
// 2.5 ms for the paper's three-pass cascade + Jacobi at 64^3 elements, although the
// gather alone is slower than the cascade alone (1.3 vs 1.0 ms). On several ranks
// the gather runs in two passes around the ghost exchange, and every global sum
// goes through MPI_Allreduce, so all ranks use the same Krylov scalars and stop on
// the same iteration.

#include "backend/distributed/unstructured/solvers/mars_cellwise_layout.hpp"

#include <cuda_runtime.h>
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace mars {
namespace cellwise {

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

template <bool SHELL>
__global__ void __launch_bounds__(kThreads)
dss_kernel(const double* __restrict__ in, double* __restrict__ out, Block blk, ElementSet set,
           GhostView ghost)
{
    const auto value = [in](long long i) { return in[i]; };
    for_each_node<SHELL>(blk, set, [&](const Node& nd) { out[nd.t] = gather_value<SHELL>(value, ghost, nd, blk); });
}

// z = P q: gather-DSS, then zero on the outer boundary (homogeneous Dirichlet) and
// divide by the assembled diagonal elsewhere. With AZ / ZZ it also forms the
// weighted dot products <a, z> / <z, z> of the result while z is in registers.
template <bool AZ, bool ZZ, bool SHELL>
__global__ void __launch_bounds__(kThreads)
precondition_kernel(const double* __restrict__ q, double* __restrict__ z,
                    const double* __restrict__ diag, const double* __restrict__ a, Block blk,
                    ElementSet set, GhostView ghost, double* __restrict__ partial)
{
    constexpr int NV = AZ + ZZ;
    [[maybe_unused]] double acc[NV > 0 ? NV : 1] = {};
    const auto value = [q](long long i) { return q[i]; };
    for_each_node<SHELL>(blk, set, [&](const Node& nd) {
        const double v = on_outer_boundary(nd, blk) ? 0.0
                                                    : gather_value<SHELL>(value, ghost, nd, blk) / diag[nd.t];
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

// The weighted dot products <a, z> (AZ) and <z, z> (ZZ) of a preconditioned vector,
// for preconditioners that cannot form them on the fly.
template <bool AZ, bool ZZ>
__global__ void __launch_bounds__(kThreads)
wdots_kernel(const double* __restrict__ z, const double* __restrict__ a, Block blk,
             double* __restrict__ partial)
{
    constexpr int NV = AZ + ZZ;
    double acc[NV] = {};
    for_each_node(blk, [&](const Node& nd) {
        const double v = z[nd.t], w = weight(nd, blk);
        if constexpr (AZ) acc[0] += a[nd.t] * v * w;
        if constexpr (ZZ) acc[NV - 1] += v * v * w;
    });
    store_partial<NV>(acc, partial);
}

// Krylov scalars live in device memory: the reductions write them and the vector
// updates read them. kT0/kT1 hold a rank's totals on their way through MPI.
enum Scalar { kRho, kAlpha, kOmega, kBeta, kRR, kT0, kT1, kScalars };

// What a finished reduction does with its global totals: store one, start the
// solve, or form alpha, omega, or beta (with the new rho) for the next update.
enum Step { kStepStore, kStepStart, kStepAlpha, kStepOmega, kStepRho };

template <Step STEP>
__device__ inline void apply_step(double* __restrict__ sc, const double* v, int slot)
{
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

// Sums the block partials of a reduction. GLOBAL (one rank): the totals are final, so
// apply the step. Otherwise leave the rank's totals in kT0/kT1 for MPI.
template <Step STEP, int NV, bool GLOBAL>
__global__ void __launch_bounds__(kThreads)
finish_kernel(const double* __restrict__ partial, int blocks, double* __restrict__ sc, int slot)
{
    double v[NV] = {};
    for (int i = threadIdx.x; i < blocks; i += kThreads)
#pragma unroll
        for (int j = 0; j < NV; ++j) v[j] += partial[i * NV + j];
    block_sum<NV>(v);
    if (threadIdx.x != 0) return;
    if constexpr (GLOBAL) {
        apply_step<STEP>(sc, v, slot);
    } else {
#pragma unroll
        for (int j = 0; j < NV; ++j) sc[kT0 + j] = v[j];
    }
}

template <Step STEP>
__global__ void step_kernel(double* __restrict__ sc, int slot)
{
    const double v[2] = {sc[kT0], sc[kT1]};
    apply_step<STEP>(sc, v, slot);
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

// Scratch for the weighted dot products: per-block partial sums (two values per
// block for two passes of any resident grid, at most 32 blocks per SM), the Krylov
// scalars, and the communicator the totals are summed over.
struct Reduction {
    double* partial = nullptr;
    double* scalars = nullptr;
    double* host = nullptr;   // pinned, for the totals on their way through MPI
    MPI_Comm comm;
    int ranks = 1;
    int grid_dot, grid_xr;
    explicit Reduction(MPI_Comm c = MPI_COMM_SELF)
        : comm(c), grid_dot(resident_grid(wdot_kernel)), grid_xr(resident_grid(update_xr_kernel))
    {
        MARS_CELLWISE_MPI(MPI_Comm_size(comm, &ranks));
        int dev = 0, sms = 0;
        MARS_CELLWISE_CK(cudaGetDevice(&dev));
        MARS_CELLWISE_CK(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev));
        MARS_CELLWISE_CK(cudaMalloc(&partial, 4 * 32 * sms * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&scalars, kScalars * sizeof(double)));
        MARS_CELLWISE_CK(cudaMallocHost(&host, 2 * sizeof(double)));
    }
    ~Reduction()
    {
        cudaFree(partial);
        cudaFree(scalars);
        cudaFreeHost(host);
    }
    Reduction(const Reduction&) = delete;
    Reduction& operator=(const Reduction&) = delete;

    // Global totals of NV values from `blocks` block partials, then STEP.
    template <Step STEP, int NV>
    void finish(int blocks, int slot, cudaStream_t stream)
    {
        if (ranks == 1) {
            finish_kernel<STEP, NV, true><<<1, kThreads, 0, stream>>>(partial, blocks, scalars, slot);
            MARS_CELLWISE_CK(cudaGetLastError());
            return;
        }
        finish_kernel<STEP, NV, false><<<1, kThreads, 0, stream>>>(partial, blocks, scalars, slot);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaMemcpyAsync(host, scalars + kT0, NV * sizeof(double),
                                         cudaMemcpyDeviceToHost, stream));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, host, NV, MPI_DOUBLE, MPI_SUM, comm));
        MARS_CELLWISE_CK(cudaMemcpyAsync(scalars + kT0, host, NV * sizeof(double),
                                         cudaMemcpyHostToDevice, stream));
        step_kernel<STEP><<<1, 1, 0, stream>>>(scalars, slot);
        MARS_CELLWISE_CK(cudaGetLastError());
    }
};

// out = DSS(in): every copy of a node gets the sum of all its copies. `halo` is null
// on one rank.
inline void dss(const double* d_in, double* d_out, const Block& b, Halo* halo,
                cudaStream_t stream = 0)
{
    run_gather_pass(halo, d_in, nullptr, stream, [&](const ElementSet& set, bool shell, int) {
        if (shell) {
            static const int grid = resident_grid(dss_kernel<true>);
            dss_kernel<true><<<grid, kThreads, 0, stream>>>(d_in, d_out, b, set, halo->view());
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        }
        static const int grid = resident_grid(dss_kernel<false>);
        dss_kernel<false><<<grid, kThreads, 0, stream>>>(d_in, d_out, b, set, GhostView{});
        MARS_CELLWISE_CK(cudaGetLastError());
        return grid;
    });
}

// DSS + Jacobi, fused with the dot products that follow it. A preconditioner maps an
// unassembled q to a continuous z = P q, writes the block partials of <a, z> (AZ) and
// <z, z> (ZZ) to `partial`, and returns how many blocks wrote them. `halo` is null
// on one rank.
struct JacobiPreconditioner {
    const double* diag;
    Block blk;
    Halo* halo;

    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial,
                   cudaStream_t stream) const
    {
        constexpr int NV = AZ + ZZ;
        return run_gather_pass(halo, q, nullptr, stream, [&](const ElementSet& set, bool shell, int first) {
            double* p = partial ? partial + (long long)first * (NV > 0 ? NV : 1) : nullptr;
            if (shell) return launch<AZ, ZZ, true>(q, z, a, set, halo->view(), p, stream);
            return launch<AZ, ZZ, false>(q, z, a, set, GhostView{}, p, stream);
        });
    }

private:
    template <bool AZ, bool ZZ, bool SHELL>
    int launch(const double* q, double* z, const double* a, const ElementSet& set,
               const GhostView& ghost, double* partial, cudaStream_t stream) const
    {
        static const int grid = resident_grid(precondition_kernel<AZ, ZZ, SHELL>);
        precondition_kernel<AZ, ZZ, SHELL><<<grid, kThreads, 0, stream>>>(q, z, diag, a, blk, set, ghost,
                                                                         partial);
        MARS_CELLWISE_CK(cudaGetLastError());
        return grid;
    }
};

// z = P q without dot products. The result is continuous, as the equivalence with the
// assembled solve requires.
inline void precondition(const double* d_q, double* d_z, const double* d_diag, const Block& b,
                         Halo* halo, cudaStream_t stream = 0)
{
    JacobiPreconditioner{d_diag, b, halo}.operator()<false, false>(d_q, d_z, nullptr, nullptr, stream);
}

// Weighted (= assembled) dot product of two continuous fields into scalar `slot`.
inline void wdot(const double* d_a, const double* d_b, const Block& b, Reduction& red, int slot,
                 cudaStream_t stream = 0)
{
    wdot_kernel<<<red.grid_dot, kThreads, 0, stream>>>(d_a, d_b, b, red.partial);
    MARS_CELLWISE_CK(cudaGetLastError());
    red.finish<kStepStore, 1>(red.grid_dot, slot, stream);
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
    Workspace(long long n, MPI_Comm comm) : red(comm)
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
// computes y = A u on this rank's element-local vectors (u continuous, y
// unassembled); `precond` is a preconditioner as JacobiPreconditioner above; `comm`
// holds the ranks that share the block. Per iteration: two operator calls, two
// preconditioner-and-dot passes, three vector updates (one also forms the two dot
// products the next iteration needs), and three global sums.
template <typename Apply, typename Precond>
SolveResult bicgstab(Apply&& apply, const Precond& precond, const Block& b, const double* d_rhs,
                     double* d_x, MPI_Comm comm, double tol, int max_iterations,
                     cudaStream_t stream = 0)
{
    const long long n = b.values();
    const unsigned grid = (unsigned)((n + kThreads - 1) / kThreads);
    Workspace ws(n, comm);
    Reduction& red = ws.red;
    double* sc = red.scalars;

    fill_kernel<<<grid, kThreads, 0, stream>>>(d_x, 0.0, n);
    fill_kernel<<<grid, kThreads, 0, stream>>>(ws.p, 0.0, n);
    fill_kernel<<<grid, kThreads, 0, stream>>>(ws.v, 0.0, n);
    int blocks = precond.template operator()<false, true>(d_rhs, ws.r, nullptr, red.partial,
                                                          stream);   // r = P (b - A 0)
    red.finish<kStepStart, 1>(blocks, 0, stream);
    MARS_CELLWISE_CK(cudaMemcpyAsync(ws.rh, ws.r, n * sizeof(double), cudaMemcpyDeviceToDevice, stream));
    MARS_CELLWISE_CK(cudaGetLastError());

    SolveResult res;
    res.residual0 = read_norm(red, stream);
    res.residual = res.residual0;
    for (int k = 0; k < max_iterations; ++k) {
        update_p_kernel<<<grid, kThreads, 0, stream>>>(ws.p, ws.r, ws.v, sc, n);
        apply(ws.p, ws.q);
        blocks = precond.template operator()<true, false>(ws.q, ws.v, ws.rh, red.partial, stream);
        red.finish<kStepAlpha, 1>(blocks, 0, stream);
        update_s_kernel<<<grid, kThreads, 0, stream>>>(ws.s, ws.r, ws.v, sc, n);
        apply(ws.s, ws.q);
        blocks = precond.template operator()<true, true>(ws.q, ws.t, ws.s, red.partial, stream);
        red.finish<kStepOmega, 2>(blocks, 0, stream);
        update_xr_kernel<<<red.grid_xr, kThreads, 0, stream>>>(d_x, ws.r, ws.p, ws.s, ws.t, ws.rh,
                                                               sc, b, red.partial);
        red.finish<kStepRho, 2>(red.grid_xr, 0, stream);
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
