#pragma once
// Communication-free flexible conjugate gradients on element-local (cell-wise) vectors:
// Algorithm 1 of Wichrowski, "Coalesced Matrix-Free Finite Elements in Cell-Wise Storage"
// (internal-notes/FlexibleCG.pdf), for a symmetric operator (e.g. the Galerkin BP1/BP3
// operators MARSIR generates from marsir-compiler/specs/bp*.op).
//
// The iterate x, the search direction p and the preconditioned residual z are
// continuous (every copy of a node holds the same value); the residual r and q = A p
// are unassembled. Every scalar pairs one of each, and for a continuous u the plain sum
// over all copies <f, u> equals the assembled pairing (the paper's Lemma 1), so the
// iteration needs no communication outside the preconditioner. The flexible beta =
// -<q, z_new> / <q, p> keeps the search directions A-conjugate when the preconditioner
// changes between iterations (e.g. a reduced-precision V-cycle).
//
// The Krylov scalars stay on the device: the reductions write them and the vector
// updates read them. The host reads one value per iteration, for the stopping test
// sqrt(<r, z> / <r0, z0>) <= tol. Reference: marsir-mlir/test/cellwise_bp_ref.py.

#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <cmath>
#include <vector>

namespace mars {
namespace cellwise {

namespace fcg_detail {

enum Scalar { kRz, kQp, kBeta, kRz0, kT0, kT1, kCount };
// What a finished reduction does with its totals: start (rz0 = rz), store <q, p>, or
// form beta from <q, z> and store the new <r, z>.
enum Step { kStart, kStoreQp, kNext };

template <Step STEP>
__device__ inline void apply_step(double* __restrict__ sc, const double* v)
{
    if constexpr (STEP == kStart) {
        sc[kRz] = v[0];
        sc[kRz0] = v[0];
        sc[kBeta] = 0.0;
    } else if constexpr (STEP == kStoreQp) {
        sc[kQp] = v[0];
    } else {
        sc[kBeta] = -v[1] / sc[kQp];
        sc[kRz] = v[0];
    }
}

// Sums the block partials. GLOBAL (one rank): apply the step; otherwise leave the
// rank's totals in kT0 / kT1 for MPI.
template <Step STEP, int NV, bool GLOBAL>
__global__ void __launch_bounds__(kThreads)
finish_kernel(const double* __restrict__ partial, int blocks, double* __restrict__ sc)
{
    double v[2] = {0.0, 0.0};
    for (int i = threadIdx.x; i < blocks; i += kThreads)
#pragma unroll
        for (int j = 0; j < NV; ++j) v[j] += partial[i * NV + j];
    double w[NV];
#pragma unroll
    for (int j = 0; j < NV; ++j) w[j] = v[j];
    block_sum<NV>(w);
    if (threadIdx.x != 0) return;
#pragma unroll
    for (int j = 0; j < NV; ++j) v[j] = w[j];
    if constexpr (GLOBAL) {
        apply_step<STEP>(sc, v);
    } else {
        sc[kT0] = v[0];
        sc[kT1] = v[1];
    }
}

template <Step STEP>
__global__ void step_kernel(double* __restrict__ sc)
{
    const double v[2] = {sc[kT0], sc[kT1]};
    apply_step<STEP>(sc, v);
}

// Block partials of <a, b> (a plain sum over all copies).
__global__ void __launch_bounds__(kThreads)
pair_kernel(const double* __restrict__ a, const double* __restrict__ b, long long n, double* __restrict__ partial)
{
    double acc[1] = {0.0};
    for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += (long long)gridDim.x * blockDim.x)
        acc[0] += a[i] * b[i];
    store_partial<1>(acc, partial);
}

// x += alpha p, r -= alpha q, alpha = <r, z> / <q, p>.
__global__ void __launch_bounds__(kThreads)
update_kernel(double* __restrict__ x, double* __restrict__ r, const double* __restrict__ p,
              const double* __restrict__ q, const double* __restrict__ sc, long long n)
{
    const double alpha = sc[kRz] / sc[kQp];
    for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += (long long)gridDim.x * blockDim.x) {
        x[i] += alpha * p[i];
        r[i] -= alpha * q[i];
    }
}

// p = z + beta p.
__global__ void __launch_bounds__(kThreads)
direction_kernel(double* __restrict__ p, const double* __restrict__ z, const double* __restrict__ sc, long long n)
{
    const double beta = sc[kBeta];
    for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += (long long)gridDim.x * blockDim.x)
        p[i] = z[i] + beta * p[i];
}

// z = s * inv_diag (inv_diag: 1 / the assembled diagonal, 0 on Dirichlet copies), with
// the block partials of <r, z> and <q, z>.
__global__ void __launch_bounds__(kThreads)
jacobi_kernel(const double* __restrict__ s, const double* __restrict__ inv_diag, const double* __restrict__ r,
              const double* __restrict__ q, double* __restrict__ z, long long n, double* __restrict__ partial)
{
    double acc[2] = {0.0, 0.0};
    for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += (long long)gridDim.x * blockDim.x) {
        const double zi = s[i] * inv_diag[i];
        z[i] = zi;
        acc[0] += r[i] * zi;
        acc[1] += q[i] * zi;
    }
    store_partial<2>(acc, partial);
}

}  // namespace fcg_detail

// The scalars, the block partials of the reductions, and the communicator.
struct FcgReduction {
    Reduction base;   // partial buffer, pinned host values, communicator
    double* sc = nullptr;
    explicit FcgReduction(MPI_Comm comm) : base(comm)
    {
        MARS_CELLWISE_CK(cudaMalloc(&sc, fcg_detail::kCount * sizeof(double)));
    }
    ~FcgReduction() { cudaFree(sc); }
    FcgReduction(const FcgReduction&) = delete;
    FcgReduction& operator=(const FcgReduction&) = delete;

    template <fcg_detail::Step STEP, int NV>
    void finish(int blocks, cudaStream_t stream)
    {
        using namespace fcg_detail;
        if (base.ranks == 1) {
            finish_kernel<STEP, NV, true><<<1, kThreads, 0, stream>>>(base.partial, blocks, sc);
            MARS_CELLWISE_CK(cudaGetLastError());
            return;
        }
        finish_kernel<STEP, NV, false><<<1, kThreads, 0, stream>>>(base.partial, blocks, sc);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(
            cudaMemcpyAsync(base.host, sc + fcg_detail::kT0, 2 * sizeof(double), cudaMemcpyDeviceToHost, stream));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, base.host, NV, MPI_DOUBLE, MPI_SUM, base.comm));
        MARS_CELLWISE_CK(
            cudaMemcpyAsync(sc + fcg_detail::kT0, base.host, 2 * sizeof(double), cudaMemcpyHostToDevice, stream));
        step_kernel<STEP><<<1, 1, 0, stream>>>(sc);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

    // sqrt(<r, z> / <r0, z0>), the stopping test.
    double relative(cudaStream_t stream) const
    {
        double h[fcg_detail::kCount];
        MARS_CELLWISE_CK(cudaMemcpyAsync(h, sc, sizeof(h), cudaMemcpyDeviceToHost, stream));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        return std::sqrt(std::fmax(h[fcg_detail::kRz], 0.0) / h[fcg_detail::kRz0]);
    }
};

// The DSS-Jacobi preconditioner: z = inv_diag * DSS(r), with the block partials of
// <r, z> and <q, z> (q is zero at the first call). dss(in, out, stream) is any DSS of
// the n local copies (the hand dss<P>, or the MARSIR-generated one).
template <class Dss>
struct DssJacobi {
    Dss dss;
    const double* inv_diag;
    double* s;   // scratch: DSS(r)
    long long n;

    int operator()(const double* r, const double* q, double* z, double* partial, cudaStream_t stream) const
    {
        dss(r, s, stream);
        static const int grid = resident_grid(fcg_detail::jacobi_kernel);
        fcg_detail::jacobi_kernel<<<grid, kThreads, 0, stream>>>(s, inv_diag, r, q, z, n, partial);
        MARS_CELLWISE_CK(cudaGetLastError());
        return grid;
    }
};

struct FcgResult {
    int iterations = 0;
    double relative = 0.0;
    std::vector<double> history;   // sqrt(<r, z> / <r0, z0>) after each iteration
};

// Solves A x = b from x = 0: apply(u, y, stream) is y = A u (u continuous, y
// unassembled); precond(r, q, z, partial, stream) as DssJacobi. n: local copies. Work
// vectors r, q, z, p of n values each.
template <class Apply, class Precond>
FcgResult fcg(Apply&& apply, const Precond& precond, const double* b, double* x, long long n, double* r, double* q,
              double* z, double* p, FcgReduction& red, double tol, int max_iterations, cudaStream_t stream)
{
    using namespace fcg_detail;
    const long long bytes = n * (long long)sizeof(double);
    MARS_CELLWISE_CK(cudaMemsetAsync(x, 0, bytes, stream));
    MARS_CELLWISE_CK(cudaMemsetAsync(q, 0, bytes, stream));
    MARS_CELLWISE_CK(cudaMemsetAsync(p, 0, bytes, stream));
    MARS_CELLWISE_CK(cudaMemcpyAsync(r, b, bytes, cudaMemcpyDeviceToDevice, stream));
    static const int pair_grid = resident_grid(pair_kernel);
    static const int update_grid = resident_grid(update_kernel);
    static const int dir_grid = resident_grid(direction_kernel);

    red.finish<kStart, 2>(precond(r, q, z, red.base.partial, stream), stream);   // z0 = P r0
    direction_kernel<<<dir_grid, kThreads, 0, stream>>>(p, z, red.sc, n);       // p = z0 (beta = 0)
    FcgResult res;
    for (int k = 0; k < max_iterations; ++k) {
        apply(p, q, stream);   // q = A p
        pair_kernel<<<pair_grid, kThreads, 0, stream>>>(q, p, n, red.base.partial);
        red.finish<kStoreQp, 1>(pair_grid, stream);
        update_kernel<<<update_grid, kThreads, 0, stream>>>(x, r, p, q, red.sc, n);
        red.finish<kNext, 2>(precond(r, q, z, red.base.partial, stream), stream);
        direction_kernel<<<dir_grid, kThreads, 0, stream>>>(p, z, red.sc, n);
        MARS_CELLWISE_CK(cudaGetLastError());
        res.relative = red.relative(stream);
        res.history.push_back(res.relative);
        res.iterations = k + 1;
        if (res.relative <= tol) break;
    }
    return res;
}

}  // namespace cellwise
}  // namespace mars
