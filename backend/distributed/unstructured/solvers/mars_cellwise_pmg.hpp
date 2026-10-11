#pragma once
// p-multigrid preconditioner of the unstructured cell-wise solve: levels p = 7 -> 3 -> 1
// on the same mesh and the same tables, the Knaus CVFEM operator rediscretized at each
// degree (the MARSIR kernel at p = 7, mars_cellwise_operator.hpp below), Chebyshev-Jacobi
// smoothing through the DSS, element-local GLL interpolation between degrees
// (restriction is its transpose on the raw residual), and Chebyshev steps on the p = 1
// level as its coarse solve. One V-cycle per BiCGStab preconditioner call.
// Design: docs/design/cellwise_unstructured_dss.md ("Multigrid on unstructured meshes");
// reference: marsir-mlir/test/cellwise_pmg_ref.py, which it must match.

#include "backend/distributed/unstructured/solvers/mars_cellwise_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <cuda_runtime.h>

#include <array>
#include <cmath>
#include <vector>

namespace mars {
namespace cellwise {

// I[a][j]: coarse GLL basis function j at fine GLL node a.
template <int PF, int PC>
struct Interpolation {
    double I[PF + 1][PC + 1];
};

template <int PF, int PC>
Interpolation<PF, PC> interpolation()
{
    const fem::HoCvfemOperators f = fem::buildHoCvfemOperators(PF), c = fem::buildHoCvfemOperators(PC);
    Interpolation<PF, PC> T;
    for (int a = 0; a <= PF; ++a)
        for (int j = 0; j <= PC; ++j) {
            double v = 1.0;
            for (int k = 0; k <= PC; ++k)
                if (k != j) v *= (f.zeta[a] - c.zeta[k]) / (c.zeta[j] - c.zeta[k]);
            T.I[a][j] = v;
        }
    return T;
}

namespace pmg_detail {

// uf (+)= I (x) I (x) I uc, element by element.
template <int PF, int PC, bool ADD>
__global__ void prolong_kernel(const double* __restrict__ uc, double* __restrict__ uf, long long E,
                               Interpolation<PF, PC> T)
{
    constexpr int nf = PF + 1, nc = PC + 1;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * nf * nf * nf) return;
    const int l = (int)(t % (nf * nf * nf)), a = l / (nf * nf), b = (l / nf) % nf, c = l % nf;
    const double* u = uc + (t / (nf * nf * nf)) * nc * nc * nc;
    double s = 0.0;
    for (int j = 0; j < nc; ++j)
        for (int k = 0; k < nc; ++k) {
            const double w = T.I[a][j] * T.I[b][k];
            for (int m = 0; m < nc; ++m) s += w * T.I[c][m] * u[(j * nc + k) * nc + m];
        }
    uf[t] = ADD ? uf[t] + s : s;
}

// rc = (I (x) I (x) I)^T rf, element by element.
template <int PF, int PC>
__global__ void restrict_kernel(const double* __restrict__ rf, double* __restrict__ rc, long long E,
                                Interpolation<PF, PC> T)
{
    constexpr int nf = PF + 1, nc = PC + 1;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * nc * nc * nc) return;
    const int l = (int)(t % (nc * nc * nc)), j = l / (nc * nc), k = (l / nc) % nc, m = l % nc;
    const double* r = rf + (t / (nc * nc * nc)) * nf * nf * nf;
    double s = 0.0;
    for (int a = 0; a < nf; ++a)
        for (int b = 0; b < nf; ++b) {
            const double w = T.I[a][j] * T.I[b][k];
            for (int c = 0; c < nf; ++c) s += w * T.I[c][m] * r[(a * nf + b) * nf + c];
        }
    rc[t] = s;
}

__global__ void axpby_kernel(double* __restrict__ out, const double* __restrict__ a, const double* __restrict__ b,
                             double beta, long long n)   // out = a + beta b (b may be null: out = a)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] = b ? a[t] + beta * b[t] : a[t];
}

__global__ void scale_kernel(double* __restrict__ v, long long n, double s)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) v[t] *= s;
}

// Deterministic start of the power iteration (the same hash as the structured multigrid).
__global__ void hash_kernel(double* __restrict__ v, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) v[t] = (double)(((unsigned long long)t * 2654435761ULL) & 0xffffffffULL) / 4294967296.0 - 0.5;
}

inline unsigned blocks(long long n) { return (unsigned)((n + kThreads - 1) / kThreads); }

}  // namespace pmg_detail

// apply7(u, y, stream): y = A u at p = 7 (the MARSIR kernel). G3, G1: the p = 3 and p = 1
// metrics (mars_cellwise_operator.hpp layout); halos[0..2]: the exchanges at p = 7, 3, 1
// (all null on one rank).
template <class Apply7>
class PMultigrid {
public:
    PMultigrid(Apply7 apply7, const UnstructuredTopology& T, const double* diag7, const double* G3, const double* G1,
               std::array<UnstructuredHalo*, 3> halos, int pre, int post, int coarse_steps, MPI_Comm comm,
               cudaStream_t stream)
        : apply7_(apply7), T_(view(T)), E_(T.local), diag7_(diag7), G3_(G3), G1_(G1), halos_(halos), pre_(pre),
          post_(post), coarse_steps_(coarse_steps), stream_(stream), M3_(cvfem_matrices<3>()),
          M1_(cvfem_matrices<1>()), I73_(interpolation<7, 3>()), I31_(interpolation<3, 1>())
    {
        r7_ = alloc(E_ * 512), d7_ = alloc(E_ * 512), y7_ = alloc(E_ * 512);
        b3_ = alloc(E_ * 64), x3_ = alloc(E_ * 64), r3_ = alloc(E_ * 64), d3_ = alloc(E_ * 64), y3_ = alloc(E_ * 64);
        diag3_ = alloc(E_ * 64);
        b1_ = alloc(E_ * 8), x1_ = alloc(E_ * 8), d1_ = alloc(E_ * 8), y1_ = alloc(E_ * 8), diag1_ = alloc(E_ * 8);
        // Assembled diagonals of the coarse operators.
        if (E_ > 0) {
            cvfem_diagonal_kernel<3><<<pmg_detail::blocks(E_ * 64), kThreads, 0, stream_>>>(G3_, y3_, E_, M3_);
            cvfem_diagonal_kernel<1><<<pmg_detail::blocks(E_ * 8), kThreads, 0, stream_>>>(G1_, y1_, E_, M1_);
        }
        MARS_CELLWISE_CK(cudaGetLastError());
        dss<3>(y3_, diag3_, T, halos_[1], stream_);
        dss<1>(y1_, diag1_, T, halos_[2], stream_);
        Reduction red(comm);
        lmax_[0] = power_iteration<7>(red);
        lmax_[1] = power_iteration<3>(red);
        lmax_[2] = power_iteration<1>(red);
    }
    ~PMultigrid()
    {
        for (double* p : {r7_, d7_, y7_, b3_, x3_, r3_, d3_, y3_, diag3_, b1_, x1_, d1_, y1_, diag1_}) cudaFree(p);
    }
    PMultigrid(const PMultigrid&) = delete;
    PMultigrid& operator=(const PMultigrid&) = delete;

    const std::array<double, 3>& lambda_max() const { return lmax_; }

    // z = one V-cycle applied to the unassembled q, then the dots BiCGStab needs.
    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial, cudaStream_t stream) const
    {
        vcycle<7>(q, z);
        if constexpr (AZ || ZZ) {
            const Block flat((int)E_, 1, 1);
            static const int grid = resident_grid(wdots_kernel<AZ, ZZ, CountingWeight>);
            wdots_kernel<AZ, ZZ, CountingWeight><<<grid, kThreads, 0, stream>>>(z, a, flat, CountingWeight{T_.counting},
                                                                                partial);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        }
        return 0;
    }

    // x = V-cycle(b) on level P (x is overwritten).
    template <int P>
    void vcycle(const double* b, double* x) const
    {
        if constexpr (P == 1) {
            chebyshev<1>(b, x, coarse_steps_, true);
        } else {
            constexpr int PC = P == 7 ? 3 : 1, n3 = (P + 1) * (P + 1) * (P + 1);
            double* r = P == 7 ? r7_ : r3_;
            double* bc = PC == 3 ? b3_ : b1_;
            double* xc = PC == 3 ? x3_ : x1_;
            if (pre_ > 0) {
                chebyshev<P>(b, x, pre_, true);
                apply<P>(x, P == 7 ? y7_ : y3_);
                pmg_detail::axpby_kernel<<<pmg_detail::blocks(E_ * n3), kThreads, 0, stream_>>>(
                    r, b, P == 7 ? y7_ : y3_, -1.0, E_ * n3);
            } else {
                pmg_detail::axpby_kernel<<<pmg_detail::blocks(E_ * n3), kThreads, 0, stream_>>>(r, b, nullptr, 0.0,
                                                                                               E_ * n3);
            }
            MARS_CELLWISE_CK(cudaGetLastError());
            zero_dirichlet<P>(r, T_, stream_);
            transfer<P, PC>(r, bc, false);
            vcycle<PC>(bc, xc);
            transfer<P, PC>(xc, x, true);
            if (post_ > 0) chebyshev<P>(b, x, post_, false);
        }
    }

private:
    static double* alloc(long long n)
    {
        double* p = nullptr;
        if (n > 0) MARS_CELLWISE_CK(cudaMalloc(&p, n * sizeof(double)));
        return p;
    }

    template <int P>
    void apply(const double* u, double* y) const
    {
        if constexpr (P == 7) apply7_(u, y, stream_);
        else if constexpr (P == 3) cvfem_apply<3>(u, y, G3_, E_, M3_, stream_);
        else cvfem_apply<1>(u, y, G1_, E_, M1_, stream_);
    }

    // PROLONG = false: restrict the fine raw residual f into c. PROLONG = true: add the
    // interpolated coarse correction c to the fine x (x = interpolant when pre_ == 0).
    template <int P, int PC>
    void transfer(const double* from, double* to, bool prolong) const
    {
        constexpr int nf = P + 1, nc = PC + 1;
        if (E_ == 0) return;
        auto go = [&](const Interpolation<P, PC>& T) {
            if (!prolong)
                pmg_detail::restrict_kernel<P, PC><<<pmg_detail::blocks(E_ * nc * nc * nc), kThreads, 0, stream_>>>(
                    from, to, E_, T);
            else if (pre_ > 0)
                pmg_detail::prolong_kernel<P, PC, true>
                    <<<pmg_detail::blocks(E_ * nf * nf * nf), kThreads, 0, stream_>>>(from, to, E_, T);
            else
                pmg_detail::prolong_kernel<P, PC, false>
                    <<<pmg_detail::blocks(E_ * nf * nf * nf), kThreads, 0, stream_>>>(from, to, E_, T);
            MARS_CELLWISE_CK(cudaGetLastError());
        };
        if constexpr (P == 7) go(I73_);
        else go(I31_);
    }

    template <int P>
    const double* diag() const { return P == 7 ? diag7_ : (P == 3 ? diag3_ : diag1_); }

    // One fused Chebyshev step: d = c1 d + c2 Jacobi(DSS(b - y)), x += d (see
    // unstructured_chebyshev_kernel).
    template <int P, bool FIRST>
    void step(const double* b, const double* y, double* d, double* x, double c1, double c2) const
    {
        UnstructuredHalo* halo = halos_[P == 7 ? 0 : (P == 3 ? 1 : 2)];
        run_unstructured_pass(halo, E_, b, y, stream_, [&](const ElementList& set, const double* ghost, int) {
            static const int grid = resident_grid(unstructured_chebyshev_kernel<P, FIRST>);
            unstructured_chebyshev_kernel<P, FIRST><<<grid, kThreads, 0, stream_>>>(b, y, ghost, diag<P>(), d, x, c1,
                                                                                   c2, T_, set);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        });
    }

    // `steps` Chebyshev-Jacobi steps on A x = b at level P, from x = 0 when zero_x, on the
    // interval [lmax / 15, 1.1 lmax] of the Jacobi-preconditioned operator.
    template <int P>
    void chebyshev(const double* b, double* x, int steps, bool zero_x) const
    {
        double* d = P == 7 ? d7_ : (P == 3 ? d3_ : d1_);
        double* y = P == 7 ? y7_ : (P == 3 ? y3_ : y1_);
        const double hi = 1.1 * lmax_[P == 7 ? 0 : (P == 3 ? 1 : 2)], lo = hi / kRange;
        const double theta = 0.5 * (hi + lo), delta = 0.5 * (hi - lo), sigma = theta / delta;
        double rho = 1.0 / sigma;
        if (!zero_x) apply<P>(x, y);
        step<P, true>(b, zero_x ? nullptr : y, d, x, 0.0, 1.0 / theta);
        for (int i = 1; i < steps; ++i) {
            const double rho_new = 1.0 / (2.0 * sigma - rho);
            apply<P>(x, y);
            step<P, false>(b, y, d, x, rho_new * rho, 2.0 * rho_new / delta);
            rho = rho_new;
        }
    }

    template <int P>
    double norm(Reduction& red, const double* v) const
    {
        static const int grid = resident_grid(level_wdot_kernel<P>);
        level_wdot_kernel<P><<<grid, kThreads, 0, stream_>>>(v, v, T_.counting, E_, red.partial);
        MARS_CELLWISE_CK(cudaGetLastError());
        red.finish<kStepStore, 1>(grid, kRR, stream_);
        return read_norm(red, stream_);
    }

    // The largest eigenvalue of the Jacobi-preconditioned operator at level P: 30 power
    // steps from a hash vector (10 underestimate it, see the structured multigrid). The
    // start depends on the local numbering, so lmax differs in the last digits between
    // rank counts; the smoother's interval has a 10% margin.
    template <int P>
    double power_iteration(Reduction& red) const
    {
        constexpr int n3 = (P + 1) * (P + 1) * (P + 1);
        const long long n = E_ * n3;
        double* d = P == 7 ? d7_ : (P == 3 ? d3_ : d1_);
        double* y = P == 7 ? y7_ : (P == 3 ? y3_ : y1_);
        double* v = alloc(n);
        double* w = alloc(n);
        if (n > 0) pmg_detail::hash_kernel<<<pmg_detail::blocks(n), kThreads, 0, stream_>>>(y, n);
        MARS_CELLWISE_CK(cudaGetLastError());
        step<P, true>(y, nullptr, d, v, 0.0, 1.0);   // v = Jacobi(DSS(hash))
        double lam = 0.0;
        for (int i = 0; i < 30; ++i) {
            apply<P>(v, y);
            step<P, true>(y, nullptr, d, w, 0.0, 1.0);   // w = Jacobi(DSS(A v))
            const double nw = norm<P>(red, w), nv = norm<P>(red, v);
            lam = nv > 0.0 ? nw / nv : 0.0;
            if (nw > 0.0 && n > 0) {
                pmg_detail::axpby_kernel<<<pmg_detail::blocks(n), kThreads, 0, stream_>>>(v, w, nullptr, 0.0, n);
                pmg_detail::scale_kernel<<<pmg_detail::blocks(n), kThreads, 0, stream_>>>(v, n, 1.0 / nw);
                MARS_CELLWISE_CK(cudaGetLastError());
            }
        }
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        cudaFree(v);
        cudaFree(w);
        return lam;
    }

    static constexpr double kRange = 15.0;

    Apply7 apply7_;
    TopologyView T_;
    long long E_;
    const double *diag7_, *G3_, *G1_;
    std::array<UnstructuredHalo*, 3> halos_;
    int pre_, post_, coarse_steps_;
    cudaStream_t stream_;
    CvfemMatrices<3> M3_;
    CvfemMatrices<1> M1_;
    Interpolation<7, 3> I73_;
    Interpolation<3, 1> I31_;
    double *r7_ = nullptr, *d7_ = nullptr, *y7_ = nullptr;
    double *b3_ = nullptr, *x3_ = nullptr, *r3_ = nullptr, *d3_ = nullptr, *y3_ = nullptr, *diag3_ = nullptr;
    double *b1_ = nullptr, *x1_ = nullptr, *d1_ = nullptr, *y1_ = nullptr, *diag1_ = nullptr;
    std::array<double, 3> lmax_{};
};

}  // namespace cellwise
}  // namespace mars
