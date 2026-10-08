#pragma once
// Geometric multigrid preconditioner on element-local vectors, after Wichrowski,
// "Coalesced Matrix-Free Geometric Multigrid on Persistent Cell-Wise Storage"
// (arXiv:2607.03413), here at p = 7 and for the nonsymmetric CVFEM operator.
//
// Levels: the structured block coarsened 2:1 down to one element, all at p = 7, so
// every level runs the same operator kernel with its own metric. Transfers are
// element-local. Prolongation evaluates the coarse polynomial at the child's GLL
// nodes. Restriction is its exact transpose, applied to the raw unassembled residual:
// restricting an assembled residual would count shared nodes several times. Elements
// talk only through the DSS inside the smoother. Smoother: Chebyshev on the
// DSS-Jacobi preconditioned operator, its upper bound from a power iteration (the
// top of that spectrum is nearly real). The one-element level is solved exactly.
// The V-cycle is a fixed linear map with continuous output, as left-preconditioned
// BiCGStab requires. Reference: marsir-mlir/test/cellwise_multigrid_ref.py.

#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

namespace mars {
namespace cellwise {

constexpr int kInner = kN - 2, kInner3 = kInner * kInner * kInner;   // interior nodes of one element
constexpr int kMetric = 3 * (kN - 1) * 3 * kNN;                      // metric doubles per element

// Coarse GLL basis at the GLL nodes of the left (0) and right (1) child: [child][fine][coarse].
__constant__ double c_child_interp[2][kN][kN];

// xf += prolongation of the coarse field xc. One thread block per fine element; the
// three 1D contractions go through shared memory.
__global__ void __launch_bounds__(kThreads)
prolong_add_kernel(const double* __restrict__ xc, double* __restrict__ xf, Block fine)
{
    __shared__ double s0[kN3], s1[kN3];
    const long long e = blockIdx.x;
    const int ez = (int)(e % fine.nz);
    const long long rest = e / fine.nz;
    const int ey = (int)(rest % fine.ny), ex = (int)(rest / fine.ny);
    const int cny = fine.ny / 2, cnz = fine.nz / 2;
    const long long ec = ((long long)(ex / 2) * cny + ey / 2) * cnz + ez / 2;
    const double(*Ix)[kN] = c_child_interp[ex & 1];
    const double(*Iy)[kN] = c_child_interp[ey & 1];
    const double(*Iz)[kN] = c_child_interp[ez & 1];
    for (int l = threadIdx.x; l < kN3; l += kThreads) s0[l] = xc[ec * kN3 + l];
    __syncthreads();
    for (int l = threadIdx.x; l < kN3; l += kThreads) {   // [j][k][m] -> [j][k][c]
        const int jk = l / kN, c = l % kN;
        double s = 0.0;
        for (int m = 0; m < kN; ++m) s += Iz[c][m] * s0[jk * kN + m];
        s1[l] = s;
    }
    __syncthreads();
    for (int l = threadIdx.x; l < kN3; l += kThreads) {   // [j][m][c] -> [j][b][c]
        const int j = l / kNN, b = (l / kN) % kN, c = l % kN;
        double s = 0.0;
        for (int m = 0; m < kN; ++m) s += Iy[b][m] * s1[j * kNN + m * kN + c];
        s0[l] = s;
    }
    __syncthreads();
    for (int l = threadIdx.x; l < kN3; l += kThreads) {   // [m][b][c] -> [a][b][c]
        const int a = l / kNN, bc = l % kNN;
        double s = 0.0;
        for (int m = 0; m < kN; ++m) s += Ix[a][m] * s0[m * kNN + bc];
        xf[e * kN3 + l] += s;
    }
}

// bc = restriction of the fine residual bf - qf, its Dirichlet copies zeroed: the
// transpose of prolongation, summed over the 8 children. One thread block per coarse
// element.
__global__ void __launch_bounds__(kThreads)
restrict_kernel(const double* __restrict__ bf, const double* __restrict__ qf,
                double* __restrict__ bc, Block fine)
{
    __shared__ double s0[kN3], s1[kN3];
    const long long ec = blockIdx.x;
    const int cny = fine.ny / 2, cnz = fine.nz / 2;
    const int EZ = (int)(ec % cnz);
    const long long rest = ec / cnz;
    const int EY = (int)(rest % cny), EX = (int)(rest / cny);
    double acc[kN3 / kThreads] = {};
    for (int child = 0; child < 8; ++child) {
        Node nd;
        nd.ex = 2 * EX + (child >> 2);
        nd.ey = 2 * EY + ((child >> 1) & 1);
        nd.ez = 2 * EZ + (child & 1);
        nd.e = ((long long)nd.ex * fine.ny + nd.ey) * fine.nz + nd.ez;
        const double(*Ix)[kN] = c_child_interp[nd.ex & 1];
        const double(*Iy)[kN] = c_child_interp[nd.ey & 1];
        const double(*Iz)[kN] = c_child_interp[nd.ez & 1];
        for (int l = threadIdx.x; l < kN3; l += kThreads) {
            nd.t = nd.e * kN3 + l;
            nd.a = l / kNN;
            nd.b = (l / kN) % kN;
            nd.c = l % kN;
            s0[l] = on_outer_boundary(nd, fine) ? 0.0 : bf[nd.t] - qf[nd.t];
        }
        __syncthreads();
        for (int l = threadIdx.x; l < kN3; l += kThreads) {   // [a][b][c] -> [j][b][c]
            const int j = l / kNN, bcl = l % kNN;
            double s = 0.0;
            for (int a = 0; a < kN; ++a) s += Ix[a][j] * s0[a * kNN + bcl];
            s1[l] = s;
        }
        __syncthreads();
        for (int l = threadIdx.x; l < kN3; l += kThreads) {   // [j][b][c] -> [j][k][c]
            const int j = l / kNN, k = (l / kN) % kN, c = l % kN;
            double s = 0.0;
            for (int b = 0; b < kN; ++b) s += Iy[b][k] * s1[j * kNN + b * kN + c];
            s0[l] = s;
        }
        __syncthreads();
        int i = 0;
        for (int l = threadIdx.x; l < kN3; l += kThreads, ++i) {   // [j][k][c] -> [j][k][m]
            const int jk = l / kN, m = l % kN;
            double s = 0.0;
            for (int c = 0; c < kN; ++c) s += Iz[c][m] * s0[jk * kN + c];
            acc[i] += s;
        }
        __syncthreads();   // the next child overwrites s0
    }
    int i = 0;
    for (int l = threadIdx.x; l < kN3; l += kThreads, ++i) bc[ec * kN3 + l] = acc[i];
}

// One Chebyshev step on P_J A: z = P_J (b - q) with q = A x, d = c_d d + c_z z, x += d.
// FIRST: d = c_z z, the old d is not read. ZERO_X: x = 0 and q is not used (a smoother
// starting from zero), so z = P_J b and x = d.
template <bool FIRST, bool ZERO_X>
__global__ void __launch_bounds__(kThreads)
chebyshev_kernel(const double* __restrict__ b, const double* __restrict__ q,
                 const double* __restrict__ diag, double* __restrict__ d, double* __restrict__ x,
                 double c_d, double c_z, Block blk)
{
    for_each_node(blk, [&](const Node& nd) {
        double z = 0.0;
        if (!on_outer_boundary(nd, blk)) {
            const double g = ZERO_X ? gather_sum(b, nd, blk)
                                    : gather_sum_of([&](long long i) { return b[i] - q[i]; }, nd, blk);
            z = g / diag[nd.t];
        }
        const double dn = FIRST ? c_z * z : c_d * d[nd.t] + c_z * z;
        d[nd.t] = dn;
        x[nd.t] = ZERO_X ? dn : x[nd.t] + dn;
    });
}

// Local index of interior node j of one element, j = ((a-1) * 6 + (b-1)) * 6 + (c-1).
__host__ __device__ inline int interior_local(int j)
{
    return ((j / (kInner * kInner) + 1) * kN + (j / kInner) % kInner + 1) * kN + j % kInner + 1;
}

// x = A^-1 b on the one-element level: interior nodes through the stored inverse,
// boundary nodes zero.
__global__ void __launch_bounds__(kThreads)
coarse_solve_kernel(const double* __restrict__ ainv, const double* __restrict__ b,
                    double* __restrict__ x)
{
    for (int l = threadIdx.x; l < kN3; l += kThreads) {
        const int a = l / kNN, bb = (l / kN) % kN, c = l % kN;
        double s = 0.0;
        if (a > 0 && a < kN - 1 && bb > 0 && bb < kN - 1 && c > 0 && c < kN - 1) {
            const int i = ((a - 1) * kInner + bb - 1) * kInner + c - 1;
            for (int j = 0; j < kInner3; ++j) s += ainv[i * kInner3 + j] * b[interior_local(j)];
        }
        x[l] = s;
    }
}

// Setup helpers.
__global__ void replicate_kernel(const double* __restrict__ src, double* __restrict__ dst,
                                 long long len, long long copies)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < len * copies) dst[t] = src[t % len];
}

__global__ void interior_units_kernel(double* __restrict__ u)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < (long long)kInner3 * kN3)
        u[t] = (t % kN3) == interior_local((int)(t / kN3)) ? 1.0 : 0.0;
}

// A reproducible start vector for the power iteration, the same in the reference.
__global__ void hash_kernel(double* __restrict__ v, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n)
        v[t] = (double)(((unsigned long long)t * 2654435761ULL) & 0xffffffffULL) / 4294967296.0 - 0.5;
}

__global__ void scale_kernel(double* __restrict__ out, const double* __restrict__ in, double s,
                             long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] = s * in[t];
}

struct MultigridLevel {
    Block blk;
    const double* G;      // metric, as the operator reads it
    const double* diag;   // assembled diagonal on every copy
};

// `apply(u, y, G, E)` computes y = A u on E elements with metric G. Level 0 is the
// finest; each next level halves every block dimension; the last has one element.
template <typename Apply>
class Multigrid {
public:
    Multigrid(Apply apply, std::vector<MultigridLevel> levels, const double* zeta_host,
              int degree = 3, double smoothing_range = 15.0)
        : apply_(std::move(apply)), lv_(std::move(levels)), degree_(degree), range_(smoothing_range)
    {
        for (size_t l = 0; l + 1 < lv_.size(); ++l) {
            const Block f = lv_[l].blk, c = lv_[l + 1].blk;
            if (c.nx * 2 != f.nx || c.ny * 2 != f.ny || c.nz * 2 != f.nz) {
                fprintf(stderr, "multigrid: level %zu is not a 2:1 coarsening\n", l + 1);
                std::abort();
            }
        }
        if (lv_.back().blk.elements() != 1) {
            fprintf(stderr, "multigrid: the coarsest level must be one element\n");
            std::abort();
        }
        upload_child_interp(zeta_host);
        b_.assign(lv_.size(), nullptr);
        x_.assign(lv_.size(), nullptr);
        q_.assign(lv_.size(), nullptr);
        d_.assign(lv_.size(), nullptr);
        for (size_t l = 0; l < lv_.size(); ++l) {
            const long long n = lv_[l].blk.values();
            MARS_CELLWISE_CK(cudaMalloc(&q_[l], n * sizeof(double)));
            MARS_CELLWISE_CK(cudaMalloc(&d_[l], n * sizeof(double)));
            if (l > 0) {
                MARS_CELLWISE_CK(cudaMalloc(&b_[l], n * sizeof(double)));
                MARS_CELLWISE_CK(cudaMalloc(&x_[l], n * sizeof(double)));
            }
        }
        build_coarse_inverse();
        for (size_t l = 0; l + 1 < lv_.size(); ++l) lmax_.push_back(power_iteration(l));
    }
    ~Multigrid()
    {
        for (auto* v : {&b_, &x_, &q_, &d_})
            for (double* p : *v) cudaFree(p);
        cudaFree(ainv_);
    }
    Multigrid(const Multigrid&) = delete;
    Multigrid& operator=(const Multigrid&) = delete;

    // z = one V-cycle applied to the unassembled q, then the dots BiCGStab needs.
    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial,
                   cudaStream_t stream) const
    {
        vcycle(0, q, z, stream);
        if constexpr (AZ || ZZ) {
            static const int grid = resident_grid(wdots_kernel<AZ, ZZ>);
            wdots_kernel<AZ, ZZ><<<grid, kThreads, 0, stream>>>(z, a, lv_[0].blk, partial);
            return grid;
        }
        return 0;
    }

    void vcycle(size_t l, const double* b, double* x, cudaStream_t stream) const
    {
        if (l + 1 == lv_.size()) {
            coarse_solve_kernel<<<1, kThreads, 0, stream>>>(ainv_, b, x);
            return;
        }
        const MultigridLevel& L = lv_[l];
        smooth(l, b, x, true, stream);
        apply_(x, q_[l], L.G, L.blk.elements());
        restrict_kernel<<<(unsigned)lv_[l + 1].blk.elements(), kThreads, 0, stream>>>(b, q_[l],
                                                                                    b_[l + 1], L.blk);
        vcycle(l + 1, b_[l + 1], x_[l + 1], stream);
        prolong_add_kernel<<<(unsigned)L.blk.elements(), kThreads, 0, stream>>>(x_[l + 1], x, L.blk);
        smooth(l, b, x, false, stream);
    }

    double lambda_max(size_t l) const { return lmax_[l]; }

private:
    void smooth(size_t l, const double* b, double* x, bool zero_x, cudaStream_t stream) const
    {
        const MultigridLevel& L = lv_[l];
        const double lmax = 1.1 * lmax_[l], lmin = lmax / range_;
        const double theta = 0.5 * (lmax + lmin), delta = 0.5 * (lmax - lmin);
        const double sigma = theta / delta;
        double rho = 1.0 / sigma;
        if (zero_x) {
            launch<true, true>(b, nullptr, L, d_[l], x, 0.0, 1.0 / theta, stream);
        } else {
            apply_(x, q_[l], L.G, L.blk.elements());
            launch<true, false>(b, q_[l], L, d_[l], x, 0.0, 1.0 / theta, stream);
        }
        for (int k = 1; k < degree_; ++k) {
            const double rho_new = 1.0 / (2.0 * sigma - rho);
            apply_(x, q_[l], L.G, L.blk.elements());
            launch<false, false>(b, q_[l], L, d_[l], x, rho_new * rho, 2.0 * rho_new / delta, stream);
            rho = rho_new;
        }
    }

    template <bool FIRST, bool ZERO_X>
    static void launch(const double* b, const double* q, const MultigridLevel& L, double* d,
                       double* x, double c_d, double c_z, cudaStream_t stream)
    {
        static const int grid = resident_grid(chebyshev_kernel<FIRST, ZERO_X>);
        chebyshev_kernel<FIRST, ZERO_X><<<grid, kThreads, 0, stream>>>(b, q, L.diag, d, x, c_d, c_z,
                                                                       L.blk);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

    void upload_child_interp(const double* zeta)
    {
        double h[2][kN][kN];
        for (int child = 0; child < 2; ++child)
            for (int i = 0; i < kN; ++i) {
                const double x = (zeta[i] + (child ? 1.0 : -1.0)) / 2;
                for (int j = 0; j < kN; ++j) {
                    double v = 1.0;
                    for (int m = 0; m < kN; ++m)
                        if (m != j) v *= (x - zeta[m]) / (zeta[j] - zeta[m]);
                    h[child][i][j] = v;
                }
            }
        MARS_CELLWISE_CK(cudaMemcpyToSymbol(c_child_interp, h, sizeof(h)));
    }

    // The interior rows and columns of the one-element operator: all 216 interior unit
    // vectors in one operator call (216 copies of the element's metric), then inverted
    // once on the host (216^3 operations, setup only).
    void build_coarse_inverse()
    {
        const MultigridLevel& L = lv_.back();
        double *g = nullptr, *u = nullptr, *y = nullptr;
        MARS_CELLWISE_CK(cudaMalloc(&g, (long long)kInner3 * kMetric * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&u, (long long)kInner3 * kN3 * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&y, (long long)kInner3 * kN3 * sizeof(double)));
        const long long gn = (long long)kInner3 * kMetric, un = (long long)kInner3 * kN3;
        replicate_kernel<<<(unsigned)((gn + kThreads - 1) / kThreads), kThreads>>>(L.G, g, kMetric, kInner3);
        interior_units_kernel<<<(unsigned)((un + kThreads - 1) / kThreads), kThreads>>>(u);
        apply_(u, y, g, kInner3);
        std::vector<double> hy(un), a((size_t)kInner3 * kInner3), inv((size_t)kInner3 * kInner3, 0.0);
        MARS_CELLWISE_CK(cudaMemcpy(hy.data(), y, un * sizeof(double), cudaMemcpyDeviceToHost));
        for (int i = 0; i < kInner3; ++i)
            for (int j = 0; j < kInner3; ++j) a[(size_t)i * kInner3 + j] = hy[(size_t)j * kN3 + interior_local(i)];
        for (int i = 0; i < kInner3; ++i) inv[(size_t)i * kInner3 + i] = 1.0;
        for (int c = 0; c < kInner3; ++c) {   // Gauss-Jordan with partial pivoting
            int piv = c;
            for (int r = c + 1; r < kInner3; ++r)
                if (std::fabs(a[(size_t)r * kInner3 + c]) > std::fabs(a[(size_t)piv * kInner3 + c])) piv = r;
            for (int k = 0; k < kInner3; ++k) {
                std::swap(a[(size_t)c * kInner3 + k], a[(size_t)piv * kInner3 + k]);
                std::swap(inv[(size_t)c * kInner3 + k], inv[(size_t)piv * kInner3 + k]);
            }
            const double s = 1.0 / a[(size_t)c * kInner3 + c];
            for (int k = 0; k < kInner3; ++k) {
                a[(size_t)c * kInner3 + k] *= s;
                inv[(size_t)c * kInner3 + k] *= s;
            }
            for (int r = 0; r < kInner3; ++r) {
                const double f = a[(size_t)r * kInner3 + c];
                if (r == c || f == 0.0) continue;
                for (int k = 0; k < kInner3; ++k) {
                    a[(size_t)r * kInner3 + k] -= f * a[(size_t)c * kInner3 + k];
                    inv[(size_t)r * kInner3 + k] -= f * inv[(size_t)c * kInner3 + k];
                }
            }
        }
        MARS_CELLWISE_CK(cudaMalloc(&ainv_, inv.size() * sizeof(double)));
        MARS_CELLWISE_CK(cudaMemcpy(ainv_, inv.data(), inv.size() * sizeof(double), cudaMemcpyHostToDevice));
        for (double* p : {g, u, y}) cudaFree(p);
    }

    // Largest eigenvalue of P_J A on level l: 30 power steps from the hash vector,
    // estimate ||w||_w / ||v||_w with v normalized each step.
    double power_iteration(size_t l)
    {
        const MultigridLevel& L = lv_[l];
        const long long n = L.blk.values();
        const unsigned grid = (unsigned)((n + kThreads - 1) / kThreads);
        double* v = nullptr;
        MARS_CELLWISE_CK(cudaMalloc(&v, n * sizeof(double)));
        double* w = d_[l];
        hash_kernel<<<grid, kThreads>>>(w, n);
        precondition(w, v, L.diag, L.blk);
        Reduction red;
        double lam = 0.0;
        for (int it = 0; it < 30; ++it) {
            apply_(v, q_[l], L.G, L.blk.elements());
            precondition(q_[l], w, L.diag, L.blk);
            wdot(w, w, L.blk, red, kRR);
            const double nw = read_norm(red, 0);
            wdot(v, v, L.blk, red, kRR);
            lam = nw / read_norm(red, 0);
            scale_kernel<<<grid, kThreads>>>(v, w, 1.0 / nw, n);
        }
        cudaFree(v);
        return lam;
    }

    Apply apply_;
    std::vector<MultigridLevel> lv_;
    int degree_;
    double range_;
    std::vector<double*> b_, x_, q_, d_;   // per level: right-hand side, iterate, A x, Chebyshev d
    std::vector<double> lmax_;
    double* ainv_ = nullptr;
};

}  // namespace cellwise
}  // namespace mars
