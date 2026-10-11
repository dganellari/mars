#pragma once
// Geometric multigrid preconditioner on element-local vectors.
//
// Credit: the design follows M. Wichrowski, "Coalesced Matrix-Free Geometric Multigrid on Persistent
// Cell-Wise Storage", arXiv:2607.03413 (2026): 2:1 coarsening at fixed degree on
// element-local storage, element-local tensor-product transfers, restriction applied
// to the raw unassembled residual, and DSS only inside the smoother. Added here: p = 7
// (the paper tests p <= 5), the nonsymmetric CVFEM operator, Chebyshev instead of
// damped Jacobi smoothing, one-sided (0,3) cycles, and the multi-GPU levels with the
// coarse level gathered onto rank 0.
//
// Levels: the structured block coarsened 2:1, all at p = 7, so every level runs the
// same operator kernel with its own metric. Transfers are element-local.
// Prolongation evaluates the coarse polynomial at the child's GLL nodes. Restriction
// is its exact transpose, applied to the raw unassembled residual: restricting an
// assembled residual would count shared nodes several times. Elements talk only
// through the DSS inside the smoother. Smoother: Chebyshev on the DSS-Jacobi
// preconditioned operator, its upper bound from a power iteration (the top of that
// spectrum is nearly real), `pre` steps before and `post` steps after the coarse
// correction. With pre = 0 the cycle restricts the right-hand side directly and needs
// no residual operator call. The one-element level is solved exactly.
//
// On several ranks every rank coarsens its own sub-block while all its dimensions
// stay even; the transfers need no communication because rank offsets stay even.
// The level where a rank can no longer halve is gathered onto rank 0, which runs the
// rest of the V-cycle on that whole (small) block as a one-rank multigrid and
// scatters the correction back. Only rank 0 sees data from every rank. The V-cycle is
// a fixed linear map with continuous output, as left-preconditioned BiCGStab
// requires. References: marsir-mlir/test/cellwise_multigrid_ref.py (one rank) and
// marsir-mlir/test/cellwise_distributed_ref.py (several).

#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"

#include <cuda_runtime.h>
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <utility>
#include <vector>

namespace mars {
namespace cellwise {

constexpr int kInner = kN - 2, kInner3 = kInner * kInner * kInner;   // interior nodes of one element
constexpr int kMetric = 3 * (kN - 1) * 3 * kNN;                      // metric doubles per element

// Coarse GLL basis at the GLL nodes of the left (0) and right (1) child: [child][fine][coarse].
__constant__ double c_child_interp[2][kN][kN];

// xf (+)= prolongation of the coarse field xc. One thread block per fine element; the
// three 1D contractions go through shared memory.
template <bool ADD>
__global__ void __launch_bounds__(kThreads)
prolong_kernel(const double* __restrict__ xc, double* __restrict__ xf, Block fine)
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
        xf[e * kN3 + l] = ADD ? xf[e * kN3 + l] + s : s;
    }
}

// bc = restriction of the fine residual bf - qf (bf alone without HAVE_Q), its
// Dirichlet copies zeroed: the transpose of prolongation, summed over the 8 children.
// One thread block per coarse element.
template <bool HAVE_Q>
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
            s0[l] = on_outer_boundary(nd, fine) ? 0.0 : (HAVE_Q ? bf[nd.t] - qf[nd.t] : bf[nd.t]);
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
// starting from zero), so z = P_J b and x = d. LAST: d is not needed again, not stored.
// SHELL: the pass over elements that read ghost copies.
template <bool FIRST, bool ZERO_X, bool LAST, bool SHELL>
__global__ void __launch_bounds__(kThreads)
chebyshev_kernel(const double* __restrict__ b, const double* __restrict__ q,
                 const double* __restrict__ diag, double* __restrict__ d, double* __restrict__ x,
                 double c_d, double c_z, Block blk, ElementSet set, GhostView ghost)
{
    const auto value_b = [b](long long i) { return b[i]; };
    const auto value_r = [b, q](long long i) { return b[i] - q[i]; };
    for_each_node<SHELL>(blk, set, [&](const Node& nd) {
        double z = 0.0;
        if (!on_outer_boundary(nd, blk)) {
            const double g = ZERO_X ? gather_value<SHELL>(value_b, ghost, nd, blk)
                                    : gather_value<SHELL>(value_r, ghost, nd, blk);
            z = g / diag[nd.t];
        }
        const double dn = FIRST ? c_z * z : c_d * d[nd.t] + c_z * z;
        if (!LAST) d[nd.t] = dn;
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

// Between the rank-major order of a gathered level (rank r's elements in its local
// order, ranks in turn) and the order of the whole block. TO_WHOLE: whole = gathered.
template <bool TO_WHOLE>
__global__ void permute_kernel(const double* __restrict__ in, double* __restrict__ out, Block whole,
                               Block part, int py, int pz)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= whole.values()) return;
    const long long e = t / kN3;
    const int l = (int)(t % kN3);
    const int gz = (int)(e % whole.nz), gy = (int)((e / whole.nz) % whole.ny);
    const int gx = (int)(e / ((long long)whole.nz * whole.ny));
    const long long rank = ((long long)(gx / part.nx) * py + gy / part.ny) * pz + gz / part.nz;
    const long long local = ((long long)(gx % part.nx) * part.ny + gy % part.ny) * part.nz + gz % part.nz;
    const long long r = (rank * part.elements() + local) * kN3 + l;
    if (TO_WHOLE) out[t] = in[r];
    else out[r] = in[t];
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

// A reproducible start vector for the power iteration, hashed from the GLOBAL value
// index (global element * 512 + node), so it is the same on any number of ranks and
// in the reference.
__global__ void hash_kernel(double* __restrict__ v, Block blk)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= blk.values()) return;
    const long long e = t / kN3;
    const int ez = (int)(e % blk.nz), ey = (int)((e / blk.nz) % blk.ny);
    const int ex = (int)(e / ((long long)blk.nz * blk.ny));
    const unsigned long long g =
        ((((unsigned long long)(ex + blk.ox) * blk.NY + ey + blk.oy) * blk.NZ + ez + blk.oz) * kN3) + t % kN3;
    v[t] = (double)((g * 2654435761ULL) & 0xffffffffULL) / 4294967296.0 - 0.5;
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

// `apply(u, y, G, E, stream)` computes y = A u on E elements with metric G.
// `make_level(blk)` allocates the metric and the UNASSEMBLED diagonal of a block
// (a rank's part or a whole block) and returns {G, diag}; the multigrid frees them.
// `finest` belongs to the caller unless `own_finest`.
template <typename Apply, typename MakeLevel>
class Multigrid {
public:
    Multigrid(Apply apply, MakeLevel make_level, MultigridLevel finest, const Decomposition& dec,
              const double* zeta_host, int pre, int post, cudaStream_t stream,
              bool own_finest = false, double smoothing_range = 15.0)
        : apply_(apply), make_(make_level), comm_(dec.comm), rank_(dec.rank), ranks_(dec.size),
          pre_(pre), post_(post), range_(smoothing_range), stream_(stream)
    {
        if (pre_ < 0 || post_ < 0 || pre_ + post_ == 0) {
            fprintf(stderr, "multigrid: need pre, post >= 0 and at least one smoothing step\n");
            std::abort();
        }
        for (int k = 0; k < 3; ++k) P_[k] = dec.P[k];
        for (int i = 0; i < kN; ++i) zeta_[i] = zeta_host[i];
        upload_child_interp();
        if (own_finest) {
            owned_.push_back(finest.G);
            owned_.push_back(finest.diag);
        }
        // Every local dimension halves at each level. One rank goes down to one element;
        // several ranks stop at the first level some rank cannot halve, which is the
        // level gathered onto rank 0.
        std::vector<Block> blocks{finest.blk};
        while (blocks.back().halvable() && blocks.back().elements() > 1) {
            blocks.push_back(blocks.back().coarse());
            if (ranks_ > 1 && !(blocks.back().halvable() && blocks.back().elements() > 1)) break;
        }
        if (ranks_ == 1 && blocks.back().elements() != 1) {
            fprintf(stderr, "multigrid: one rank needs ne = 2^k (coarsest level is %d x %d x %d)\n",
                    blocks.back().nx, blocks.back().ny, blocks.back().nz);
            std::abort();
        }
        for (size_t l = 0; l < blocks.size(); ++l) {
            const bool gathered_level = ranks_ > 1 && l + 1 == blocks.size();
            halos_.emplace_back(ranks_ > 1 && !gathered_level ? new Halo(dec, blocks[l]) : nullptr);
            if (l == 0) lv_.push_back(finest);
            else if (gathered_level) lv_.push_back({blocks[l], nullptr, nullptr});   // rank 0 makes its own
            else lv_.push_back(make_assembled(blocks[l], halos_[l].get()));
        }
        for (size_t l = 0; l < lv_.size(); ++l) {
            const long long n = lv_[l].blk.values();
            q_.push_back(nullptr);
            d_.push_back(nullptr);
            b_.push_back(nullptr);
            x_.push_back(nullptr);
            if (l + 1 < lv_.size()) {
                MARS_CELLWISE_CK(cudaMalloc(&q_[l], n * sizeof(double)));
                MARS_CELLWISE_CK(cudaMalloc(&d_[l], n * sizeof(double)));
            }
            if (l > 0) {
                MARS_CELLWISE_CK(cudaMalloc(&b_[l], n * sizeof(double)));
                MARS_CELLWISE_CK(cudaMalloc(&x_[l], n * sizeof(double)));
            }
        }
        if (ranks_ == 1) build_coarse_inverse();
        else setup_gather();
        for (size_t l = 0; l + 1 < lv_.size(); ++l) lmax_.push_back(power_iteration(l));
    }
    ~Multigrid()
    {
        for (auto* v : {&b_, &x_, &q_, &d_})
            for (double* p : *v) cudaFree(p);
        for (const double* p : owned_) cudaFree(const_cast<double*>(p));
        for (double* p : {ainv_, d_all_, root_b_, root_x_}) cudaFree(p);
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
            static const int grid = resident_grid(wdots_kernel<AZ, ZZ, BlockWeight>);
            wdots_kernel<AZ, ZZ, BlockWeight><<<grid, kThreads, 0, stream>>>(z, a, lv_[0].blk,
                                                                              BlockWeight{lv_[0].blk}, partial);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        }
        return 0;
    }

    void vcycle(size_t l, const double* b, double* x, cudaStream_t stream) const
    {
        if (l + 1 == lv_.size()) {
            if (ranks_ > 1) gather_solve(b, x, stream);
            else coarse_solve_kernel<<<1, kThreads, 0, stream>>>(ainv_, b, x);
            return;
        }
        const MultigridLevel& L = lv_[l];
        const unsigned fine = (unsigned)L.blk.elements(), coarse = (unsigned)lv_[l + 1].blk.elements();
        if (pre_ > 0) {
            smooth(l, b, x, true, pre_, stream);
            apply_(x, q_[l], L.G, L.blk.elements(), stream);
            restrict_kernel<true><<<coarse, kThreads, 0, stream>>>(b, q_[l], b_[l + 1], L.blk);
        } else {
            restrict_kernel<false><<<coarse, kThreads, 0, stream>>>(b, nullptr, b_[l + 1], L.blk);
        }
        vcycle(l + 1, b_[l + 1], x_[l + 1], stream);
        if (pre_ > 0)
            prolong_kernel<true><<<fine, kThreads, 0, stream>>>(x_[l + 1], x, L.blk);
        else
            prolong_kernel<false><<<fine, kThreads, 0, stream>>>(x_[l + 1], x, L.blk);
        if (post_ > 0) smooth(l, b, x, false, post_, stream);
    }

    // Largest eigenvalue estimates of P_J A per smoothed level, rank 0's levels after.
    std::vector<double> lambda_max() const
    {
        std::vector<double> all = lmax_;
        if (root_) {
            const std::vector<double> r = root_->lambda_max();
            all.insert(all.end(), r.begin(), r.end());
        }
        return all;
    }

private:
    MultigridLevel make_assembled(const Block& b, Halo* halo)
    {
        const std::pair<double*, double*> made = make_(b);
        double* diag = nullptr;
        MARS_CELLWISE_CK(cudaMalloc(&diag, b.values() * sizeof(double)));
        dss(made.second, diag, b, halo, stream_);
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        cudaFree(made.second);
        owned_.push_back(made.first);
        owned_.push_back(diag);
        return {b, made.first, diag};
    }

    void smooth(size_t l, const double* b, double* x, bool zero_x, int steps,
                cudaStream_t stream) const
    {
        const double lmax = 1.1 * lmax_[l], lmin = lmax / range_;
        const double theta = 0.5 * (lmax + lmin), delta = 0.5 * (lmax - lmin);
        const double sigma = theta / delta;
        double rho = 1.0 / sigma;
        const bool last = steps == 1;
        if (zero_x) {
            if (last) launch<true, true, true>(l, b, x, 0.0, 1.0 / theta, stream);
            else launch<true, true, false>(l, b, x, 0.0, 1.0 / theta, stream);
        } else {
            apply_(x, q_[l], lv_[l].G, lv_[l].blk.elements(), stream);
            if (last) launch<true, false, true>(l, b, x, 0.0, 1.0 / theta, stream);
            else launch<true, false, false>(l, b, x, 0.0, 1.0 / theta, stream);
        }
        for (int k = 1; k < steps; ++k) {
            const double rho_new = 1.0 / (2.0 * sigma - rho);
            const double c_d = rho_new * rho, c_z = 2.0 * rho_new / delta;
            apply_(x, q_[l], lv_[l].G, lv_[l].blk.elements(), stream);
            if (k + 1 == steps) launch<false, false, true>(l, b, x, c_d, c_z, stream);
            else launch<false, false, false>(l, b, x, c_d, c_z, stream);
            rho = rho_new;
        }
    }

    template <bool FIRST, bool ZERO_X, bool LAST>
    void launch(size_t l, const double* b, double* x, double c_d, double c_z, cudaStream_t stream) const
    {
        const MultigridLevel& L = lv_[l];
        double* d = d_[l];
        const double* q = q_[l];
        Halo* halo = halos_[l].get();
        run_gather_pass(halo, b, ZERO_X ? nullptr : q, stream, [&](const ElementSet& set, bool shell, int) {
            if (shell) {
                static const int grid = resident_grid(chebyshev_kernel<FIRST, ZERO_X, LAST, true>);
                chebyshev_kernel<FIRST, ZERO_X, LAST, true><<<grid, kThreads, 0, stream>>>(
                    b, q, L.diag, d, x, c_d, c_z, L.blk, set, halo->view());
                MARS_CELLWISE_CK(cudaGetLastError());
                return grid;
            }
            static const int grid = resident_grid(chebyshev_kernel<FIRST, ZERO_X, LAST, false>);
            chebyshev_kernel<FIRST, ZERO_X, LAST, false><<<grid, kThreads, 0, stream>>>(
                b, q, L.diag, d, x, c_d, c_z, L.blk, set, GhostView{});
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        });
    }

    void upload_child_interp()
    {
        double h[2][kN][kN];
        for (int child = 0; child < 2; ++child)
            for (int i = 0; i < kN; ++i) {
                const double x = (zeta_[i] + (child ? 1.0 : -1.0)) / 2;
                for (int j = 0; j < kN; ++j) {
                    double v = 1.0;
                    for (int m = 0; m < kN; ++m)
                        if (m != j) v *= (x - zeta_[m]) / (zeta_[j] - zeta_[m]);
                    h[child][i][j] = v;
                }
            }
        MARS_CELLWISE_CK(cudaMemcpyToSymbol(c_child_interp, h, sizeof(h)));
        // A pageable upload may still be in flight when the call returns, and the
        // solver's stream does not wait for the default stream.
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
    }

    // The interior rows and columns of the one-element operator: all 216 interior unit
    // vectors in one operator call (216 copies of the element's metric), then inverted
    // once on the host (216^3 operations, setup only).
    void build_coarse_inverse()
    {
        const MultigridLevel& L = lv_.back();
        double *g = nullptr, *u = nullptr, *y = nullptr;
        const long long gn = (long long)kInner3 * kMetric, un = (long long)kInner3 * kN3;
        MARS_CELLWISE_CK(cudaMalloc(&g, gn * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&u, un * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&y, un * sizeof(double)));
        replicate_kernel<<<(unsigned)((gn + kThreads - 1) / kThreads), kThreads, 0, stream_>>>(L.G, g, kMetric,
                                                                                            kInner3);
        interior_units_kernel<<<(unsigned)((un + kThreads - 1) / kThreads), kThreads, 0, stream_>>>(u);
        apply_(u, y, g, kInner3, stream_);
        std::vector<double> hy(un), a((size_t)kInner3 * kInner3), inv((size_t)kInner3 * kInner3, 0.0);
        MARS_CELLWISE_CK(cudaMemcpyAsync(hy.data(), y, un * sizeof(double), cudaMemcpyDeviceToHost, stream_));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
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
        MARS_CELLWISE_CK(cudaDeviceSynchronize());   // see upload_child_interp
        for (double* p : {g, u, y}) cudaFree(p);
    }

    // The gathered level: every rank's part travels to rank 0, which owns a one-rank
    // multigrid of the whole block at this level (made, assembled and coarsened there).
    // Device buffers go straight to MPI, as for the ghost exchange.
    void setup_gather()
    {
        const Block part = lv_.back().blk;
        if (rank_ != 0) return;
        whole_ = Block(part.NX, part.NY, part.NZ);
        const long long all = whole_.values();
        MARS_CELLWISE_CK(cudaMalloc(&d_all_, all * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&root_b_, all * sizeof(double)));
        MARS_CELLWISE_CK(cudaMalloc(&root_x_, all * sizeof(double)));
        const Decomposition one(MPI_COMM_SELF, whole_.NX, whole_.NY, whole_.NZ);
        const std::pair<double*, double*> made = make_(whole_);
        double* diag = nullptr;
        MARS_CELLWISE_CK(cudaMalloc(&diag, all * sizeof(double)));
        dss(made.second, diag, whole_, nullptr, stream_);
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        cudaFree(made.second);
        root_.reset(new Multigrid(apply_, make_, {whole_, made.first, diag}, one, zeta_, pre_, post_, stream_,
                                  true, range_));
    }

    void gather_solve(const double* b, double* x, cudaStream_t stream) const
    {
        const Block part = lv_.back().blk;
        const int n = (int)part.values();
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));   // b complete before MPI reads it
        MARS_CELLWISE_MPI(MPI_Gather(b, n, MPI_DOUBLE, d_all_, n, MPI_DOUBLE, 0, comm_));
        if (rank_ == 0) {
            const long long all = whole_.values();
            const unsigned grid = (unsigned)((all + kThreads - 1) / kThreads);
            permute_kernel<true><<<grid, kThreads, 0, stream>>>(d_all_, root_b_, whole_, part, P_[1], P_[2]);
            root_->vcycle(0, root_b_, root_x_, stream);
            permute_kernel<false><<<grid, kThreads, 0, stream>>>(root_x_, d_all_, whole_, part, P_[1], P_[2]);
            MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        }
        MARS_CELLWISE_MPI(MPI_Scatter(d_all_, n, MPI_DOUBLE, x, n, MPI_DOUBLE, 0, comm_));
    }

    // Largest eigenvalue of P_J A on level l: 30 power steps from the hash vector,
    // estimate ||w||_w / ||v||_w with v normalized each step. Fewer steps underestimate
    // it, and Chebyshev then amplifies the modes above the bound (5 -> 29 iterations
    // with 10 steps at 8^3, deform 0.1).
    double power_iteration(size_t l)
    {
        const MultigridLevel& L = lv_[l];
        const long long n = L.blk.values();
        const unsigned grid = (unsigned)((n + kThreads - 1) / kThreads);
        Halo* halo = halos_[l].get();
        double* v = nullptr;
        MARS_CELLWISE_CK(cudaMalloc(&v, n * sizeof(double)));
        double* w = d_[l];
        hash_kernel<<<grid, kThreads, 0, stream_>>>(w, L.blk);
        precondition(w, v, L.diag, L.blk, halo, stream_);
        Reduction red(comm_);
        double lam = 0.0;
        for (int it = 0; it < 30; ++it) {
            apply_(v, q_[l], L.G, L.blk.elements(), stream_);
            precondition(q_[l], w, L.diag, L.blk, halo, stream_);
            wdot(w, w, L.blk, red, kRR, stream_);
            const double nw = read_norm(red, stream_);
            wdot(v, v, L.blk, red, kRR, stream_);
            lam = nw / read_norm(red, stream_);
            scale_kernel<<<grid, kThreads, 0, stream_>>>(v, w, 1.0 / nw, n);
        }
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        cudaFree(v);
        return lam;
    }

    Apply apply_;
    MakeLevel make_;
    MPI_Comm comm_;
    int rank_, ranks_;
    int P_[3];
    int pre_, post_;
    double range_;
    cudaStream_t stream_;
    double zeta_[kN];
    std::vector<MultigridLevel> lv_;
    std::vector<std::unique_ptr<Halo>> halos_;   // null on one rank and on the gathered level
    std::vector<double*> b_, x_, q_, d_;         // per level: right-hand side, iterate, A x, Chebyshev d
    std::vector<double> lmax_;
    std::vector<const double*> owned_;           // metrics and diagonals this object made
    double* ainv_ = nullptr;
    // The gathered level: rank 0's whole block, its rank-major copy, and its solver.
    Block whole_;
    double *d_all_ = nullptr, *root_b_ = nullptr, *root_x_ = nullptr;
    std::unique_ptr<Multigrid> root_;
};

}  // namespace cellwise
}  // namespace mars
