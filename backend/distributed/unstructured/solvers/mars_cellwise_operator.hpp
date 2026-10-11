#pragma once
// The Knaus CVFEM operator at any degree P on element-local vectors, for the coarse
// levels of the cell-wise p-multigrid (the p = 7 level uses the MARSIR kernel). The 1D
// matrices travel as a kernel argument rather than through the constant bank of
// fem/mars_cvfem_ho_matfree.hpp, which holds one degree at a time, while a V-cycle
// alternates degrees. The metric layout is the MARSIR kernel's: per element
// [dir][face l][g0 g1 g2][s][r].
//
// The apply follows the layout of ho_cvfem_apply_kernel: one thread per tangential slot
// (s, r) of an element, the element's values staged in shared memory, and per direction
// the B/D sweeps along the thread's own column, the tangential derivatives and the flux
// through shared memory, and the W smoothing back onto the column.

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"

#include <cuda_runtime.h>

namespace mars {
namespace cellwise {

template <int P>
struct CvfemMatrices {
    double Btil[P][P + 1], Dtil[P][P + 1], D[P + 1][P + 1], W[P + 1][P + 1];
    double xi[P], zeta[P + 1];
};

template <int P>
CvfemMatrices<P> cvfem_matrices()
{
    const fem::HoCvfemOperators op = fem::buildHoCvfemOperators(P);
    constexpr int n = P + 1;
    CvfemMatrices<P> M;
    for (int l = 0; l < P; ++l)
        for (int a = 0; a < n; ++a) {
            M.Btil[l][a] = op.Btil[l * n + a];
            M.Dtil[l][a] = op.Dtil[l * n + a];
        }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            M.D[i][j] = op.D[i * n + j];
            M.W[i][j] = op.W[i * n + j];
        }
    for (int l = 0; l < P; ++l) M.xi[l] = op.xi[l];
    for (int i = 0; i < n; ++i) M.zeta[i] = op.zeta[i];
    return M;
}

template <int P>
constexpr long long metric_size = 3LL * P * 3 * (P + 1) * (P + 1);   // doubles per element

// Node a of the column through tangential slot (s, r) along direction d.
template <int P>
__device__ inline int column_node(int d, int a, int s, int r)
{
    return d == 0 ? node_at<P>(a, s, r) : (d == 1 ? node_at<P>(s, a, r) : node_at<P>(s, r, a));
}

// One thread per metric point (e, dir, l, s, r). corners_of(e, c) gives the corners.
template <int P, class Corners>
__global__ void cvfem_metric_kernel(double* __restrict__ G, long long E, Corners corners_of, CvfemMatrices<P> M)
{
    constexpr int n = P + 1, nn = n * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * (3LL * P * nn)) return;
    const long long e = t / (3LL * P * nn);
    int q = (int)(t % (3LL * P * nn));
    const int r = q % n, s = (q / n) % n, l = (q / nn) % P, dir = q / (nn * P);
    double c[8][3], g[3];
    corners_of(e, c);
    fem::ho_cvfem_metric_at(c, dir, M.xi[l], M.zeta[s], M.zeta[r], g);
    double* out = G + e * metric_size<P> + (dir * P + l) * 3 * nn + s * n + r;
    out[0] = g[0];
    out[nn] = g[1];
    out[2 * nn] = g[2];
}

// One thread per node: the element operator's diagonal in closed form (see the p = 7
// diagonal_kernel of mars_marsir_cellwise_solve.cu).
template <int P>
__global__ void cvfem_diagonal_kernel(const double* __restrict__ G, double* __restrict__ diag, long long E,
                                      CvfemMatrices<P> M)
{
    constexpr int n = P + 1, nn = n * n, n3 = nn * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * n3) return;
    const int node = (int)(t % n3);
    const int a = node / nn, b = (node / n) % n, c = node % n;
    const double* Ge = G + (t / n3) * metric_size<P>;
    double sum = 0.0;
    for (int dir = 0; dir < 3; ++dir) {
        const int A = dir == 0 ? a : (dir == 1 ? b : c), B = dir == 0 ? b : a, C = dir == 2 ? b : c;
        for (int side = 0; side < 2; ++side) {
            const int l = side == 0 ? A - 1 : A;
            if (l < 0 || l >= P) continue;
            const double* g0 = Ge + (dir * P + l) * 3 * nn;
            const double* g1 = g0 + nn;
            const double* g2 = g0 + 2 * nn;
            double term = M.W[B][B] * M.W[C][C] * g2[B * n + C] * M.Dtil[l][A];
            double t0 = 0.0, t1 = 0.0;
            for (int q = 0; q < n; ++q) {
                t0 += M.W[C][q] * g0[B * n + q] * M.D[q][C];
                t1 += M.W[B][q] * g1[q * n + C] * M.D[q][B];
            }
            term += M.W[B][B] * M.Btil[l][A] * t0 + M.W[C][C] * M.Btil[l][A] * t1;
            sum += side == 0 ? term : -term;
        }
    }
    diag[t] = sum;
}

// y = A u on every element (y is written, not accumulated).
template <int P>
__global__ void __launch_bounds__(kThreads)
cvfem_apply_kernel(const double* __restrict__ u, double* __restrict__ y, const double* __restrict__ G, long long E,
                   CvfemMatrices<P> M)
{
    constexpr int n = P + 1, nn = n * n, n3 = nn * n, per = kThreads / nn;
    static_assert(per >= 1, "one element's slots must fit a thread block");
    __shared__ double su[per][n3], sy[per][n3], sa[per][P][nn], sb[per][P][nn];
    const int k = threadIdx.x / nn, slot = threadIdx.x % nn, s = slot / n, r = slot % n;
    for (long long e0 = (long long)blockIdx.x * per; e0 < E; e0 += (long long)gridDim.x * per) {
        const long long e = e0 + k;
        const bool on = k < per && e < E;
        for (int i = threadIdx.x; i < per * n3; i += blockDim.x) {
            const long long g = e0 * n3 + i;
            (&su[0][0])[i] = g < E * n3 ? u[g] : 0.0;
            (&sy[0][0])[i] = 0.0;
        }
        __syncthreads();
        const double* Ge = G + (on ? e : 0) * metric_size<P>;
        for (int d = 0; d < 3; ++d) {
            double deriv[P];
            if (on) {
                double col[n];
                for (int a = 0; a < n; ++a) col[a] = su[k][column_node<P>(d, a, s, r)];
                for (int l = 0; l < P; ++l) {
                    double bi = 0.0, di = 0.0;
                    for (int a = 0; a < n; ++a) {
                        bi += M.Btil[l][a] * col[a];
                        di += M.Dtil[l][a] * col[a];
                    }
                    sa[k][l][slot] = bi;
                    deriv[l] = di;
                }
            }
            __syncthreads();
            if (on)
                for (int l = 0; l < P; ++l) {
                    double dt1 = 0.0, dt2 = 0.0;   // derivatives of the interpolant along s and r
                    for (int q = 0; q < n; ++q) {
                        dt1 += M.D[s][q] * sa[k][l][q * n + r];
                        dt2 += M.D[r][q] * sa[k][l][s * n + q];
                    }
                    const double* g = Ge + (d * P + l) * 3 * nn + slot;
                    sb[k][l][slot] = g[2 * nn] * deriv[l] + g[0] * dt2 + g[nn] * dt1;
                }
            __syncthreads();
            if (on)
                for (int l = 0; l < P; ++l) {
                    double f = 0.0;
                    for (int s2 = 0; s2 < n; ++s2) {
                        double row = 0.0;
                        for (int r2 = 0; r2 < n; ++r2) row += M.W[r][r2] * sb[k][l][s2 * n + r2];
                        f += M.W[s][s2] * row;
                    }
                    sy[k][column_node<P>(d, l, s, r)] -= f;
                    sy[k][column_node<P>(d, l + 1, s, r)] += f;
                }
        }
        __syncthreads();
        for (int i = threadIdx.x; i < per * n3; i += blockDim.x) {
            const long long g = e0 * n3 + i;
            if (g < E * n3) y[g] = (&sy[0][0])[i];
        }
        __syncthreads();
    }
}

template <int P>
void cvfem_apply(const double* d_u, double* d_y, const double* d_G, long long E, const CvfemMatrices<P>& M,
                 cudaStream_t stream)
{
    static const int grid = resident_grid(cvfem_apply_kernel<P>);
    if (E > 0) cvfem_apply_kernel<P><<<grid, kThreads, 0, stream>>>(d_u, d_y, d_G, E, M);
    MARS_CELLWISE_CK(cudaGetLastError());
}

}  // namespace cellwise
}  // namespace mars
