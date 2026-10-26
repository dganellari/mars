// The CEED benchmark problems BP1 (mass) and BP3 (Laplacian) on element-local
// (cell-wise) vectors, with the operators MARSIR generates from marsir-compiler/specs/
// bp1.op and bp3.op (PTX from marsir-mlir/test/bp_gpu_pipeline.py), the unstructured DSS
// (solvers/mars_cellwise_topology.hpp), and the communication-free flexible CG of
// Wichrowski's "Coalesced Matrix-Free Finite Elements in Cell-Wise Storage"
// (solvers/mars_cellwise_fcg.hpp) with DSS + Jacobi.
//
// Mesh: the unit cube as ne^3 trilinear hexahedra deformed by x += a sin(pi x) sin(pi y)
// sin(pi z) (1,1,1), GLL nodes of degree p, q = p + 2 Gauss points. Geometry per Gauss
// point on the GPU: w det J (BP1), the symmetric w det J J^-1 J^-T (BP3).
//
// Checks: BP1 applied to 1 sums to the volume (1); BP3 applied to 1 is 0; the operator on
// the first elements equals direct quadrature sums on the host. Solve: b = A u with u =
// sin(pi x) sin(pi y) sin(pi z) at the nodes (zero on the boundary for BP3), so FCG
// must return u; its history must match marsir-mlir/test/cellwise_bp_ref.py
// (GPU-comparable run) for the same --spec, --p, --ne, --deform.
//
// Timing (one GPU): the operator apply, the DSS, apply + DSS (the assembled operator, as
// the paper's Table 7), an FCG iteration. Rates are per unique DoF, (ne p + 1)^3; the
// bandwidth counts what one pass must move at least: u and y, and the geometry.
//
// --dss-ptx <dss_p{p}_sm90.ptx>: the DSS MARSIR generates from mir.dss
// (marsir-mlir/test/dss_gpu_pipeline.py p; same bits as the hand gather, so the same
// history) instead of the hand dss<P>.
//
// Run: mars_marsir_bp --ptx-dir <marsir-mlir/generated> --spec bp3 --p 3 [--ne 32]
//                     [--deform 0.05] [--tol 1e-10] [--maxit 2000] [--reps 20] [--dss-ptx file]

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_cellwise.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_galerkin.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_fcg.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <mpi.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace mars;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kMaxQ = 9;   // Gauss points per axis at p = 7
__constant__ double c_gauss[kMaxQ], c_weight[kMaxQ], c_gll[kMaxQ];

unsigned blocks_for(long long n) { return (unsigned)((n + cellwise::kThreads - 1) / cellwise::kThreads); }

__device__ void lattice_vertex(int i, int j, int k, int ne, double deform, double x[3])
{
    x[0] = (double)i / ne;
    x[1] = (double)j / ne;
    x[2] = (double)k / ne;
    const double bump = deform * sin(kPi * x[0]) * sin(kPi * x[1]) * sin(kPi * x[2]);
    for (int d = 0; d < 3; ++d) x[d] += bump;
}

// Element e = (ex ne + ey) ne + ez, corners in VTK order (corner_bits).
__device__ void element_corners(long long e, int ne, double deform, double c[8][3])
{
    const int ez = (int)(e % ne), ey = (int)((e / ne) % ne), ex = (int)(e / ((long long)ne * ne));
    for (int k = 0; k < 8; ++k) {
        int bx, by, bz;
        cellwise::corner_bits(k, bx, by, bz);
        lattice_vertex(ex + bx, ey + by, ez + bz, ne, deform, c[k]);
    }
}

// The trilinear map and its Jacobian J[a][b] = dx_a / dxi_b at reference point r.
__device__ void trilinear(const double c[8][3], const double r[3], double x[3], double J[3][3])
{
    for (int a = 0; a < 3; ++a) {
        x[a] = 0.0;
        for (int b = 0; b < 3; ++b) J[a][b] = 0.0;
    }
    for (int k = 0; k < 8; ++k) {
        int bits[3];
        cellwise::corner_bits(k, bits[0], bits[1], bits[2]);
        double f[3], s[3];
        for (int d = 0; d < 3; ++d) {
            s[d] = 2.0 * bits[d] - 1.0;
            f[d] = 1.0 + s[d] * r[d];
        }
        const double w = 0.125 * f[0] * f[1] * f[2];
        const double dw[3] = {0.125 * s[0] * f[1] * f[2], 0.125 * f[0] * s[1] * f[2], 0.125 * f[0] * f[1] * s[2]};
        for (int a = 0; a < 3; ++a) {
            x[a] += w * c[k][a];
            for (int b = 0; b < 3; ++b) J[a][b] += dw[b] * c[k][a];
        }
    }
}

// One thread per (element, Gauss point): G[e][k][point], K = 1 (BP1: w det J) or 6 (BP3:
// w det J J^-1 J^-T, entries 00 01 02 11 12 22).
__global__ void geometry_kernel(double* __restrict__ G, long long E, int Q, int K, int ne, double deform)
{
    const long long q3 = (long long)Q * Q * Q;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * q3) return;
    const long long e = t / q3;
    const int pt = (int)(t % q3), a = pt / (Q * Q), b = (pt / Q) % Q, c = pt % Q;
    double corners[8][3], x[3], J[3][3];
    element_corners(e, ne, deform, corners);
    const double r[3] = {c_gauss[a], c_gauss[b], c_gauss[c]};
    trilinear(corners, r, x, J);
    const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                       J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                       J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    const double wdet = det * c_weight[a] * c_weight[b] * c_weight[c];
    double* g = G + e * K * q3 + pt;
    if (K == 1) {
        g[0] = wdet;
        return;
    }
    double Ji[3][3];   // J^-1 = adj(J) / det
    Ji[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
    Ji[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det;
    Ji[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
    Ji[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / det;
    Ji[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
    Ji[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det;
    Ji[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
    Ji[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det;
    Ji[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
    const int sym[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
    for (int k = 0; k < 6; ++k) {
        const int d = sym[k][0], f = sym[k][1];
        g[k * q3] = wdet * (Ji[d][0] * Ji[f][0] + Ji[d][1] * Ji[f][1] + Ji[d][2] * Ji[f][2]);
    }
}

// u = sin(pi x) sin(pi y) sin(pi z) at every copy; zero on the cube's boundary copies
// when `dirichlet` (sin vanishes there up to rounding; the copies are set exactly).
template <int P>
__global__ void exact_kernel(double* __restrict__ u, long long E, int ne, double deform, bool dirichlet)
{
    constexpr int n = P + 1, n3 = n * n * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * n3) return;
    const long long e = t / n3;
    const int node = (int)(t % n3), a = node / (n * n), b = (node / n) % n, c = node % n;
    const int ez = (int)(e % ne), ey = (int)((e / ne) % ne), ex = (int)(e / ((long long)ne * ne));
    const int I = ex * P + a, Jx = ey * P + b, K = ez * P + c, last = ne * P;
    if (dirichlet && (I == 0 || Jx == 0 || K == 0 || I == last || Jx == last || K == last)) {
        u[t] = 0.0;
        return;
    }
    double corners[8][3], x[3], J[3][3];
    element_corners(e, ne, deform, corners);
    const double r[3] = {c_gll[a], c_gll[b], c_gll[c]};
    trilinear(corners, r, x, J);
    u[t] = sin(kPi * x[0]) * sin(kPi * x[1]) * sin(kPi * x[2]);
}

// The element operator's diagonal, closed form (as cellwise_bp_ref.diagonal): BP1
// sum_q g0 phi^2, BP3 sum_q sum_rs M_rs d_r phi d_s phi, phi_ijk = B_ai B_bj B_ck, DB =
// D B the basis' derivative at the Gauss points.
template <int P>
__global__ void diagonal_kernel(const double* __restrict__ G, const double* __restrict__ B,
                                const double* __restrict__ DB, double* __restrict__ diag, long long E, int K)
{
    constexpr int n = P + 1, n3 = n * n * n, Q = P + 2, q3 = Q * Q * Q;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * n3) return;
    const long long e = t / n3;
    const int node = (int)(t % n3), i = node / (n * n), j = (node / n) % n, k = node % n;
    const double* g = G + e * K * q3;
    double s = 0.0;
    for (int a = 0; a < Q; ++a)
        for (int b = 0; b < Q; ++b)
            for (int c = 0; c < Q; ++c) {
                const int pt = (a * Q + b) * Q + c;
                const double bi = B[a * n + i], bj = B[b * n + j], bk = B[c * n + k];
                if (K == 1) {
                    const double phi = bi * bj * bk;
                    s += g[pt] * phi * phi;
                    continue;
                }
                const double d0 = DB[a * n + i] * bj * bk, d1 = bi * DB[b * n + j] * bk, d2 = bi * bj * DB[c * n + k];
                s += g[pt] * d0 * d0 + g[3 * q3 + pt] * d1 * d1 + g[5 * q3 + pt] * d2 * d2 +
                     2.0 * (g[q3 + pt] * d0 * d1 + g[2 * q3 + pt] * d0 * d2 + g[4 * q3 + pt] * d1 * d2);
            }
    diag[t] = s;
}

// inv_diag = mask / diag, mask = 0 on the cube's boundary copies when `dirichlet`.
template <int P>
__global__ void inverse_kernel(const double* __restrict__ diag, double* __restrict__ inv, long long E, int ne,
                               bool dirichlet)
{
    constexpr int n = P + 1, n3 = n * n * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * n3) return;
    const long long e = t / n3;
    const int node = (int)(t % n3), a = node / (n * n), b = (node / n) % n, c = node % n;
    const int ez = (int)(e % ne), ey = (int)((e / ne) % ne), ex = (int)(e / ((long long)ne * ne));
    const int I = ex * P + a, Jx = ey * P + b, K = ez * P + c, last = ne * P;
    const bool boundary = I == 0 || Jx == 0 || K == 0 || I == last || Jx == last || K == last;
    inv[t] = dirichlet && boundary ? 0.0 : 1.0 / diag[t];
}

__global__ void mesh_kernel(long long E, int ne, unsigned long long* const* key, int* const* lid)
{
    const long long e = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= E) return;
    const int ez = (int)(e % ne), ey = (int)((e / ne) % ne), ex = (int)(e / ((long long)ne * ne));
    for (int c = 0; c < 8; ++c) {
        int bx, by, bz;
        cellwise::corner_bits(c, bx, by, bz);
        const long long v = ((long long)(ex + bx) * (ne + 1) + ey + by) * (ne + 1) + ez + bz;
        key[c][e] = (unsigned long long)v;
        lid[c][e] = (int)v;
    }
}

__global__ void difference_kernel(double* __restrict__ out, const double* __restrict__ a,
                                  const double* __restrict__ b, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] = fabs(a[t] - b[t]);
}

#define BP_CK(call) MARS_CELLWISE_CK(call)

struct Options {
    std::string ptx_dir, spec = "bp3", dss_ptx;
    int p = 3, ne = 32, maxit = 2000, reps = 20;
    double deform = 0.05, tol = 1e-10;
};

template <int P>
int run(const Options& o)
{
    constexpr int n = P + 1, n3 = n * n * n, Q = P + 2, q3 = Q * Q * Q;
    const bool bp3 = o.spec == "bp3";
    const int K = bp3 ? 6 : 1;
    const long long E = (long long)o.ne * o.ne * o.ne, N = E * n3;
    const long long unique = (long long)(o.ne * P + 1) * (o.ne * P + 1) * (o.ne * P + 1);

    // The 1D matrices: B (GLL basis at the Gauss points), D (Gauss-point derivative),
    // their transposes, and DB = D B for the diagonal.
    std::vector<double> zeta, wgll, Dgll, gq, wq;
    fem::gllBasis(P, zeta, wgll, Dgll);
    fem::gaussLegendre(Q, gq, wq);
    std::vector<double> B(Q * n), Bt(n * Q), D(Q * Q), Dt(Q * Q), DB(Q * n, 0.0), val, der;
    for (int a = 0; a < Q; ++a) {
        fem::lagrangeValDeriv(zeta, gq[a], val, der);
        for (int i = 0; i < n; ++i) B[a * n + i] = Bt[i * Q + a] = val[i];
        fem::lagrangeValDeriv(gq, gq[a], val, der);
        for (int b = 0; b < Q; ++b) D[a * Q + b] = Dt[b * Q + a] = der[b];
    }
    for (int a = 0; a < Q; ++a)
        for (int i = 0; i < n; ++i)
            for (int m = 0; m < Q; ++m) DB[a * n + i] += D[a * Q + m] * B[m * n + i];
    BP_CK(cudaMemcpyToSymbol(c_gauss, gq.data(), Q * sizeof(double)));
    BP_CK(cudaMemcpyToSymbol(c_weight, wq.data(), Q * sizeof(double)));
    BP_CK(cudaMemcpyToSymbol(c_gll, zeta.data(), n * sizeof(double)));
    thrust::device_vector<double> dB(B), dBt(Bt), dD(D), dDt(Dt), dDB(DB);
    auto raw = [](thrust::device_vector<double>& v) { return thrust::raw_pointer_cast(v.data()); };

    cudaStream_t stream;
    BP_CK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    thrust::device_vector<double> G(E * K * q3), u(N), y(N), b(N), x(N), diag(N), inv(N), s(N), r(N), q(N), z(N),
        pdir(N), ones(N, 1.0);
    geometry_kernel<<<blocks_for(E * q3), cellwise::kThreads, 0, stream>>>(raw(G), E, Q, K, o.ne, o.deform);
    BP_CK(cudaGetLastError());

    const std::string kernel = o.spec + "_apply_p" + std::to_string(P);
    const std::string path = o.ptx_dir + "/" + o.spec + "_p" + std::to_string(P) + "_sm90.ptx";
    marsir::GalerkinPtx op(path, kernel, P, Q, bp3, K);
    auto apply = [&](const double* in, double* out, cudaStream_t st) {
        op.apply(in, raw(dB), raw(dBt), raw(dD), raw(dDt), raw(G), out, E, st);
    };

    // The unstructured tables of the cube (one rank: every element local).
    thrust::device_vector<unsigned long long> key[8];
    thrust::device_vector<int> lid[8];
    std::vector<unsigned long long*> kp(8);
    std::vector<int*> lp(8);
    for (int c = 0; c < 8; ++c) {
        key[c].resize(E);
        lid[c].resize(E);
        kp[c] = thrust::raw_pointer_cast(key[c].data());
        lp[c] = thrust::raw_pointer_cast(lid[c].data());
    }
    thrust::device_vector<unsigned long long*> d_kp(kp.begin(), kp.end());
    thrust::device_vector<int*> d_lp(lp.begin(), lp.end());
    mesh_kernel<<<blocks_for(E), cellwise::kThreads, 0, stream>>>(E, o.ne, thrust::raw_pointer_cast(d_kp.data()),
                                                                  thrust::raw_pointer_cast(d_lp.data()));
    BP_CK(cudaStreamSynchronize(stream));
    const std::vector<const unsigned long long*> ckp(kp.begin(), kp.end());
    const std::vector<const int*> clp(lp.begin(), lp.end());
    thrust::device_vector<const unsigned long long*> d_ckp(ckp.begin(), ckp.end());
    thrust::device_vector<const int*> d_clp(clp.begin(), clp.end());
    const cellwise::UnstructuredTopology topo =
        cellwise::build_topology(thrust::raw_pointer_cast(d_ckp.data()), thrust::raw_pointer_cast(d_clp.data()), E, E);
    std::unique_ptr<marsir::DssPtx> dss_gen(o.dss_ptx.empty() ? nullptr : new marsir::DssPtx(o.dss_ptx, P));
    auto dss = [&topo, gen = dss_gen.get(), E](const double* in, double* out, cudaStream_t st) {
        if (gen) gen->run(in, nullptr, out, topo, cellwise::ElementList{nullptr, E}, st);
        else cellwise::dss<P>(in, out, topo, nullptr, st);
    };

    // Checks of the operator.
    apply(raw(ones), raw(y), stream);
    BP_CK(cudaStreamSynchronize(stream));
    const double sum1 = thrust::reduce(y.begin(), y.end(), 0.0);
    exact_kernel<P><<<blocks_for(N), cellwise::kThreads, 0, stream>>>(raw(u), E, o.ne, o.deform, bp3);
    BP_CK(cudaGetLastError());
    apply(raw(u), raw(b), stream);   // b = A u
    BP_CK(cudaStreamSynchronize(stream));
    double max1 = 0.0, maxu = 0.0;
    {
        std::vector<double> hy(N), hb(N);
        thrust::copy(y.begin(), y.end(), hy.begin());
        thrust::copy(b.begin(), b.end(), hb.begin());
        for (long long i = 0; i < N; ++i) {
            max1 = std::fmax(max1, std::fabs(hy[i]));
            maxu = std::fmax(maxu, std::fabs(hb[i]));
        }
    }
    // The first elements against direct quadrature sums on the host.
    const long long Ec = std::min<long long>(E, 4);
    double oracle_err = 0.0, oracle_scale = 0.0;
    {
        std::vector<double> hu(Ec * n3), hb(Ec * n3), hG(Ec * K * q3);
        thrust::copy(u.begin(), u.begin() + Ec * n3, hu.begin());
        thrust::copy(b.begin(), b.begin() + Ec * n3, hb.begin());
        thrust::copy(G.begin(), G.begin() + Ec * K * q3, hG.begin());
        const int sym[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        for (long long e = 0; e < Ec; ++e) {
            std::vector<double> uq(q3, 0.0), d[3], w(q3, 0.0);
            for (int pt = 0; pt < q3; ++pt) {
                const int a = pt / (Q * Q), bb = (pt / Q) % Q, c = pt % Q;
                for (int node = 0; node < n3; ++node)
                    uq[pt] += B[a * n + node / (n * n)] * B[bb * n + (node / n) % n] * B[c * n + node % n] *
                              hu[e * n3 + node];
            }
            const double* g = &hG[e * K * q3];
            if (!bp3) {
                for (int pt = 0; pt < q3; ++pt) w[pt] = g[pt] * uq[pt];
            } else {
                for (int rr = 0; rr < 3; ++rr) d[rr].assign(q3, 0.0);
                for (int pt = 0; pt < q3; ++pt) {
                    const int a = pt / (Q * Q), bb = (pt / Q) % Q, c = pt % Q;
                    for (int m = 0; m < Q; ++m) {
                        d[0][pt] += D[a * Q + m] * uq[(m * Q + bb) * Q + c];
                        d[1][pt] += D[bb * Q + m] * uq[(a * Q + m) * Q + c];
                        d[2][pt] += D[c * Q + m] * uq[(a * Q + bb) * Q + m];
                    }
                }
                std::vector<double> dv[3];
                for (int rr = 0; rr < 3; ++rr) {
                    dv[rr].assign(q3, 0.0);
                    for (int pt = 0; pt < q3; ++pt)
                        for (int ss = 0; ss < 3; ++ss) dv[rr][pt] += g[sym[rr][ss] * q3 + pt] * d[ss][pt];
                }
                for (int pt = 0; pt < q3; ++pt) {
                    const int a = pt / (Q * Q), bb = (pt / Q) % Q, c = pt % Q;
                    for (int m = 0; m < Q; ++m)
                        w[pt] += D[m * Q + a] * dv[0][(m * Q + bb) * Q + c] + D[m * Q + bb] * dv[1][(a * Q + m) * Q + c] +
                                 D[m * Q + c] * dv[2][(a * Q + bb) * Q + m];
                }
            }
            for (int node = 0; node < n3; ++node) {
                double sum = 0.0;
                for (int pt = 0; pt < q3; ++pt) {
                    const int a = pt / (Q * Q), bb = (pt / Q) % Q, c = pt % Q;
                    sum += B[a * n + node / (n * n)] * B[bb * n + (node / n) % n] * B[c * n + node % n] * w[pt];
                }
                oracle_err = std::fmax(oracle_err, std::fabs(sum - hb[e * n3 + node]));
                oracle_scale = std::fmax(oracle_scale, std::fabs(sum));
            }
        }
    }

    // The preconditioner: the assembled diagonal on every copy, inverted, Dirichlet zeroed.
    diagonal_kernel<P><<<blocks_for(N), cellwise::kThreads, 0, stream>>>(raw(G), raw(dB), raw(dDB), raw(s), E, K);
    BP_CK(cudaGetLastError());
    dss(raw(s), raw(diag), stream);
    inverse_kernel<P><<<blocks_for(N), cellwise::kThreads, 0, stream>>>(raw(diag), raw(inv), E, o.ne, bp3);
    BP_CK(cudaGetLastError());
    const cellwise::DssJacobi<decltype(dss)> precond{dss, raw(inv), raw(s), N};

    cellwise::FcgReduction red(MPI_COMM_SELF);
    cudaEvent_t t0, t1;
    BP_CK(cudaEventCreate(&t0));
    BP_CK(cudaEventCreate(&t1));
    BP_CK(cudaEventRecord(t0, stream));
    const cellwise::FcgResult res =
        cellwise::fcg(apply, precond, raw(b), raw(x), N, raw(r), raw(q), raw(z), raw(pdir), red, o.tol, o.maxit, stream);
    BP_CK(cudaEventRecord(t1, stream));
    BP_CK(cudaEventSynchronize(t1));
    float solve_ms = 0.0f;
    BP_CK(cudaEventElapsedTime(&solve_ms, t0, t1));
    difference_kernel<<<blocks_for(N), cellwise::kThreads, 0, stream>>>(raw(s), raw(x), raw(u), N);
    BP_CK(cudaStreamSynchronize(stream));
    const double err = thrust::reduce(s.begin(), s.end(), 0.0, thrust::maximum<double>());
    double umax = 0.0;
    {
        std::vector<double> hu(N);
        thrust::copy(u.begin(), u.end(), hu.begin());
        for (double v : hu) umax = std::fmax(umax, std::fabs(v));
    }

    auto time_ms = [&](auto&& f) {
        f();
        BP_CK(cudaEventRecord(t0, stream));
        for (int i = 0; i < o.reps; ++i) f();
        BP_CK(cudaEventRecord(t1, stream));
        BP_CK(cudaEventSynchronize(t1));
        float ms = 0.0f;
        BP_CK(cudaEventElapsedTime(&ms, t0, t1));
        return ms / o.reps;
    };
    const float op_ms = time_ms([&] { apply(raw(u), raw(y), stream); });
    const float dss_ms = time_ms([&] { dss(raw(y), raw(s), stream); });
    const float both_ms = time_ms([&] {
        apply(raw(u), raw(y), stream);
        dss(raw(y), raw(s), stream);
    });
    const double op_bytes = (2.0 * N + (double)E * K * q3) * 8.0, dss_bytes = 2.0 * N * 8.0;
    const bool pass = std::fabs(oracle_err) <= 1e-12 * oracle_scale &&
                      (bp3 ? max1 <= 1e-12 * maxu : std::fabs(sum1 - 1.0) <= 1e-12) && err <= 1e-8 * umax;

    printf("%s, p=%d, q=%d, %d^3 elements (%lld), %lld unique DoFs, %lld cell-wise, deform %.3f\n", o.spec.c_str(), P,
           Q, o.ne, E, unique, N, o.deform);
    if (bp3)
        printf("  BP3 applied to 1: max %.1e (of max %.1e for u)\n", max1, maxu);
    else
        printf("  BP1 applied to 1, summed: %.15f (the volume: 1)\n", sum1);
    printf("  first %lld elements vs direct quadrature sums: %.1e of max %.1e\n", Ec, oracle_err, oracle_scale);
    printf("  FCG + DSS-Jacobi: %d iterations, sqrt(<r,z>/<r0,z0>) = %.3e\n", res.iterations, res.relative);
    for (int k : {1, 5, 10, 20, 50, 100})
        if (k <= (int)res.history.size()) printf("  it %3d  %.6e\n", k, res.history[k - 1]);
    printf("  max |x - u| / max |u| = %.3e  %s\n", err / umax, pass ? "PASS" : "FAIL");
    printf("  operator     %8.3f ms  %6.2f GDoF/s unique  %6.2f cell-wise  %6.0f GB/s (u, y, geometry)\n", op_ms,
           unique / (op_ms * 1e6), N / (op_ms * 1e6), op_bytes / (op_ms * 1e6));
    printf("  DSS          %8.3f ms  %6.2f GDoF/s unique  %6.0f GB/s  (%s)\n", dss_ms, unique / (dss_ms * 1e6),
           dss_bytes / (dss_ms * 1e6), dss_gen ? "MARSIR-generated" : "hand gather");
    printf("  operator+DSS %8.3f ms  %6.2f GDoF/s unique  (the assembled operator)\n", both_ms,
           unique / (both_ms * 1e6));
    printf("  FCG iteration %7.3f ms  (solve %.1f ms)\n", res.iterations ? solve_ms / res.iterations : 0.0f, solve_ms);
    BP_CK(cudaStreamDestroy(stream));
    return pass ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Options o;
    const char* error = nullptr;
    for (int i = 1; i < argc && !error; i += 2) {
        if (i + 1 >= argc) error = "an option needs a value";
        else if (!strcmp(argv[i], "--ptx-dir")) o.ptx_dir = argv[i + 1];
        else if (!strcmp(argv[i], "--spec")) o.spec = argv[i + 1];
        else if (!strcmp(argv[i], "--p")) o.p = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--ne")) o.ne = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--deform")) o.deform = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--tol")) o.tol = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--maxit")) o.maxit = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--reps")) o.reps = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--dss-ptx")) o.dss_ptx = argv[i + 1];
        else error = "unknown option";
    }
    int ranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (!error && o.ptx_dir.empty()) error = "--ptx-dir is required";
    if (!error && o.spec != "bp1" && o.spec != "bp3") error = "--spec is bp1 or bp3";
    if (!error && o.p != 1 && o.p != 3 && o.p != 7) error = "--p is 1, 3 or 7 (the generated degrees)";
    if (!error && ranks != 1) error = "one rank for now";
    if (error) {
        fprintf(stderr, "%s\nusage: %s --ptx-dir <dir> --spec bp1|bp3 --p 1|3|7 [--ne N] [--deform a] [--tol t] "
                        "[--maxit k] [--reps r] [--dss-ptx file]\n", error, argv[0]);
        MPI_Finalize();
        return 1;
    }
    const int status = o.p == 1 ? run<1>(o) : o.p == 3 ? run<3>(o) : run<7>(o);
    MPI_Finalize();
    return status;
}
