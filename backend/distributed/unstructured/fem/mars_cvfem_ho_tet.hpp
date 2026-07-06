#ifndef MARS_CVFEM_HO_TET_HPP
#define MARS_CVFEM_HO_TET_HPP

// High-order CVFEM on tetrahedra via collapsed-coordinate (Duffy) sum-factorization with a
// BOX-PARTITION of control volumes. Tet counterpart of the hex Knaus CVFEM (mars_cvfem_ho_apply).
//
// DOF are nodal on the collapsed GLL tensor grid (degenerate-hex view of the tet; coincident
// collapsed nodes are merged by the DOF handler). The control volume of node (i,j,k) is the mapped
// collapsed box [e_i,e_i+1]^3, with edges e = {0, the P interior Gauss stations, 1}. Sub-control
// surfaces are the mapped box faces at the stations; each panel is integrated over its own node's
// tangential box with m-point Gauss. The panel flux integrand (Nanson area-normal form) is
// POLYNOMIAL in the tangential collapsed coordinates, so the quadrature is EXACT and the scheme is
// exactly conservative AND exactly consistent (A.1=0 and A.linear=0 to machine precision) --
// validated in tet_ho_ref rung 10c. Non-symmetric (Petrov-Galerkin): use GMRES, not CG.
//
// Geometry per affine element is only MMt = B^-1 B^-T and detB (10 doubles) -- the metric is
// recomputed from the collapse Jacobian at apply time (matrix-free style, no stored d_G).

#include "mars_ho_tet_basis.hpp"
#include <vector>
#include <array>
#include <cmath>

namespace mars {
namespace fem {

// GLL nodes on [0,1]: endpoints + roots of P'_{N-1} via Newton.
inline std::vector<double> gllNodes01(int N) {
    std::vector<double> z(N);
    z[0] = -1.0; z[N-1] = 1.0;
    const int deg = N - 1;
    for (int i = 1; i < N - 1; ++i) {
        double x = -std::cos(M_PI * i / (double)deg);
        for (int it = 0; it < 100; ++it) {
            const double f  = djacobiP(deg, 0.0, 0.0, x);
            const double fp = 0.25 * (deg + 1) * (deg + 2) * jacobiP(deg - 2, 2.0, 2.0, x);
            const double dx = f / fp;
            x -= dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        z[i] = x;
    }
    for (int i = 0; i < N; ++i) z[i] = 0.5 * (1.0 + z[i]);
    return z;
}

inline double lagr01(int j, double x, const std::vector<double>& z) {
    const int N = (int)z.size(); double v = 1;
    for (int m = 0; m < N; ++m) if (m != j) v *= (x - z[m]) / (z[j] - z[m]);
    return v;
}
inline double dlagr01(int j, double x, const std::vector<double>& z) {
    const int N = (int)z.size();
    double denom = 1; for (int m = 0; m < N; ++m) if (m != j) denom *= (z[j] - z[m]);
    double num = 0;
    for (int i = 0; i < N; ++i) if (i != j) { double pr = 1;
        for (int m = 0; m < N; ++m) if (m != j && m != i) pr *= (x - z[m]); num += pr; }
    return num / denom;
}

// Per-order tables. All coordinates on [0,1].
template <typename RealType>
struct HoCvfemTetOps {
    int P = 0, n = 0, m = 0, NN = 0;         // n=P+1 nodes/axis, m Gauss pts per box, NN=n^3
    std::vector<RealType> Z;                  // n GLL nodes
    std::vector<RealType> xs;                 // P stations (interior Gauss)
    std::vector<RealType> yq, wqv;            // box-quad pts/weights [n*m] (box-major)
    std::vector<RealType> Ls, Ds;             // [P*n]  lagr/dlagr at stations
    std::vector<RealType> Lq, Dq;             // [n*m*n] lagr/dlagr at box-quad pts
};

// stationType selects the sub-control-surface placement (box edges / control-volume boundaries):
//   0 = P interior Gauss-Legendre (default; Lin-Yang-Zou-optimal on tensor meshes),
//   1 = midpoints between the P+1 GLL DOF nodes (median-dual style),
//   2 = P uniform interior points,  3 = interior GLL nodes of order P+1.
// All keep the scheme exactly conservative + consistent (per-box quadrature stays exact regardless
// of where the boundaries sit); only the L2 rate changes -- p=2 order-condition study, see paper.
template <typename RealType>
HoCvfemTetOps<RealType> buildHoCvfemTetOps(int P, int mOverride = -1, int stationType = 0) {
    HoCvfemTetOps<RealType> o;
    // panel integrand degree <= P+1 per tangential variable (verified), so 2m-1 >= P+1 is exact
    o.P = P; o.n = P + 1; o.m = (mOverride > 0) ? mOverride : (P + 3) / 2; o.NN = o.n * o.n * o.n;
    const int n = o.n, m = o.m;
    auto Zd = gllNodes01(n);
    o.Z.assign(Zd.begin(), Zd.end());
    o.xs.resize(P);
    if (stationType == 1) {
        for (int k = 0; k < P; ++k) o.xs[k] = 0.5 * (Zd[k] + Zd[k + 1]);
    } else if (stationType == 2) {
        for (int k = 0; k < P; ++k) o.xs[k] = (k + 1.0) / (P + 1.0);
    } else if (stationType == 3) {
        auto Zi = gllNodes01(P + 2);
        for (int k = 0; k < P; ++k) o.xs[k] = Zi[k + 1];
    } else {
        GaussJacobi gs = gaussJacobi(P, 0.0, 0.0);
        for (int k = 0; k < P; ++k) o.xs[k] = 0.5 * (1.0 + gs.x[k]);
    }
    std::vector<double> e(n + 1); e[0] = 0.0; for (int k = 0; k < P; ++k) e[k+1] = o.xs[k]; e[n] = 1.0;
    GaussJacobi gq = gaussJacobi(m, 0.0, 0.0);
    o.yq.resize((size_t)n * m); o.wqv.resize((size_t)n * m);
    for (int i = 0; i < n; ++i) for (int q = 0; q < m; ++q) {
        const double len = e[i+1] - e[i];
        o.yq[(size_t)i*m+q] = e[i] + len * 0.5 * (1.0 + gq.x[q]);
        o.wqv[(size_t)i*m+q] = len * 0.5 * gq.w[q];
    }
    o.Ls.resize((size_t)P * n); o.Ds.resize((size_t)P * n);
    for (int k = 0; k < P; ++k) for (int j = 0; j < n; ++j) {
        o.Ls[(size_t)k*n+j] = lagr01(j, o.xs[k], Zd); o.Ds[(size_t)k*n+j] = dlagr01(j, o.xs[k], Zd); }
    o.Lq.resize((size_t)n * m * n); o.Dq.resize((size_t)n * m * n);
    for (int t = 0; t < n * m; ++t) for (int j = 0; j < n; ++j) {
        o.Lq[(size_t)t*n+j] = lagr01(j, o.yq[t], Zd); o.Dq[(size_t)t*n+j] = dlagr01(j, o.yq[t], Zd); }
    return o;
}

// Affine element geometry: MMt = B^-1 B^-T (row-major 9) and |detB|, from the 4 corners.
template <typename RealType>
void computeTetGeom(const RealType corners[4][3], RealType MMt[9], RealType& detB) {
    double B[9];
    for (int d = 0; d < 3; ++d) {
        B[d*3+0] = corners[1][d] - corners[0][d];
        B[d*3+1] = corners[2][d] - corners[0][d];
        B[d*3+2] = corners[3][d] - corners[0][d];
    }
    const double det = B[0]*(B[4]*B[8]-B[5]*B[7]) - B[1]*(B[3]*B[8]-B[5]*B[6]) + B[2]*(B[3]*B[7]-B[4]*B[6]);
    detB = (RealType)std::fabs(det);
    const double inv = 1.0 / det;
    double M[9];
    M[0]=(B[4]*B[8]-B[5]*B[7])*inv; M[1]=(B[2]*B[7]-B[1]*B[8])*inv; M[2]=(B[1]*B[5]-B[2]*B[4])*inv;
    M[3]=(B[5]*B[6]-B[3]*B[8])*inv; M[4]=(B[0]*B[8]-B[2]*B[6])*inv; M[5]=(B[2]*B[3]-B[0]*B[5])*inv;
    M[6]=(B[3]*B[7]-B[4]*B[6])*inv; M[7]=(B[1]*B[6]-B[0]*B[7])*inv; M[8]=(B[0]*B[4]-B[1]*B[3])*inv;
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
        double v = 0; for (int c = 0; c < 3; ++c) v += M[a*3+c]*M[b*3+c]; MMt[a*3+b] = (RealType)v; }
}

#if defined(USE_CUDA)
#define MARS_TET_HD __host__ __device__
#else
#define MARS_TET_HD
#endif

// Directional flux factor at a collapsed point: row `dir` of Ghat = J_c MMt J_c^T times detB*detJd.
// Pointwise identical to the Nanson form grad_x(u).cof(B J_d)e_dir; all evaluations interior.
template <typename RealType>
MARS_TET_HD inline void cvfemTetFluxRow(int dir, RealType a, RealType b, RealType c,
                                        const RealType* MMt, RealType detB, RealType row[3]) {
    const RealType r = a*(1-b)*(1-c), s = b*(1-c);
    const RealType D1 = (RealType)1 - s - c, D2 = (RealType)1 - c;
    const RealType Jc[9] = { (RealType)1/D1, r/(D1*D1), r/(D1*D1),
                             0, (RealType)1/D2, s/(D2*D2),
                             0, 0, (RealType)1 };
    // row_t = sum_u (Jc MMt Jc^T)[dir][t] * detB * (1-b)(1-c)^2
    RealType JM[3];
    for (int t = 0; t < 3; ++t) { RealType v = 0;
        for (int u = 0; u < 3; ++u) v += Jc[dir*3+u]*MMt[u*3+t]; JM[t] = v; }
    const RealType w = detB * ((RealType)1-b) * ((RealType)1-c) * ((RealType)1-c);
    for (int t = 0; t < 3; ++t) { RealType v = 0;
        for (int u = 0; u < 3; ++u) v += JM[u]*Jc[t*3+u]; row[t] = w * v; }
}

// Reference apply (dense per-panel-point gradient). Oracle for the sum-fac version and the GPU.
template <typename RealType>
void applyHoCvfemTetElementRef(const HoCvfemTetOps<RealType>& o, const RealType* MMt, RealType detB,
                               const RealType* u, RealType* y) {
    const int P = o.P, n = o.n, m = o.m;
    for (int I = 0; I < o.NN; ++I) y[I] = 0;
    for (int dir = 0; dir < 3; ++dir) {
        const int t0 = (dir+1)%3, t1 = (dir+2)%3;
        for (int k = 0; k < P; ++k)
        for (int tq0 = 0; tq0 < n*m; ++tq0)
        for (int tq1 = 0; tq1 < n*m; ++tq1) {
            RealType g[3] = {0,0,0};
            for (int ia = 0; ia < n; ++ia) for (int ib = 0; ib < n; ++ib) for (int ic = 0; ic < n; ++ic) {
                int idx[3] = {ia, ib, ic};
                RealType li[3], di[3];
                li[dir] = o.Ls[(size_t)k*n+idx[dir]];   di[dir] = o.Ds[(size_t)k*n+idx[dir]];
                li[t0]  = o.Lq[(size_t)tq0*n+idx[t0]];  di[t0]  = o.Dq[(size_t)tq0*n+idx[t0]];
                li[t1]  = o.Lq[(size_t)tq1*n+idx[t1]];  di[t1]  = o.Dq[(size_t)tq1*n+idx[t1]];
                const RealType ua = u[((size_t)ia*n+ib)*n+ic];
                g[0] += ua*(di[0]*li[1]*li[2]);
                g[1] += ua*(li[0]*di[1]*li[2]);
                g[2] += ua*(li[0]*li[1]*di[2]);
            }
            RealType crd[3]; crd[dir] = o.xs[k]; crd[t0] = o.yq[tq0]; crd[t1] = o.yq[tq1];
            RealType row[3];
            cvfemTetFluxRow<RealType>(dir, crd[0], crd[1], crd[2], MMt, detB, row);
            const RealType piece = (row[0]*g[0]+row[1]*g[1]+row[2]*g[2]) * o.wqv[tq0] * o.wqv[tq1];
            const int b0 = tq0 / m, b1 = tq1 / m;
            for (int in_ = 0; in_ < n; ++in_) {
                const RealType wn = (k == in_-1) ? (RealType)1 : (k == in_) ? (RealType)-1 : (RealType)0;
                if (wn == (RealType)0) continue;
                int nid[3]; nid[dir] = in_; nid[t0] = b0; nid[t1] = b1;
                y[((size_t)nid[0]*n+nid[1])*n+nid[2]] += wn * piece;
            }
        }
    }
}

// Sum-factorized apply: per direction, the gradient at the (station x box-quad x box-quad) tensor
// grid is built by three successive 1D contractions (O(p^4 m^2) instead of O(p^6 m^2) dense).
template <typename RealType>
void applyHoCvfemTetElement(const HoCvfemTetOps<RealType>& o, const RealType* MMt, RealType detB,
                            const RealType* u, RealType* y) {
    const int P = o.P, n = o.n, m = o.m, nm = n*m;
    for (int I = 0; I < o.NN; ++I) y[I] = 0;
    // scratch: s1 [n_out0 * n * n], s2 [n_out0 * nm * n], g [3][P * nm * nm]
    std::vector<RealType> s1((size_t)nm*n*n), s2((size_t)nm*nm*n), g((size_t)3*P*nm*nm);

    for (int dir = 0; dir < 3; ++dir) {
        const int t0 = (dir+1)%3, t1 = (dir+2)%3;
        // for each gradient component: pick per-axis table (station table on axis dir, box-quad
        // table on tangential axes; derivative table on the component's axis).
        for (int comp = 0; comp < 3; ++comp) {
            auto tab = [&](int axis, int station) -> const RealType* {
                // returns row pointer for this axis at output index `station`
                if (axis == dir) return ((axis == comp) ? &o.Ds[(size_t)station*n] : &o.Ls[(size_t)station*n]);
                return ((axis == comp) ? &o.Dq[(size_t)station*n] : &o.Lq[(size_t)station*n]);
            };
            const int out0 = (0 == dir) ? P : nm, out1 = (1 == dir) ? P : nm, out2 = (2 == dir) ? P : nm;
            // contract axis 0 (ia)
            for (int oa = 0; oa < out0; ++oa) { const RealType* T = tab(0, oa);
                for (int ib = 0; ib < n; ++ib) for (int ic = 0; ic < n; ++ic) {
                    RealType v = 0;
                    for (int ia = 0; ia < n; ++ia) v += T[ia] * u[((size_t)ia*n+ib)*n+ic];
                    s1[((size_t)oa*n+ib)*n+ic] = v; } }
            // contract axis 1 (ib)
            for (int oa = 0; oa < out0; ++oa) for (int ob = 0; ob < out1; ++ob) { const RealType* T = tab(1, ob);
                for (int ic = 0; ic < n; ++ic) {
                    RealType v = 0;
                    for (int ib = 0; ib < n; ++ib) v += T[ib] * s1[((size_t)oa*n+ib)*n+ic];
                    s2[((size_t)oa*out1+ob)*n+ic] = v; } }
            // contract axis 2 (ic) -> g[comp] on the panel grid (k, tq0, tq1) in axis order
            for (int oa = 0; oa < out0; ++oa) for (int ob = 0; ob < out1; ++ob)
            for (int oc = 0; oc < out2; ++oc) { const RealType* T = tab(2, oc);
                RealType v = 0;
                for (int ic = 0; ic < n; ++ic) v += T[ic] * s2[((size_t)oa*out1+ob)*n+ic];
                // map (oa,ob,oc) to (k, tq0, tq1): k is the dir axis, tq0/tq1 the tangential axes
                int idx[3] = {oa, ob, oc};
                const int k = idx[dir], q0 = idx[t0], q1 = idx[t1];
                g[(((size_t)comp*P+k)*nm+q0)*nm+q1] = v; }
        }
        // flux + incidence scatter (tangential node = box index; no cross-node coupling)
        for (int k = 0; k < P; ++k)
        for (int tq0 = 0; tq0 < nm; ++tq0)
        for (int tq1 = 0; tq1 < nm; ++tq1) {
            RealType crd[3]; crd[dir] = o.xs[k]; crd[t0] = o.yq[tq0]; crd[t1] = o.yq[tq1];
            RealType row[3];
            cvfemTetFluxRow<RealType>(dir, crd[0], crd[1], crd[2], MMt, detB, row);
            const RealType g0 = g[(((size_t)0*P+k)*nm+tq0)*nm+tq1];
            const RealType g1 = g[(((size_t)1*P+k)*nm+tq0)*nm+tq1];
            const RealType g2 = g[(((size_t)2*P+k)*nm+tq0)*nm+tq1];
            const RealType piece = (row[0]*g0+row[1]*g1+row[2]*g2) * o.wqv[tq0] * o.wqv[tq1];
            const int b0 = tq0 / m, b1 = tq1 / m;
            for (int in_ = 0; in_ < n; ++in_) {
                const RealType wn = (k == in_-1) ? (RealType)1 : (k == in_) ? (RealType)-1 : (RealType)0;
                if (wn == (RealType)0) continue;
                int nid[3]; nid[dir] = in_; nid[t0] = b0; nid[t1] = b1;
                y[((size_t)nid[0]*n+nid[1])*n+nid[2]] += wn * piece;
            }
        }
    }
}

// Universal reference matrices for AFFINE elements: the operator is LINEAR in the entries of
// MMt = B^-1 B^-T and proportional to detB, so A_e = detB * sum_c coef_c(e) * K_c with SIX
// element-independent NN x NN matrices (c over the symmetric index pairs of MMt). Built once per
// order with exact quadrature -- exactness/conservation/consistency are inherited -- and shared by
// every element. The apply then needs no quadrature and no rational metric at runtime; per-element
// geometry stays 10 doubles. This is the affine fast path; the panel kernels remain for curved.
// Layout: K[(c*NN + i)*NN + j], coef order {MMt00, MMt01, MMt02, MMt11, MMt12, MMt22}.
template <typename RealType>
std::vector<RealType> buildCvfemTetRefMatrices(const HoCvfemTetOps<RealType>& o) {
    const int NN = o.NN;
    static const int uv[6][2] = {{0,0},{0,1},{0,2},{1,1},{1,2},{2,2}};
    std::vector<RealType> K((size_t)6*NN*NN);
    std::vector<RealType> u(NN), y(NN);
    for (int c = 0; c < 6; ++c) {
        RealType E[9] = {0,0,0, 0,0,0, 0,0,0};
        E[uv[c][0]*3 + uv[c][1]] = 1;
        E[uv[c][1]*3 + uv[c][0]] = 1;   // no-op on the diagonal pairs
        for (int j = 0; j < NN; ++j) {
            std::fill(u.begin(), u.end(), (RealType)0); u[j] = 1;
            applyHoCvfemTetElement<RealType>(o, E, (RealType)1, u.data(), y.data());
            for (int i = 0; i < NN; ++i) K[((size_t)c*NN + i)*NN + j] = y[i];
        }
    }
    return K;
}

#if defined(USE_CUDA)
// One thread per element, dense per-panel-point gradient (correctness-gate kernel; the staged
// shared-memory sum-fac kernel is the perf follow-up). Gathers nodal DOF via d_elemDof (repeated
// indices at merged collapsed nodes are fine), scatters with atomicAdd.
// Tables laid out exactly as HoCvfemTetOps; geometry per element: d_geom[e*10] = {MMt[9], detB}.
template <typename RealType, int P, int M>
__global__ void ho_cvfem_tet_apply_cg(
        const RealType* __restrict__ d_u_dof, RealType* __restrict__ d_y_dof,
        const int* __restrict__ d_elemDof,
        const RealType* __restrict__ d_xs, const RealType* __restrict__ d_yq, const RealType* __restrict__ d_wq,
        const RealType* __restrict__ d_Ls, const RealType* __restrict__ d_Ds,
        const RealType* __restrict__ d_Lq, const RealType* __restrict__ d_Dq,
        const RealType* __restrict__ d_geom, size_t numElements) {
    constexpr int n = P+1, NN = n*n*n, nm = n*M;
    const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    const int* dof = d_elemDof + e*NN;
    const RealType* MMt = d_geom + e*10;
    const RealType detB = d_geom[e*10+9];

    RealType ue[NN], ye[NN];
    for (int i = 0; i < NN; ++i) { ue[i] = d_u_dof[dof[i]]; ye[i] = 0; }

    for (int dir = 0; dir < 3; ++dir) {
        const int t0 = (dir+1)%3, t1 = (dir+2)%3;
        for (int k = 0; k < P; ++k)
        for (int tq0 = 0; tq0 < nm; ++tq0)
        for (int tq1 = 0; tq1 < nm; ++tq1) {
            RealType g[3] = {0,0,0};
            for (int ia = 0; ia < n; ++ia) for (int ib = 0; ib < n; ++ib) for (int ic = 0; ic < n; ++ic) {
                int idx[3] = {ia, ib, ic};
                RealType li[3], di[3];
                li[dir] = d_Ls[k*n+idx[dir]];   di[dir] = d_Ds[k*n+idx[dir]];
                li[t0]  = d_Lq[tq0*n+idx[t0]];  di[t0]  = d_Dq[tq0*n+idx[t0]];
                li[t1]  = d_Lq[tq1*n+idx[t1]];  di[t1]  = d_Dq[tq1*n+idx[t1]];
                const RealType ua = ue[(ia*n+ib)*n+ic];
                g[0] += ua*(di[0]*li[1]*li[2]);
                g[1] += ua*(li[0]*di[1]*li[2]);
                g[2] += ua*(li[0]*li[1]*di[2]);
            }
            RealType crd[3]; crd[dir] = d_xs[k]; crd[t0] = d_yq[tq0]; crd[t1] = d_yq[tq1];
            RealType row[3];
            cvfemTetFluxRow<RealType>(dir, crd[0], crd[1], crd[2], MMt, detB, row);
            const RealType piece = (row[0]*g[0]+row[1]*g[1]+row[2]*g[2]) * d_wq[tq0] * d_wq[tq1];
            const int b0 = tq0 / M, b1 = tq1 / M;
            for (int in_ = 0; in_ < n; ++in_) {
                const RealType wn = (k == in_-1) ? (RealType)1 : (k == in_) ? (RealType)-1 : (RealType)0;
                if (wn == (RealType)0) continue;
                int nid[3]; nid[dir] = in_; nid[t0] = b0; nid[t1] = b1;
                ye[(nid[0]*n+nid[1])*n+nid[2]] += wn * piece;
            }
        }
    }
    for (int i = 0; i < NN; ++i) atomicAdd(&d_y_dof[dof[i]], ye[i]);
}
// Performance kernel: one BLOCK per element, staged shared-memory sum-factorization.
// Per direction: (A) contract the normal axis to the P stations (two variants: interp Ls and
// derivative Ds -- the tangential components share the interp one), then per station (B),(C)
// contract the two tangential axes to the box-quad grid, (D) flux + incidence scatter into the
// element accumulator in shared memory (smem atomics), finally one global atomicAdd per node.
// Static smem fits p<=5; higher orders need the dynamic-smem variant (follow-up).
template <typename RealType, int P, int M, int BlockSize = 256>
__global__ void __launch_bounds__(BlockSize)
ho_cvfem_tet_apply_cg_smem(
        const RealType* __restrict__ d_u_dof, RealType* __restrict__ d_y_dof,
        const int* __restrict__ d_elemDof,
        const RealType* __restrict__ d_xs, const RealType* __restrict__ d_yq, const RealType* __restrict__ d_wq,
        const RealType* __restrict__ d_Ls, const RealType* __restrict__ d_Ds,
        const RealType* __restrict__ d_Lq, const RealType* __restrict__ d_Dq,
        const RealType* __restrict__ d_geom, size_t numElements) {
    constexpr int n = P+1, NN = n*n*n, nm = n*M;
    const size_t e = blockIdx.x;
    if (e >= numElements) return;
    const int* dof = d_elemDof + e*NN;
    const RealType* MMt = d_geom + e*10;
    const RealType detB = d_geom[e*10+9];
    const int tid = threadIdx.x;

    __shared__ RealType su[NN], sy[NN];
    __shared__ RealType f1L[P*n*n], f1D[P*n*n];      // normal-axis contraction (interp / deriv)
    __shared__ RealType tg[3][nm*n];                  // after tangential-axis-0 contraction
    __shared__ RealType g[3][nm*nm];                  // gradient components on the panel grid

    for (int i = tid; i < NN; i += BlockSize) { su[i] = d_u_dof[dof[i]]; sy[i] = 0; }
    __syncthreads();

    for (int dir = 0; dir < 3; ++dir) {
        const int t0 = (dir+1)%3, t1 = (dir+2)%3;

        // (A) contract the normal axis: f1{L,D}[k][jt0][jt1], jt0/jt1 = GLL indices on t0/t1
        for (int w = tid; w < P*n*n; w += BlockSize) {
            const int k = w/(n*n), jt0 = (w/n)%n, jt1 = w%n;
            RealType vL = 0, vD = 0;
            for (int jn = 0; jn < n; ++jn) {
                int idx[3]; idx[dir] = jn; idx[t0] = jt0; idx[t1] = jt1;
                const RealType ua = su[(idx[0]*n+idx[1])*n+idx[2]];
                vL += d_Ls[k*n+jn]*ua; vD += d_Ds[k*n+jn]*ua;
            }
            f1L[w] = vL; f1D[w] = vD;
        }
        __syncthreads();

        for (int k = 0; k < P; ++k) {
            // (B) contract tangential axis t0 to the box-quad grid:
            //   comp==dir : f1D through Lq ; comp==t0 : f1L through Dq ; comp==t1 : f1L through Lq
            for (int w = tid; w < 3*nm*n; w += BlockSize) {
                const int comp = w/(nm*n), q0 = (w/n)%nm, jt1 = w%n;
                const RealType* src = (comp == 0) ? &f1D[k*n*n] : &f1L[k*n*n];   // comp index: 0=dir,1=t0,2=t1
                const RealType* T   = (comp == 1) ? &d_Dq[q0*n] : &d_Lq[q0*n];
                RealType v = 0;
                for (int jt0 = 0; jt0 < n; ++jt0) v += T[jt0]*src[jt0*n+jt1];
                tg[comp][q0*n+jt1] = v;
            }
            __syncthreads();
            // (C) contract tangential axis t1
            for (int w = tid; w < 3*nm*nm; w += BlockSize) {
                const int comp = w/(nm*nm), q0 = (w/nm)%nm, q1 = w%nm;
                const RealType* T = (comp == 2) ? &d_Dq[q1*n] : &d_Lq[q1*n];
                RealType v = 0;
                for (int jt1 = 0; jt1 < n; ++jt1) v += T[jt1]*tg[comp][q0*n+jt1];
                g[comp][q0*nm+q1] = v;
            }
            __syncthreads();
            // (D) flux + incidence: node k+1 gets +piece, node k gets -piece (Del convention)
            for (int w = tid; w < nm*nm; w += BlockSize) {
                const int q0 = w/nm, q1 = w%nm;
                RealType crd[3]; crd[dir] = d_xs[k]; crd[t0] = d_yq[q0]; crd[t1] = d_yq[q1];
                RealType row[3];
                cvfemTetFluxRow<RealType>(dir, crd[0], crd[1], crd[2], MMt, detB, row);
                // g comp order is (dir,t0,t1); row is (a,b,c): remap
                RealType gx[3]; gx[dir] = g[0][w]; gx[t0] = g[1][w]; gx[t1] = g[2][w];
                const RealType piece = (row[0]*gx[0]+row[1]*gx[1]+row[2]*gx[2]) * d_wq[q0]*d_wq[q1];
                const int b0 = q0/M, b1 = q1/M;
                int nid[3]; nid[t0] = b0; nid[t1] = b1;
                nid[dir] = k+1; atomicAdd(&sy[(nid[0]*n+nid[1])*n+nid[2]],  piece);
                nid[dir] = k;   atomicAdd(&sy[(nid[0]*n+nid[1])*n+nid[2]], -piece);
            }
            __syncthreads();
        }
    }
    for (int i = tid; i < NN; i += BlockSize) atomicAdd(&d_y_dof[dof[i]], sy[i]);
}

// Affine fast path: y_e = detB * sum_c coef_c (K_c u_e) with the universal reference matrices.
// One block per element; u_e staged in smem; the K reads are shared across all blocks (L2).
template <typename RealType, int P, int BlockSize = 128>
__global__ void __launch_bounds__(BlockSize)
ho_cvfem_tet_apply_cg_univ(
        const RealType* __restrict__ d_u_dof, RealType* __restrict__ d_y_dof,
        const int* __restrict__ d_elemDof,
        const RealType* __restrict__ d_K,      // [6*NN*NN]
        const RealType* __restrict__ d_geom, size_t numElements) {
    constexpr int n = P+1, NN = n*n*n;
    const size_t e = blockIdx.x;
    if (e >= numElements) return;
    const int* dof = d_elemDof + e*NN;
    const RealType* MMt = d_geom + e*10;
    const RealType detB = d_geom[e*10+9];
    const RealType c6[6] = { detB*MMt[0], detB*MMt[1], detB*MMt[2],
                             detB*MMt[4], detB*MMt[5], detB*MMt[8] };
    __shared__ RealType su[NN];
    for (int i = threadIdx.x; i < NN; i += BlockSize) su[i] = d_u_dof[dof[i]];
    __syncthreads();
    for (int i = threadIdx.x; i < NN; i += BlockSize) {
        RealType acc = 0;
        const RealType *K0 = &d_K[(size_t)0*NN*NN + (size_t)i*NN], *K1 = &d_K[(size_t)1*NN*NN + (size_t)i*NN],
                       *K2 = &d_K[(size_t)2*NN*NN + (size_t)i*NN], *K3 = &d_K[(size_t)3*NN*NN + (size_t)i*NN],
                       *K4 = &d_K[(size_t)4*NN*NN + (size_t)i*NN], *K5 = &d_K[(size_t)5*NN*NN + (size_t)i*NN];
        for (int j = 0; j < NN; ++j) {
            const RealType uj = su[j];
            acc += (c6[0]*K0[j] + c6[1]*K1[j] + c6[2]*K2[j] +
                    c6[3]*K3[j] + c6[4]*K4[j] + c6[5]*K5[j]) * uj;
        }
        atomicAdd(&d_y_dof[dof[i]], acc);
    }
}

#endif  // USE_CUDA

}  // namespace fem
}  // namespace mars
#endif
