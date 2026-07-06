#ifndef MARS_HO_TET_BASIS_HPP
#define MARS_HO_TET_BASIS_HPP

// Collapsed-coordinate (Duffy / PKD) basis + quadrature for high-order tetrahedra.
// Tet twin of mars_cvfem_ho_basis.hpp (which is GLL tensor product for hexes). Here the reference
// element is the unit tet {r,s,t>=0, r+s+t<=1}, expanded in the Proriol-Koornwinder-Dubiner modal
// basis and integrated on the Gauss-Jacobi Duffy grid. Host-only (no device code); the tables it
// builds are uploaded once per order and consumed by mars_ho_laplacian_tet.hpp.
//
// The math is the validated tet_ho_ref host reference (rung 9). Modal index (p,q,r), p+q+r<=P.
//   psi_pqr = A_p(a) B_pq(b) C_pqr(c),   a=r/(1-s-t), b=s/(1-t), c=t   (collapse)
//   A_p(a)   = P_p^{(0,0)}(2a-1)
//   B_pq(b)  = (1-b)^p       P_q^{(2p+1,0)}(2b-1)
//   C_pqr(c) = (1-c)^{p+q}   P_r^{(2p+2q+2,0)}(2c-1)

#include <cmath>
#include <cstddef>
#include <vector>
#include <array>

namespace mars {
namespace fem {

// Jacobi polynomial P_n^{(a,b)}(x) via the 3-term recurrence (a=b=0 -> Legendre).
inline double jacobiP(int n, double a, double b, double x) {
    if (n == 0) return 1.0;
    double p0 = 1.0, p1 = 0.5 * (a - b) + 0.5 * (a + b + 2.0) * x;
    if (n == 1) return p1;
    for (int k = 2; k <= n; ++k) {
        const double k2ab = 2.0 * k + a + b;
        const double c1 = 2.0 * k * (k + a + b) * (k2ab - 2.0);
        const double c2 = (k2ab - 1.0) * (a * a - b * b);
        const double c3 = (k2ab - 2.0) * (k2ab - 1.0) * k2ab;
        const double c4 = 2.0 * (k + a - 1.0) * (k + b - 1.0) * k2ab;
        const double p2 = ((c2 + c3 * x) * p1 - c4 * p0) / c1;
        p0 = p1; p1 = p2;
    }
    return p1;
}
inline double djacobiP(int n, double a, double b, double x) {
    if (n == 0) return 0.0;
    return 0.5 * (n + a + b + 1.0) * jacobiP(n - 1, a + 1.0, b + 1.0, x);
}

struct GaussJacobi { std::vector<double> x, w; };

// n-point Gauss-Jacobi rule on [-1,1] for weight (1-x)^a (1+x)^b. Newton with Maehly deflation so the
// asymmetric (a>0) roots do not collapse -- the c-axis rule (a=2) fails without it.
inline GaussJacobi gaussJacobi(int n, double a, double b) {
    GaussJacobi q; q.x.resize(n); q.w.resize(n);
    for (int i = 0; i < n; ++i) {
        double xi = -std::cos(M_PI * (2 * i + 1) / (double)(2 * n));
        for (int it = 0; it < 200; ++it) {
            const double f = jacobiP(n, a, b, xi), fp = djacobiP(n, a, b, xi);
            double defl = 0.0;
            for (int j = 0; j < i; ++j) defl += 1.0 / (xi - q.x[j]);
            const double dx = f / (fp - f * defl);
            xi -= dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        q.x[i] = xi;
    }
    const double mu0 = std::pow(2.0, a + b + 1.0) * std::tgamma(a + 1.0) * std::tgamma(b + 1.0)
                     / std::tgamma(a + b + 2.0);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double dp = djacobiP(n, a, b, q.x[i]);
        q.w[i] = 1.0 / ((1.0 - q.x[i] * q.x[i]) * dp * dp);
        sum += q.w[i];
    }
    for (int i = 0; i < n; ++i) q.w[i] *= mu0 / sum;
    return q;
}

// Per-order tables consumed by the element apply. n = P+1 quad points per collapsed direction.
// Factor tables are indexed exactly as in tet_ho_ref rung 9:
//   A [p*n + i]                        (p<=P, i<n)      at a-nodes (Gauss-Legendre)
//   B [(p*(P+1)+q)*n + j]              at b-nodes (Gauss-Jacobi 1,0)
//   C [((p*(P+1)+q)*(P+1)+r)*n + k]    at c-nodes (Gauss-Jacobi 2,0)
// Ad/Bd/Cd are the coordinate derivatives (d/da etc = 2 d/dx on [-1,1]).
template <typename RealType>
struct HoTetOps {
    int P = 0, n = 0, W = 0;                 // W = P+1 (mode-index stride)
    int nq = 0, nModes = 0;                  // nq = n^3 quad points, nModes = (P+1)(P+2)(P+3)/6
    GaussJacobi ga, gb, gc;                  // a: GJ(0,0), b: GJ(1,0), c: GJ(2,0)
    std::vector<RealType> A, Ad, B, Bd, C, Cd;
    std::vector<std::array<int,3>> modes;    // (p,q,r), p+q+r<=P
};

template <typename RealType>
HoTetOps<RealType> buildHoTetOps(int P) {
    HoTetOps<RealType> o;
    o.P = P; o.n = P + 1; o.W = P + 1;
    const int n = o.n, W = o.W;
    o.nq = n * n * n;
    o.ga = gaussJacobi(n, 0.0, 0.0);
    o.gb = gaussJacobi(n, 1.0, 0.0);
    o.gc = gaussJacobi(n, 2.0, 0.0);
    for (int p = 0; p <= P; ++p) for (int q = 0; q <= P - p; ++q) for (int r = 0; r <= P - p - q; ++r)
        o.modes.push_back({p, q, r});
    o.nModes = (int)o.modes.size();

    o.A.assign((size_t)W * n, 0);  o.Ad.assign((size_t)W * n, 0);
    o.B.assign((size_t)W * W * n, 0);  o.Bd.assign((size_t)W * W * n, 0);
    o.C.assign((size_t)W * W * W * n, 0);  o.Cd.assign((size_t)W * W * W * n, 0);
    for (int p = 0; p <= P; ++p) for (int i = 0; i < n; ++i) {
        o.A[(size_t)p * n + i]  = jacobiP(p, 0.0, 0.0, o.ga.x[i]);
        o.Ad[(size_t)p * n + i] = 2.0 * djacobiP(p, 0.0, 0.0, o.ga.x[i]);
    }
    for (int p = 0; p <= P; ++p) for (int q = 0; q <= P - p; ++q) for (int j = 0; j < n; ++j) {
        const double x = o.gb.x[j], w = 0.5 * (1 - x);
        const double J = jacobiP(q, 2.0*p+1, 0.0, x), dJ = djacobiP(q, 2.0*p+1, 0.0, x);
        o.B[((size_t)p*W+q)*n+j]  = std::pow(w, p) * J;
        o.Bd[((size_t)p*W+q)*n+j] = ((p>0) ? -(double)p*std::pow(w,p-1)*J : 0.0) + std::pow(w,p)*2.0*dJ;
    }
    for (int p = 0; p <= P; ++p) for (int q = 0; q <= P - p; ++q) for (int r = 0; r <= P - p - q; ++r)
    for (int k = 0; k < n; ++k) {
        const int e = p + q; const double x = o.gc.x[k], w = 0.5 * (1 - x);
        const double J = jacobiP(r, 2.0*e+2, 0.0, x), dJ = djacobiP(r, 2.0*e+2, 0.0, x);
        o.C[(((size_t)p*W+q)*W+r)*n+k]  = std::pow(w, e) * J;
        o.Cd[(((size_t)p*W+q)*W+r)*n+k] = ((e>0) ? -(double)e*std::pow(w,e-1)*J : 0.0) + std::pow(w,e)*2.0*dJ;
    }
    return o;
}

// PKD basis value at a reference-tet point, with collapse-limit guards: at s+t=1 the a-collapse is
// 0/0 (and b at t=1); the (1-b)^p / (1-c)^{p+q} prefactors kill those modes, so any finite a,b give
// the correct limit -- we just must not divide by zero.
inline double pkdEval(int p, int q, int r, double rr, double ss, double tt) {
    const double d1 = 1.0 - ss - tt, d2 = 1.0 - tt;
    const double a = (d1 > 1e-12) ? rr / d1 : 0.0;
    const double b = (d2 > 1e-12) ? ss / d2 : 0.0;
    const double Pa = jacobiP(p, 0.0, 0.0, 2*a - 1);
    const double Pb = jacobiP(q, 2.0*p + 1.0, 0.0, 2*b - 1);
    const double Pc = jacobiP(r, 2.0*p + 2.0*q + 2.0, 0.0, 2*tt - 1);
    return Pa * std::pow(1.0 - b, p) * Pb * std::pow(1.0 - tt, p + q) * Pc;
}

// Dense NxN inverse, Gauss-Jordan with partial pivoting (host, once per order).
inline bool invertDense(std::vector<double> A, int N, std::vector<double>& Ainv) {
    Ainv.assign((size_t)N * N, 0.0);
    for (int i = 0; i < N; ++i) Ainv[(size_t)i * N + i] = 1.0;
    for (int c = 0; c < N; ++c) {
        int piv = c;
        for (int rr = c + 1; rr < N; ++rr)
            if (std::fabs(A[(size_t)rr * N + c]) > std::fabs(A[(size_t)piv * N + c])) piv = rr;
        if (std::fabs(A[(size_t)piv * N + c]) < 1e-14) return false;
        if (piv != c) for (int j = 0; j < N; ++j) {
            std::swap(A[(size_t)c * N + j], A[(size_t)piv * N + j]);
            std::swap(Ainv[(size_t)c * N + j], Ainv[(size_t)piv * N + j]);
        }
        const double d = 1.0 / A[(size_t)c * N + c];
        for (int j = 0; j < N; ++j) { A[(size_t)c * N + j] *= d; Ainv[(size_t)c * N + j] *= d; }
        for (int rr = 0; rr < N; ++rr) if (rr != c) {
            const double f = A[(size_t)rr * N + c];
            if (f == 0.0) continue;
            for (int j = 0; j < N; ++j) { A[(size_t)rr * N + j] -= f * A[(size_t)c * N + j]; Ainv[(size_t)rr * N + j] -= f * Ainv[(size_t)c * N + j]; }
        }
    }
    return true;
}

// Nodal layer for C0 assembly: equispaced barycentric nodes (vertex/edge/face/interior structure)
// + the PKD Vandermonde. Element apply becomes V^-T K_modal V^-1 -- the sum-fac core is untouched,
// and cross-element continuity reduces to plain nodal DOF matching. Equispaced conditioning is fine
// for p<=6; warp-blend nodes are a later refinement.
template <typename RealType>
struct HoTetNodal {
    int P = 0, Np = 0;
    std::vector<std::array<int,3>> bary;     // (i,j,k), i+j+k<=P; node (r,s,t)=(i,j,k)/P
    std::vector<RealType> V, Vinv;           // Np x Np, V[node*Np + mode], modes in HoTetOps order
};

template <typename RealType>
HoTetNodal<RealType> buildHoTetNodal(const HoTetOps<RealType>& o) {
    HoTetNodal<RealType> nd;
    nd.P = o.P; nd.Np = o.nModes;
    const int P = o.P, Np = nd.Np;
    for (int i = 0; i <= P; ++i) for (int j = 0; j <= P - i; ++j) for (int k = 0; k <= P - i - j; ++k)
        nd.bary.push_back({i, j, k});
    std::vector<double> V((size_t)Np * Np);
    for (int nn = 0; nn < Np; ++nn) {
        const double r = nd.bary[nn][0] / (double)P, s = nd.bary[nn][1] / (double)P, t = nd.bary[nn][2] / (double)P;
        for (int m = 0; m < Np; ++m)
            V[(size_t)nn * Np + m] = pkdEval(o.modes[m][0], o.modes[m][1], o.modes[m][2], r, s, t);
    }
    std::vector<double> Vi;
    invertDense(V, Np, Vi);
    nd.V.assign(V.begin(), V.end());
    nd.Vinv.assign(Vi.begin(), Vi.end());
    return nd;
}

}  // namespace fem
}  // namespace mars
#endif
