#pragma once

// Stage A of the high-order matrix-free CVFEM (Knaus, SAND2022-3366J): the 1D
// reference operators for a tensor-product hex of order p. Host-side, arbitrary
// p. These are the building blocks of the sum-factorized stiffness apply
// (mars_cvfem_ho_apply.hpp, next stage) via Kronecker products.
//
// Node sets (reference interval [-1,1]):
//   zeta : (p+1) Gauss-Lobatto-Legendre "solution" nodes  (the DOFs)
//   xi   : p Gauss-Legendre quadrature points (the subcontrol-surface points)
// Lagrange basis l_j is through the (p+1) GLL nodes zeta.
//
// Operators (Knaus Eq 6-8):
//   Btil[i][j] = l_j(xi_i)      interp GLL nodes -> Gauss/SCS points   (p x (p+1))
//   Dtil[i][j] = l'_j(xi_i)     derivative GLL -> Gauss/SCS points     (p x (p+1))
//   D[i][j]    = l'_j(zeta_i)   GLL differentiation matrix             ((p+1)^2)
//   W          integration operator over the (p+1) subcontrol intervals
//              defined by W^{-1}[i][j] = h_j(zeta_i), h = histopolation edge
//              functions through the padded Gauss points {-1,xi,+1}     ((p+1)^2)
//   Deltatil   (p+1) x p subcontrol-face incidence (bidiagonal -1/+1)   (Eq 7)
//
// p=1 reduces to standard linear CVFEM: zeta={-1,1}, xi={0},
//   Btil=[1/2,1/2], Dtil=[-1/2,1/2], Deltatil=[-1;+1].

#include "mars_spectral_basis.hpp"
#include <vector>
#include <cmath>

namespace mars {
namespace fem {

// p Gauss-Legendre nodes + weights on [-1,1] (roots of Legendre P_p).
inline void gaussLegendre(int p, std::vector<double>& x, std::vector<double>& w)
{
    x.assign(p, 0.0); w.assign(p, 0.0);
    for (int i = 0; i < p; ++i)
    {
        double xi = std::cos(M_PI * (i + 0.75) / (p + 0.5));   // Chebyshev guess
        for (int it = 0; it < 100; ++it)
        {
            double L, dL; legendreP(p, xi, L, dL);
            double dx = -L / dL;
            xi += dx;
            if (std::abs(dx) < 1e-15) break;
        }
        x[i] = xi;
        double L, dL; legendreP(p, xi, L, dL);
        w[i] = 2.0 / ((1.0 - xi * xi) * dL * dL);
    }
    // sort ascending (Newton from cos guesses lands them descending)
    for (int i = 0; i < p; ++i)
        for (int j = i + 1; j < p; ++j)
            if (x[j] < x[i]) { std::swap(x[i], x[j]); std::swap(w[i], w[j]); }
}

// Lagrange basis through `t` (m nodes), evaluated at x: value[k], deriv[k].
// Handles x exactly at a node (nodal-derivative formula).
inline void lagrangeValDeriv(const std::vector<double>& t, double x,
                             std::vector<double>& val, std::vector<double>& der)
{
    int m = (int)t.size();
    val.assign(m, 0.0); der.assign(m, 0.0);
    int atNode = -1;
    for (int j = 0; j < m; ++j) if (std::abs(x - t[j]) < 1e-13) { atNode = j; break; }

    if (atNode < 0) {
        for (int k = 0; k < m; ++k) {
            double prod = 1.0, s = 0.0;
            for (int i = 0; i < m; ++i) if (i != k) { prod *= (x - t[i]) / (t[k] - t[i]); s += 1.0 / (x - t[i]); }
            val[k] = prod;
            der[k] = prod * s;
        }
        return;
    }
    // x == t[atNode]: value is delta; derivative via nodal formulas.
    int mnode = atNode;
    for (int k = 0; k < m; ++k) val[k] = (k == mnode) ? 1.0 : 0.0;
    for (int k = 0; k < m; ++k) {
        if (k != mnode) {
            double prod = 1.0 / (t[k] - t[mnode]);
            for (int i = 0; i < m; ++i) if (i != k && i != mnode) prod *= (t[mnode] - t[i]) / (t[k] - t[i]);
            der[k] = prod;
        } else {
            double s = 0.0;
            for (int i = 0; i < m; ++i) if (i != mnode) s += 1.0 / (t[mnode] - t[i]);
            der[k] = s;
        }
    }
}

// Invert a small square matrix A[n*n] (row-major) via Gauss-Jordan -> Ainv.
inline void invertSmall(const std::vector<double>& A, int n, std::vector<double>& Ainv)
{
    std::vector<double> M(A);
    Ainv.assign(n * n, 0.0);
    for (int i = 0; i < n; ++i) Ainv[i * n + i] = 1.0;
    for (int c = 0; c < n; ++c) {
        int piv = c; double best = std::abs(M[c * n + c]);
        for (int r = c + 1; r < n; ++r) if (std::abs(M[r * n + c]) > best) { best = std::abs(M[r * n + c]); piv = r; }
        for (int j = 0; j < n; ++j) { std::swap(M[c*n+j], M[piv*n+j]); std::swap(Ainv[c*n+j], Ainv[piv*n+j]); }
        double d = M[c * n + c];
        for (int j = 0; j < n; ++j) { M[c*n+j] /= d; Ainv[c*n+j] /= d; }
        for (int r = 0; r < n; ++r) if (r != c) {
            double f = M[r * n + c];
            for (int j = 0; j < n; ++j) { M[r*n+j] -= f * M[c*n+j]; Ainv[r*n+j] -= f * Ainv[c*n+j]; }
        }
    }
}

struct HoCvfemOperators {
    int p = 0;
    std::vector<double> zeta;       // (p+1) GLL nodes
    std::vector<double> xi;         // p Gauss points
    std::vector<double> Btil;       // p x (p+1)
    std::vector<double> Dtil;       // p x (p+1)
    std::vector<double> D;          // (p+1)^2 GLL derivative
    std::vector<double> W;          // (p+1)^2 integration
    std::vector<double> Deltatil;   // (p+1) x p
};

inline HoCvfemOperators buildHoCvfemOperators(int P)
{
    HoCvfemOperators op; op.p = P;
    const int n = P + 1;

    std::vector<double> w_unused, Dgll;
    gllBasis(P, op.zeta, w_unused, Dgll);
    op.D = Dgll;

    std::vector<double> wg;
    gaussLegendre(P, op.xi, wg);

    // Btil, Dtil: Lagrange (through zeta) interp/deriv at the Gauss points xi.
    op.Btil.assign(P * n, 0.0);
    op.Dtil.assign(P * n, 0.0);
    for (int i = 0; i < P; ++i) {
        std::vector<double> v, d;
        lagrangeValDeriv(op.zeta, op.xi[i], v, d);
        for (int j = 0; j < n; ++j) { op.Btil[i * n + j] = v[j]; op.Dtil[i * n + j] = d[j]; }
    }

    // Deltatil: (p+1) x p bidiagonal subcontrol-face incidence (Eq 7).
    op.Deltatil.assign(n * P, 0.0);
    for (int i = 0; i < n; ++i) {
        if (i < P)     op.Deltatil[i * P + i]       = -1.0;
        if (i - 1 >= 0) op.Deltatil[i * P + (i - 1)] = +1.0;
    }

    // W: padded Gauss points {-1, xi_0..xi_{p-1}, +1}; edge (histopolation)
    // functions e_i(x) = -sum_{k<=i} d'_k(x) over the padded-node Lagrange basis
    // d_k. W^{-1}[i][j] = e_j(zeta_i); W = inverse. (Knaus Eq 5-6.)
    std::vector<double> pad; pad.reserve(P + 2);
    pad.push_back(-1.0);
    for (double x : op.xi) pad.push_back(x);
    pad.push_back(1.0);

    std::vector<double> Winv(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        std::vector<double> v, d;
        lagrangeValDeriv(pad, op.zeta[i], v, d);   // d = d'_k(zeta_i)
        // edge function e_j = -sum_{k=0..j} d'_k ; j = 0..p
        double cum = 0.0;
        for (int j = 0; j < n; ++j) { cum += d[j]; Winv[i * n + j] = -cum; }
    }
    invertSmall(Winv, n, op.W);

    return op;
}

} // namespace fem
} // namespace mars
