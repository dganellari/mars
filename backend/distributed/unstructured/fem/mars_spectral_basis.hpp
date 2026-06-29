#pragma once

// Gauss-Lobatto-Legendre (GLL) nodal basis for tensor-product spectral hex
// elements: nodes, quadrature weights, and the 1D differentiation matrix.
// Host-side, arbitrary order. These feed the matrix-free Laplacian kernel
// (mars_ho_laplacian.hpp): D is the per-direction contraction operator, the
// weights build the per-element metric, the nodes place the high-order DOFs.
//
// Order P  ->  n = P+1 nodes on the reference interval [-1, 1].
// GLL nodes are the endpoints +-1 plus the roots of L'_P (extrema of the
// degree-P Legendre polynomial). Newton iteration; the update simplifies to
//   x <- x + ((1-x^2) L'_P) / (P(P+1) L_P)
// because d/dx[(1-x^2)L'_P] = -P(P+1) L_P via the Legendre ODE.

#include <vector>
#include <cmath>

namespace mars {
namespace fem {

// Legendre L_P and derivative L'_P at x, via the standard recurrence.
inline void legendreP(int P, double x, double& L, double& dL)
{
    if (P == 0) { L = 1.0; dL = 0.0; return; }
    if (P == 1) { L = x;   dL = 1.0; return; }
    double Lkm1 = 1.0, Lk = x;
    for (int k = 1; k < P; ++k)
    {
        double Lkp1 = ((2 * k + 1) * x * Lk - k * Lkm1) / (k + 1);
        Lkm1 = Lk;
        Lk   = Lkp1;
    }
    L = Lk;                                   // L_P
    // L'_P = P (x L_P - L_{P-1}) / (x^2 - 1); interior nodes avoid x=+-1.
    double denom = x * x - 1.0;
    dL = (std::abs(denom) < 1e-14) ? 0.0 : P * (x * Lk - Lkm1) / denom;
}

// Build GLL nodes x[n], weights w[n], and differentiation matrix D[n*n]
// (row-major, D[i*n+j]) for order P. n = P+1.
inline void gllBasis(int P,
                     std::vector<double>& x,
                     std::vector<double>& w,
                     std::vector<double>& D)
{
    int n = P + 1;
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    D.assign(n * n, 0.0);

    if (P == 0) { x[0] = 0.0; w[0] = 2.0; return; }

    // Endpoints are exact GLL nodes.
    x[0]     = -1.0;
    x[P]     =  1.0;
    // Interior nodes: Newton from the Chebyshev-Gauss-Lobatto guess.
    for (int j = 1; j < P; ++j)
    {
        double xj = -std::cos(M_PI * j / P);
        for (int it = 0; it < 100; ++it)
        {
            double L, dL;
            legendreP(P, xj, L, dL);
            double dx = ((1.0 - xj * xj) * dL) / (double(P) * (P + 1) * L);
            xj += dx;
            if (std::abs(dx) < 1e-15) break;
        }
        x[j] = xj;
    }

    // Weights: w_j = 2 / (P(P+1) [L_P(x_j)]^2), valid at all nodes incl. ends.
    for (int j = 0; j < n; ++j)
    {
        double L, dL;
        legendreP(P, x[j], L, dL);
        w[j] = 2.0 / (double(P) * (P + 1) * L * L);
    }

    // Differentiation matrix (Canuto et al.):
    //   D[i][j] = L_P(x_i) / (L_P(x_j) (x_i - x_j))   i != j
    //   D[0][0]   = -P(P+1)/4
    //   D[P][P]   =  P(P+1)/4
    //   D[i][i]   = 0 otherwise
    std::vector<double> LP(n);
    for (int i = 0; i < n; ++i) { double dl; legendreP(P, x[i], LP[i], dl); }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
                D[i * n + j] = LP[i] / (LP[j] * (x[i] - x[j]));
            else if (i == 0)
                D[i * n + j] = -0.25 * P * (P + 1);
            else if (i == P)
                D[i * n + j] =  0.25 * P * (P + 1);
            else
                D[i * n + j] = 0.0;
        }
}

} // namespace fem
} // namespace mars
