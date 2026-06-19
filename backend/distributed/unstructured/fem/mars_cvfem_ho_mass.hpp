#pragma once

// Consistent CVFEM control-volume mass/load apply (Knaus Alg 1 volume term).
// (M f)_node = integral over the node's sub-control volume of the nodal
// interpolant of f = detJ * (W (x) W (x) W) f, applied sum-factorized (W along
// each axis). Used to form the consistent load b = M f for the manufactured-
// solution convergence study. Affine cube: detJ = (h/2)^3 is constant.

#include "mars_cvfem_ho_basis.hpp"
#include <vector>

namespace mars {
namespace fem {

inline void applyMassCube(const HoCvfemOperators& op, double detJ,
                          const std::vector<double>& f, std::vector<double>& y)
{
    const int p = op.p, n = p + 1, n3 = n * n * n;
    y.assign(n3, 0.0);
    auto id = [&](int a, int b, int c) { return (a * n + b) * n + c; };
    std::vector<double> t1(n3, 0.0), t2(n3, 0.0);
    // contract W along axis 0 (output i <- input a)
    for (int i = 0; i < n; ++i) for (int b = 0; b < n; ++b) for (int c = 0; c < n; ++c) {
        double s = 0; for (int a = 0; a < n; ++a) s += op.W[i * n + a] * f[id(a, b, c)];
        t1[id(i, b, c)] = s;
    }
    // axis 1 (j <- b)
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int c = 0; c < n; ++c) {
        double s = 0; for (int b = 0; b < n; ++b) s += op.W[j * n + b] * t1[id(i, b, c)];
        t2[id(i, j, c)] = s;
    }
    // axis 2 (k <- c) + the physical volume factor
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
        double s = 0; for (int c = 0; c < n; ++c) s += op.W[k * n + c] * t2[id(i, j, c)];
        y[id(i, j, k)] = detJ * s;
    }
}

} // namespace fem
} // namespace mars
