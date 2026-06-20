// Stage C validation (host-only): the element diffusion stiffness K_e built by
// the sum-factorization apply must satisfy the structural properties of a
// Laplacian on an orthogonal cube element:
//   nullspace  K_e * 1 = 0          (constant in the kernel)
//   symmetry   K_e[I][J] = K_e[J][I] (orthogonal cube -> symmetric)
//   PSD        u^T K_e u >= 0
// Scaling-invariant, so coeff=1. Absolute scaling pinned later vs --kernel=tensor.

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"

#include <cstdio>
#include <cmath>
#include <random>
#include <vector>

using namespace mars::fem;

static bool checkElement(int P)
{
    auto op = buildHoCvfemOperators(P);
    const int n = P + 1, n3 = n * n * n;

    // Build K_e by applying to unit vectors: column J = K_e e_J.
    std::vector<double> K((size_t)n3 * n3, 0.0);
    for (int J = 0; J < n3; ++J) {
        std::vector<double> e(n3, 0.0), y; e[J] = 1.0;
        applyHoCvfemElementCube(op, 1.0, e, y);
        for (int I = 0; I < n3; ++I) K[(size_t)I * n3 + J] = y[I];
    }

    // nullspace: K * 1
    double nullMax = 0;
    for (int I = 0; I < n3; ++I) { double s = 0; for (int J = 0; J < n3; ++J) s += K[(size_t)I*n3+J]; nullMax = std::max(nullMax, std::abs(s)); }

    // symmetry
    double symMax = 0, scale = 0;
    for (size_t e = 0; e < (size_t)n3*n3; ++e) scale = std::max(scale, std::abs(K[e]));
    for (int I = 0; I < n3; ++I) for (int J = 0; J < n3; ++J)
        symMax = std::max(symMax, std::abs(K[(size_t)I*n3+J] - K[(size_t)J*n3+I]));

    // PSD: random u, u^T K u >= 0
    std::mt19937 rng(3); std::uniform_real_distribution<double> uni(-1, 1);
    double minQ = 1e300;
    for (int t = 0; t < 20; ++t) {
        std::vector<double> u(n3); for (auto& z : u) z = uni(rng);
        double q = 0; for (int I = 0; I < n3; ++I) { double Ku = 0; for (int J = 0; J < n3; ++J) Ku += K[(size_t)I*n3+J]*u[J]; q += u[I]*Ku; }
        minQ = std::min(minQ, q);
    }

    bool ok = (nullMax < 1e-10*std::max(scale,1.0)) && (symMax < 1e-10*std::max(scale,1.0)) && (minQ > -1e-9*std::max(scale,1.0));
    printf("p=%d (n3=%d) | null=%.2e | sym=%.2e | min uTKu=%.3e | %s\n",
           P, n3, nullMax, symMax, minQ, ok ? "PASS" : "FAIL");
    return ok;
}

int main()
{
    bool ok = true;
    ok &= checkElement(1);
    ok &= checkElement(2);
    ok &= checkElement(4);
    printf("Stage C (HO-CVFEM element stiffness, structural): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
