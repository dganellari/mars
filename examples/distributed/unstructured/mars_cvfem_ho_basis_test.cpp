// Stage A validation (host-only, runs anywhere): the high-order CVFEM reference
// operators must reproduce polynomials exactly up to degree p.
//   Btil : interp GLL->Gauss        Btil * f(zeta) == f(xi)        for deg<=p
//   Dtil : derivative GLL->Gauss    Dtil * f(zeta) == f'(xi)       for deg<=p
//   W    : subcontrol integration   sum_m (W f(zeta))_m == int_-1^1 f  for deg<=p
// plus p=1 sanity (Btil=[1/2,1/2], Dtil=[-1/2,1/2], Deltatil=[-1;+1]).

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"

#include <cstdio>
#include <cmath>
#include <vector>

using namespace mars::fem;

static double exactIntMonomial(int k) { return (k % 2 == 0) ? 2.0 / (k + 1) : 0.0; }

static bool checkOrder(int P)
{
    auto op = buildHoCvfemOperators(P);
    const int n = P + 1;
    double maxB = 0, maxD = 0, maxW = 0;

    for (int k = 0; k <= P; ++k)               // monomial x^k, k = 0..p
    {
        std::vector<double> f(n);
        for (int i = 0; i < n; ++i) f[i] = std::pow(op.zeta[i], k);

        // Btil * f  vs  xi^k ;  Dtil * f  vs  k*xi^(k-1)
        for (int i = 0; i < P; ++i) {
            double b = 0, d = 0;
            for (int j = 0; j < n; ++j) { b += op.Btil[i*n+j]*f[j]; d += op.Dtil[i*n+j]*f[j]; }
            maxB = std::max(maxB, std::abs(b - std::pow(op.xi[i], k)));
            double dref = (k == 0) ? 0.0 : k * std::pow(op.xi[i], k - 1);
            maxD = std::max(maxD, std::abs(d - dref));
        }

        // sum_m (W f)_m  vs  int_-1^1 x^k
        double s = 0;
        for (int m = 0; m < n; ++m) for (int j = 0; j < n; ++j) s += op.W[m*n+j]*f[j];
        maxW = std::max(maxW, std::abs(s - exactIntMonomial(k)));
    }

    bool ok = (maxB < 1e-10) && (maxD < 1e-10) && (maxW < 1e-10);
    printf("p=%d | Btil=%.1e Dtil=%.1e W=%.1e | %s\n", P, maxB, maxD, maxW, ok ? "PASS" : "FAIL");
    return ok;
}

int main()
{
    bool ok = true;

    // p=1 explicit sanity
    {
        auto op = buildHoCvfemOperators(1);
        bool s = std::abs(op.Btil[0]-0.5)<1e-12 && std::abs(op.Btil[1]-0.5)<1e-12
              && std::abs(op.Dtil[0]+0.5)<1e-12 && std::abs(op.Dtil[1]-0.5)<1e-12
              && std::abs(op.Deltatil[0]+1.0)<1e-12 && std::abs(op.Deltatil[1]-1.0)<1e-12;
        printf("p=1 sanity (Btil=[.5,.5] Dtil=[-.5,.5] Deltatil=[-1;1]): %s\n", s ? "PASS" : "FAIL");
        ok &= s;
    }

    ok &= checkOrder(1);
    ok &= checkOrder(2);
    ok &= checkOrder(4);
    ok &= checkOrder(7);

    printf("Stage A (HO-CVFEM reference operators): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
