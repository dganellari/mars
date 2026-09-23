// Host reference shared by the GPU gate (run_hl_mma.cpp) and the CPU warp
// emulator gate (emu_hl_mma.cpp): the Knaus Alg-2 apply for one element, and the
// per-element data hash both harnesses feed the kernel.
#pragma once

#include <vector>

// CPU port of the Knaus Alg-2 apply (applyHoCvfemElement) for one element.
static void oracle(int p, const double* u, const double* Bt, const double* Dt,
                   const double* Dm, const double* W, const double* G,
                   double* y)
{
    const int n = p + 1, nn = n * n;
    auto idx = [&](int d, int nrm, int t1, int t2) {
        if (d == 0) return nrm * nn + t1 * n + t2;
        if (d == 1) return t1 * nn + nrm * n + t2;
        return t1 * nn + t2 * n + nrm;
    };
    for (int i = 0; i < n * nn; ++i) y[i] = 0.0;
    std::vector<double> interp(nn), deriv(nn), flux(nn), tmp(nn), intf(nn);
    for (int d = 0; d < 3; ++d)
        for (int l = 0; l < p; ++l) {
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double bi = 0, di = 0;
                    for (int q = 0; q < n; ++q) {
                        bi += Bt[l * n + q] * u[idx(d, q, s, r)];
                        di += Dt[l * n + q] * u[idx(d, q, s, r)];
                    }
                    interp[s * n + r] = bi; deriv[s * n + r] = di;
                }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double dt2 = 0, dt1 = 0;
                    for (int q = 0; q < n; ++q) dt2 += Dm[r * n + q] * interp[s * n + q];
                    for (int q = 0; q < n; ++q) dt1 += Dm[s * n + q] * interp[q * n + r];
                    const double* g = G + (((d * p + l) * n + s) * n + r) * 3;
                    flux[s * n + r] = g[2] * deriv[s * n + r] + g[0] * dt2 + g[1] * dt1;
                }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) { double v = 0;
                    for (int q = 0; q < n; ++q) v += W[r * n + q] * flux[s * n + q];
                    tmp[s * n + r] = v; }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) { double v = 0;
                    for (int q = 0; q < n; ++q) v += W[s * n + q] * tmp[q * n + r];
                    intf[s * n + r] = v; }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    y[idx(d, l, s, r)]     -= intf[s * n + r];
                    y[idx(d, l + 1, s, r)] += intf[s * n + r];
                }
        }
}

// Per-element U and metric are a deterministic hash of (element, slot, salt), so
// any element can be regenerated on the host without a full mirror. Data stays
// DISTINCT per element, so a gate can catch a kernel that reads the wrong one.
static double hashVal(long long e, size_t i, unsigned salt)
{
    unsigned h = (unsigned)((unsigned long long)e * 2654435761ull +
                            (unsigned long long)i * 40503ull + salt);
    h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
    return 2.0 * (h / 4294967296.0) - 1.0;
}
