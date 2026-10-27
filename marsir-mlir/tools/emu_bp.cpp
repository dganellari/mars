// CPU gate of a Galerkin operator generated from its spec (test/bp_gpu_pipeline.py
// --emulate): the emulated kernel on random u, 1D matrices and geometry, against the
// operator written out as direct sums over all nodes and Gauss points (no sum
// factorization, no shared code with the compiler):
//   uq[a,b,c]  = sum_ijk B[a,i] B[b,j] B[c,k] u[i,j,k]
//   du_d       = D along axis d of uq
//   outputs    = the spec's pointwise action (BP1: v = g0 u; BP3: dv = G du, G symmetric)
//   w          = v + sum_d Dt along axis d of dv_d
//   y[i,j,k]   = sum_abc Bt[i,a] Bt[j,b] Bt[k,c] w[a,b,c]
// Defines: BP_NAME (kernel), BP_SPEC_BP1 | BP_SPEC_BP3, BP_Q, BP_K (geometry components),
// BP_GRAD (D and Dt are arguments), MIR_EMU_REAL (float for the f32 kernel).
// BP_TRILINEAR (geometry = trilinear): the kernel takes gq (random points in [-1, 1]), wq
// (random positive weights) and the corners X (a perturbed unit cube per element) instead
// of G; the oracle forms G from them by direct sums over the corners (w det J, or the
// symmetric w det J J^-1 J^-T) and checks the same action.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <thread>
#include <vector>

#include "emu_warp_runtime.h"

#ifndef MIR_EMU_REAL
#define MIR_EMU_REAL double
#endif
using Real = MIR_EMU_REAL;
#define BP_CAT_(a, b) a##b
#define BP_CAT(a, b) BP_CAT_(a, b)
#define BP_STR_(x) #x
#define BP_STR(x) BP_STR_(x)

using D1 = Desc<Real, 1>;
using D2 = Desc<Real, 2>;
using D4 = Desc<Real, 4>;
using D5 = Desc<Real, 5>;
extern "C" {
void BP_CAT(_mlir_ciface_, BP_NAME)(D4*, D2*, D2*,
#ifdef BP_GRAD
                                    D2*, D2*,
#endif
#ifdef BP_TRILINEAR
                                    D1*, D1*, D5*,
#elif BP_K > 0
                                    D5*,
#endif
                                    D4*);
int64_t BP_CAT(BP_NAME, __threads)();
}

int main(int argc, char** argv)
{
    const long long E = argc > 1 ? atoll(argv[1]) : 3;
    const int p = argc > 2 ? atoi(argv[2]) : 3, n = p + 1, Q = BP_Q, K = BP_K;
    const long long n3 = (long long)n * n * n, q3 = (long long)Q * Q * Q;
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    auto rand = [&](long long size) {
        std::vector<double> v(size);
        for (auto& x : v) x = val(rng);
        return v;
    };
    std::vector<double> u = rand(E * n3), B = rand((long long)Q * n), D = rand((long long)Q * Q);
    std::vector<double> G = rand(E * (K > 0 ? K : 1) * q3), Bt(n * Q), Dt(Q * Q);
    for (int i = 0; i < Q; ++i)
        for (int j = 0; j < n; ++j) Bt[j * Q + i] = B[i * n + j];
    for (int i = 0; i < Q; ++i)
        for (int j = 0; j < Q; ++j) Dt[j * Q + i] = D[i * Q + j];
#ifdef BP_TRILINEAR
    // Random points, positive weights, corners of a perturbed unit cube.
    std::vector<double> gq = rand(Q), wq1(Q), X(E * 24);
    for (auto& w : wq1) w = 1.0 + 0.5 * val(rng);
    for (long long e = 0; e < E; ++e)
        for (int c = 0; c < 3; ++c)
            for (int corner = 0; corner < 8; ++corner) {
                const int bit = (corner >> (2 - c)) & 1;   // corner = (i, j, k) bits, X[c][i][j][k]
                X[e * 24 + c * 8 + corner] = bit + 0.15 * val(rng);
            }
    // G from the corners: J[a][b] = sum_ijk X[a][i][j][k] M0(i) M1(j) M2(k) with M = the
    // linear basis (1 -+ xi) / 2 at the point's coordinate, or its derivative -+ 1/2 for
    // axis b.
    for (long long e = 0; e < E; ++e)
        for (long long pt = 0; pt < q3; ++pt) {
            const int qa[3] = {(int)(pt / (Q * Q)), (int)((pt / Q) % Q), (int)(pt % Q)};
            double J[3][3];
            for (int a = 0; a < 3; ++a)
                for (int b2 = 0; b2 < 3; ++b2) {
                    double s = 0;
                    for (int corner = 0; corner < 8; ++corner) {
                        const int bits[3] = {corner >> 2, (corner >> 1) & 1, corner & 1};
                        double f = X[e * 24 + a * 8 + corner];
                        for (int ax = 0; ax < 3; ++ax) {
                            const double sign = bits[ax] ? 0.5 : -0.5;
                            f *= ax == b2 ? sign : 0.5 + sign * gq[qa[ax]];
                        }
                        s += f;
                    }
                    J[a][b2] = s;
                }
            const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                               J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                               J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
            double inv[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) {
                    const int r1 = (c + 1) % 3, r2 = (c + 2) % 3, c1 = (r + 1) % 3, c2 = (r + 2) % 3;
                    inv[r][c] = (J[r1][c1] * J[r2][c2] - J[r1][c2] * J[r2][c1]) / det;
                }
            const double wdet = wq1[qa[0]] * wq1[qa[1]] * wq1[qa[2]] * det;
            if (K == 1) {
                G[e * q3 + pt] = wdet;
                continue;
            }
            const int sym[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
            for (int k = 0; k < 6; ++k) {
                const int d = sym[k][0], f = sym[k][1];
                G[(e * K + k) * q3 + pt] = wdet * (inv[d][0] * inv[f][0] + inv[d][1] * inv[f][1] + inv[d][2] * inv[f][2]);
            }
        }
#endif

    // The kernel's inputs in its precision; the oracle sees the same rounded values.
    auto narrow = [](std::vector<double>& v) {
        std::vector<Real> r(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            r[i] = (Real)v[i];
            v[i] = (double)r[i];
        }
        return r;
    };
    std::vector<Real> ru = narrow(u), rB = narrow(B), rBt = narrow(Bt), rD = narrow(D), rDt = narrow(Dt),
                      rG = narrow(G), Y(E * n3, (Real)std::nan(""));
#ifdef BP_TRILINEAR
    std::vector<Real> rgq = narrow(gq), rwq = narrow(wq1), rX = narrow(X);
    D1 dgq{rgq.data(), rgq.data(), 0, {Q}, {1}}, dwq{rwq.data(), rwq.data(), 0, {Q}, {1}};
    D5 dX{rX.data(), rX.data(), 0, {E, 3, 2, 2, 2}, {24, 8, 4, 2, 1}};
#endif

    D4 du{ru.data(), ru.data(), 0, {E, n, n, n}, {n3, (int64_t)n * n, n, 1}};
    D4 dy{Y.data(), Y.data(), 0, {E, n, n, n}, {n3, (int64_t)n * n, n, 1}};
    D2 dB{rB.data(), rB.data(), 0, {Q, n}, {n, 1}}, dBt{rBt.data(), rBt.data(), 0, {n, Q}, {Q, 1}};
    D2 dD{rD.data(), rD.data(), 0, {Q, Q}, {Q, 1}}, dDt{rDt.data(), rDt.data(), 0, {Q, Q}, {Q, 1}};
    D5 dG{rG.data(), rG.data(), 0, {E, K, Q, Q, Q}, {K * q3, q3, (int64_t)Q * Q, Q, 1}};
    (void)dD, (void)dDt, (void)dG;

    Warp w;
    w.bar.threads = (int)BP_CAT(BP_NAME, __threads)();
    warp = &w;
    for (long long b = 0; b < E; ++b) {   // a block per element covers any launch shape
        std::vector<std::thread> lanes;
        for (int L = 0; L < w.bar.threads; ++L)
            lanes.emplace_back([&, L, b] {
                tlLane = L;
                tlBlock = b;
                tlGrid = E;
                BP_CAT(_mlir_ciface_, BP_NAME)(&du, &dB, &dBt,
#ifdef BP_GRAD
                                               &dD, &dDt,
#endif
#ifdef BP_TRILINEAR
                                               &dgq, &dwq, &dX,
#elif BP_K > 0
                                               &dG,
#endif
                                               &dy);
            });
        for (auto& t : lanes) t.join();
    }

    double err = 0, scale = 0;
    for (long long e = 0; e < E; ++e) {
        const double* ue = &u[e * n3];
        std::vector<double> uq(q3, 0.0), wq(q3, 0.0);
        for (int a = 0; a < Q; ++a)
            for (int b2 = 0; b2 < Q; ++b2)
                for (int c = 0; c < Q; ++c) {
                    double s = 0;
                    for (int i = 0; i < n; ++i)
                        for (int j = 0; j < n; ++j)
                            for (int k = 0; k < n; ++k)
                                s += B[a * n + i] * B[b2 * n + j] * B[c * n + k] * ue[(i * n + j) * n + k];
                    uq[(a * Q + b2) * Q + c] = s;
                }
        auto at = [&](int a, int b2, int c) { return (a * Q + b2) * Q + c; };
        auto g = [&](int comp, long long pt) { return G[(e * K + comp) * q3 + pt]; };
#if defined(BP_SPEC_BP1)
        for (long long pt = 0; pt < q3; ++pt) wq[pt] = g(0, pt) * uq[pt];
#elif defined(BP_SPEC_BP3)
        std::vector<double> d[3], dv[3];
        for (int r = 0; r < 3; ++r) d[r].assign(q3, 0.0), dv[r].assign(q3, 0.0);
        for (int a = 0; a < Q; ++a)
            for (int b2 = 0; b2 < Q; ++b2)
                for (int c = 0; c < Q; ++c)
                    for (int m = 0; m < Q; ++m) {
                        d[0][at(a, b2, c)] += D[a * Q + m] * uq[at(m, b2, c)];
                        d[1][at(a, b2, c)] += D[b2 * Q + m] * uq[at(a, m, c)];
                        d[2][at(a, b2, c)] += D[c * Q + m] * uq[at(a, b2, m)];
                    }
        const int sym[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        for (long long pt = 0; pt < q3; ++pt)
            for (int r = 0; r < 3; ++r)
                for (int s2 = 0; s2 < 3; ++s2) dv[r][pt] += g(sym[r][s2], pt) * d[s2][pt];
        for (int a = 0; a < Q; ++a)
            for (int b2 = 0; b2 < Q; ++b2)
                for (int c = 0; c < Q; ++c)
                    for (int m = 0; m < Q; ++m)
                        wq[at(a, b2, c)] += Dt[a * Q + m] * dv[0][at(m, b2, c)] + Dt[b2 * Q + m] * dv[1][at(a, m, c)] +
                                            Dt[c * Q + m] * dv[2][at(a, b2, m)];
#endif
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    double s = 0;
                    for (int a = 0; a < Q; ++a)
                        for (int b2 = 0; b2 < Q; ++b2)
                            for (int c = 0; c < Q; ++c)
                                s += Bt[i * Q + a] * Bt[j * Q + b2] * Bt[k * Q + c] * wq[at(a, b2, c)];
                    const long long t = e * n3 + (i * n + j) * n + k;
                    err = fmax(err, fabs((double)Y[t] - s));
                    if (std::isnan((double)Y[t])) err = INFINITY;
                    scale = fmax(scale, fabs(s));
                }
    }
    const bool single = sizeof(Real) == sizeof(float);
    const bool ok = err <= (single ? 1e-5 : 1e-12) * scale;
    printf("GALERKIN WARP EMULATION (%s, p=%d, q=%d, E=%lld, %d threads per block, %s): max|err| = %.3e of max|y| "
           "%.3e  %s\n",
           BP_STR(BP_NAME), p, Q, E, w.bar.threads, single ? "f32" : "f64", err, scale, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
