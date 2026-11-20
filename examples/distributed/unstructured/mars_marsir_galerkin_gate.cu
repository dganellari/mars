// The generated Galerkin BP3 operator (MARSIR, bp3_apply_p{p}) on one GPU, against direct
// quadrature sums on the host: the same sums as mars_marsir_bp's oracle, on every element.
// The kernel file is PTX on NVIDIA and a gfx942 code object on AMD, so this checks a new
// backend without the rest of the cell-wise solver. Inputs are random: element vectors, and
// a symmetric geometry tensor with a dominant diagonal (any symmetric G defines the same
// contractions, so it need not come from a mesh). Also checks A 1 = 0 and the symmetry
// <v, A u> = <u, A v>. Then times the apply: median of `reps` warm launches, host wall time
// around a device synchronization (the same measure on CUDA and HIP). For large runs only the
// first 512 elements are compared on the host.
// Usage: mars_marsir_galerkin_gate <kernel file> [p = 3] [elements = 512] [reps = 0]
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_galerkin.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using namespace mars;

namespace {

struct Tables {
    std::vector<double> B, Bt, D, Dt;
};

// B: GLL basis at the q = p + 2 Gauss points (q x n); D: Gauss-point derivative (q x q).
Tables tables(int p)
{
    const int n = p + 1, Q = p + 2;
    std::vector<double> zeta, wgll, Dgll, gq, wq, val, der;
    fem::gllBasis(p, zeta, wgll, Dgll);
    fem::gaussLegendre(Q, gq, wq);
    Tables t;
    t.B.resize(Q * n), t.Bt.resize(n * Q), t.D.resize(Q * Q), t.Dt.resize(Q * Q);
    for (int a = 0; a < Q; ++a) {
        fem::lagrangeValDeriv(zeta, gq[a], val, der);
        for (int i = 0; i < n; ++i) t.B[a * n + i] = t.Bt[i * Q + a] = val[i];
        fem::lagrangeValDeriv(gq, gq[a], val, der);
        for (int k = 0; k < Q; ++k) t.D[a * Q + k] = t.Dt[k * Q + a] = der[k];
    }
    return t;
}

// y_e = A_e u_e by direct quadrature sums (mars_marsir_bp.cu's oracle, BP3 branch).
void reference(const Tables& t, int p, const double* u, const double* g, double* y)
{
    const int n = p + 1, Q = p + 2, n3 = n * n * n, q3 = Q * Q * Q;
    const int sym[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    const std::vector<double>& B = t.B;
    const std::vector<double>& D = t.D;
    std::vector<double> uq(q3, 0.0), w(q3, 0.0), d[3], dv[3];
    for (int pt = 0; pt < q3; ++pt) {
        const int a = pt / (Q * Q), b = (pt / Q) % Q, c = pt % Q;
        for (int node = 0; node < n3; ++node)
            uq[pt] += B[a * n + node / (n * n)] * B[b * n + (node / n) % n] * B[c * n + node % n] * u[node];
    }
    for (int r = 0; r < 3; ++r) d[r].assign(q3, 0.0), dv[r].assign(q3, 0.0);
    for (int pt = 0; pt < q3; ++pt) {
        const int a = pt / (Q * Q), b = (pt / Q) % Q, c = pt % Q;
        for (int m = 0; m < Q; ++m) {
            d[0][pt] += D[a * Q + m] * uq[(m * Q + b) * Q + c];
            d[1][pt] += D[b * Q + m] * uq[(a * Q + m) * Q + c];
            d[2][pt] += D[c * Q + m] * uq[(a * Q + b) * Q + m];
        }
    }
    for (int r = 0; r < 3; ++r)
        for (int pt = 0; pt < q3; ++pt)
            for (int s = 0; s < 3; ++s) dv[r][pt] += g[sym[r][s] * q3 + pt] * d[s][pt];
    for (int pt = 0; pt < q3; ++pt) {
        const int a = pt / (Q * Q), b = (pt / Q) % Q, c = pt % Q;
        for (int m = 0; m < Q; ++m)
            w[pt] += D[m * Q + a] * dv[0][(m * Q + b) * Q + c] + D[m * Q + b] * dv[1][(a * Q + m) * Q + c] +
                     D[m * Q + c] * dv[2][(a * Q + b) * Q + m];
    }
    for (int node = 0; node < n3; ++node) {
        double sum = 0.0;
        for (int pt = 0; pt < q3; ++pt) {
            const int a = pt / (Q * Q), b = (pt / Q) % Q, c = pt % Q;
            sum += B[a * n + node / (n * n)] * B[b * n + (node / n) % n] * B[c * n + node % n] * w[pt];
        }
        y[node] = sum;
    }
}

#define GATE_CK(call)                                                                       \
    do {                                                                                    \
        cudaError_t e_ = (call);                                                            \
        if (e_ != cudaSuccess) {                                                            \
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_));     \
            std::exit(1);                                                                   \
        }                                                                                   \
    } while (0)

double* to_device(const std::vector<double>& h)
{
    double* d = nullptr;
    GATE_CK(cudaMalloc(&d, h.size() * sizeof(double)));
    GATE_CK(cudaMemcpy(d, h.data(), h.size() * sizeof(double), cudaMemcpyHostToDevice));
    return d;
}

std::vector<double> apply(const marsir::GalerkinPtx& op, const Tables& t, const double* dB, const double* dBt,
                          const double* dD, const double* dDt, const double* dG, const std::vector<double>& u,
                          long long E)
{
    double* du = to_device(u);
    double* dy = nullptr;
    GATE_CK(cudaMalloc(&dy, u.size() * sizeof(double)));
    op.apply(du, dB, dBt, dD, dDt, dG, dy, E);
    GATE_CK(cudaGetLastError());
    GATE_CK(cudaDeviceSynchronize());
    std::vector<double> y(u.size());
    GATE_CK(cudaMemcpy(y.data(), dy, y.size() * sizeof(double), cudaMemcpyDeviceToHost));
    GATE_CK(cudaFree(du));
    GATE_CK(cudaFree(dy));
    return y;
}

}   // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <kernel file> [p = 3] [elements = 512]\n", argv[0]);
        return 2;
    }
    const std::string file = argv[1];
    const int p = argc > 2 ? std::atoi(argv[2]) : 3;
    const long long E = argc > 3 ? std::atoll(argv[3]) : 512;
    const int reps = argc > 4 ? std::atoi(argv[4]) : 0;
    const long long Echeck = std::min<long long>(E, 512);
    const int n = p + 1, Q = p + 2, n3 = n * n * n, q3 = Q * Q * Q, K = 6;

    const Tables t = tables(p);
    std::mt19937_64 rng(1234);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::vector<double> u(E * n3), v(E * n3), ones(E * n3, 1.0), G(E * K * q3);
    for (auto& x : u) x = unif(rng);
    for (auto& x : v) x = unif(rng);
    for (long long e = 0; e < E; ++e)
        for (int k = 0; k < K; ++k)
            for (int pt = 0; pt < q3; ++pt)   // diagonal entries 0, 3, 5 dominate: symmetric positive definite
                G[(e * K + k) * q3 + pt] = (k == 0 || k == 3 || k == 5) ? 2.0 + unif(rng) : 0.3 * unif(rng);

    double* dB = to_device(t.B);
    double* dBt = to_device(t.Bt);
    double* dD = to_device(t.D);
    double* dDt = to_device(t.Dt);
    double* dG = to_device(G);
    marsir::GalerkinPtx op(file, "bp3_apply_p" + std::to_string(p), p, Q, true, K);

    const std::vector<double> Au = apply(op, t, dB, dBt, dD, dDt, dG, u, E);
    const std::vector<double> Av = apply(op, t, dB, dBt, dD, dDt, dG, v, E);
    const std::vector<double> A1 = apply(op, t, dB, dBt, dD, dDt, dG, ones, E);

    double err = 0.0, scale = 0.0, null = 0.0, vAu = 0.0, uAv = 0.0;
    bool finite = true;
    std::vector<double> ref(n3);
    for (long long e = 0; e < Echeck; ++e) {
        reference(t, p, &u[e * n3], &G[e * K * q3], ref.data());
        for (int i = 0; i < n3; ++i) {
            finite = finite && std::isfinite(Au[e * n3 + i]) && std::isfinite(A1[e * n3 + i]);
            err = std::fmax(err, std::fabs(Au[e * n3 + i] - ref[i]));
            scale = std::fmax(scale, std::fabs(ref[i]));
            null = std::fmax(null, std::fabs(A1[e * n3 + i]));
        }
    }
    for (long long i = 0; i < E * n3; ++i) vAu += v[i] * Au[i], uAv += u[i] * Av[i];
    const double asym = std::fabs(vAu - uAv) / std::fmax(std::fabs(vAu), 1e-300);

    const bool pass = finite && err <= 1e-12 * scale && null <= 1e-12 * scale && asym <= 1e-12;
    printf("bp3_apply_p%d, %lld elements, kernel %s\n", p, E, file.c_str());
    printf("  vs direct quadrature sums: %.3e of max %.3e\n", err, scale);
    printf("  A 1: max %.3e; symmetry |<v,Au>-<u,Av>|/|<v,Au>|: %.3e\n", null, asym);
    if (reps > 0) {
        double* du = to_device(u);
        double* dy = nullptr;
        GATE_CK(cudaMalloc(&dy, u.size() * sizeof(double)));
        for (int w = 0; w < 3; ++w) op.apply(du, dB, dBt, dD, dDt, dG, dy, E);
        GATE_CK(cudaDeviceSynchronize());
        std::vector<double> ms(reps);
        for (int r = 0; r < reps; ++r) {
            const auto t0 = std::chrono::steady_clock::now();
            op.apply(du, dB, dBt, dD, dDt, dG, dy, E);
            GATE_CK(cudaDeviceSynchronize());
            ms[r] = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        }
        GATE_CK(cudaGetLastError());
        std::sort(ms.begin(), ms.end());
        const double med = ms[reps / 2];
        printf("  timing: %d warm applies, median %.4f ms (min %.4f, max %.4f), %.3f GDoF/s element-local\n", reps,
               med, ms.front(), ms.back(), (double)E * n3 / (med * 1e6));
    }
    printf("compared %lld of %lld elements on the host\n", Echeck, E);
    printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
