// Stage 0: element-local parity gate for the matrix-free high-order Laplacian.
// One affine reference-cube hex element, order p. Validates the sum-factorized
// kernel (mars_ho_laplacian.hpp) + GLL basis (mars_spectral_basis.hpp) against
// operator properties that need no external reference matrix:
//   1. nullspace:  A * 1 = 0     (constant field has zero Laplacian)
//   2. symmetry:   u^T A v = v^T A u
//   3. PSD:        u^T A u >= 0
//
// Single GPU, one thread per element (thread-per-element is fine for a 1-element
// correctness check; the production p=7 kernel will be block-per-element).
//
// Affine cube of edge h, reference [-1,1]^3: J=(h/2)I, detJ=(h/2)^3,
// J^{-T}=(2/h)I. The Laplacian weak term w_q detJ (J^{-T} grad)^2 collapses to
// a diagonal metric G = (h/2) * w_i w_j w_k per GLL point (off-diagonals zero).

#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian.hpp"

#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>

using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return false; } } while(0)

template<int P>
bool runStage0()
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    using Real = double;

    std::vector<double> x, w, D;
    gllBasis(P, x, w, D);   // D is n*n row-major

    // Affine-cube metric: diagonal, G = (h/2) * w_i w_j w_k. h = 1 -> h/2 = 0.5.
    std::vector<Real> G(N3 * 6, 0.0);
    const Real hhalf = 0.5;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                int p   = i * n * n + j * n + k;
                Real wp = w[i] * w[j] * w[k];
                G[p * 6 + 0] = hhalf * wp;   // Gxx
                G[p * 6 + 3] = hhalf * wp;   // Gyy
                G[p * 6 + 5] = hhalf * wp;   // Gzz
            }

    std::vector<Real> Dr(D.begin(), D.end());

    Real *dD, *dG, *du, *dy;
    CK(cudaMalloc(&dD, n * n * sizeof(Real)));
    CK(cudaMalloc(&dG, N3 * 6 * sizeof(Real)));
    CK(cudaMalloc(&du, N3 * sizeof(Real)));
    CK(cudaMalloc(&dy, N3 * sizeof(Real)));
    CK(cudaMemcpy(dD, Dr.data(), n * n * sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(),  N3 * 6 * sizeof(Real), cudaMemcpyHostToDevice));

    auto apply = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpy(du, u.data(), N3 * sizeof(Real), cudaMemcpyHostToDevice));
        ho_laplacian_apply_scalar<Real, P><<<1, 1>>>(du, dy, dD, dG, 1);
        CK(cudaDeviceSynchronize());
        y.resize(N3);
        CK(cudaMemcpy(y.data(), dy, N3 * sizeof(Real), cudaMemcpyDeviceToHost));
        return true;
    };

    // (1) nullspace
    std::vector<Real> ones(N3, 1.0), y;
    if (!apply(ones, y)) return false;
    Real nullMax = 0;
    for (Real v : y) nullMax = std::max(nullMax, std::abs(v));

    // (2) symmetry + (3) PSD, random u, v
    std::mt19937 rng(12345);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);
    std::vector<Real> uu(N3), vv(N3), Au, Av;
    for (auto& z : uu) z = uni(rng);
    for (auto& z : vv) z = uni(rng);
    if (!apply(uu, Au)) return false;
    if (!apply(vv, Av)) return false;
    Real uAv = 0, vAu = 0, uAu = 0;
    for (int p = 0; p < N3; ++p) { uAv += uu[p] * Av[p]; vAu += vv[p] * Au[p]; uAu += uu[p] * Au[p]; }

    Real symErr = std::abs(uAv - vAu) / (std::abs(uAv) + 1.0);
    bool ok = (nullMax < 1e-9) && (symErr < 1e-10) && (uAu > -1e-12);

    printf("P=%d n=%d N3=%d | null max|A*1|=%.3e | sym rel|uAv-vAu|=%.3e | uAu=%.6e | %s\n",
           P, n, N3, nullMax, symErr, uAu, ok ? "PASS" : "FAIL");

    cudaFree(dD); cudaFree(dG); cudaFree(du); cudaFree(dy);
    return ok;
}

int main()
{
    bool ok = true;
    ok &= runStage0<2>();
    ok &= runStage0<4>();
    ok &= runStage0<7>();
    printf("Stage 0 (element-local parity): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
