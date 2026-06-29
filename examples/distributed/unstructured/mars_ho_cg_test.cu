// Stage 1: single-rank continuous-Galerkin matrix-free Laplacian.
// Validates the element->global DOF map + additive scatter (direct stiffness
// summation) — the bridge from the DG block-diagonal operator to a CG operator,
// identified as the biggest risk in the Phase B design.
//
// DOF map here is a STRUCTURED LATTICE (Ex*Ey*Ez axis-aligned hexes, order p):
// element (ex,ey,ez) local node (i,j,k) -> global GLL lattice node
// (ex*p+i, ey*p+j, ez*p+k) on the (E*p+1)^3 grid. Continuity is guaranteed by
// construction (shared faces/edges/corners map to one global DOF), so this
// isolates the gather/scatter-add machinery from the topological entity
// enumeration (next increment, which is a drop-in replacement for elemDof).
//
// Gates (no external reference matrix needed):
//   1. DOF count == (E*p+1)^3
//   2. nullspace:  A * 1 = 0           (all DOFs)
//   3. symmetry:   u^T A v = v^T A u
//   4. continuity: A * (linear field) = 0 at INTERIOR DOFs   <- the real test:
//      a broken seam (missing shared DOF / bad scatter) leaves nonzero residual
//      for a linear field exactly at element interfaces.

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
bool runStage1(int E)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    using Real = double;

    const int    G1   = E * P + 1;             // global nodes per dimension
    const long   nDof = long(G1) * G1 * G1;    // total global DOFs
    const size_t nEl  = size_t(E) * E * E;
    const Real   h    = 1.0 / E;               // element edge (domain [0,1]^3)

    std::vector<double> x, w, D;
    gllBasis(P, x, w, D);
    std::vector<Real> Dr(D.begin(), D.end());

    // elemDof[e*N3 + (i*n*n + j*n + k)] = global lattice index.
    std::vector<int> elemDof(nEl * N3);
    auto gidx = [&](int gx, int gy, int gz) { return (long(gx) * G1 + gy) * G1 + gz; };
    for (int ex = 0; ex < E; ++ex)
    for (int ey = 0; ey < E; ++ey)
    for (int ez = 0; ez < E; ++ez)
    {
        size_t e = (size_t(ex) * E + ey) * E + ez;
        for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k)
        {
            int l = i * n * n + j * n + k;
            elemDof[e * N3 + l] = int(gidx(ex * P + i, ey * P + j, ez * P + k));
        }
    }

    // Affine-cube metric: diagonal, G = (h/2) * w_i w_j w_k, same for all elements.
    std::vector<Real> Gelem(N3 * 6, 0.0);
    const Real hhalf = 0.5 * h;
    for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
    for (int k = 0; k < n; ++k)
    {
        int p   = i * n * n + j * n + k;
        Real wp = w[i] * w[j] * w[k];
        Gelem[p * 6 + 0] = hhalf * wp;
        Gelem[p * 6 + 3] = hhalf * wp;
        Gelem[p * 6 + 5] = hhalf * wp;
    }
    std::vector<Real> G(nEl * N3 * 6);
    for (size_t e = 0; e < nEl; ++e)
        std::copy(Gelem.begin(), Gelem.end(), G.begin() + e * N3 * 6);

    // Device
    int   *dElemDof;
    Real  *dD, *dG, *dU, *dY;
    CK(cudaMalloc(&dElemDof, elemDof.size() * sizeof(int)));
    CK(cudaMalloc(&dD, Dr.size() * sizeof(Real)));
    CK(cudaMalloc(&dG, G.size() * sizeof(Real)));
    CK(cudaMalloc(&dU, nDof * sizeof(Real)));
    CK(cudaMalloc(&dY, nDof * sizeof(Real)));
    CK(cudaMemcpy(dElemDof, elemDof.data(), elemDof.size() * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dD, Dr.data(), Dr.size() * sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(),  G.size() * sizeof(Real), cudaMemcpyHostToDevice));

    int block = 128;
    int grid  = int((nEl + block - 1) / block);
    auto apply = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpy(dU, u.data(), nDof * sizeof(Real), cudaMemcpyHostToDevice));
        CK(cudaMemset(dY, 0, nDof * sizeof(Real)));
        ho_laplacian_apply_cg<Real, P><<<grid, block>>>(dU, dY, dElemDof, dD, dG, nEl);
        CK(cudaDeviceSynchronize());
        y.resize(nDof);
        CK(cudaMemcpy(y.data(), dY, nDof * sizeof(Real), cudaMemcpyDeviceToHost));
        return true;
    };

    // (2) nullspace
    std::vector<Real> ones(nDof, 1.0), y;
    if (!apply(ones, y)) return false;
    Real nullMax = 0;
    for (Real v : y) nullMax = std::max(nullMax, std::abs(v));

    // (3) symmetry + PSD
    std::mt19937 rng(2024);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);
    std::vector<Real> uu(nDof), vv(nDof), Au, Av;
    for (auto& z : uu) z = uni(rng);
    for (auto& z : vv) z = uni(rng);
    if (!apply(uu, Au)) return false;
    if (!apply(vv, Av)) return false;
    Real uAv = 0, vAu = 0, uAu = 0;
    for (long p = 0; p < nDof; ++p) { uAv += uu[p] * Av[p]; vAu += vv[p] * Au[p]; uAu += uu[p] * Au[p]; }
    Real symErr = std::abs(uAv - vAu) / (std::abs(uAv) + 1.0);

    // (4) continuity: u = physical x-coordinate; A*u must be 0 at interior DOFs.
    // The global node gx maps to element ex=gx/P, local GLL index i=gx-ex*P;
    // physical x = (ex + (x_gll[i]+1)/2) * h. NOTE: lattice index is NOT
    // physical position because GLL nodes are clustered, not equispaced
    // (only coincides for p=2). Using the index here would inject a kinked
    // (non-linear) field and produce a spurious interior residual at p>=4.
    auto physCoord = [&](int g) -> Real {
        int ex = (g == G1 - 1) ? (E - 1) : (g / P);
        int i  = g - ex * P;
        return (ex + 0.5 * (x[i] + 1.0)) * h;
    };
    std::vector<Real> ux(nDof), Ax;
    for (int gx = 0; gx < G1; ++gx)
    for (int gy = 0; gy < G1; ++gy)
    for (int gz = 0; gz < G1; ++gz)
        ux[gidx(gx, gy, gz)] = physCoord(gx);
    if (!apply(ux, Ax)) return false;
    Real interMax = 0;
    for (int gx = 1; gx < G1 - 1; ++gx)
    for (int gy = 1; gy < G1 - 1; ++gy)
    for (int gz = 1; gz < G1 - 1; ++gz)
        interMax = std::max(interMax, std::abs(Ax[gidx(gx, gy, gz)]));

    bool ok = (nDof == long(G1) * G1 * G1) && (nullMax < 1e-9) && (symErr < 1e-10)
              && (uAu > -1e-12) && (interMax < 1e-9);

    printf("P=%d E=%d | nDof=%ld | null=%.2e | sym=%.2e | uAu=%.4e | continuity A*x interior=%.2e | %s\n",
           P, E, nDof, nullMax, symErr, uAu, interMax, ok ? "PASS" : "FAIL");

    cudaFree(dElemDof); cudaFree(dD); cudaFree(dG); cudaFree(dU); cudaFree(dY);
    return ok;
}

int main()
{
    bool ok = true;
    ok &= runStage1<2>(4);
    ok &= runStage1<4>(3);
    ok &= runStage1<7>(2);
    printf("Stage 1 (single-rank CG, lattice DOF map): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
