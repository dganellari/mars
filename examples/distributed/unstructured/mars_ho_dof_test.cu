// Stage 1b: validate the topological HO DOF handler (mars_ho_dof_handler.hpp)
// against the known-good lattice result on a structured cube. The handler keys
// shared edges/faces on global corner ids (the real Phase B machinery), with no
// lattice assumption -- so this proves the entity enumeration + continuity match
// the structured ground truth before we feed it SFC-ordered ElementDomain meshes.
//
// Gates:
//   1. numDof == (E*p+1)^3        (entity enumeration completeness, no double-count)
//   2. nullspace  A*1 = 0
//   3. symmetry   u^T A v = v^T A u
//   4. continuity A*(physical x) = 0 at interior DOFs

#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian.hpp"

#include <cuda_runtime.h>
#include <array>
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>

using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return false; } } while(0)

template<int P>
bool runStage1b(int E)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    using Real = double;
    const Real h = 1.0 / E;

    // Structured cube corner connectivity, hex order matching the handler.
    auto cgid = [&](int cx, int cy, int cz) { return (cx * (E + 1) + cy) * (E + 1) + cz; };
    std::vector<std::array<int,8>> elemCorners;
    elemCorners.reserve(size_t(E) * E * E);
    std::vector<std::array<int,3>> elemIJK;   // remember (ex,ey,ez) for coords
    for (int ex = 0; ex < E; ++ex)
    for (int ey = 0; ey < E; ++ey)
    for (int ez = 0; ez < E; ++ez)
    {
        elemCorners.push_back({
            cgid(ex,ey,ez),     cgid(ex+1,ey,ez),   cgid(ex+1,ey+1,ez),   cgid(ex,ey+1,ez),
            cgid(ex,ey,ez+1),   cgid(ex+1,ey,ez+1), cgid(ex+1,ey+1,ez+1), cgid(ex,ey+1,ez+1) });
        elemIJK.push_back({ex, ey, ez});
    }
    long nCorner = long(E + 1) * (E + 1) * (E + 1);

    HODofHandler dh;
    dh.build(elemCorners, nCorner, P);
    long nDof = dh.numDof;
    long expect = long(E * P + 1) * (E * P + 1) * (E * P + 1);

    std::vector<double> x, w, D;
    gllBasis(P, x, w, D);
    std::vector<Real> Dr(D.begin(), D.end());

    // Per-element affine metric (diagonal), and per-DOF physical x + boundary flag.
    const Real hhalf = 0.5 * h;
    std::vector<Real> Gelem(N3 * 6, 0.0);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
        int p = i*n*n + j*n + k; Real wp = w[i]*w[j]*w[k];
        Gelem[p*6+0] = hhalf*wp; Gelem[p*6+3] = hhalf*wp; Gelem[p*6+5] = hhalf*wp;
    }
    size_t nEl = elemCorners.size();
    std::vector<Real> G(nEl * N3 * 6);
    for (size_t e = 0; e < nEl; ++e) std::copy(Gelem.begin(), Gelem.end(), G.begin() + e*N3*6);

    std::vector<Real>    dofX(nDof, 0.0);
    std::vector<uint8_t> isBdry(nDof, 0);
    auto coord1d = [&](int e_ijk, int loc) { return (e_ijk + 0.5*(x[loc]+1.0)) * h; };
    for (size_t e = 0; e < nEl; ++e)
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
            int l = i*n*n + j*n + k;
            int dof = dh.elemDof[e*N3 + l];
            Real px = coord1d(elemIJK[e][0], i), py = coord1d(elemIJK[e][1], j), pz = coord1d(elemIJK[e][2], k);
            dofX[dof] = px;
            bool b = (px<1e-12||px>1-1e-12||py<1e-12||py>1-1e-12||pz<1e-12||pz>1-1e-12);
            isBdry[dof] = b ? 1 : 0;
        }

    int *dElemDof; Real *dD, *dG, *dU, *dY;
    CK(cudaMalloc(&dElemDof, dh.elemDof.size()*sizeof(int)));
    CK(cudaMalloc(&dD, Dr.size()*sizeof(Real)));
    CK(cudaMalloc(&dG, G.size()*sizeof(Real)));
    CK(cudaMalloc(&dU, nDof*sizeof(Real)));
    CK(cudaMalloc(&dY, nDof*sizeof(Real)));
    CK(cudaMemcpy(dElemDof, dh.elemDof.data(), dh.elemDof.size()*sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dD, Dr.data(), Dr.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(), G.size()*sizeof(Real), cudaMemcpyHostToDevice));

    int block = 128, grid = int((nEl + block - 1)/block);
    auto apply = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpy(dU, u.data(), nDof*sizeof(Real), cudaMemcpyHostToDevice));
        CK(cudaMemset(dY, 0, nDof*sizeof(Real)));
        ho_laplacian_apply_cg<Real, P><<<grid, block>>>(dU, dY, dElemDof, dD, dG, nEl);
        CK(cudaDeviceSynchronize());
        y.resize(nDof);
        CK(cudaMemcpy(y.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
        return true;
    };

    std::vector<Real> ones(nDof, 1.0), y;
    if (!apply(ones, y)) return false;
    Real nullMax = 0; for (Real v : y) nullMax = std::max(nullMax, std::abs(v));

    std::mt19937 rng(7); std::uniform_real_distribution<Real> uni(-1,1);
    std::vector<Real> uu(nDof), vv(nDof), Au, Av;
    for (auto& z : uu) z = uni(rng);
    for (auto& z : vv) z = uni(rng);
    if (!apply(uu, Au)) return false;
    if (!apply(vv, Av)) return false;
    Real uAv=0, vAu=0; for (long p=0;p<nDof;++p){ uAv+=uu[p]*Av[p]; vAu+=vv[p]*Au[p]; }
    Real symErr = std::abs(uAv-vAu)/(std::abs(uAv)+1.0);

    std::vector<Real> Ax;
    if (!apply(dofX, Ax)) return false;
    Real interMax = 0;
    for (long p=0;p<nDof;++p) if (!isBdry[p]) interMax = std::max(interMax, std::abs(Ax[p]));

    bool ok = (nDof==expect) && (nullMax<1e-9) && (symErr<1e-10) && (interMax<1e-9);
    printf("P=%d E=%d | nDof=%ld (expect %ld) | nEdge=%ld nFace=%ld | null=%.2e | sym=%.2e | continuity=%.2e | %s\n",
           P, E, nDof, expect, dh.nEdge, dh.nFace, nullMax, symErr, interMax, ok?"PASS":"FAIL");

    cudaFree(dElemDof); cudaFree(dD); cudaFree(dG); cudaFree(dU); cudaFree(dY);
    return ok;
}

int main()
{
    bool ok = true;
    ok &= runStage1b<2>(4);
    ok &= runStage1b<4>(3);
    ok &= runStage1b<7>(2);
    printf("Stage 1b (topological DOF handler, structured validation): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
