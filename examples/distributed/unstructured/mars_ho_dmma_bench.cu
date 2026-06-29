// DMMA throughput benchmark: FP64 tensor-core Laplacian apply
// (ho_laplacian_apply_cg_dmma) vs the scalar CUDA-core reference
// (ho_laplacian_apply_cg), over a cube-scale element set. Reports per-apply
// time, achieved useful FP64 GFLOP/s, and the DMMA-over-scalar speedup.
//
// "Useful" FLOPs = the math the operator must do (not padded tile FLOPs):
//   6 contractions * (2 * n^4) + metric (N3 * 15), per element. The tensor
//   cores may move more than this; we credit only the work the algorithm needs.
//
// Both paths scatter additively, so each timed apply includes the d_y memset
// (as a real CG matvec would). P=7 (n=8) is the native m8n8k4 tile.
//
//   --E=<elems/dim>   structured cube E^3 elements (default 32)
//   --iters=<N>       timed applies per kernel (default 30)

#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian_dmma.hpp"

#include <cuda_runtime.h>
#include <array>
#include <cstdio>
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return false; } } while(0)

template<int P>
bool runBench(int E, int iters)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    using Real = double;
    const Real h = 1.0 / E;

    auto cgid = [&](int cx, int cy, int cz) { return (cx * (E + 1) + cy) * (E + 1) + cz; };
    std::vector<std::array<int,8>> elemCorners;
    elemCorners.reserve(size_t(E) * E * E);
    for (int ex = 0; ex < E; ++ex)
    for (int ey = 0; ey < E; ++ey)
    for (int ez = 0; ez < E; ++ez)
        elemCorners.push_back({
            cgid(ex,ey,ez),     cgid(ex+1,ey,ez),   cgid(ex+1,ey+1,ez),   cgid(ex,ey+1,ez),
            cgid(ex,ey,ez+1),   cgid(ex+1,ey,ez+1), cgid(ex+1,ey+1,ez+1), cgid(ex,ey+1,ez+1) });
    long nCorner = long(E + 1) * (E + 1) * (E + 1);

    HODofHandler dh;
    dh.build(elemCorners, nCorner, P);
    long nDof = dh.numDof;
    size_t nEl = elemCorners.size();

    std::vector<double> x, w, D;
    gllBasis(P, x, w, D);
    std::vector<Real> Dr(D.begin(), D.end());

    // Affine-cube diagonal metric (value irrelevant to timing).
    const Real hhalf = 0.5 * h;
    std::vector<Real> Gelem(N3 * 6, 0.0);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
        int p = i*n*n + j*n + k; Real d = hhalf * w[i]*w[j]*w[k];
        Gelem[p*6+0] = d; Gelem[p*6+3] = d; Gelem[p*6+5] = d;
    }
    std::vector<Real> G(nEl * N3 * 6);
    for (size_t e = 0; e < nEl; ++e) std::copy(Gelem.begin(), Gelem.end(), G.begin() + e*N3*6);

    // Affine path: constant geometric metric (6 doubles/element, weightless), with
    // the GLL weight applied in-kernel from the 1D weights. G = w_ijk * Ghat, so
    // Ghat is the per-point G with the weight divided out -- here diag(h/2).
    std::vector<Real> Ghat(nEl * 6, 0.0);
    for (size_t e = 0; e < nEl; ++e) { Ghat[e*6+0] = hhalf; Ghat[e*6+3] = hhalf; Ghat[e*6+5] = hhalf; }
    std::vector<Real> wgt(w.begin(), w.end());   // n GLL quadrature weights

    std::mt19937 rng(7); std::uniform_real_distribution<Real> uni(-1,1);
    std::vector<Real> uu(nDof);
    for (auto& z : uu) z = uni(rng);

    int *dElemDof; Real *dD, *dG, *dGhat, *dW, *dU, *dY;
    CK(cudaMalloc(&dElemDof, dh.elemDof.size()*sizeof(int)));
    CK(cudaMalloc(&dD, Dr.size()*sizeof(Real)));
    CK(cudaMalloc(&dG, G.size()*sizeof(Real)));
    CK(cudaMalloc(&dGhat, Ghat.size()*sizeof(Real)));
    CK(cudaMalloc(&dW, wgt.size()*sizeof(Real)));
    CK(cudaMalloc(&dU, nDof*sizeof(Real)));
    CK(cudaMalloc(&dY, nDof*sizeof(Real)));
    CK(cudaMemcpy(dElemDof, dh.elemDof.data(), dh.elemDof.size()*sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dD, Dr.data(), Dr.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(), G.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dGhat, Ghat.data(), Ghat.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dW, wgt.data(), wgt.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dU, uu.data(), nDof*sizeof(Real), cudaMemcpyHostToDevice));

    const int blockScalar = 128, gridScalar = int((nEl + blockScalar - 1)/blockScalar);
    auto scalarApply = [&]() {
        cudaMemset(dY, 0, nDof*sizeof(Real));
        ho_laplacian_apply_cg<Real, P><<<gridScalar, blockScalar>>>(dU, dY, dElemDof, dD, dG, nEl);
    };
    auto dmmaApply = [&]() {
        cudaMemset(dY, 0, nDof*sizeof(Real));
        ho_laplacian_dmma_launch<Real, P, 64>(dU, dY, dElemDof, dD, dG, nEl);
    };
    auto affineApply = [&]() {
        cudaMemset(dY, 0, nDof*sizeof(Real));
        ho_laplacian_dmma_launch_affine<Real, P, 64>(dU, dY, dElemDof, dD, dGhat, dW, nEl);
    };
    // Matched-structure scalar: same warp+smem kernel as affine DMMA, scalar tile.
    auto smemScalarApply = [&]() {
        cudaMemset(dY, 0, nDof*sizeof(Real));
        ho_laplacian_smem_scalar_launch_affine<Real, P, 64>(dU, dY, dElemDof, dD, dGhat, dW, nEl);
    };

    // Warmup (also surfaces launch errors before timing).
    scalarApply(); dmmaApply(); affineApply(); smemScalarApply();
    CK(cudaGetLastError());
    CK(cudaDeviceSynchronize());

    // Parity: affine must reproduce per-point dmma (G = w_ijk * Ghat is exact).
    std::vector<Real> Yd(nDof), Ya(nDof);
    dmmaApply();   CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(Yd.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
    affineApply(); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(Ya.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
    Real ynorm = 0, maxd = 0;
    for (long p = 0; p < nDof; ++p) ynorm = std::max(ynorm, std::abs(Yd[p]));
    for (long p = 0; p < nDof; ++p) maxd = std::max(maxd, std::abs(Ya[p] - Yd[p]));
    const Real parity = maxd / (ynorm > 0 ? ynorm : 1.0);

    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    auto timeIt = [&](auto fn) -> float {
        cudaEventRecord(t0);
        for (int it = 0; it < iters; ++it) fn();
        cudaEventRecord(t1); cudaEventSynchronize(t1);
        float ms = 0; cudaEventElapsedTime(&ms, t0, t1); return ms / iters;
    };

    float msScalar     = timeIt(scalarApply);
    float msDmma       = timeIt(dmmaApply);
    float msAffine     = timeIt(affineApply);
    float msSmemScalar = timeIt(smemScalarApply);
    CK(cudaGetLastError());

    const double flopPerElem = 6.0 * 2.0 * std::pow((double)n, 4) + (double)N3 * 15.0;
    const double totalFlop   = flopPerElem * (double)nEl;
    auto gflops = [&](float ms) { return totalFlop / (ms * 1e-3) / 1e9; };

    // Two speedups: vs the naive thread-per-element scalar (structure+TC bundled),
    // and the ISOLATED tensor-core speedup vs the matched warp+smem scalar twin.
    printf("P=%d E=%d | elems=%zu nDof=%ld | affine-vs-perpoint parity=%.2e\n"
           "  scalar (thread/elem) : %.3f ms/apply  %.1f GFLOP/s  (%.2fx)\n"
           "  scalar (warp+smem)   : %.3f ms/apply  %.1f GFLOP/s  (%.2fx)\n"
           "  dmma   (G/point)     : %.3f ms/apply  %.1f GFLOP/s  (%.2fx)\n"
           "  dmma   (affine)      : %.3f ms/apply  %.1f GFLOP/s  (%.2fx)\n"
           "  -> isolated tensor-core speedup (affine dmma vs matched scalar) = %.2fx\n",
           P, E, nEl, nDof, parity,
           msScalar,     gflops(msScalar),     1.0,
           msSmemScalar, gflops(msSmemScalar), msScalar / msSmemScalar,
           msDmma,       gflops(msDmma),       msScalar / msDmma,
           msAffine,     gflops(msAffine),     msScalar / msAffine,
           msSmemScalar / msAffine);

    cudaEventDestroy(t0); cudaEventDestroy(t1);
    cudaFree(dElemDof); cudaFree(dD); cudaFree(dG); cudaFree(dGhat); cudaFree(dW);
    cudaFree(dU); cudaFree(dY);
    return true;
}

int main(int argc, char** argv)
{
    int E = 32, iters = 30;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.find("--E=")==0) E = std::stoi(a.substr(4));
        else if (a.find("--iters=")==0) iters = std::stoi(a.substr(8)); }
    bool ok = runBench<7>(E, iters);
    return ok ? 0 : 1;
}
