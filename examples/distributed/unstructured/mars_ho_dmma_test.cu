// DMMA parity test: the FP64 tensor-core Laplacian apply
// (ho_laplacian_apply_cg_dmma) must reproduce the scalar reference
// (ho_laplacian_apply_cg) bit-for-bit up to FP64 reduction order. Both kernels
// run on the SAME elemDof / G / D and the SAME random u, on a structured cube
// HO DOF map (mirrors mars_ho_dof_test.cu). P=7 (n=8) is the native m8n8k4 tile.
//
// Gates:
//   1. parity      max relative diff (dmma vs scalar) < 1e-10
//   2. nullspace   A*1 = 0                 (both paths)
//   3. symmetry    u^T A v = v^T A u       (dmma path)
//   4. continuity  A*(physical x) = 0 at interior DOFs (dmma path)

#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian_dmma.hpp"

#include "backend/distributed/unstructured/utils/mars_gpu_runtime.hpp"
#include <array>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>

using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return false; } } while(0)

template<int P, bool Affine = false, bool Varying = false>
bool runDmmaParity(int E)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    static_assert(!Affine || !Varying);
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

    std::vector<double> x, w, D;
    gllBasis(P, x, w, D);
    std::vector<Real> Dr(D.begin(), D.end());

    // General (non-diagonal) symmetric SPD metric per quad point. The off-
    // diagonal Gxy/Gxz/Gyz terms exercise the dir0/dir1/dir2 coupling in the
    // metric step -- a diagonal G would hide an x/y/z mix-up. Built once per
    // element (affine cube), copied to every element.
    const Real hhalf = 0.5 * h;
    std::vector<Real> Gelem(N3 * 6, 0.0);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
        int p = i*n*n + j*n + k; Real wp = w[i]*w[j]*w[k];
        Real d = hhalf*wp;
        Gelem[p*6+0] = d;            // Gxx
        Gelem[p*6+1] = 0.15 * d;     // Gxy
        Gelem[p*6+2] = 0.10 * d;     // Gxz
        Gelem[p*6+3] = d;            // Gyy
        Gelem[p*6+4] = 0.05 * d;     // Gyz
        Gelem[p*6+5] = d;            // Gzz
    }
    size_t nEl = elemCorners.size();
    std::vector<Real> G(nEl * N3 * 6);
    for (size_t e = 0; e < nEl; ++e) std::copy(Gelem.begin(), Gelem.end(), G.begin() + e*N3*6);

    if constexpr (Varying)
        for (size_t e = 0; e < nEl; ++e)
            for (int p = 0; p < N3; ++p)
                for (int c = 0; c < 6; ++c)
                    G[(e * N3 + p) * 6 + c] *= 1 + .1 * std::sin(double(e * N3 + p));

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

    int *dElemDof; Real *dD, *dG, *dU, *dY, *dGhat, *dW;
    CK(cudaMalloc(&dElemDof, dh.elemDof.size()*sizeof(int)));
    CK(cudaMalloc(&dD, Dr.size()*sizeof(Real)));
    CK(cudaMalloc(&dG, G.size()*sizeof(Real)));
    std::vector<Real> Ghat(nEl * 6);
    const Real metric[6] = {hhalf, .15*hhalf, .10*hhalf, hhalf, .05*hhalf, hhalf};
    for(size_t e=0; e<nEl; ++e) std::copy(metric, metric+6, Ghat.data()+6*e);
    CK(cudaMalloc(&dGhat, Ghat.size()*sizeof(Real)));
    CK(cudaMalloc(&dW, w.size()*sizeof(Real)));
    CK(cudaMemcpy(dGhat, Ghat.data(), Ghat.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dW, w.data(), w.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&dU, nDof*sizeof(Real)));
    CK(cudaMalloc(&dY, nDof*sizeof(Real)));
    CK(cudaMemcpy(dElemDof, dh.elemDof.data(), dh.elemDof.size()*sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dD, Dr.data(), Dr.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(), G.size()*sizeof(Real), cudaMemcpyHostToDevice));

    int block = 128, grid = int((nEl + block - 1)/block);

    auto applyScalar = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpy(dU, u.data(), nDof*sizeof(Real), cudaMemcpyHostToDevice));
        CK(cudaMemset(dY, 0, nDof*sizeof(Real)));
        ho_laplacian_apply_cg<Real, P><<<grid, block>>>(dU, dY, dElemDof, dD, dG, nEl);
        CK(cudaGetLastError());
        CK(cudaDeviceSynchronize());
        y.resize(nDof);
        CK(cudaMemcpy(y.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
        return std::all_of(y.begin(), y.end(), [](Real v) { return std::isfinite(v); });
    };

    // A non-default stream checks that the public launcher preserves caller ordering.
    cudaStream_t stream;
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
    CK(hipStreamCreateWithFlags(&stream, hipStreamNonBlocking));
#else
    CK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
#endif
    auto applyDmma = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpyAsync(dU, u.data(), nDof*sizeof(Real), cudaMemcpyHostToDevice, stream));
        CK(cudaMemsetAsync(dY, 0, nDof*sizeof(Real), stream));
        if constexpr (Affine)
            CK((ho_laplacian_dmma_launch_affine<Real, P>(dU, dY, dElemDof, dD, dGhat, dW, nEl, stream)));
        else
            CK((ho_laplacian_dmma_launch<Real, P>(dU, dY, dElemDof, dD, dG, nEl, stream)));
        CK(cudaStreamSynchronize(stream));
        y.resize(nDof);
        CK(cudaMemcpy(y.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
        return std::all_of(y.begin(), y.end(), [](Real v) { return std::isfinite(v); });
    };

    // --- Gate 1: parity on a random u ---
    std::mt19937 rng(7); std::uniform_real_distribution<Real> uni(-1,1);
    std::vector<Real> uu(nDof), Ys, Yd;
    for (auto& z : uu) z = uni(rng);
    if (!applyScalar(uu, Ys)) return false;
    if (!applyDmma(uu, Yd)) return false;

    // The scatter must add to existing output, not overwrite it.
    if constexpr (Affine)
        CK((ho_laplacian_dmma_launch_affine<Real, P>(dU, dY, dElemDof, dD, dGhat, dW, nEl, stream)));
    else
        CK((ho_laplacian_dmma_launch<Real, P>(dU, dY, dElemDof, dD, dG, nEl, stream)));
    CK(cudaStreamSynchronize(stream));
    std::vector<Real> doubled(nDof);
    CK(cudaMemcpy(doubled.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
    Real additiveError = 0, additiveScale = 1;
    for (long p = 0; p < nDof; ++p) {
        if (!std::isfinite(doubled[p])) return false;
        additiveError = std::max(additiveError, std::abs(doubled[p] - 2 * Ys[p]));
        additiveScale = std::max(additiveScale, std::abs(2 * Ys[p]));
    }

    Real maxRel = 0, ynorm = 0;
    for (long p = 0; p < nDof; ++p) ynorm = std::max(ynorm, std::abs(Ys[p]));
    const Real scale = ynorm > 0 ? ynorm : 1.0;
    for (long p = 0; p < nDof; ++p)
        maxRel = std::max(maxRel, std::abs(Yd[p] - Ys[p]) / scale);

    // --- Gate 2: nullspace A*1 = 0 (dmma) ---
    std::vector<Real> ones(nDof, 1.0), Yn;
    if (!applyDmma(ones, Yn)) return false;
    Real nullMax = 0; for (Real v : Yn) nullMax = std::max(nullMax, std::abs(v));

    // --- Gate 3: symmetry u^T A v = v^T A u (dmma) ---
    std::vector<Real> vv(nDof), Au, Av;
    for (auto& z : vv) z = uni(rng);
    if (!applyDmma(uu, Au)) return false;
    if (!applyDmma(vv, Av)) return false;
    Real uAv=0, vAu=0; for (long p=0;p<nDof;++p){ uAv+=uu[p]*Av[p]; vAu+=vv[p]*Au[p]; }
    Real symErr = std::abs(uAv-vAu)/(std::abs(uAv)+1.0);

    // A linear field cancels only for the constant metric fixtures.
    // The varying metric fixture is checked by scalar parity instead.
    std::vector<Real> Ax;
    if (!applyDmma(dofX, Ax)) return false;
    Real interMax = 0;
    for (long p=0;p<nDof;++p) if (!Varying && !isBdry[p]) interMax = std::max(interMax, std::abs(Ax[p]));

    bool ok = (additiveError/additiveScale < 1e-10) && (maxRel < 1e-10) && (nullMax < 1e-9) && (symErr < 1e-10) && (interMax < 1e-9);
    printf("path=public-launcher metric=%s additive=%.2e ", Varying ? "varying" : (Affine ? "affine" : "per-point"), additiveError/additiveScale);
    printf("P=%d E=%d | nDof=%ld | parity(dmma vs scalar)=%.2e | null=%.2e | sym=%.2e | continuity=%.2e | %s\n",
           P, E, nDof, maxRel, nullMax, symErr, interMax, ok?"PASS":"FAIL");

    if constexpr (Varying) printf("Linear-field cancellation skipped for varying metric.\n");
    CK(cudaStreamDestroy(stream));
    cudaFree(dGhat); cudaFree(dW); cudaFree(dElemDof); cudaFree(dD); cudaFree(dG); cudaFree(dU); cudaFree(dY);
    return ok;
}

int main()
{
    bool ok = true;
    for (int E : {1, 2, 3})
    {
        ok &= runDmmaParity<7, false>(E);
        ok &= runDmmaParity<7, true>(E);
    }
    ok &= runDmmaParity<7, false, true>(3);
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
    // Empty partitions must not launch a zero-sized grid or touch any pointer.
    ok &= ho_laplacian_dmma_launch<double, 7>(nullptr, nullptr, nullptr, nullptr, nullptr, 0) == cudaSuccess;
    ok &= ho_laplacian_dmma_launch_affine<double, 7>(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0) == cudaSuccess;
#endif
    printf("DMMA FP64 tensor-core Laplacian parity (vs scalar reference): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
