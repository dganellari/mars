#ifndef MARS_GPU_ACM_PRECONDITIONER_HPP
#define MARS_GPU_ACM_PRECONDITIONER_HPP

// ACM Stage 4b: the GPU agglomeration V-cycle (Stage 3 coarse-op + Stage 4a cycle) wrapped as a
// pluggable Preconditioner for FlexGMRES. setup() builds the multilevel hierarchy ONCE; apply(r,z)
// runs ONE V-cycle (z = M^-1 r): load r into the finest bvec, zero xvec, cycle, read xvec back.
//
// The coarsest damped-Jacobi is non-contractive on the indefinite saddle operator (Stage 4a note),
// so M is approximate / non-stationary -> it MUST be driven by FLEXIBLE GMRES (the Z-basis), not
// stationary GMRES. apply() re-zeros xvec and reloads bvec every call, so z is output-only with no
// state carried between Krylov iterations.

#include "mars_cg_solver_with_preconditioner.hpp"   // Preconditioner base (mars::fem)
#include "mars_acm_vcycle.hpp"                       // AcmLevel, acmBuildHierarchy, acmVcycleGpu (mars)
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <type_traits>
#include <cusolverDn.h>   // dense LU for the factor-once coarsest solve

namespace mars
{
namespace fem
{

template<typename RealType, typename IndexType, typename AcceleratorTag>
class GpuAcmPreconditioner : public Preconditioner<RealType, IndexType, AcceleratorTag>
{
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    // maxCoarseND=256: directional aggregation coarsens fast, so a moderate coarsest is reached in few
    // levels and is solved directly by QR. beta defaults to 0.5 (the paper's directional factor); env
    // MARS_ACM_BETA overrides it (set 0 for isotropic, to A/B against directional in one binary).
    GpuAcmPreconditioner(int kmax = 4, int maxCoarseND = 256,
                         int pre = 2, int post = 1, int coarse = 10, RealType omega = RealType(0.7))
        : kmax_(kmax), maxCoarseND_(maxCoarseND), pre_(pre), post_(post), coarse_(coarse),
          omega_(omega), beta_(RealType(0.5))
    {
        if (const char* e = std::getenv("MARS_ACM_BETA")) beta_ = static_cast<RealType>(std::atof(e));
        if (const char* e = std::getenv("MARS_ACM_PRE"))  pre_  = std::atoi(e);   // V-cycle pre-smooth sweeps (default 2)
        if (const char* e = std::getenv("MARS_ACM_POST")) post_ = std::atoi(e);   // post-smooth (default 1); fewer sweeps = cheaper smoother (66% of solve) vs more Krylov iters
        cusolverSpCreate(&cs_);
        cusparseCreateMatDescr(&descr_);
        cusparseSetMatType(descr_, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr_, CUSPARSE_INDEX_BASE_ZERO);
        cudaStreamCreate(&capStream_);              // dedicated capturable stream for the graphed V-cycle
        cusolverDnCreate(&cd_);                     // dense LU: factor the coarsest ONCE/Picard, solve per V-cycle
    }

    ~GpuAcmPreconditioner()
    {
        destroyGraphs();
        if (cd_) cusolverDnDestroy(cd_);
        if (capStream_) cudaStreamDestroy(capStream_);
        if (descr_) cusparseDestroyMatDescr(descr_);
        if (cs_) cusolverSpDestroy(cs_);
    }

    GpuAcmPreconditioner(const GpuAcmPreconditioner&) = delete;             // owns a cusolver handle
    GpuAcmPreconditioner& operator=(const GpuAcmPreconditioner&) = delete;

    // build the multilevel hierarchy once from the finest coupled CSR (interleaved 4*node+comp)
    void setup(const Matrix& A) override
    {
        if (reuse_ && !levels_.empty()) { lastSetupMs_ = 0; return; }   // frozen preconditioner: reuse hierarchy, no build
        cudaDeviceSynchronize();
        auto tSetup0 = std::chrono::high_resolution_clock::now();
        const int ND = static_cast<int>(A.numRows());
        const int nz = static_cast<int>(A.nnz());
        thrust::device_vector<IndexType> ro(ND + 1), ci(nz);
        thrust::device_vector<RealType>  va(nz);
        thrust::copy(thrust::device_pointer_cast(A.rowOffsetsPtr()),
                     thrust::device_pointer_cast(A.rowOffsetsPtr() + ND + 1), ro.begin());
        thrust::copy(thrust::device_pointer_cast(A.colIndicesPtr()),
                     thrust::device_pointer_cast(A.colIndicesPtr() + nz), ci.begin());
        thrust::copy(thrust::device_pointer_cast(A.valuesPtr()),
                     thrust::device_pointer_cast(A.valuesPtr() + nz), va.begin());
        // acmBuildHierarchy wants int CSR; IndexType is int in this driver's instantiation.
        mars::acmBuildHierarchy<RealType>(ro, ci, va, ND / 4, kmax_, maxCoarseND_, levels_, beta_);
        nd_ = ND;
        factorCoarseLU();   // dense LU of the small coarsest op, ONCE (vs the old per-V-cycle cusolverSp QR re-factor)
        cudaDeviceSynchronize();
        lastSetupMs_ = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - tSetup0).count();   // hierarchy build + coarse factor
        captureGraphs();   // (re)capture the V-cycle into CUDA graphs on the new hierarchy (DILU only; eager fallback)
    }

    // z = M^-1 r  (one V-cycle, fresh zero initial guess; direct QR at the coarsest)
    void apply(const Vector& r, Vector& z) override
    {
        if (static_cast<int>(z.size()) != nd_) z.resize(nd_);
        if (haveGraph_) {
            // captured-graph replay: down-graph -> eager QR coarse -> up-graph, all on capStream_.
            cudaStreamSynchronize(0);                                   // r was produced on the default stream
            cudaMemcpyAsync(thrust::raw_pointer_cast(levels_[0].bvec.data()), r.data(),
                            sizeof(RealType) * nd_, cudaMemcpyDeviceToDevice, capStream_);
            cudaMemsetAsync(thrust::raw_pointer_cast(levels_[0].xvec.data()), 0,
                            sizeof(RealType) * nd_, capStream_);        // fresh zero initial guess
            cudaGraphLaunch(gDownExec_, capStream_);
            solveCoarseLU(capStream_);                                  // factor-once dense LU (replaces per-V-cycle QR)
            cudaGraphLaunch(gUpExec_, capStream_);
            cudaMemcpyAsync(z.data(), thrust::raw_pointer_cast(levels_[0].xvec.data()),
                            sizeof(RealType) * nd_, cudaMemcpyDeviceToDevice, capStream_);
            cudaStreamSynchronize(capStream_);                         // z ready before FlexGMRES consumes it
            return;
        }
        // eager fallback (non-DILU smoother, or capture unavailable): explicit down -> dense-LU coarse -> up
        thrust::copy(thrust::device_pointer_cast(r.data()),
                     thrust::device_pointer_cast(r.data() + nd_), levels_[0].bvec.begin());
        thrust::fill(levels_[0].xvec.begin(), levels_[0].xvec.end(), RealType(0));
        applyCycleEager(0);
        thrust::copy(levels_[0].xvec.begin(), levels_[0].xvec.end(),
                     thrust::device_pointer_cast(z.data()));
    }

    int numLevels() const { return static_cast<int>(levels_.size()); }
    RealType beta() const { return beta_; }
    void setReuse(bool b) { reuse_ = b; }   // true -> skip the next setup() rebuild, reuse the current hierarchy
    double getLastSetupMs() const { return lastSetupMs_; }   // hierarchy-build wall-time of the last setup() (0 if reused)
    bool usingGraph() const { return haveGraph_; }           // true if the V-cycle is replayed from a captured graph

private:
    void destroyGraphs()
    {
        if (gDownExec_) { cudaGraphExecDestroy(gDownExec_); gDownExec_ = nullptr; }
        if (gUpExec_)   { cudaGraphExecDestroy(gUpExec_);   gUpExec_   = nullptr; }
        haveGraph_ = false;
    }

    // Capture the down/up V-cycle halves into CUDA graphs (the cusolver QR coarse solve stays eager
    // between them -- it is non-capturable). Only the DILU smoother is graph-safe (no host buffer swap).
    // ANY capture failure leaves haveGraph_=false -> apply() uses the correct eager path. cudaGraphInstantiate
    // address-binds the graph, so this must re-run on every hierarchy rebuild (called at the end of setup()).
    void captureGraphs()
    {
        destroyGraphs();
        if (!mars::acmUseDilu()) return;                     // Jacobi/ILU paths have a host swap -> not capturable
        // warmup one eager cycle: lazy CUDA/cusolver module+handle init is itself uncapturable on first touch
        applyCycleEager(capStream_);
        cudaStreamSynchronize(capStream_);
        cudaGraph_t g = nullptr;
        bool ok = false;
        do {
            if (cudaStreamBeginCapture(capStream_, cudaStreamCaptureModeThreadLocal) != cudaSuccess) break;
            mars::acmVcycleDownGpu<RealType>(levels_, 0, pre_, omega_, true, capStream_);
            if (cudaStreamEndCapture(capStream_, &g) != cudaSuccess) break;
            if (cudaGraphInstantiateWithFlags(&gDownExec_, g, 0) != cudaSuccess) break;
            cudaGraphDestroy(g); g = nullptr;
            if (cudaStreamBeginCapture(capStream_, cudaStreamCaptureModeThreadLocal) != cudaSuccess) break;
            mars::acmVcycleUpGpu<RealType>(levels_, 0, post_, omega_, true, capStream_);
            if (cudaStreamEndCapture(capStream_, &g) != cudaSuccess) break;
            if (cudaGraphInstantiateWithFlags(&gUpExec_, g, 0) != cudaSuccess) break;
            cudaGraphDestroy(g); g = nullptr;
            ok = true;
        } while (false);
        if (ok) {
            haveGraph_ = true;
            std::fprintf(stderr, "[acm-graph] V-cycle captured (down+up graphs); dense-LU coarse eager between\n");
        } else {
            if (g) cudaGraphDestroy(g);
            destroyGraphs();
            cudaGetLastError();   // clear the sticky capture error so a later driver cudaGetLastError() won't misattribute it
            std::fprintf(stderr, "[acm-graph] V-cycle capture unavailable -> eager path (correct, slower)\n");
        }
    }

    // dense LU of the coarsest operator: factor ONCE per Picard (setup), solve per V-cycle. Replaces the old
    // per-V-cycle cusolverSp QR which re-analyzed/re-factored/re-allocated the same small op ~35x/Picard.
    void factorCoarseLU()
    {
        coarseLUok_ = false;
        mars::AcmLevel<RealType>& C = levels_.back();
        coarseN_ = C.ND;
        coarseLU_.assign((size_t)coarseN_ * coarseN_, RealType(0));      // dense col-major, zeroed
        coarseIpiv_.resize(coarseN_);
        coarseInfo_.resize(1);
        const int blk = 256, g = (coarseN_ + blk - 1) / blk;
        mars::acmDenseFromCsrKernel<RealType><<<g, blk>>>(mars::acmRaw(C.rowOff), mars::acmRaw(C.colInd),
            mars::acmRaw(C.vals), thrust::raw_pointer_cast(coarseLU_.data()), coarseN_);
        RealType* A = thrust::raw_pointer_cast(coarseLU_.data());
        int lwork = 0;
        if constexpr (std::is_same_v<RealType, double>) cusolverDnDgetrf_bufferSize(cd_, coarseN_, coarseN_, A, coarseN_, &lwork);
        else                                            cusolverDnSgetrf_bufferSize(cd_, coarseN_, coarseN_, A, coarseN_, &lwork);
        coarseWork_.resize(lwork);
        RealType* W  = thrust::raw_pointer_cast(coarseWork_.data());
        int* ip = thrust::raw_pointer_cast(coarseIpiv_.data());
        int* inf = thrust::raw_pointer_cast(coarseInfo_.data());
        if constexpr (std::is_same_v<RealType, double>) cusolverDnDgetrf(cd_, coarseN_, coarseN_, A, coarseN_, W, ip, inf);
        else                                            cusolverDnSgetrf(cd_, coarseN_, coarseN_, A, coarseN_, W, ip, inf);
        int info = 0;
        cudaMemcpy(&info, inf, sizeof(int), cudaMemcpyDeviceToHost);     // one-int control readback (setup only)
        coarseLUok_ = (info == 0);
        std::fprintf(stderr, "[acm-coarse] dense LU n=%d -> %s\n", coarseN_,
                     coarseLUok_ ? "factored once (solve per V-cycle)" : "singular -> Jacobi fallback");
    }

    // z = A_c^-1 b on the coarsest via the cached LU (getrs overwrites the RHS). Stream-ordered, no malloc,
    // no re-factor. Singular -> reset + smooth (the old QR fallback). Stays eager between the captured graphs.
    void solveCoarseLU(cudaStream_t s)
    {
        mars::AcmLevel<RealType>& C = levels_.back();
        RealType* x = thrust::raw_pointer_cast(C.xvec.data());
        if (!coarseLUok_) {
            cudaMemsetAsync(x, 0, sizeof(RealType) * coarseN_, s);
            mars::acmSmoothGpu(C, coarse_, omega_, /*useBlock=*/true, s);
            return;
        }
        cudaMemcpyAsync(x, thrust::raw_pointer_cast(C.bvec.data()), sizeof(RealType) * coarseN_,
                        cudaMemcpyDeviceToDevice, s);
        cusolverDnSetStream(cd_, s);
        RealType* A = thrust::raw_pointer_cast(coarseLU_.data());
        int* ip = thrust::raw_pointer_cast(coarseIpiv_.data());
        int* inf = thrust::raw_pointer_cast(coarseInfo_.data());
        if constexpr (std::is_same_v<RealType, double>) cusolverDnDgetrs(cd_, CUBLAS_OP_N, coarseN_, 1, A, coarseN_, ip, x, coarseN_, inf);
        else                                            cusolverDnSgetrs(cd_, CUBLAS_OP_N, coarseN_, 1, A, coarseN_, ip, x, coarseN_, inf);
    }

    // one eager V-cycle on levels_[0].bvec/xvec on stream s: down -> dense-LU coarse -> up.
    void applyCycleEager(cudaStream_t s)
    {
        mars::acmVcycleDownGpu<RealType>(levels_, 0, pre_, omega_, /*useBlock=*/true, s);
        solveCoarseLU(s);
        mars::acmVcycleUpGpu<RealType>(levels_, 0, post_, omega_, /*useBlock=*/true, s);
    }

    int kmax_, maxCoarseND_, pre_, post_, coarse_, nd_ = 0;
    RealType omega_, beta_;
    bool reuse_ = false;
    double lastSetupMs_ = 0;
    std::vector<mars::AcmLevel<RealType>> levels_;
    cusolverSpHandle_t cs_ = nullptr;
    cusparseMatDescr_t descr_ = nullptr;
    cudaStream_t capStream_ = nullptr;
    cudaGraphExec_t gDownExec_ = nullptr, gUpExec_ = nullptr;   // captured V-cycle halves (Step 1)
    bool haveGraph_ = false;
    cusolverDnHandle_t cd_ = nullptr;                          // dense LU for the coarsest (factor once/Picard)
    int coarseN_ = 0;
    bool coarseLUok_ = false;
    thrust::device_vector<RealType> coarseLU_, coarseWork_;    // col-major LU factors + getrf workspace
    thrust::device_vector<int> coarseIpiv_, coarseInfo_;
};

}  // namespace fem
}  // namespace mars

#endif  // MARS_GPU_ACM_PRECONDITIONER_HPP
