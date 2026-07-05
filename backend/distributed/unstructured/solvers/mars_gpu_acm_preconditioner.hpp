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
        if (dist_) {
            mars::acmBuildHierarchyDist<RealType>(ro, ci, va, ND / 4, nOwned0_, ownNode0_, exch0_,
                                                  rank_, nRanks_, maxCoarseND_, levels_, G_, beta_);
            nd_ = ND;
            factorGlobalCoarse();   // replicated dense LU of the GLOBAL coarsest (every rank holds it)
        } else {
            mars::acmBuildHierarchy<RealType>(ro, ci, va, ND / 4, kmax_, maxCoarseND_, levels_, beta_);
            nd_ = ND;
            factorCoarseLU();   // dense LU of the small coarsest op, ONCE (vs the old per-V-cycle cusolverSp QR re-factor)
        }
        cudaDeviceSynchronize();
        lastSetupMs_ = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - tSetup0).count();   // hierarchy build + coarse factor
        if (dist_) destroyGraphs();   // MPI inside the V-cycle is not capturable -> eager path (perf TODO)
        else captureGraphs();   // (re)capture the V-cycle into CUDA graphs on the new hierarchy (DILU only; eager fallback)
    }

    // z = M^-1 r  (one V-cycle, fresh zero initial guess; direct QR at the coarsest)
    void apply(const Vector& r, Vector& z) override
    {
        if (static_cast<int>(z.size()) != nd_) z.resize(nd_);
        if (haveGraph_) {
            // Replay the V-cycle on the NULL stream -- the SAME stream FlexGMRES runs its BLAS/SpMV on. r->bvec, the
            // two captured halves, the eager getrs coarse, and xvec->z are then all stream-ordered with the Krylov
            // work, so the two cudaStreamSynchronize that drained the device every iteration (the host-latency floor)
            // are gone. The graphs were captured on capStream_ but launch on any stream (stream-agnostic at launch).
            cudaMemcpyAsync(thrust::raw_pointer_cast(levels_[0].bvec.data()), r.data(),
                            sizeof(RealType) * nd_, cudaMemcpyDeviceToDevice, 0);
            cudaMemsetAsync(thrust::raw_pointer_cast(levels_[0].xvec.data()), 0, sizeof(RealType) * nd_, 0);
            cudaGraphLaunch(gDownExec_, 0);
            solveCoarseLU(0);                                          // cusolver getrs on stream 0 (eager, no cross-stream drain)
            cudaGraphLaunch(gUpExec_, 0);
            cudaMemcpyAsync(z.data(), thrust::raw_pointer_cast(levels_[0].xvec.data()),
                            sizeof(RealType) * nd_, cudaMemcpyDeviceToDevice, 0);
            return;                                                    // no sync: z is produced + consumed on stream 0
        }
        // eager fallback (non-DILU smoother, or capture unavailable): explicit down -> dense-LU coarse -> up
        thrust::copy(thrust::device_pointer_cast(r.data()),
                     thrust::device_pointer_cast(r.data() + nd_), levels_[0].bvec.begin());
        thrust::fill(levels_[0].xvec.begin(), levels_[0].xvec.end(), RealType(0));
        applyCycleEager(0);
        thrust::copy(levels_[0].xvec.begin(), levels_[0].xvec.end(),
                     thrust::device_pointer_cast(z.data()));
    }

    // multi-rank: seam-consistent distributed hierarchy + replicated global coarsest (mars_acm_dist.hpp).
    // exch = level-0 halo (vec, ncomp): ncomp==1 -> node scalar, ncomp==4 -> interleaved DOF vector.
    void setDistributed(int rank, int numRanks, int nOwnedNodes, const uint8_t* d_ownNode,
                        std::function<void(thrust::device_vector<RealType>&, int)> exch)
    {
        dist_ = (numRanks > 1); rank_ = rank; nRanks_ = numRanks;
        nOwned0_ = nOwnedNodes; ownNode0_ = d_ownNode; exch0_ = std::move(exch);
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
        if (!mars::acmUseDilu()) return;                     // ILU has host-side level loops (not capturable); the
                                                             // Stage-4a jacobi gate wants the eager path. The BJ sweep
                                                             // itself is capturable now (ping-pong) -> hybrid captures.
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

    // Replicated GLOBAL coarsest (dist): the dense matrix was Allgather-assembled identically on every rank
    // (mars_acm_dist.hpp); factor it once per Picard. Same getrf machinery as the local path.
    void factorGlobalCoarse()
    {
        static_assert(std::is_same_v<RealType, double>, "distributed ACM coarsest assembles in double");
        coarseLUok_ = false;
        coarseN_ = G_.nGlobalDof;
        if (coarseN_ == 0) return;
        coarseLU_ = G_.dense;                                            // col-major, replicated
        coarseIpiv_.resize(coarseN_);
        coarseInfo_.resize(1);
        RealType* A = thrust::raw_pointer_cast(coarseLU_.data());
        int lwork = 0;
        cusolverDnDgetrf_bufferSize(cd_, coarseN_, coarseN_, A, coarseN_, &lwork);
        coarseWork_.resize(lwork);
        cusolverDnDgetrf(cd_, coarseN_, coarseN_, A, coarseN_,
                         thrust::raw_pointer_cast(coarseWork_.data()),
                         thrust::raw_pointer_cast(coarseIpiv_.data()),
                         thrust::raw_pointer_cast(coarseInfo_.data()));
        int info = 0;
        cudaMemcpy(&info, thrust::raw_pointer_cast(coarseInfo_.data()), sizeof(int), cudaMemcpyDeviceToHost);
        coarseLUok_ = (info == 0);
        if (rank_ == 0)
            std::fprintf(stderr, "[acm-coarse] global dense LU n=%d -> %s\n", coarseN_,
                         coarseLUok_ ? "factored once (replicated)" : "singular -> Jacobi fallback");
        gRhs_.resize(coarseN_);
        h_rhsLoc_.resize(G_.myDofCount); h_rhsAll_.resize(coarseN_);
    }

    // dist coarsest solve: Allgatherv the owned RHS slices (global-DOF order == rank order by construction),
    // every rank solves the same replicated system, then takes its owned + ghost slices. Host staging is fine
    // here: the coarsest is tiny and this path is eager (no graph at numRanks>1).
    void solveCoarseGlobal(cudaStream_t s)
    {
        mars::AcmLevel<RealType>& C = levels_.back();
        RealType* x = thrust::raw_pointer_cast(C.xvec.data());
        static const bool noCoarse = std::getenv("MARS_ACM_DIST_NOCOARSE") != nullptr;   // debug: pure-Schwarz A/B
        if (noCoarse) {
            cudaMemsetAsync(x, 0, sizeof(RealType) * C.ND, s);
            mars::acmSmoothGpu(C, coarse_, omega_, /*useBlock=*/true, s);
            return;
        }
        if (!coarseLUok_) {
            cudaMemsetAsync(x, 0, sizeof(RealType) * C.ND, s);
            mars::acmSmoothGpu(C, coarse_, omega_, /*useBlock=*/true, s);
            return;
        }
        cudaStreamSynchronize(s);
        cudaMemcpy(h_rhsLoc_.data(), thrust::raw_pointer_cast(C.bvec.data()),
                   sizeof(RealType) * G_.myDofCount, cudaMemcpyDeviceToHost);   // owned coarse rows are [0, 4*nOwned)
        MPI_Allgatherv(h_rhsLoc_.data(), G_.myDofCount, MPI_DOUBLE,
                       h_rhsAll_.data(), G_.dofCounts.data(), G_.dofDispls.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        thrust::copy(h_rhsAll_.begin(), h_rhsAll_.end(), gRhs_.begin());
        cusolverDnSetStream(cd_, s);
        cusolverDnDgetrs(cd_, CUBLAS_OP_N, coarseN_, 1,
                         thrust::raw_pointer_cast(coarseLU_.data()), coarseN_,
                         thrust::raw_pointer_cast(coarseIpiv_.data()),
                         thrust::raw_pointer_cast(gRhs_.data()), coarseN_,
                         thrust::raw_pointer_cast(coarseInfo_.data()));
        cudaStreamSynchronize(s);
        // owned slice (contiguous global range), then the ghost ring from the replicated solution
        cudaMemcpy(x, thrust::raw_pointer_cast(gRhs_.data()) + G_.myDofBegin,
                   sizeof(RealType) * G_.myDofCount, cudaMemcpyDeviceToDevice);
        if (C.nNodes > C.nOwned) {
            thrust::copy(gRhs_.begin(), gRhs_.end(), h_rhsAll_.begin());
            std::vector<RealType> h_g(4 * (C.nNodes - C.nOwned));
            for (int k = 0; k < C.nNodes - C.nOwned; ++k)
                for (int c = 0; c < 4; ++c) h_g[4 * k + c] = h_rhsAll_[4 * C.ghostAggGid[k] + c];
            cudaMemcpy(x + 4 * C.nOwned, h_g.data(), sizeof(RealType) * h_g.size(), cudaMemcpyHostToDevice);
        }
        // MARS_ACM_DIST_PROBE=1: does the replicated solution actually solve the LEVEL operator?
        // ||b - A_level x|| / ||b|| over OWNED coarse rows, global. ~1e-12 -> the gathered dense matches the
        // level CSR and the bug is downstream (inject); large -> assembly/ordering bug in the gather.
        static const bool probe = std::getenv("MARS_ACM_DIST_PROBE") != nullptr;
        if (probe) {
            const int blk = 256, grd = (C.ND + blk - 1) / blk;
            mars::acmResidualKernel<RealType><<<grd, blk>>>(mars::acmRaw(C.rowOff), mars::acmRaw(C.colInd),
                mars::acmRaw(C.vals), x, thrust::raw_pointer_cast(C.bvec.data()),
                thrust::raw_pointer_cast(C.rtmp.data()), C.ND);
            cudaDeviceSynchronize();
            std::vector<RealType> h_r(C.ND), h_b(C.ND);
            thrust::copy(C.rtmp.begin(), C.rtmp.end(), h_r.begin());
            thrust::copy(C.bvec.begin(), C.bvec.end(), h_b.begin());
            double rr = 0, bb = 0;
            for (int d = 0; d < 4 * C.nOwned; ++d) { rr += (double)h_r[d] * h_r[d]; bb += (double)h_b[d] * h_b[d]; }
            double g[2] = {rr, bb};
            MPI_Allreduce(MPI_IN_PLACE, g, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (rank_ == 0)
                std::fprintf(stderr, "[dist-probe] coarse |b-Ax|/|b| = %.3e (n=%d)\n",
                             std::sqrt(g[0] / (g[1] > 0 ? g[1] : 1.0)), coarseN_);
        }
    }

    // z = A_c^-1 b on the coarsest via the cached LU (getrs overwrites the RHS). Stream-ordered, no malloc,
    // no re-factor. Singular -> reset + smooth (the old QR fallback). Stays eager between the captured graphs.
    void solveCoarseLU(cudaStream_t s)
    {
        if (dist_) { solveCoarseGlobal(s); return; }
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
    // ---- distributed (numRanks > 1) ----
    bool dist_ = false;
    int rank_ = 0, nRanks_ = 1, nOwned0_ = 0;
    const uint8_t* ownNode0_ = nullptr;                        // level-0 node ownership (device)
    std::function<void(thrust::device_vector<RealType>&, int)> exch0_;   // level-0 halo (vec, ncomp)
    mars::AcmGlobalCoarse G_;                                  // replicated global coarsest
    thrust::device_vector<RealType> gRhs_;                     // global coarse RHS/solution (device)
    std::vector<RealType> h_rhsLoc_, h_rhsAll_;                // host staging for the Allgatherv
};

}  // namespace fem
}  // namespace mars

#endif  // MARS_GPU_ACM_PRECONDITIONER_HPP
