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
#include <chrono>

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
        cusolverSpCreate(&cs_);
        cusparseCreateMatDescr(&descr_);
        cusparseSetMatType(descr_, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr_, CUSPARSE_INDEX_BASE_ZERO);
    }

    ~GpuAcmPreconditioner()
    {
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
        cudaDeviceSynchronize();
        lastSetupMs_ = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - tSetup0).count();   // hierarchy build (host-heavy for ILU/DILU factor)
    }

    // z = M^-1 r  (one V-cycle, fresh zero initial guess; direct QR at the coarsest)
    void apply(const Vector& r, Vector& z) override
    {
        if (static_cast<int>(z.size()) != nd_) z.resize(nd_);
        thrust::copy(thrust::device_pointer_cast(r.data()),
                     thrust::device_pointer_cast(r.data() + nd_), levels_[0].bvec.begin());
        thrust::fill(levels_[0].xvec.begin(), levels_[0].xvec.end(), RealType(0));
        mars::acmVcycleGpu<RealType>(levels_, 0, pre_, post_, coarse_, omega_, /*useBlock=*/true, cs_, descr_);
        thrust::copy(levels_[0].xvec.begin(), levels_[0].xvec.end(),
                     thrust::device_pointer_cast(z.data()));
    }

    int numLevels() const { return static_cast<int>(levels_.size()); }
    RealType beta() const { return beta_; }
    void setReuse(bool b) { reuse_ = b; }   // true -> skip the next setup() rebuild, reuse the current hierarchy
    double getLastSetupMs() const { return lastSetupMs_; }   // hierarchy-build wall-time of the last setup() (0 if reused)

private:
    int kmax_, maxCoarseND_, pre_, post_, coarse_, nd_ = 0;
    RealType omega_, beta_;
    bool reuse_ = false;
    double lastSetupMs_ = 0;
    std::vector<mars::AcmLevel<RealType>> levels_;
    cusolverSpHandle_t cs_ = nullptr;
    cusparseMatDescr_t descr_ = nullptr;
};

}  // namespace fem
}  // namespace mars

#endif  // MARS_GPU_ACM_PRECONDITIONER_HPP
