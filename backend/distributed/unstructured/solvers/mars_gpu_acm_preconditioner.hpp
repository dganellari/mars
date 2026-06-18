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

    GpuAcmPreconditioner(int kmax = 4, int maxCoarseND = 64,
                         int pre = 2, int post = 1, int coarse = 10, RealType omega = RealType(0.7))
        : kmax_(kmax), maxCoarseND_(maxCoarseND), pre_(pre), post_(post), coarse_(coarse), omega_(omega) {}

    // build the multilevel hierarchy once from the finest coupled CSR (interleaved 4*node+comp)
    void setup(const Matrix& A) override
    {
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
        mars::acmBuildHierarchy<RealType>(ro, ci, va, ND / 4, kmax_, maxCoarseND_, levels_);
        nd_ = ND;
    }

    // z = M^-1 r  (one V-cycle, fresh zero initial guess)
    void apply(const Vector& r, Vector& z) override
    {
        if (static_cast<int>(z.size()) != nd_) z.resize(nd_);
        thrust::copy(thrust::device_pointer_cast(r.data()),
                     thrust::device_pointer_cast(r.data() + nd_), levels_[0].bvec.begin());
        thrust::fill(levels_[0].xvec.begin(), levels_[0].xvec.end(), RealType(0));
        mars::acmVcycleGpu<RealType>(levels_, 0, pre_, post_, coarse_, omega_);
        thrust::copy(levels_[0].xvec.begin(), levels_[0].xvec.end(),
                     thrust::device_pointer_cast(z.data()));
    }

    int numLevels() const { return static_cast<int>(levels_.size()); }

private:
    int kmax_, maxCoarseND_, pre_, post_, coarse_, nd_ = 0;
    RealType omega_;
    std::vector<mars::AcmLevel<RealType>> levels_;
};

}  // namespace fem
}  // namespace mars

#endif  // MARS_GPU_ACM_PRECONDITIONER_HPP
