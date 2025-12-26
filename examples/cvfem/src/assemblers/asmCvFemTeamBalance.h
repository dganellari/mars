#ifndef ASM_CVFEM_TEAM_BALANCE_H
#define ASM_CVFEM_TEAM_BALANCE_H

#include "asmBase.h"
#include "FEhex.h"
// #define USE_INV_JACOBIAN

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          typename ElemType,
          bool includeAdv = true,
          bool isShifted = true>
class AsmCvFemTeamBalance : public AsmBase<scalar,
                                         label,
                                         BLOCKSIZE,
                                         SPATIAL_DIM,
                                         includeAdv,
                                         isShifted>
{
public:
    // Type aliases from base class
    using Context = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Context;
    using Matrix = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Matrix;
    using Vector = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Vector;
    using ExecSpace = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ExecSpace;
    
    // Scratch memory view types
    using ScratchViewScalar = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewScalar;
    using ScratchViewScalar2 = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewScalar2;
    using TeamScratchViewScalar = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewScalar2;
    using TeamScratchViewScalar3 = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewScalar3;
    
    using ScratchViewLabel = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewLabel;
    using ScratchViewEntity = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewEntity;
    using ScratchViewIndex = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchViewIndex;
    using ScratchView2D = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView2D;
    using ScratchView3D = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView3D;
    using TeamHandleType = typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::TeamHandleType;

    // Shape function and derivative view types
    using HexDerivView = typename ElemType::HexDerivView;
    using HexDerivViewHost = typename ElemType::HexDerivViewHost;
    using ConstHexDerivView = typename ElemType::ConstHexDerivView;
    using SFView = typename ElemType::SFView;
    using SFViewHost = typename ElemType::SFViewHost;
    using ConstSFView = typename ElemType::ConstSFView;

    // Constants
    static constexpr unsigned SCRATCH_SPACE_LEVEL = 1;

    // Constructor
    AsmCvFemTeamBalance(const stk::mesh::BulkData& bulk)
        : AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>(bulk)
    {
        // Initialize pre-computed shape function and derivative views
        hexDerivView_ = HexDerivView("hexDerivView");
        hexDerivViewHost_ = Kokkos::create_mirror_view(hexDerivView_);

        sfView_ = SFView("sfView");
        sfViewHost_ = Kokkos::create_mirror_view(sfView_);
        
        // Pre-compute shape functions and derivatives
        ElemType::preCalcShapeFunction(sfViewHost_);
        Kokkos::deep_copy(sfView_, sfViewHost_);
        sfViewConst_ = sfView_;

        ElemType::preCalcDeriv(hexDerivViewHost_);
        Kokkos::deep_copy(hexDerivView_, hexDerivViewHost_);
        hexDerivViewConst_ = hexDerivView_;
    }

    // Interface methods
    virtual std::string name() const override
    {
        return std::string("CvFemTeamBalanceOptimized");
    }

    void assemble(Matrix& A,
                  Vector& b,
                  const stk::mesh::MetaData& metaData,
                  const stk::mesh::BulkData& bulkData) const override;

    // helper methods
    void getScratchSizes(size_t& perTeam, size_t& perThread) const override
    {
        constexpr int maxBucketSize = 512;
        // stk::mesh::impl::default_maximum_bucket_capacity;
        size_t systemSize = ElemType::nodesPerElement * BLOCKSIZE;

        // scalar array sizes
        const label lhsSize = systemSize * systemSize;
        const label rhsSize = systemSize;

        const label phiSize = systemSize;
        const label betaSize = systemSize;
        const label gradPhiSize = systemSize * SPATIAL_DIM;
        const label gammaSize = systemSize;

        const size_t lhsScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, lhsSize);
        const size_t rhsScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, rhsSize);

        const size_t phiScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, phiSize);
        const size_t betaScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, betaSize);
        const size_t gradPhiScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, gradPhiSize);
        const size_t gammaScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, gammaSize);

        // DoubleType array sizes (for master element data)
        const size_t coordinatesScratchSize =
            TeamScratchViewScalar3::shmem_size(
                maxBucketSize, ElemType::nodesPerElement, SPATIAL_DIM);
        
        // scratch size required for each thread
        const size_t totalThreadScratchSize = 0;

        // scratch size required for each team
        const size_t totalTeamScratchSize =
            lhsScratchSize + rhsScratchSize + phiScratchSize + betaScratchSize +
            gradPhiScratchSize + gammaScratchSize + coordinatesScratchSize;


        perTeam = totalTeamScratchSize;
        perThread = totalThreadScratchSize;
    }              
    
    // Utility function for sorting matrix column indices
    template <size_t numVals>
    KOKKOS_INLINE_FUNCTION void performSortPermutations(
        label arrayA[numVals],
        const label arrayB[numVals]) const
    {
        for (size_t i = 0; i < numVals; ++i) {
            arrayA[i] = i;
        }

        for (size_t i = 1; i < numVals; ++i) {
            size_t key = arrayA[i];
            int keyVal = arrayB[key];
            size_t j = i;

            while (j > 0 && arrayB[arrayA[j - 1]] > keyVal) {
                arrayA[j] = arrayA[j - 1];
                --j;
            }
            arrayA[j] = key;
        }
    }

private:
    // Matrix/vector assembly function
    KOKKOS_FUNCTION
    void applyCoeff_(
        const Matrix& A,
        const Vector& b,
        const stk::mesh::Entity connectedNodes[ElemType::nodesPerElement],
        const scalar* rhs,
        const scalar* lhs,
        bool deductUnfound = true) const;

    // Pre-computed shape function and derivative views
    HexDerivView hexDerivView_;
    HexDerivViewHost hexDerivViewHost_;
    ConstHexDerivView hexDerivViewConst_;

    SFView sfView_;
    SFViewHost sfViewHost_;
    ConstSFView sfViewConst_;
};

#endif // ASM_CVFEM_TEAM_BALANCE_H