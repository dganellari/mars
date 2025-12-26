#ifndef ASM_CVFEM_INLINE_SF_H
#define ASM_CVFEM_INLINE_SF_H

#include "asmBase.h"
#include <Eigen/Core>
#include <Eigen/Dense>

#include "FEhex.h"
// if defined, the gradient weights are calculated per ip, per node
// if not defined, weights are calculated per ip, for all nodes.
// #define USE_INV_JACOBIAN

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          typename ElemType,
          bool includeAdv = true,
          bool isShifted = true>
class AsmCvFemInlineSF : public AsmBase<scalar,
                                        label,
                                        BLOCKSIZE,
                                        SPATIAL_DIM,
                                        includeAdv,
                                        isShifted>
{
public:
    using Context = AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Context;
    using Matrix =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Matrix;
    using Vector =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Vector;
    using ExecSpace =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ExecSpace;
    using ScratchViewScalar =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewScalar;
    using ScratchViewScalar2 =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewScalar2;
    using TeamScratchViewScalar =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewScalar2;
    using TeamScratchViewScalar3 =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewScalar3;

    using ScratchViewLabel =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewLabel;
    using ScratchViewEntity =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewEntity;
    using ScratchViewIndex =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewIndex;
    using ScratchView2D =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView2D;
    using ScratchView3D =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView3D;
    using TeamHandleType =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::TeamHandleType;

    // views for shape function data
    using HexDerivView = typename ElemType::HexDerivView;
    using HexDerivViewHost = typename ElemType::HexDerivViewHost;
    using ConstHexDerivView = typename ElemType::ConstHexDerivView;
    using SFView = typename ElemType::SFView;
    using SFViewHost = typename ElemType::SFViewHost;
    using ConstSFView = typename ElemType::ConstSFView;

    static constexpr unsigned SCRATCH_SPACE_LEVEL = 1;

    // using AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv,
    // isShifted>::AsmBase; // inherit all
    //                                                      // constructors

    // make sure constructor of base class is called
    AsmCvFemInlineSF(const stk::mesh::BulkData& bulk)
        : AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>(
              bulk)
    {
        hexDerivView_ = HexDerivView("hexDerivView");
        hexDerivViewHost_ = Kokkos::create_mirror_view(hexDerivView_);

        sfView_ = SFView("sfView");
        sfViewHost_ = Kokkos::create_mirror_view(sfView_);
        ElemType::preCalcShapeFunction(sfViewHost_);
        Kokkos::deep_copy(sfView_, sfViewHost_);
        sfViewConst_ = sfView_;

        ElemType::preCalcDeriv(hexDerivViewHost_);
        Kokkos::deep_copy(hexDerivView_, hexDerivViewHost_);
        hexDerivViewConst_ = hexDerivView_;
    }

    virtual std::string name() const override
    {
        return std::string("CvFemInlineShapeFunction");
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

    template <size_t numVals>
    KOKKOS_INLINE_FUNCTION void performSortPermutations(
        label arrayA[numVals],             // sortPermutations
        const label arrayB[numVals]) const // matrixColumnIds
    {
        for (size_t i = 0; i < numVals; ++i)
        {
            arrayA[i] = i;
        }

        // Insertion sort arrayA based on values in arrayB
        for (size_t i = 1; i < numVals; ++i)
        {
            size_t key = arrayA[i];
            int keyVal = arrayB[key];
            size_t j = i;

            while (j > 0 && arrayB[arrayA[j - 1]] > keyVal)
            {
                arrayA[j] = arrayA[j - 1];
                --j;
            }
            arrayA[j] = key;
        }
    }

    // helper methods
    KOKKOS_FUNCTION
    void applyCoeff_(
        const Matrix& A,
        const Vector& b,
        const stk::mesh::Entity connectedNodes[ElemType::nodesPerElement],
        const TeamScratchViewScalar& rhs,
        const TeamScratchViewScalar& lhs,
        const int& iElement,
        bool deductUnfound = true) const;

    HexDerivView hexDerivView_;
    HexDerivViewHost hexDerivViewHost_;
    ConstHexDerivView hexDerivViewConst_;

    SFView sfView_;
    SFViewHost sfViewHost_;
    ConstSFView sfViewConst_;
};

#endif // ASM_CVFEM_INLINE_SF_H
