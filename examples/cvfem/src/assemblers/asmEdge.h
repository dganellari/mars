#ifndef ASM_EDGE_H
#define ASM_EDGE_H

#include "asmBase.h"

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
class AsmEdge : public AsmBase<scalar,
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

    static constexpr unsigned SCRATCH_SPACE_LEVEL = 1;

    // maximum sizes for edge data
    static constexpr unsigned nodesPerEdge = 2;
    static constexpr unsigned numScsIpMax = 1;

    // block matrices not currently supported for edge
    static_assert(BLOCKSIZE == 1);

    // make sure constructor of base class is called
    AsmEdge(const stk::mesh::BulkData& bulk)
        : AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>(
              bulk)
    {
    }

    virtual std::string name() const override
    {
        return std::string("Edge");
    }

    void assemble(Matrix& A,
                  Vector& b,
                  const stk::mesh::MetaData& metaData,
                  const stk::mesh::BulkData& bulkData) const override;

    // helper methods
    void getScratchSizes(size_t& perTeam, size_t& perThread) const override
    {
        constexpr int maxBucketSize = 512;

        // scalar array sizes
        const label lhsSize =
            nodesPerEdge * BLOCKSIZE * nodesPerEdge * BLOCKSIZE;
        const label rhsSize = nodesPerEdge * BLOCKSIZE;

        const size_t lhsScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, lhsSize);
        const size_t rhsScratchSize =
            TeamScratchViewScalar::shmem_size(maxBucketSize, rhsSize);

        perTeam = lhsScratchSize + rhsScratchSize;
        perThread = 0;
    }

    KOKKOS_FUNCTION
    void applyCoeff_(const Matrix& A,
                     const Vector& b,
                     const stk::mesh::Entity connectedNodes[2],
                     const TeamScratchViewScalar& rhs,
                     const TeamScratchViewScalar& lhs,
                     const int& iEdge,
                     bool deductUnfound = true) const;
};

#endif // ASM_EDGE_H
