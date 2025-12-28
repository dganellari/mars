#ifndef ASM_CVFEM_TEAM_SCRATCH_H
#define ASM_CVFEM_TEAM_SCRATCH_H

#include "asmBase.h"

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
class AsmCvFemTeamScratch : public AsmBase<scalar,
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
    // maximum sizes for element data (hex)
    static const int nodesPerElementMax = 8;
    static const int numScsIpMax = 12;
    static const int numScvIpMax = 8;

    AsmCvFemTeamScratch(const stk::mesh::BulkData& bulk)
        : AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>(
              bulk)
    {
        const stk::mesh::NgpMesh& ngpMesh =
            stk::mesh::get_updated_ngp_mesh(bulk);

        const stk::mesh::MetaData& metaData = bulk.mesh_meta_data();

        const stk::mesh::Selector selAllElements = metaData.universal_part();

        stk::NgpVector<unsigned> bucketIds =
            ngpMesh.get_bucket_ids(stk::topology::ELEMENT_RANK, selAllElements);
        unsigned numBuckets = bucketIds.size();

        std::set<stk::topology> all_topos;
        for (unsigned i = 0; i < numBuckets; ++i)
        {
            const stk::mesh::NgpMesh::BucketType& myelementBucket =
                ngpMesh.get_bucket(stk::topology::ELEM_RANK, i);
            all_topos.insert(myelementBucket.topology());
        }

        for (const stk::topology& topo : all_topos)
        {
            label topo_index = topo.value();
            meSCS_ptrs_[topo_index] =
                accel::MasterElementRepo::get_surface_master_element_on_dev(
                    topo);
        }
    }

    virtual std::string name() const override
    {
        return std::string("CvFemTeamScratch");
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
        size_t systemSize = nodesPerElementMax * BLOCKSIZE;

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
            ScratchView2D::shmem_size(nodesPerElementMax, SPATIAL_DIM);
        const size_t scsAreaVScratchSize =
            ScratchView2D::shmem_size(numScsIpMax, SPATIAL_DIM);
        const size_t dndxScratchSize = ScratchView3D::shmem_size(
            numScsIpMax, nodesPerElementMax, SPATIAL_DIM);
        const size_t derivScratchSize = ScratchView3D::shmem_size(
            numScsIpMax, nodesPerElementMax, SPATIAL_DIM);
        const size_t shapeFunctionScratchSize =
            ScratchView2D::shmem_size(numScsIpMax, nodesPerElementMax);

        // scratch size required for each thread
        const size_t totalThreadScratchSize =
            (coordinatesScratchSize + scsAreaVScratchSize + dndxScratchSize +
             derivScratchSize);

        // scratch size required for each team
        const size_t totalTeamScratchSize =
            lhsScratchSize + rhsScratchSize + phiScratchSize + betaScratchSize +
            gradPhiScratchSize + gammaScratchSize + shapeFunctionScratchSize;

        perTeam = totalTeamScratchSize;
        perThread = totalThreadScratchSize;
    }

    KOKKOS_INLINE_FUNCTION
    void performSortPermutations(
        label arrayA[nodesPerElementMax], // sortPermutations
        label arrayB[nodesPerElementMax], // matrixColumnIds
        size_t numVals) const
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

    KOKKOS_FUNCTION
    void applyCoeff_(const Matrix& A,
                     const Vector& b,
                     const stk::mesh::Entity connectedNodes[nodesPerElementMax],
                     const TeamScratchViewScalar& rhs,
                     const TeamScratchViewScalar& lhs,
                     const int& iElement,
                     label nodesPerElement,
                     bool deductUnfound = true) const;

private:
    accel::MasterElement* meSCS_ptrs_[stk::topology::NUM_TOPOLOGIES] = {
        nullptr};
};

#endif // ASM_CVFEM_TEAM_SCRATCH_H
