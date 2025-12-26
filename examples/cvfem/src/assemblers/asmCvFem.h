#ifndef ASM_CVFEM_H
#define ASM_CVFEM_H

#include "asmBase.h"

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
class AsmCvFem : public AsmBase<scalar,
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

    AsmCvFem(const stk::mesh::BulkData& bulk)
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
        return std::string("CvFem");
    }

    void assemble(Matrix& A,
                  Vector& b,
                  const stk::mesh::MetaData& metaData,
                  const stk::mesh::BulkData& bulkData) const override;

    // helper methods
    void getScratchSizes(size_t& perTeam, size_t& perThread) const override
    {
        // scalar array sizes
        const label lhsSize =
            nodesPerElementMax * BLOCKSIZE * nodesPerElementMax * BLOCKSIZE;
        const label rhsSize = nodesPerElementMax * BLOCKSIZE;
        const label scratchValsSize = rhsSize;

        const label phiSize = nodesPerElementMax * BLOCKSIZE;
        const label betaSize = nodesPerElementMax * BLOCKSIZE;
        const label gradPhiSize = nodesPerElementMax * BLOCKSIZE * SPATIAL_DIM;
        const label gammaSize = nodesPerElementMax;

        // label array sizes
        const label scratchIdsSize = rhsSize;

        // entity array sizes
        const label connectedNodesSize = nodesPerElementMax;

        // totals
        const label totalScalarSize = lhsSize + rhsSize + scratchValsSize +
                                      phiSize + betaSize + gradPhiSize +
                                      gammaSize;

        const label totalLabelSize = scratchIdsSize;
        const label totalEntitySize = connectedNodesSize;

        const label totalIndexSize = 2 * connectedNodesSize;

        const size_t scalarScratchSize =
            ScratchViewScalar::shmem_size(totalScalarSize);
        const size_t labelScratchSize =
            ScratchViewLabel::shmem_size(totalLabelSize);
        const size_t entityScratchSize =
            ScratchViewEntity::shmem_size(totalEntitySize);
        const size_t indexScratchSize =
            ScratchViewIndex::shmem_size(totalIndexSize);

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
            (scalarScratchSize + labelScratchSize + entityScratchSize +
             indexScratchSize + coordinatesScratchSize + scsAreaVScratchSize +
             dndxScratchSize + derivScratchSize);

        // scratch size required for each team
        const size_t totalTeamScratchSize = shapeFunctionScratchSize;

        perTeam = totalTeamScratchSize;
        perThread = totalThreadScratchSize;
    }

    KOKKOS_INLINE_FUNCTION
    void
    performSortPermutations(const ScratchViewIndex& arrayA, // sortPermutations
                            const ScratchViewIndex& arrayB, // matrixColumnIds
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
            int keyVal = arrayB(key);
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
                     const ScratchViewEntity& connectedNodes,
                     const ScratchViewLabel& scratchIds,
                     const ScratchViewIndex& matrixColumnIds,
                     const ScratchViewIndex& sortPermutation,
                     const ScratchViewScalar& rhs,
                     const ScratchViewScalar& lhs,
                     label nodesPerElement,
                     bool deductUnfound = true) const;

private:
    accel::MasterElement* meSCS_ptrs_[stk::topology::NUM_TOPOLOGIES] = {
        nullptr};
};

#endif // ASM_CVFEM_H
