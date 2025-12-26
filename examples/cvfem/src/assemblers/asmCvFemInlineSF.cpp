#include "asmCvFemInlineSF.h"
#include <iostream>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          typename ElemType,
          bool includeAdv,
          bool isShifted>
void AsmCvFemInlineSF<scalar,
                      label,
                      BLOCKSIZE,
                      SPATIAL_DIM,
                      ElemType,
                      includeAdv,
                      isShifted>::assemble(Matrix& A,
                                           Vector& b,
                                           const stk::mesh::MetaData& metaData,
                                           const stk::mesh::BulkData& bulkData)
    const
{
    // get fields references
    const stk::mesh::Field<scalar>& GammaSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "gamma");
    const stk::mesh::Field<scalar>& mdotSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::ELEMENT_RANK, "mdot");
    const stk::mesh::Field<scalar>& phiSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phi");
    const stk::mesh::Field<scalar>& gradPhiSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phiGrad");
    const stk::mesh::Field<scalar>& betaSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "beta");
    const stk::mesh::Field<scalar>& areaVSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::ELEMENT_RANK, "areaV");
    const stk::mesh::Field<scalar>& coordinatesRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "coordinates");

    // init ngp data
    const stk::mesh::NgpMesh& ngpMesh =
        stk::mesh::get_updated_ngp_mesh(bulkData);
    stk::mesh::NgpField<scalar>& ngpGammaSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(GammaSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpMdotSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(mdotSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpPhiSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(phiSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpGradPhiSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(gradPhiSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpBetaSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(betaSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpAreaVSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(areaVSTKFieldRef);
    stk::mesh::NgpField<double>& ngpCoordinatesRef =
        stk::mesh::get_updated_ngp_field<double>(coordinatesRef);

    const stk::mesh::Selector selAllElements = metaData.universal_part();

    stk::NgpVector<unsigned> bucketIds =
        ngpMesh.get_bucket_ids(stk::topology::ELEMENT_RANK, selAllElements);
    unsigned numBuckets = bucketIds.size();

    auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(numBuckets, Kokkos::AUTO);

    size_t totalTeamScratchSize;
    size_t totalThreadScratchSize;

    getScratchSizes(totalTeamScratchSize, totalThreadScratchSize);

    size_t systemSize = ElemType::nodesPerElement * BLOCKSIZE;

    Kokkos::parallel_for(
        "asm_cvfem_inline",
        teamPolicy.set_scratch_size(SCRATCH_SPACE_LEVEL,
                                    Kokkos::PerTeam(totalTeamScratchSize),
                                    Kokkos::PerThread(totalThreadScratchSize)),
        KOKKOS_CLASS_LAMBDA( // TODO: class lambda here for numNodes_
            const TeamHandleType& teamMember) {
            const int bucketIndex =
                bucketIds.get<ExecSpace>(teamMember.league_rank());
            const stk::mesh::NgpMesh::BucketType& elementBucket =
                ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
            unsigned nElementsPerBucket = elementBucket.size();

            TeamScratchViewScalar lhs(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize * systemSize);
            TeamScratchViewScalar rhs(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize);

            // nodal fields to gather
            TeamScratchViewScalar ws_phi(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize);
            TeamScratchViewScalar ws_beta(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize);
            TeamScratchViewScalar ws_gradPhi(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize * SPATIAL_DIM);
            TeamScratchViewScalar ws_Gamma(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize);

            // geometry related to populate
            TeamScratchViewScalar3 ws_coordinates(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                ElemType::nodesPerElement,
                SPATIAL_DIM);

            const label lhsSize = systemSize *
                                  systemSize;
            const label rhsSize = systemSize;

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(teamMember, 0u, nElementsPerBucket),
                [&](const int& iElement)
            {
                // ip values
                scalar phiIp[BLOCKSIZE];
                scalar coordIp[SPATIAL_DIM];

                // extrapolated value from the L/R direction
                scalar phiIpL[BLOCKSIZE];
                scalar phiIpR[BLOCKSIZE];

                // get elem
                stk::mesh::Entity elem = elementBucket[iElement];
                stk::mesh::FastMeshIndex elemFastIndex =
                    ngpMesh.fast_mesh_index(elem);

                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                {
                    lhs(iElement, p) = 0.0;
                }
                for (label p = 0; p < rhsSize; ++p)
                {
                    rhs(iElement, p) = 0.0;
                }

                //===================
                // gather nodal data
                //===================
                stk::mesh::NgpMesh::ConnectedNodes nodeRels =
                    ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemFastIndex);
                label numNodes = nodeRels.size();

                // sanity check on num nodes
                STK_NGP_ThrowAssert(numNodes == ElemType::nodesPerElement);

                stk::mesh::Entity connectedNodes[ElemType::nodesPerElement];

                for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    stk::mesh::FastMeshIndex nodeFastIndex =
                        ngpMesh.fast_mesh_index(node);

                    // gather scalars
                    ws_Gamma(iElement, ni) =
                        ngpGammaSTKFieldRef(nodeFastIndex, 0);

                    // gather vectors
                    for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                    {
                        ws_coordinates(iElement, ni, i) =
                            ngpCoordinatesRef(nodeFastIndex, i);
                    }

                    // gather BLOCKSIZE-dim fields
                    for (unsigned i = 0; i < BLOCKSIZE; ++i)
                    {
                        ws_phi(iElement, ni * BLOCKSIZE + i) =
                            ngpPhiSTKFieldRef(nodeFastIndex, i);
                        ws_beta(iElement, ni * BLOCKSIZE + i) =
                            ngpBetaSTKFieldRef(nodeFastIndex, i);

                        for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                        {
                            ws_gradPhi(iElement,
                                       ni * BLOCKSIZE * SPATIAL_DIM +
                                           i * SPATIAL_DIM + j) =
                                ngpGradPhiSTKFieldRef(nodeFastIndex,
                                                      i * SPATIAL_DIM + j);
                        }
                    }
                }

                // TODO: move to ElemType
                constexpr int lrscv[24] = {0, 1, 1, 2, 2, 3, 0, 3, 4, 5, 5, 6,
                                           6, 7, 4, 7, 0, 4, 1, 5, 2, 6, 3, 7};

                for (unsigned ip = 0; ip < ElemType::numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    // save off mdot
                    const scalar tmdot =
                        includeAdv ? ngpMdotSTKFieldRef(elemFastIndex, ip)
                                   : 0.0;

                    // zero out values of interest for this ip
                    for (unsigned j = 0; j < BLOCKSIZE; ++j)
                    {
                        phiIp[j] = 0.0;
                    }
                    for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                    {
                        coordIp[j] = 0.0;
                    }

                    // save off ip values; offset to Shape Function
                    scalar GammaIp = 0.0;

                    for (unsigned ic = 0; ic < ElemType::nodesPerElement; ++ic)
                    {
                        const scalar r = sfViewConst_(ip, ic);

                        GammaIp += r * ws_Gamma(iElement, ic);

                        // compute scs ip value
                        for (unsigned i = 0; i < BLOCKSIZE; ++i)
                        {
                            phiIp[i] +=
                                r * ws_phi(iElement, ic * BLOCKSIZE + i);
                        }
                        for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                        {
                            coordIp[i] += r * ws_coordinates(iElement, ic, i);
                        }
                    }

                    //===========
                    // Advection
                    //===========

                    // final upwind extrapolation;
                    for (unsigned i = 0; i < BLOCKSIZE; ++i)
                    {
                        phiIpL[i] = ws_phi(iElement, il * BLOCKSIZE + i);
                        phiIpR[i] = ws_phi(iElement, ir * BLOCKSIZE + i);
                    }

                    // assemble advection; rhs and upwind contributions
                    for (unsigned i = 0; i < BLOCKSIZE; ++i)
                    {
                        scalar phiUpwind;
                        scalar dcorr = 0;
                        if (tmdot > 0)
                        {
                            phiUpwind = phiIpL[i];

                            // deferred correction
                            for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar dxj =
                                    coordIp[j] -
                                    ws_coordinates(iElement, il, j);
                                dcorr +=
                                    ws_beta(iElement, il * BLOCKSIZE + i) *
                                    dxj *
                                    ws_gradPhi(iElement,
                                               il * BLOCKSIZE * SPATIAL_DIM +
                                                   i * SPATIAL_DIM + j);
                            }
                        }
                        else
                        {
                            phiUpwind = phiIpR[i];

                            // deferred correction
                            for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar dxj =
                                    coordIp[j] -
                                    ws_coordinates(iElement, ir, j);
                                dcorr +=
                                    ws_beta(iElement, ir * BLOCKSIZE + i) *
                                    dxj *
                                    ws_gradPhi(iElement,
                                               ir * BLOCKSIZE * SPATIAL_DIM +
                                                   i * SPATIAL_DIM + j);
                            }
                        }

                        // total upwind advection
                        const scalar aflux = tmdot * (phiUpwind + dcorr);

                        const label indexL = il * BLOCKSIZE + i;
                        const label indexR = ir * BLOCKSIZE + i;

                        const label rowL =
                            indexL * ElemType::nodesPerElement * BLOCKSIZE;
                        const label rowR =
                            indexR * ElemType::nodesPerElement * BLOCKSIZE;

                        const label rLiL_i = rowL + il * BLOCKSIZE + i;
                        const label rLiR_i = rowL + ir * BLOCKSIZE + i;
                        const label rRiR_i = rowR + ir * BLOCKSIZE + i;
                        const label rRiL_i = rowR + il * BLOCKSIZE + i;

                        // right hand side; L and R
                        rhs(iElement, indexL) -= aflux;
                        rhs(iElement, indexR) += aflux;

                        // upwind advection left node
                        const scalar alhsfacL = 0.5 * (tmdot + std::abs(tmdot));
                        lhs(iElement, rLiL_i) += alhsfacL;
                        lhs(iElement, rRiL_i) -= alhsfacL;

                        // upwind advection right node
                        const scalar alhsfacR = 0.5 * (tmdot - std::abs(tmdot));
                        lhs(iElement, rRiR_i) -= alhsfacR;
                        lhs(iElement, rLiR_i) += alhsfacR;
                    }

                    //===========
                    // Diffusion
                    //===========
                    auto elemNodeCoords = Kokkos::subview(
                        ws_coordinates, iElement, Kokkos::ALL, Kokkos::ALL);

#ifdef USE_INV_JACOBIAN
                    Eigen::Matrix3d invJac = ElemType::calcInvJac(
                        hexDerivViewConst_, elemNodeCoords, ip);
#else
                    scalar dndx[ElemType::nodesPerElement][SPATIAL_DIM];
                    ElemType::gradOpInline(
                        hexDerivViewConst_, elemNodeCoords, dndx, ip);
#endif
                    // // with eigenGradOp
                    // Eigen::Matrix<scalar, ElemType::nodesPerElement, 3> dndx;
                    // eigenGradOP(
                    //     hexDerivViewConst_, ws_coordinates, dndx, ip);

                    for (unsigned ic = 0; ic < ElemType::nodesPerElement; ++ic)
                    {
                        const label icNdim = ic * BLOCKSIZE;
#ifdef USE_INV_JACOBIAN
                        // with invJac
                        Eigen::Vector3d dnds(hexDerivViewConst_(ip, ic, 0),
                                             hexDerivViewConst_(ip, ic, 1),
                                             hexDerivViewConst_(ip, ic, 2));
                        Eigen::Vector3d dndx = invJac * dnds;
#endif
                        for (unsigned i = 0; i < BLOCKSIZE; ++i)
                        {
                            const label indexL = il * BLOCKSIZE + i;
                            const label indexR = ir * BLOCKSIZE + i;

                            const label rowL =
                                indexL * ElemType::nodesPerElement * BLOCKSIZE;
                            const label rowR =
                                indexR * ElemType::nodesPerElement * BLOCKSIZE;

                            scalar lhs_riC_i = 0.0;
                            for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj = ngpAreaVSTKFieldRef(
                                    elemFastIndex, ip * SPATIAL_DIM + j);

#ifdef USE_INV_JACOBIAN
                                const scalar lhsfacDiff_i =
                                    -GammaIp * dndx[j] * axj;
#else
                                const scalar lhsfacDiff_i =
                                    -GammaIp * dndx[ic][j] * axj;
#endif

                                lhs_riC_i += lhsfacDiff_i;
                            }

                            lhs(iElement, rowL + icNdim + i) += lhs_riC_i;
                            lhs(iElement, rowR + icNdim + i) -= lhs_riC_i;

                            const scalar phi = ws_phi(iElement, icNdim + i);
                            rhs(iElement, indexL) -= lhs_riC_i * phi;
                            rhs(iElement, indexR) += lhs_riC_i * phi;
                        }
                    }
                }

                applyCoeff_(A, b, connectedNodes, rhs, lhs, iElement);
            });
        });

    A.modifyDevice();
    b.modifyDevice();
}

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          typename ElemType,
          bool includeAdv,
          bool isShifted>
KOKKOS_FUNCTION void AsmCvFemInlineSF<
    scalar,
    label,
    BLOCKSIZE,
    SPATIAL_DIM,
    ElemType,
    includeAdv,
    isShifted>::applyCoeff_(const Matrix& A,
                            const Vector& b,
                            const stk::mesh::Entity
                                connectedNodes[ElemType::nodesPerElement],
                            const TeamScratchViewScalar& rhs,
                            const TeamScratchViewScalar& lhs,
                            const int& iElement,
                            bool deductUnfound) const
{
    const size_t nConnectedNodes = ElemType::nodesPerElement;
    const unsigned numRows = nConnectedNodes;
    const label nOwnedRows = this->numNodes_;
    STK_NGP_ThrowAssert(BLOCKSIZE * numRows <=
                        static_cast<unsigned>(rhs.extent(1)));
    STK_NGP_ThrowAssert(BLOCKSIZE * BLOCKSIZE * numRows * numRows <=
                        static_cast<unsigned>(lhs.extent(1)));

    label scratchIds[ElemType::nodesPerElement];
    label matrixColumnIds[ElemType::nodesPerElement];
    label sortPermutation[ElemType::nodesPerElement];

    // assume local column ordering, i.e. local id = global id = local_offset-1
    for (size_t iNode = 0; iNode < nConnectedNodes; iNode++)
    {
        stk::mesh::Entity node = connectedNodes[iNode];
        matrixColumnIds[iNode] = node.local_offset() - 1;
        scratchIds[iNode] = node.local_offset() - 1;
    }

    performSortPermutations<ElemType::nodesPerElement>(sortPermutation,
                                                       matrixColumnIds);

    for (unsigned r = 0; r < numRows; r++)
    {
        const label cur_perm_index = sortPermutation[r];
        const label rowID = scratchIds[cur_perm_index];

        const label lhsCurrentIndex =
            BLOCKSIZE * BLOCKSIZE * cur_perm_index * numRows;
        const label rhsCurrentIndex = BLOCKSIZE * cur_perm_index;

#ifndef NDEBUG
        for (unsigned k = 0; k < BLOCKSIZE; k++)
        {
            STK_NGP_ThrowAssertMsg(
                std::isfinite(rhs(iElement, rhsCurrentIndex + k)),
                "Invalid rhs");
        }
#endif /* NDEBUG */

        // assemble only owned rows (this is always true here)
        if (rowID < nOwnedRows)
        {
            auto rowVals = A.rowVals(rowID);
            const auto rowCols = A.rowCols(rowID);

            const label length = rowCols.size();
            const label numCols = nConnectedNodes;

            label offset = 0;
            label offsetOld = 0;
            for (label j = 0; j < numCols; j++)
            {
                const label permIndex = sortPermutation[j];
                const label curColumnIndex = matrixColumnIds[permIndex];

                // must be able to short-circuit!
                while (offset < length && rowCols[offset] != curColumnIndex)
                {
                    ++offset;
                }

                if (offset < length)
                {
                    for (unsigned m = 0; m < BLOCKSIZE; m++)
                    {
                        for (unsigned n = 0; n < BLOCKSIZE; n++)
                        {
                            const label column_idx =
                                BLOCKSIZE * BLOCKSIZE * offset + BLOCKSIZE * m +
                                n;
                            const label local_idx = BLOCKSIZE * numRows * m +
                                                    BLOCKSIZE * permIndex + n;
#ifndef NDEBUG
                            STK_NGP_ThrowAssertMsg(
                                std::isfinite(
                                    lhs(iElement, lhsCurrentIndex + local_idx)),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            Kokkos::atomic_add(
                                &rowVals[column_idx],
                                lhs(iElement, lhsCurrentIndex + local_idx));
                        }
                    }

                    // store this offset
                    offsetOld = offset;
                }
                else if (deductUnfound)
                {
                    // reduced stencil (remove reduced nodes contribution to
                    // diagonal)
                    scalar* diag = A.diag(rowID);
                    for (unsigned m = 0; m < BLOCKSIZE; m++)
                    {
                        for (unsigned n = 0; n < BLOCKSIZE; n++)
                        {
                            const label local_idx = BLOCKSIZE * numRows * m +
                                                    BLOCKSIZE * permIndex + n;
#ifndef NDEBUG
                            STK_NGP_ThrowAssertMsg(
                                std::isfinite(
                                    lhs(iElement, lhsCurrentIndex + local_idx)),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            Kokkos::atomic_add(
                                &diag[BLOCKSIZE * m + n],
                                lhs(iElement, lhsCurrentIndex + local_idx));
                        }
                    }

                    // reset offset
                    offset = offsetOld;
                }
            }
            for (unsigned k = 0; k < BLOCKSIZE; k++)
            {
                Kokkos::atomic_add(&b[BLOCKSIZE * rowID + k],
                                   rhs(iElement, rhsCurrentIndex + k));
            }
        }
    }
}

// explicit instantiation
template class AsmCvFemInlineSF<
    double,
    int,
    1,
    3,
    FEhex<AsmBase<double, int, 1, 3, true, true>::DataType>,
    true,
    true>;

template class AsmCvFemInlineSF<
    double,
    int,
    3,
    3,
    FEhex<AsmBase<double, int, 3, 3, true, true>::DataType>,
    true,
    true>;
