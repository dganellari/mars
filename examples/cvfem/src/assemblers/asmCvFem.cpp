#include "asmCvFem.h"
#include <iostream>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv,
          bool isShifted>
void AsmCvFem<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>::
    assemble(Matrix& A,
             Vector& b,
             const stk::mesh::MetaData& metaData,
             const stk::mesh::BulkData& bulkData) const
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

    size_t systemSize = nodesPerElementMax * BLOCKSIZE;

    Kokkos::
        parallel_for(
            "asm_cvfem",
            teamPolicy.set_scratch_size(
                SCRATCH_SPACE_LEVEL,
                Kokkos::PerTeam(totalTeamScratchSize),
                Kokkos::PerThread(totalThreadScratchSize)),
            KOKKOS_CLASS_LAMBDA( // class lambda for this->numNodes_
                const TeamHandleType& teamMember) {
                const int bucketIndex =
                    bucketIds.get<ExecSpace>(teamMember.league_rank());
                const stk::mesh::NgpMesh::BucketType& elementBucket =
                    ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
                unsigned nElementsPerBucket = elementBucket.size();

                ScratchViewScalar lhs(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    systemSize * systemSize);
                ScratchViewScalar rhs(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL), systemSize);
                ScratchViewEntity connectedNodes(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    nodesPerElementMax);

                // used in applyCoeff()
                ScratchViewLabel scratchIds(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    nodesPerElementMax);
                ScratchViewIndex matrixColumnIds(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    nodesPerElementMax);
                ScratchViewIndex sortPermutation(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    nodesPerElementMax);

                // nodal fields to gather
                ScratchViewScalar ws_phi(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL), systemSize);
                ScratchViewScalar ws_beta(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL), systemSize);
                ScratchViewScalar ws_gradPhi(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    systemSize * SPATIAL_DIM);
                ScratchViewScalar ws_Gamma(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL), systemSize);

                // geometry related to populate
                ScratchView2D ws_coordinates(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    nodesPerElementMax,
                    SPATIAL_DIM);
                ScratchView2D ws_scs_areav(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    numScsIpMax,
                    SPATIAL_DIM);
                ScratchView3D ws_dndx(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    numScsIpMax,
                    nodesPerElementMax,
                    SPATIAL_DIM);
                ScratchView3D ws_deriv(
                    teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                    numScsIpMax,
                    nodesPerElementMax,
                    SPATIAL_DIM);

                // allocated for entire team
                ScratchView2D ws_shape_function(
                    teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                    numScsIpMax,
                    nodesPerElementMax);

                // extract master element for current topology type
                label topology_index = elementBucket.topology();
                accel::MasterElement* meSCS = meSCS_ptrs_[topology_index];

                // extract master element specifics
                const label nodesPerElement = meSCS->nodesPerElement_;
                const label numScsIp = meSCS->num_integration_points();
                const label* lrscv = meSCS->adjacentNodes();

                const label lhsSize = systemSize * systemSize;
                const label rhsSize = systemSize;

                // extract shape function
                if (teamMember.team_rank() == 0)
                {
                    if (isShifted)
                    {
                        meSCS->shifted_shape_fcn<>(ws_shape_function);
                    }
                    else
                    {
                        meSCS->shape_fcn<>(ws_shape_function);
                    }
                }
                teamMember.team_barrier();

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
                        lhs[p] = 0.0;
                    }
                    for (label p = 0; p < rhsSize; ++p)
                    {
                        rhs[p] = 0.0;
                    }

                    //===================
                    // gather nodal data
                    //===================
                    stk::mesh::NgpMesh::ConnectedNodes nodeRels =
                        ngpMesh.get_nodes(stk::topology::ELEM_RANK,
                                          elemFastIndex);
                    label numNodes =
                        nodeRels.size(); // bulkData.num_nodes(elem);

                    // sanity check on num nodes
                    STK_NGP_ThrowAssert(numNodes == nodesPerElement);

                    for (label ni = 0; ni < numNodes; ++ni)
                    {
                        stk::mesh::Entity node = nodeRels[ni];

                        // set connected nodes
                        connectedNodes[ni] = node;

                        stk::mesh::FastMeshIndex nodeFastIndex =
                            ngpMesh.fast_mesh_index(node);

                        // gather scalars
                        ws_Gamma[ni] = ngpGammaSTKFieldRef(nodeFastIndex, 0);

                        // gather vectors
                        for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                        {
                            ws_coordinates(ni, i) =
                                ngpCoordinatesRef(nodeFastIndex, i);
                        }

                        // gather BLOCKSIZE-dim fields
                        for (unsigned i = 0; i < BLOCKSIZE; ++i)
                        {
                            ws_phi[ni * BLOCKSIZE + i] =
                                ngpPhiSTKFieldRef(nodeFastIndex, i);
                            ws_beta[ni * BLOCKSIZE + i] =
                                ngpBetaSTKFieldRef(nodeFastIndex, i);

                            for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                            {
                                ws_gradPhi[ni * BLOCKSIZE * SPATIAL_DIM +
                                           i * SPATIAL_DIM + j] =
                                    ngpGradPhiSTKFieldRef(nodeFastIndex,
                                                          i * SPATIAL_DIM + j);
                            }
                        }
                    }

                    // compute geometry
                    meSCS->determinant(ws_coordinates, ws_scs_areav);
                    // compute dndx
                    if (isShifted) //(phi_->isShifted())
                    {
                        meSCS->shifted_grad_op(
                            ws_coordinates, ws_dndx, ws_deriv);
                    }
                    else
                    {
                        meSCS->grad_op(ws_coordinates, ws_dndx, ws_deriv);
                    }

                    for (label ip = 0; ip < numScsIp; ++ip)
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

                        for (label ic = 0; ic < nodesPerElement; ++ic)
                        {
                            const scalar r =
                                ws_shape_function(ip, ic)._data.get();

                            GammaIp += r * ws_Gamma[ic];

                            // compute scs ip value
                            for (unsigned i = 0; i < BLOCKSIZE; ++i)
                            {
                                phiIp[i] += r * ws_phi[ic * BLOCKSIZE + i];
                            }
                            for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                            {
                                coordIp[i] +=
                                    r * ws_coordinates(ic, i)._data.get();
                            }
                        }

                        //===========
                        // Advection
                        //===========

                        // final upwind extrapolation;
                        for (unsigned i = 0; i < BLOCKSIZE; ++i)
                        {
                            phiIpL[i] = ws_phi[il * BLOCKSIZE + i];
                            phiIpR[i] = ws_phi[ir * BLOCKSIZE + i];
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
                                        ws_coordinates(il, j)._data.get();
                                    dcorr += ws_beta[il * BLOCKSIZE + i] * dxj *
                                             ws_gradPhi[il * BLOCKSIZE *
                                                            SPATIAL_DIM +
                                                        i * SPATIAL_DIM + j];
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
                                        ws_coordinates(ir, j)._data.get();
                                    dcorr += ws_beta[ir * BLOCKSIZE + i] * dxj *
                                             ws_gradPhi[ir * BLOCKSIZE *
                                                            SPATIAL_DIM +
                                                        i * SPATIAL_DIM + j];
                                }
                            }

                            // total upwind advection
                            const scalar aflux = tmdot * (phiUpwind + dcorr);

                            const label indexL = il * BLOCKSIZE + i;
                            const label indexR = ir * BLOCKSIZE + i;

                            const label rowL =
                                indexL * nodesPerElement * BLOCKSIZE;
                            const label rowR =
                                indexR * nodesPerElement * BLOCKSIZE;

                            const label rLiL_i = rowL + il * BLOCKSIZE + i;
                            const label rLiR_i = rowL + ir * BLOCKSIZE + i;
                            const label rRiR_i = rowR + ir * BLOCKSIZE + i;
                            const label rRiL_i = rowR + il * BLOCKSIZE + i;

                            // right hand side; L and R
                            rhs[indexL] -= aflux;
                            rhs[indexR] += aflux;

                            // upwind advection left node
                            const scalar alhsfacL =
                                0.5 * (tmdot + std::abs(tmdot));
                            lhs[rLiL_i] += alhsfacL;
                            lhs[rRiL_i] -= alhsfacL;

                            // upwind advection right node
                            const scalar alhsfacR =
                                0.5 * (tmdot - std::abs(tmdot));
                            lhs[rRiR_i] -= alhsfacR;
                            lhs[rLiR_i] += alhsfacR;
                        }

                        //===========
                        // Diffusion
                        //===========

                        for (label ic = 0; ic < nodesPerElement; ++ic)
                        {
                            const label icNdim = ic * BLOCKSIZE;

                            for (unsigned i = 0; i < BLOCKSIZE; ++i)
                            {
                                const label indexL = il * BLOCKSIZE + i;
                                const label indexR = ir * BLOCKSIZE + i;

                                const label rowL =
                                    indexL * nodesPerElement * BLOCKSIZE;
                                const label rowR =
                                    indexR * nodesPerElement * BLOCKSIZE;

                                scalar lhs_riC_i = 0.0;
                                for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    const scalar axj =
                                        ws_scs_areav(ip, j)._data.get();

                                    const scalar lhsfacDiff_i =
                                        -GammaIp *
                                        ws_dndx(ip, ic, j)._data.get() * axj;

                                    lhs_riC_i += lhsfacDiff_i;
                                }

                                lhs[rowL + icNdim + i] += lhs_riC_i;
                                lhs[rowR + icNdim + i] -= lhs_riC_i;

                                const scalar phi = ws_phi[icNdim + i];
                                rhs[indexL] -= lhs_riC_i * phi;
                                rhs[indexR] += lhs_riC_i * phi;
                            }
                        }
                    }

                    applyCoeff_(A,
                                b,
                                connectedNodes,
                                scratchIds,
                                matrixColumnIds,
                                sortPermutation,
                                rhs,
                                lhs,
                                nodesPerElement);
                });
            });

    A.modifyDevice();
    b.modifyDevice();
}

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv,
          bool isShifted>
KOKKOS_FUNCTION void
AsmCvFem<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>::
    applyCoeff_(const Matrix& A,
                const Vector& b,
                const ScratchViewEntity& connectedNodes,
                const ScratchViewLabel& scratchIds,
                const ScratchViewIndex& matrixColumnIds,
                const ScratchViewIndex& sortPermutation,
                const ScratchViewScalar& rhs,
                const ScratchViewScalar& lhs,
                label nodesPerElement,
                bool deductUnfound) const
{
    const size_t nConnectedNodes = nodesPerElement;
    const unsigned numRows = nConnectedNodes;
    const label nOwnedRows = this->numNodes_;

    STK_NGP_ThrowAssert(BLOCKSIZE * numRows <=
                        static_cast<unsigned>(rhs.size()));
    STK_NGP_ThrowAssert(BLOCKSIZE * BLOCKSIZE * numRows * numRows <=
                        static_cast<unsigned>(lhs.size()));
    // assume local column ordering, i.e. local id = global id = local_offset-1
    for (size_t iNode = 0; iNode < nConnectedNodes; iNode++)
    {
        stk::mesh::Entity node = connectedNodes[iNode];
        matrixColumnIds[iNode] = node.local_offset() - 1;
        scratchIds[iNode] = node.local_offset() - 1;
    }

    performSortPermutations(sortPermutation, matrixColumnIds, nConnectedNodes);

    for (unsigned r = 0; r < numRows; r++)
    {
        const label cur_perm_index = sortPermutation[r];
        const label rowID = scratchIds[cur_perm_index];

        const scalar* const lhsCurrent =
            &lhs[BLOCKSIZE * BLOCKSIZE * cur_perm_index * numRows];
        const scalar* const rhsCurrent = &rhs[BLOCKSIZE * cur_perm_index];

#ifndef NDEBUG
        for (unsigned k = 0; k < BLOCKSIZE; k++)
        {
            STK_NGP_ThrowAssertMsg(std::isfinite(rhsCurrent[k]), "Invalid rhs");
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
                                std::isfinite(lhsCurrent[local_idx]),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            Kokkos::atomic_add(&rowVals[column_idx],
                                               lhsCurrent[local_idx]);
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
                                std::isfinite(lhsCurrent[local_idx]),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            Kokkos::atomic_add(&diag[BLOCKSIZE * m + n],
                                               lhsCurrent[local_idx]);
                        }
                    }

                    // reset offset
                    offset = offsetOld;
                }
            }
            for (unsigned k = 0; k < BLOCKSIZE; k++)
            {
                Kokkos::atomic_add(&b[BLOCKSIZE * rowID + k], rhsCurrent[k]);
            }
        }
    }
}

// explicit instantiation
template class AsmCvFem<double, int, 1, 3, true, true>;
template class AsmCvFem<double, int, 3, 3, true, true>;
