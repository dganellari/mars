#include "asmEdge.h"
#include <iostream>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv,
          bool isShifted>
void AsmEdge<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>::
    assemble(Matrix& A,
             Vector& b,
             const stk::mesh::MetaData& metaData,
             const stk::mesh::BulkData& bulkData) const
{
    // get fields references
    // TODO: choose one version (ptr/ref)
    const stk::mesh::Field<scalar>& GammaSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "gamma");
    const stk::mesh::Field<scalar>& mdotSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::EDGE_RANK, "mdotEdge");
    const stk::mesh::Field<scalar>& phiSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phi");
    const stk::mesh::Field<scalar>& gradPhiSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phiGrad");
    const stk::mesh::Field<scalar>& betaSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "beta");
    const stk::mesh::Field<scalar>& coordinatesRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "coordinates");
    const stk::mesh::Field<scalar>& areaVEdgeRef =
        *metaData.get_field<scalar>(stk::topology::EDGE_RANK, "areaVEdge");

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
    stk::mesh::NgpField<double>& ngpAreaVEdgeRef =
        stk::mesh::get_updated_ngp_field<double>(areaVEdgeRef);

    const stk::mesh::Selector selAllEdges = metaData.universal_part();

    stk::NgpVector<unsigned> bucketIds =
        ngpMesh.get_bucket_ids(stk::topology::EDGE_RANK, selAllEdges);
    unsigned numBuckets = bucketIds.size();

    auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(numBuckets, Kokkos::AUTO);

    size_t totalTeamScratchSize;
    size_t totalThreadScratchSize;

    getScratchSizes(totalTeamScratchSize, totalThreadScratchSize);

    size_t systemSize = nodesPerEdge * BLOCKSIZE;

    Kokkos::parallel_for(
        "asm_edge",
        teamPolicy.set_scratch_size(SCRATCH_SPACE_LEVEL,
                                    Kokkos::PerTeam(totalTeamScratchSize),
                                    Kokkos::PerThread(totalThreadScratchSize)),
        KOKKOS_CLASS_LAMBDA( // TODO: class lambda here for numNodes_
            const TeamHandleType& teamMember) {
            const int bucketIndex =
                bucketIds.get<ExecSpace>(teamMember.league_rank());
            const stk::mesh::NgpMesh::BucketType& edgeBucket =
                ngpMesh.get_bucket(stk::topology::EDGE_RANK, bucketIndex);
            unsigned nEdgesPerBucket = edgeBucket.size();

            TeamScratchViewScalar lhs(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nEdgesPerBucket,
                systemSize * systemSize);
            TeamScratchViewScalar rhs(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nEdgesPerBucket,
                systemSize);

            const label lhsSize = systemSize * systemSize;
            const label rhsSize = systemSize;

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(teamMember, 0u, nEdgesPerBucket),
                [&](const int& iEdge)
            {
                constexpr unsigned numNodes = nodesPerEdge;

                stk::mesh::Entity connectedNodes[numNodes];
                scalar ws_Gamma[numNodes];
                scalar ws_coordinates[numNodes][SPATIAL_DIM];
                scalar ws_phi[numNodes * BLOCKSIZE];
                scalar ws_beta[numNodes * BLOCKSIZE];
                scalar ws_gradPhi[numNodes * BLOCKSIZE * SPATIAL_DIM];
                scalar phiIp[BLOCKSIZE];
                scalar coordIp[SPATIAL_DIM];

                // extrapolated value from the L/R direction
                scalar phiIpL[BLOCKSIZE];
                scalar phiIpR[BLOCKSIZE];

                scalar ws_scs_areav[SPATIAL_DIM];

                // get edge
                auto edge = edgeBucket[iEdge];
                const auto edgeIndex = ngpMesh.fast_mesh_index(edge);

                // zero lhs/rhs
                for (label p = 0; p < lhsSize; ++p)
                {
                    lhs(iEdge, p) = 0.0;
                }
                for (label p = 0; p < rhsSize; ++p)
                {
                    rhs(iEdge, p) = 0.0;
                }

                //===================
                // gather nodal data
                //===================

                for (unsigned ni = 0; ni < numNodes; ++ni)
                {
                    const auto node = ngpMesh.get_nodes(
                        stk::topology::EDGE_RANK, edgeIndex)[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    stk::mesh::FastMeshIndex nodeFastIndex =
                        ngpMesh.fast_mesh_index(node);

                    // gather scalars
                    ws_Gamma[ni] = ngpGammaSTKFieldRef(nodeFastIndex, 0);

                    for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                    {
                        ws_coordinates[ni][i] =
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

                // left and right nodes for this ip
                const label il = 0;
                const label ir = 1;

                // save off mdot
                // EDGE: need edge field for mdot
                const scalar tmdot =
                    includeAdv ? ngpMdotSTKFieldRef(edgeIndex, 0) : 0.0;

                // gather ip values
                const scalar GammaIp =
                    scalar(0.5) * (ws_Gamma[il] + ws_Gamma[ir]);
                for (unsigned j = 0; j < BLOCKSIZE; ++j)
                {
                    phiIp[j] = scalar(0.5) * (ws_phi[il] + ws_phi[ir]);
                }
                for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                {
                    coordIp[j] = scalar(0.5) * (ws_coordinates[il][j] +
                                                ws_coordinates[ir][j]);
                }
                for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_scs_areav[j] = ngpAreaVEdgeRef(edgeIndex, j);
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
                                coordIp[j] - ws_coordinates[il][j];
                            dcorr += ws_beta[il * BLOCKSIZE + i] * dxj *
                                     ws_gradPhi[il * BLOCKSIZE * SPATIAL_DIM +
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
                                coordIp[j] - ws_coordinates[ir][j];
                            dcorr += ws_beta[ir * BLOCKSIZE + i] * dxj *
                                     ws_gradPhi[ir * BLOCKSIZE * SPATIAL_DIM +
                                                i * SPATIAL_DIM + j];
                        }
                    }

                    // total upwind advection
                    const scalar aflux = tmdot * (phiUpwind + dcorr);

                    const label indexL = il * BLOCKSIZE + i;
                    const label indexR = ir * BLOCKSIZE + i;

                    const label rowL = indexL * nodesPerEdge * BLOCKSIZE;
                    const label rowR = indexR * nodesPerEdge * BLOCKSIZE;

                    const label rLiL_i = rowL + il * BLOCKSIZE + i;
                    const label rLiR_i = rowL + ir * BLOCKSIZE + i;
                    const label rRiR_i = rowR + ir * BLOCKSIZE + i;
                    const label rRiL_i = rowR + il * BLOCKSIZE + i;

                    // right hand side; L and R
                    rhs(iEdge, indexL) -= aflux;
                    rhs(iEdge, indexR) += aflux;

                    // upwind advection left node
                    const scalar alhsfacL = 0.5 * (tmdot + std::abs(tmdot));
                    lhs(iEdge, rLiL_i) += alhsfacL;
                    lhs(iEdge, rRiL_i) -= alhsfacL;

                    // upwind advection right node
                    const scalar alhsfacR = 0.5 * (tmdot - std::abs(tmdot));
                    lhs(iEdge, rRiR_i) -= alhsfacR;
                    lhs(iEdge, rLiR_i) += alhsfacR;
                }

                //===========
                // Diffusion
                //===========
                // TODO: non-orthogonality correction is OFF

                // // Compute area vector related quantities and (U dot
                // areaVec)
                scalar axdx = 0.0;
                scalar asq = 0.0;
                //   scalar udotx = 0.0;
                for (unsigned d = 0; d < SPATIAL_DIM; ++d)
                {
                    const scalar dxj =
                        ws_coordinates[il][d] - ws_coordinates[ir][d];
                    asq += ws_scs_areav[d] * ws_scs_areav[d];
                    axdx += ws_scs_areav[d] * dxj;
                    // udotx += 0.5 * dxj * (vrtm.get(nodeR, d) +
                    // vrtm.get(nodeL, d));
                }
                const scalar inv_axdx = scalar(1.0) / axdx;
                // const auto nodeL =
                // ngpMesh.get_nodes(stk::topology::EDGE_RANK,
                // edgeIndex)[0]; const auto nodeR =
                // ngpMesh.get_nodes(stk::topology::EDGE_RANK,
                // edgeIndex)[1];


                // Diffusive flux
                const scalar phiL = ws_phi[il];
                const scalar phiR = ws_phi[ir];
                const scalar lhsfac = -GammaIp * asq * inv_axdx;
                const scalar diffFlux = lhsfac * (phiR - phiL); // + nonOrth;

                // LHS
                // Left node
                lhs(iEdge, 0) = lhsfac;  //(0, 0)
                lhs(iEdge, 1) = -lhsfac; //(0, 1)
                // Right node
                lhs(iEdge, 2) = -lhsfac; //(1, 0)
                lhs(iEdge, 3) = lhsfac;  //(1, 1)

                // RHS
                // Left node
                rhs(iEdge, 0) = diffFlux;
                // Right node
                rhs(iEdge, 1) = -diffFlux;

                applyCoeff_(A, b, connectedNodes, rhs, lhs, iEdge);
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
AsmEdge<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>::
    applyCoeff_(const Matrix& A,
                const Vector& b,
                const stk::mesh::Entity connectedNodes[2],
                const TeamScratchViewScalar& rhs,
                const TeamScratchViewScalar& lhs,
                const int& iEdge,
                bool deductUnfound) const
{
    constexpr size_t nConnectedNodes = 2;
    const unsigned numRows = nConnectedNodes;
    const label nOwnedRows = this->numNodes_;
    STK_NGP_ThrowAssert(BLOCKSIZE * numRows <=
                        static_cast<unsigned>(rhs.extent(1)));
    STK_NGP_ThrowAssert(BLOCKSIZE * BLOCKSIZE * numRows * numRows <=
                        static_cast<unsigned>(lhs.extent(1)));

    size_t matrixColumnIds[2];
    size_t scratchIds[2];

    // assume local column ordering, i.e. local id = global id = local_offset-1
    for (size_t iNode = 0; iNode < nConnectedNodes; iNode++)
    {
        stk::mesh::Entity node = connectedNodes[iNode];
        matrixColumnIds[iNode] = node.local_offset() - 1;
        scratchIds[iNode] = node.local_offset() - 1;
    }

    size_t sortPermutation[2] = {0, 1};
    if (matrixColumnIds[1] < matrixColumnIds[0])
    {
        sortPermutation[0] = 1;
        sortPermutation[1] = 0;
    }

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
                std::isfinite(rhs(iEdge, rhsCurrentIndex + k)), "Invalid rhs");
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
                                    lhs(iEdge, lhsCurrentIndex + local_idx)),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            Kokkos::atomic_add(
                                &rowVals[column_idx],
                                lhs(iEdge, lhsCurrentIndex + local_idx));
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
                                    lhs(iEdge, lhsCurrentIndex + local_idx)),
                                "Inf or NAN lhs");
#endif /* NDEBUG */
                            Kokkos::atomic_add(
                                &diag[BLOCKSIZE * m + n],
                                lhs(iEdge, lhsCurrentIndex + local_idx));
                        }
                    }

                    // reset offset
                    offset = offsetOld;
                }
            }
            for (unsigned k = 0; k < BLOCKSIZE; k++)
            {
                Kokkos::atomic_add(&b[BLOCKSIZE * rowID + k],
                                   rhs(iEdge, rhsCurrentIndex + k));
            }
        }
    }
}

// explicit instantiation
template class AsmEdge<double, int, 1, 3, true, true>;
// template class AsmEdge<double, int, 3, 3, true, true>; // not yet supported
