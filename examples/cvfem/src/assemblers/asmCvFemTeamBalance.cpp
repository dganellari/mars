#include "asmCvFemTeamBalance.h"
#include <iostream>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <Kokkos_Core.hpp>

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          typename ElemType,
          bool includeAdv,
          bool isShifted>
void AsmCvFemTeamBalance<scalar, label, BLOCKSIZE, SPATIAL_DIM, ElemType, includeAdv, isShifted>::assemble(
    Matrix& A, Vector& b,
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData) const
{
    // Get field references
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

    // Initialize NGP data
    const stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulkData);
    stk::mesh::NgpField<scalar>& ngpGammaSTKFieldRef = stk::mesh::get_updated_ngp_field<scalar>(GammaSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpMdotSTKFieldRef = stk::mesh::get_updated_ngp_field<scalar>(mdotSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpPhiSTKFieldRef = stk::mesh::get_updated_ngp_field<scalar>(phiSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpGradPhiSTKFieldRef = stk::mesh::get_updated_ngp_field<scalar>(gradPhiSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpBetaSTKFieldRef = stk::mesh::get_updated_ngp_field<scalar>(betaSTKFieldRef);
    stk::mesh::NgpField<scalar>& ngpAreaVSTKFieldRef = stk::mesh::get_updated_ngp_field<scalar>(areaVSTKFieldRef);
    stk::mesh::NgpField<double>& ngpCoordinatesRef = stk::mesh::get_updated_ngp_field<double>(coordinatesRef);

    const stk::mesh::Selector selAllElements = metaData.universal_part();

    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEMENT_RANK, selAllElements);
    unsigned numBuckets = bucketIds.size();

    const int teamsPerBucket = 8;  // 8x teams per bucket
    const int elementsPerThread = 1;  // Multiple elements per thread
    const int teamSize = 256;  // Explicit team size
    // Increase PARALLELISM and improve load balancing - More teams, fewer elements per team
    auto teamPolicy = stk::ngp::TeamPolicy<ExecSpace>(numBuckets * teamsPerBucket, teamSize); // teamsPerBucket x teams

    size_t totalTeamScratchSize;
    size_t totalThreadScratchSize;
    getScratchSizes(totalTeamScratchSize, totalThreadScratchSize);

    constexpr size_t systemSize = ElemType::nodesPerElement * BLOCKSIZE;
    printf("System size team balance: %zu\n", systemSize);

    Kokkos::parallel_for(
        "asm_cvfem_ultra_optimized",
        teamPolicy.set_scratch_size(SCRATCH_SPACE_LEVEL,
                                    Kokkos::PerTeam(2 * 1024 * 1024),  // Large scratch per team (4MB, adjust if needed)
                                    Kokkos::PerThread(0)),  // No thread scratch
        KOKKOS_CLASS_LAMBDA(const TeamHandleType& teamMember) {
            const int virtualBucketIndex = teamMember.league_rank();
            const int bucketIndex = bucketIds.get<ExecSpace>(virtualBucketIndex / teamsPerBucket);
            const int bucketOffset = virtualBucketIndex % teamsPerBucket;
            const stk::mesh::NgpMesh::BucketType& elementBucket = 
                ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
            unsigned totalElementsInBucket = elementBucket.size();
            
            // WORK DISTRIBUTION - Each team handles subset of elements
            unsigned elementsPerTeam = teamSize * elementsPerThread;
            unsigned startElement = bucketOffset * elementsPerTeam;
            unsigned endElement = Kokkos::min(startElement + elementsPerTeam, totalElementsInBucket);
            unsigned nElementsPerBucket = endElement - startElement;

            if (nElementsPerBucket == 0) return;  // No work for this team

            // MINIMAL SCRATCH MEMORY - Only essential arrays
            // TeamScratchViewScalar lhs(teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
            //                          nElementsPerBucket, systemSize * systemSize);
            // TeamScratchViewScalar rhs(teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
            //                          nElementsPerBucket, systemSize);

            TeamScratchViewScalar3 ws_coordinates(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                ElemType::nodesPerElement,
                SPATIAL_DIM);

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

            TeamScratchViewScalar ws_gamma(
                teamMember.team_scratch(SCRATCH_SPACE_LEVEL),
                nElementsPerBucket,
                systemSize);

            constexpr int lrscv[24] = {0, 1, 1, 2, 2, 3, 0, 3, 4, 5, 5, 6,
                                       6, 7, 4, 7, 0, 4, 1, 5, 2, 6, 3, 7};

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(teamMember, 0u, teamSize),
                [&](const int& threadIndex) {
                    const int startLocal = threadIndex * elementsPerThread;
                    const int endLocal = Kokkos::min(startLocal + elementsPerThread, static_cast<int>(nElementsPerBucket));
                    
                    for (int localElementIndex = startLocal; localElementIndex < endLocal; ++localElementIndex) {
                        const int iElement = startElement + localElementIndex;
                        
                        scalar phiIp[BLOCKSIZE];
                        scalar phiL[BLOCKSIZE], phiR[BLOCKSIZE];
                        scalar coordIp[SPATIAL_DIM];

                        scalar lhs[systemSize * systemSize] = {0.0};
                        scalar rhs[systemSize] = {0.0};

                        // get elem
                        stk::mesh::Entity elem = elementBucket[iElement];
                        stk::mesh::FastMeshIndex elemFastIndex =
                            ngpMesh.fast_mesh_index(elem);

                        // for (label p = 0; p < lhsSize; ++p) { lhs(localElementIndex, p) = 0.0; }
                        // for (label p = 0; p < rhsSize; ++p) { rhs(localElementIndex, p) = 0.0; }

                        // COALESCED DATA GATHERING (FIX L1TEX STALLS)
                        stk::mesh::NgpMesh::ConnectedNodes nodeRels = 
                            ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemFastIndex);
                        label numNodes = nodeRels.size();

                        stk::mesh::Entity connectedNodes[ElemType::nodesPerElement];
                        STK_NGP_ThrowAssert(ElemType::nodesPerElement == numNodes);
                        stk::mesh::FastMeshIndex nodeFastIndices[ElemType::nodesPerElement];

                        for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni) {
                            connectedNodes[ni] = nodeRels[ni];
                            nodeFastIndices[ni] = ngpMesh.fast_mesh_index(nodeRels[ni]);
                        }

                        for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni) {
                            ws_gamma(localElementIndex, ni) = ngpGammaSTKFieldRef(nodeFastIndices[ni], 0);
                        }

                        for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni) {
                            for(unsigned i = 0; i < BLOCKSIZE; ++i) {
                                ws_phi(localElementIndex, ni * BLOCKSIZE + i) = ngpPhiSTKFieldRef(nodeFastIndices[ni], i);
                            }
                        }

                        for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni) {
                            for(unsigned i = 0; i < BLOCKSIZE; ++i) {
                                ws_beta(localElementIndex, ni * BLOCKSIZE + i) = ngpBetaSTKFieldRef(nodeFastIndices[ni], i);
                            }
                        }
                        
                        for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni) {
                            for(unsigned i = 0; i < SPATIAL_DIM; ++i) {
                                ws_coordinates(localElementIndex, ni, i) = ngpCoordinatesRef(nodeFastIndices[ni], i);
                            }
                        }
                        
                        for (unsigned ni = 0; ni < ElemType::nodesPerElement; ++ni) {
                            for (unsigned i= 0; i < BLOCKSIZE; ++i) {
                                for(unsigned j = 0; j < SPATIAL_DIM; ++j) {
                                    ws_gradPhi(localElementIndex, ni * BLOCKSIZE * SPATIAL_DIM + i * SPATIAL_DIM + j) = ngpGradPhiSTKFieldRef(nodeFastIndices[ni], i * SPATIAL_DIM + j);
                                }
                            }
                        }

                        for (unsigned ip = 0; ip < ElemType::numScsIp; ++ip) {
                            const label il = lrscv[2 * ip];
                            const label ir = lrscv[2 * ip + 1];
                            
                            const scalar tmdot = includeAdv ? ngpMdotSTKFieldRef(elemFastIndex, ip) : 0.0;

                            // VECTORIZED INTERPOLATION - Access shared memory
                            for(unsigned i = 0; i < BLOCKSIZE; ++i) {
                                phiIp[i] = 0.0;
                            }
                            for(unsigned i = 0; i < SPATIAL_DIM; ++i) {
                                coordIp[i] = 0.0;
                            }
                            scalar GammaIp = 0.0;
                            
                            #pragma unroll
                            for (unsigned ic = 0; ic < ElemType::nodesPerElement; ++ic) {
                                const scalar r = sfViewConst_(ip, ic); // Using precomputed shape functions
                            
                                GammaIp += r * ws_gamma(localElementIndex, ic);
                                for(unsigned i = 0; i < BLOCKSIZE; ++i) {
                                    phiIp[i] += r * ws_phi(localElementIndex, ic * BLOCKSIZE + i);
                                }
                                for(unsigned i = 0; i < SPATIAL_DIM; ++i) {
                                    coordIp[i] += r * ws_coordinates(localElementIndex, ic, i);
                                }
                            }

                            // VECTORIZED ADVECTION - Access shared memory
                            for(unsigned i = 0; i < BLOCKSIZE; ++i) {
                                phiL[i] = ws_phi(localElementIndex, il * BLOCKSIZE + i);
                                phiR[i] = ws_phi(localElementIndex, ir * BLOCKSIZE + i);
                            }

                            for (unsigned i = 0; i < BLOCKSIZE; ++i) {
                                scalar phiUpwind, dcorr;
                                dcorr = 0.0;
                                if (tmdot > 0) {
                                    phiUpwind = phiL[i];
                                    const scalar beta = ws_beta(localElementIndex, il * BLOCKSIZE + i);

                                    for (unsigned j = 0; j < SPATIAL_DIM; ++j) {
                                        const scalar dxj = coordIp[j] - ws_coordinates(localElementIndex, il, j);
                                        const scalar grad = ws_gradPhi(localElementIndex, il * BLOCKSIZE * SPATIAL_DIM + i * SPATIAL_DIM + j);
                                        dcorr += dxj * grad;
                                    }
                                    dcorr *= beta;
                                } else {
                                    phiUpwind = phiR[i];
                                        
                                    const scalar beta = ws_beta(localElementIndex, ir * BLOCKSIZE + i);

                                    for (unsigned j = 0; j < SPATIAL_DIM; ++j) {
                                        const scalar dxj = coordIp[j] - ws_coordinates(localElementIndex, ir, j);
                                        const scalar grad = ws_gradPhi(localElementIndex, ir * BLOCKSIZE * SPATIAL_DIM + i * SPATIAL_DIM + j);
                                        dcorr += dxj * grad;
                                    }
                                    dcorr *= beta;
                                }

                                const scalar aflux = tmdot * (phiUpwind + dcorr);
                                const label indexL = il * BLOCKSIZE + i;
                                const label indexR = ir * BLOCKSIZE + i;
                                const label rowL = indexL * systemSize;
                                const label rowR = indexR * systemSize;
                                const label rLiL_i = rowL + il * BLOCKSIZE + i;
                                const label rRiL_i = rowR + il * BLOCKSIZE + i;
                                const label rRiR_i = rowR + ir * BLOCKSIZE + i;
                                const label rLiR_i = rowL + ir * BLOCKSIZE + i;

                                rhs[indexL] -= aflux;
                                rhs[indexR] += aflux;

                                const scalar alhsfacL = 0.5 * (tmdot + std::abs(tmdot));
                                const scalar alhsfacR = 0.5 * (tmdot - std::abs(tmdot));

                                lhs[rLiL_i] += alhsfacL;
                                lhs[rRiL_i] -= alhsfacL;
                                lhs[rLiR_i] += alhsfacR;
                                lhs[rRiR_i] -= alhsfacR;
                            }
                            
                            scalar dndx[ElemType::nodesPerElement][SPATIAL_DIM];
                            auto elemNodeCoords = Kokkos::subview(ws_coordinates, localElementIndex, Kokkos::ALL(), Kokkos::ALL()); 

                            #ifdef USE_INV_JACOBIAN
                                // Device-compatible inverse Jacobian (3x3 array)
                                // scalar invJac[3][3];
                                // ElemType::calcInvJac(hexDerivViewConst_, elemNodeCoords, ip, invJac); // You must implement this
                            
                                for (unsigned ic = 0; ic < ElemType::nodesPerElement; ++ic) {
                                    for (unsigned j = 0; j < SPATIAL_DIM; ++j) {
                                        dndx[ic][j] = 0.0;
                                    }
                                }
                            #else
                                ElemType::gradOpInline(hexDerivViewConst_, elemNodeCoords, dndx, ip);
                            #endif

                            scalar areaV_local[SPATIAL_DIM];
                            for (unsigned j = 0; j < SPATIAL_DIM; ++j) {
                                areaV_local[j] = ngpAreaVSTKFieldRef(elemFastIndex, ip * SPATIAL_DIM + j);
                            }

                            #pragma unroll
                            for (unsigned ic = 0; ic < ElemType::nodesPerElement; ++ic) {
                                const label icNdim = ic * BLOCKSIZE;
                                for(unsigned i = 0; i < BLOCKSIZE; ++i) {
                                    const label indexL = il * BLOCKSIZE + i;
                                    const label indexR = ir * BLOCKSIZE + i;
                                    const label rowL = indexL * systemSize;
                                    const label rowR = indexR * systemSize;

                                    scalar lhs_riC = 0.0;
                                    #pragma unroll
                                    for (unsigned j = 0; j < SPATIAL_DIM; ++j) {
                                        // const scalar axj = ngpAreaVSTKFieldRef(elemFastIndex, ip * SPATIAL_DIM + j);
                                        const scalar axj = areaV_local[j];
                                        lhs_riC -= GammaIp * dndx[ic][j] * axj;
                                    }

                                    lhs[rowL + icNdim + i] += lhs_riC;
                                    lhs[rowR + icNdim + i] -= lhs_riC;
                                    rhs[indexL] -= lhs_riC * ws_phi(localElementIndex, icNdim + i);
                                    rhs[indexR] += lhs_riC * ws_phi(localElementIndex, icNdim + i);
                                }
                            }
                        }

                        applyCoeff_(A, b, connectedNodes, rhs, lhs);
                    }
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
KOKKOS_FUNCTION void AsmCvFemTeamBalance<scalar, label, BLOCKSIZE, SPATIAL_DIM, ElemType, includeAdv, isShifted>::applyCoeff_(
    const Matrix& A, const Vector& b,
    const stk::mesh::Entity connectedNodes[ElemType::nodesPerElement],
    const scalar* rhs,
    const scalar* lhs,
    bool deductUnfound) const
{
    const size_t nConnectedNodes = ElemType::nodesPerElement;
    const unsigned numRows = nConnectedNodes;
    const label nOwnedRows = this->numNodes_;
    // STK_NGP_ThrowAssert(BLOCKSIZE * numRows <= static_cast<unsigned>(rhs.extent(1)));
    // STK_NGP_ThrowAssert(BLOCKSIZE * BLOCKSIZE * numRows * numRows <= static_cast<unsigned>(lhs.extent(1)));

    label scratchIds[ElemType::nodesPerElement];
    label matrixColumnIds[ElemType::nodesPerElement];
    label sortPermutation[ElemType::nodesPerElement];

    // Assume local column ordering, i.e. local id = global id = local_offset-1
    for (size_t iNode = 0; iNode < nConnectedNodes; iNode++) {
        stk::mesh::Entity node = connectedNodes[iNode];
        matrixColumnIds[iNode] = node.local_offset() - 1;
        scratchIds[iNode] = node.local_offset() - 1;
    }

    performSortPermutations<ElemType::nodesPerElement>(sortPermutation, matrixColumnIds);

        for (unsigned r = 0; r < numRows; r++) {
        const label cur_perm_index = sortPermutation[r];
        const label rowID = scratchIds[cur_perm_index];
    
        const label lhsCurrentIndex = BLOCKSIZE * BLOCKSIZE * cur_perm_index * numRows;
        const label rhsCurrentIndex = BLOCKSIZE * cur_perm_index;
    
        // Assemble only owned rows
        if (rowID < nOwnedRows) {
            auto rowVals = A.rowVals(rowID);
            const auto rowCols = A.rowCols(rowID);
    
            const label length = rowCols.size();
            const label numCols = nConnectedNodes;
    
            // Local accumulation for diagonal to reduce atomic calls
            scalar localDiag[BLOCKSIZE * BLOCKSIZE] = {0.0};

                        // Collect sorted targets
            label targets[ElemType::nodesPerElement];
            for (label jj = 0; jj < numCols; jj++) {
                const label permIndex = sortPermutation[jj];
                targets[jj] = matrixColumnIds[permIndex];
            }
            
            // Vectorized lower_bound for all targets
            label ks[ElemType::nodesPerElement] = {0};
            #pragma unroll
            for (label i = 0; i < length; ++i) {
                #pragma unroll
                for (label d = 0; d < numCols; ++d) {
                    ks[d] += (rowCols[i] < targets[d]);
                }
            }
    
            for (label j = 0; j < numCols; j++) {
                const label offset = ks[j];
                const label curColumnIndex = targets[j];
                const label permIndex = sortPermutation[j];
    
                if (offset < length && rowCols[offset] == curColumnIndex) {
                    for (unsigned m = 0; m < BLOCKSIZE; m++) {
                        for (unsigned n = 0; n < BLOCKSIZE; n++) {
                            const label column_idx = BLOCKSIZE * BLOCKSIZE * offset + BLOCKSIZE * m + n;
                            const label local_idx = BLOCKSIZE * numRows * m + BLOCKSIZE * permIndex + n;
                            Kokkos::atomic_add(&rowVals[column_idx], lhs[lhsCurrentIndex + local_idx]);
                        }
                    }
                } else if (deductUnfound) {
                    // Accumulate to localDiag instead of direct atomic_add
                    for (unsigned m = 0; m < BLOCKSIZE; m++) {
                        for (unsigned n = 0; n < BLOCKSIZE; n++) {
                            const label local_idx = BLOCKSIZE * numRows * m + BLOCKSIZE * permIndex + n;
                            localDiag[BLOCKSIZE * m + n] += lhs[lhsCurrentIndex + local_idx];
                        }
                    }
                }
            }
            
            // Apply accumulated diagonal adds in one go
            scalar* diag = A.diag(rowID);
            for (unsigned mn = 0; mn < BLOCKSIZE * BLOCKSIZE; mn++) {
                Kokkos::atomic_add(&diag[mn], localDiag[mn]);
            }
            
            // RHS adds
            for (unsigned k = 0; k < BLOCKSIZE; k++) {
                Kokkos::atomic_add(&b[BLOCKSIZE * rowID + k], rhs[rhsCurrentIndex + k]);
            }
        }
    }
}

// explicit instantiation
template class AsmCvFemTeamBalance<
    double,
    int,
    1,
    3,
    FEhex<AsmBase<double, int, 1, 3, true, true>::DataType>,
    true,
    true>;

template class AsmCvFemTeamBalance<
    double,
    int,
    3,
    3,
    FEhex<AsmBase<double, int, 3, 3, true, true>::DataType>,
    true,
    true>;