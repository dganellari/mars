#pragma once

#include "mars_sparse_matrix.hpp"
#include "mars_h1_fe_space.hpp"
#include "mars_reference_element.hpp"
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <vector>
#include <algorithm>

namespace mars {
namespace fem {

// Kernel to check mesh quality
template<typename ElementTag, typename RealType, typename IndexType>
__global__ void checkMeshQualityKernel(
    const RealType* node_x,
    const RealType* node_y,
    const RealType* node_z,
    const IndexType* conn0,
    const IndexType* conn1,
    const IndexType* conn2,
    const IndexType* conn3,
    IndexType numElements,
    int* numInverted,
    int* numDegenerate,
    RealType* minDetJ,
    RealType* maxDetJ)
{
    using RefElem = ReferenceElement<ElementTag, RealType>;
    
    IndexType elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;
    
    IndexType dofs[4];
    dofs[0] = conn0[elemIdx];
    dofs[1] = conn1[elemIdx];
    dofs[2] = conn2[elemIdx];
    dofs[3] = conn3[elemIdx];
    
    RealType elem_x[4], elem_y[4], elem_z[4];
    for (int i = 0; i < 4; ++i) {
        elem_x[i] = node_x[dofs[i]];
        elem_y[i] = node_y[dofs[i]];
        elem_z[i] = node_z[dofs[i]];
    }
    
    RealType J[3][3];
    RefElem::computeJacobian(elem_x, elem_y, elem_z, J);
    RealType detJ = RefElem::computeJacobianDeterminant(J);
    
    if (detJ < 0.0) {
        atomicAdd(numInverted, 1);
    }
    if (fabs(detJ) < 1e-12) {
        atomicAdd(numDegenerate, 1);
    }
    
    atomicMin((int*)minDetJ, __float_as_int(detJ));
    atomicMax((int*)maxDetJ, __float_as_int(detJ));
}

// CUDA kernel for stiffness matrix assembly
// Computes local element stiffness matrices and atomically adds to global matrix
template<typename ElementTag, typename RealType, typename IndexType>
__global__ void assembleStiffnessKernel(
    const RealType* node_x,
    const RealType* node_y,
    const RealType* node_z,
    const IndexType* conn0,  // Connectivity node 0
    const IndexType* conn1,  // Connectivity node 1
    const IndexType* conn2,  // Connectivity node 2
    const IndexType* conn3,  // Connectivity node 3
    IndexType numElements,
    const IndexType* rowOffsets,   // CSR row pointers
    const IndexType* colIndices,   // CSR column indices
    RealType* values,              // CSR values (output)
    const IndexType* nodeToDof     // Node-to-DOF mapping (-1 for ghosts)
)
{
    using RefElem = ReferenceElement<ElementTag, RealType>;
    
    IndexType elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Get element nodes from tuple connectivity (already local IDs)
    IndexType dofs[4];
    dofs[0] = conn0[elemIdx];
    dofs[1] = conn1[elemIdx];
    dofs[2] = conn2[elemIdx];
    dofs[3] = conn3[elemIdx];

    // Check if ALL nodes are owned (no ghosts) and map to DOF indices
    bool allOwned = true;
    IndexType ownedDofs[4];
    for (int i = 0; i < 4; ++i) {
        ownedDofs[i] = nodeToDof[dofs[i]];
        if (ownedDofs[i] == static_cast<IndexType>(-1)) {
            allOwned = false;
            break;
        }
    }
    
    if (!allOwned) return;  // Skip elements with ghost nodes

    RealType elem_x[4], elem_y[4], elem_z[4];
    for (int i = 0; i < 4; ++i) {
        elem_x[i] = node_x[dofs[i]];
        elem_y[i] = node_y[dofs[i]];
        elem_z[i] = node_z[dofs[i]];
    }    // Compute Jacobian (constant for tetrahedra)
    RealType J[3][3];
    RefElem::computeJacobian(elem_x, elem_y, elem_z, J);
    
    RealType detJ = RefElem::computeJacobianDeterminant(J);
    
    // Check for degenerate elements
    if (fabs(detJ) < 1e-12) {
        // printf("Warning: Degenerate element %lu with detJ = %f\n", (unsigned long)elemIdx, detJ);
        return;
    }
    
    // Use absolute value to handle inverted elements
    detJ = fabs(detJ);
    
    // Compute inverse Jacobian
    RealType Jinv[3][3];
    RefElem::computeJacobianInverse(J, Jinv);
    
    // Compute local stiffness matrix: K_ij = ∫ ∇φ_i · ∇φ_j dV
    // For linear elements with constant Jacobian, this simplifies to:
    // K_ij = detJ * (∇φ_i_ref^T * Jinv^T * Jinv * ∇φ_j_ref) * (sum of quadrature weights)
    
    RealType localK[RefElem::numNodes][RefElem::numNodes];
    
    // Initialize local matrix
    for (int i = 0; i < RefElem::numNodes; ++i) {
        for (int j = 0; j < RefElem::numNodes; ++j) {
            localK[i][j] = 0.0;
        }
    }
    
    // Numerical integration over quadrature points
    for (int q = 0; q < RefElem::numQuadPoints; ++q) {
        auto qpt = RefElem::getQuadraturePoint(q);
        RealType weight = qpt.weight;
        
        // For each basis function pair
        for (int i = 0; i < RefElem::numNodes; ++i) {
            // Gradient in reference space
            RealType grad_i_ref[3];
            RefElem::evaluateGradient(i, grad_i_ref);
            
            // Transform to physical space
            RealType grad_i[3];
            RefElem::transformGradient(Jinv, grad_i_ref, grad_i);
            
            for (int j = 0; j < RefElem::numNodes; ++j) {
                RealType grad_j_ref[3];
                RefElem::evaluateGradient(j, grad_j_ref);
                
                RealType grad_j[3];
                RefElem::transformGradient(Jinv, grad_j_ref, grad_j);
                
                // Compute dot product of gradients
                RealType gradDot = grad_i[0] * grad_j[0] + 
                                  grad_i[1] * grad_j[1] + 
                                  grad_i[2] * grad_j[2];
                
                localK[i][j] += weight * detJ * gradDot;
            }
        }
    }
    
    // Assemble into global matrix (atomic operations)
    for (int i = 0; i < RefElem::numNodes; ++i) {
        IndexType globalRow = ownedDofs[i];
        IndexType rowStart = rowOffsets[globalRow];
        IndexType rowEnd = rowOffsets[globalRow + 1];
        
        for (int j = 0; j < RefElem::numNodes; ++j) {
            IndexType globalCol = ownedDofs[j];
            
            // Find column index and add atomically
            atomicAddSparseEntry(values, colIndices, rowStart, rowEnd, 
                               globalCol, localK[i][j]);
        }
    }
}

// Stiffness matrix assembler class
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class StiffnessAssembler {
public:
    using FESpace = H1FESpace<ElementTag, RealType, KeyType, AcceleratorTag>;
    using Matrix = SparseMatrix<KeyType, RealType, AcceleratorTag>;
    using Domain = typename FESpace::Domain;
    
    void assemble(FESpace& fes, Matrix& K) {
        auto& domain = fes.domain();
        
        // Build sparsity pattern if not already done
        if (K.nnz() == 0) {
            buildSparsityPattern(fes, K);
        }
        
        // Zero matrix values
        K.zero();
        
        // Get mesh data
        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();
        
        const auto& conn_tuple = domain.getElementToNodeConnectivity();
        const auto& conn0 = std::get<0>(conn_tuple);
        const auto& conn1 = std::get<1>(conn_tuple);
        const auto& conn2 = std::get<2>(conn_tuple);
        const auto& conn3 = std::get<3>(conn_tuple);
        
        size_t numElements = domain.localElementCount();
        int nodesPerElem = FESpace::dofsPerElement();
        
        // Build node-to-DOF mapping
        const auto& ownership = domain.getNodeOwnershipMap();
        thrust::host_vector<uint8_t> h_ownership(domain.getNodeCount());
        thrust::copy(thrust::device_pointer_cast(ownership.data()),
                     thrust::device_pointer_cast(ownership.data() + domain.getNodeCount()),
                     h_ownership.begin());
        
        std::vector<KeyType> nodeToDof(domain.getNodeCount(), static_cast<KeyType>(-1));
        KeyType dofIdx = 0;
        for (size_t nodeIdx = 0; nodeIdx < domain.getNodeCount(); ++nodeIdx) {
            if (h_ownership[nodeIdx] == 1) {
                nodeToDof[nodeIdx] = dofIdx++;
            }
        }
        
        cstone::DeviceVector<KeyType> d_nodeToDof = nodeToDof;
        
        // Launch assembly kernel
        const int blockSize = 256;
        const int gridSize = (numElements + blockSize - 1) / blockSize;
        
        assembleStiffnessKernel<ElementTag, RealType, KeyType>
            <<<gridSize, blockSize>>>(
                d_x.data(),
                d_y.data(),
                d_z.data(),
                conn0.data(),
                conn1.data(),
                conn2.data(),
                conn3.data(),
                numElements,
                K.rowOffsetsPtr(),
                K.colIndicesPtr(),
                K.valuesPtr(),
                d_nodeToDof.data()
            );
        
        cudaDeviceSynchronize();
        
        // Check for errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA error in stiffness assembly: ") + 
                                   cudaGetErrorString(err));
        }
    }
    
    // Check mesh quality
    void checkMeshQuality(FESpace& fes) {
        auto& domain = fes.domain();
        
        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();
        
        const auto& conn_tuple = domain.getElementToNodeConnectivity();
        const auto& conn0 = std::get<0>(conn_tuple);
        const auto& conn1 = std::get<1>(conn_tuple);
        const auto& conn2 = std::get<2>(conn_tuple);
        const auto& conn3 = std::get<3>(conn_tuple);
        
        size_t numElements = domain.localElementCount();
        
        // Allocate device counters
        cstone::DeviceVector<int> d_counters(2, 0);  // [numInverted, numDegenerate]
        cstone::DeviceVector<RealType> d_minmax(2);
        thrust::fill(thrust::device_pointer_cast(d_minmax.data()), 
                    thrust::device_pointer_cast(d_minmax.data() + 1), RealType(1e30));
        thrust::fill(thrust::device_pointer_cast(d_minmax.data() + 1), 
                    thrust::device_pointer_cast(d_minmax.data() + 2), RealType(-1e30));
        
        const int blockSize = 256;
        const int gridSize = (numElements + blockSize - 1) / blockSize;
        
        checkMeshQualityKernel<ElementTag, RealType, KeyType>
            <<<gridSize, blockSize>>>(
                d_x.data(),
                d_y.data(),
                d_z.data(),
                conn0.data(),
                conn1.data(),
                conn2.data(),
                conn3.data(),
                numElements,
                thrust::raw_pointer_cast(d_counters.data()),
                thrust::raw_pointer_cast(d_counters.data() + 1),
                thrust::raw_pointer_cast(d_minmax.data()),
                thrust::raw_pointer_cast(d_minmax.data() + 1)
            );
        
        cudaDeviceSynchronize();
        
        // Copy results to host
        std::vector<int> h_counters(2);
        std::vector<RealType> h_minmax(2);
        thrust::copy(thrust::device_pointer_cast(d_counters.data()),
                    thrust::device_pointer_cast(d_counters.data() + 2),
                    h_counters.begin());
        thrust::copy(thrust::device_pointer_cast(d_minmax.data()),
                    thrust::device_pointer_cast(d_minmax.data() + 2),
                    h_minmax.begin());
        
        std::cout << "   Mesh Quality Check:\n"
                  << "     Inverted elements: " << h_counters[0] << "\n"
                  << "     Degenerate elements: " << h_counters[1] << "\n"
                  << "     Min detJ: " << h_minmax[0] << "\n"
                  << "     Max detJ: " << h_minmax[1] << "\n";
    }
    
private:
    void buildSparsityPattern(FESpace& fes, Matrix& K);
};

// Build sparsity pattern for finite element matrix
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void StiffnessAssembler<ElementTag, RealType, KeyType, AcceleratorTag>::buildSparsityPattern(
    FESpace& fes, Matrix& K) 
{
    auto& domain = fes.domain();
    
    size_t numDofs = fes.numDofs();
    size_t numElements = domain.localElementCount();
    
    // Get connectivity
    const auto& conn_tuple = domain.getElementToNodeConnectivity();
    const auto& conn0 = std::get<0>(conn_tuple);
    const auto& conn1 = std::get<1>(conn_tuple);
    const auto& conn2 = std::get<2>(conn_tuple);
    const auto& conn3 = std::get<3>(conn_tuple);
    
    // Copy connectivity to host for sparsity pattern construction
    thrust::host_vector<KeyType> h_conn0(conn0.size());
    thrust::host_vector<KeyType> h_conn1(conn1.size());
    thrust::host_vector<KeyType> h_conn2(conn2.size());
    thrust::host_vector<KeyType> h_conn3(conn3.size());
    
    thrust::copy(thrust::device_pointer_cast(conn0.data()),
                 thrust::device_pointer_cast(conn0.data() + conn0.size()),
                 h_conn0.begin());
    thrust::copy(thrust::device_pointer_cast(conn1.data()),
                 thrust::device_pointer_cast(conn1.data() + conn1.size()),
                 h_conn1.begin());
    thrust::copy(thrust::device_pointer_cast(conn2.data()),
                 thrust::device_pointer_cast(conn2.data() + conn2.size()),
                 h_conn2.begin());
    thrust::copy(thrust::device_pointer_cast(conn3.data()),
                 thrust::device_pointer_cast(conn3.data() + conn3.size()),
                 h_conn3.begin());
    
    // Get ownership map to filter owned DOFs
    size_t numNodes = domain.getNodeCount();
    const auto& ownership = domain.getNodeOwnershipMap();
    thrust::host_vector<uint8_t> h_ownership(numNodes);
    thrust::copy(thrust::device_pointer_cast(ownership.data()),
                 thrust::device_pointer_cast(ownership.data() + numNodes),
                 h_ownership.begin());
    
    // Build node-to-DOF mapping (owned nodes only)
    std::vector<KeyType> nodeToLocalDof(numNodes, static_cast<KeyType>(-1));
    KeyType dofIdx = 0;
    for (size_t i = 0; i < numNodes; ++i) {
        if (h_ownership[i] == 1) {
            nodeToLocalDof[i] = dofIdx++;
        }
    }
    
    // Build sparsity pattern by collecting all DOF pairs that appear in elements
    std::vector<std::vector<KeyType>> rowCols(numDofs);
    
    // For each element, add all DOF pairs (skip elements with ghost nodes)
    for (size_t elem = 0; elem < numElements; ++elem) {
        KeyType nodes[4] = {h_conn0[elem], h_conn1[elem], h_conn2[elem], h_conn3[elem]};
        
        // Check if element has any ghost nodes
        bool hasGhost = false;
        for (int i = 0; i < 4; ++i) {
            if (h_ownership[nodes[i]] != 1) {
                hasGhost = true;
                break;
            }
        }
        if (hasGhost) continue;
        
        // Get local DOF indices for this element
        KeyType localDofs[4];
        for (int i = 0; i < 4; ++i) {
            localDofs[i] = nodeToLocalDof[nodes[i]];
        }
        
        // Add all pairs within this element
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                rowCols[localDofs[i]].push_back(localDofs[j]);
            }
        }
    }
    
    // Remove duplicates and sort each row
    size_t totalNnz = 0;
    for (size_t i = 0; i < numDofs; ++i) {
        auto& row = rowCols[i];
        std::sort(row.begin(), row.end());
        auto last = std::unique(row.begin(), row.end());
        row.erase(last, row.end());
        totalNnz += row.size();
    }
    
    // Allocate matrix
    K.allocate(numDofs, numDofs, totalNnz);
    
    // Build CSR format
    auto rowOffsets = K.rowOffsetsPtr();
    auto colIndices = K.colIndicesPtr();
    
    // Copy to device
    thrust::host_vector<KeyType> h_rowOffsets(numDofs + 1, 0);
    thrust::host_vector<KeyType> h_colIndices(totalNnz);
    
    size_t offset = 0;
    for (size_t i = 0; i < numDofs; ++i) {
        h_rowOffsets[i] = offset;
        const auto& row = rowCols[i];
        for (size_t j = 0; j < row.size(); ++j) {
            h_colIndices[offset + j] = row[j];
        }
        offset += row.size();
    }
    h_rowOffsets[numDofs] = offset;
    
    // Copy to device
    thrust::copy(h_rowOffsets.begin(), h_rowOffsets.end(), 
                 thrust::device_pointer_cast(rowOffsets));
    thrust::copy(h_colIndices.begin(), h_colIndices.end(),
                 thrust::device_pointer_cast(colIndices));
}

} // namespace fem
} // namespace mars
