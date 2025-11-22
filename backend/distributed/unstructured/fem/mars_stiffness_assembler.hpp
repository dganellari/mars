#pragma once

#include "mars_sparse_matrix.hpp"
#include "mars_h1_fe_space.hpp"
#include "mars_reference_element.hpp"
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <vector>
#include <unordered_set>
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
    const IndexType* rowOffsets,      // CSR row pointers (local DOF indexing)
    const IndexType* colIndices,      // CSR column indices (local DOF indexing with ghosts)
    RealType* values,                 // CSR values (output)
    const IndexType* nodeToLocalDof,  // Node to local DOF mapping (includes ghosts)
    const uint8_t* nodeOwnership      // Node ownership (1 = owned, 0 = ghost)
)
{
    using RefElem = ReferenceElement<ElementTag, RealType>;
    
    IndexType elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Get element nodes from tuple connectivity (local node indices)
    IndexType nodes[4];
    nodes[0] = conn0[elemIdx];
    nodes[1] = conn1[elemIdx];
    nodes[2] = conn2[elemIdx];
    nodes[3] = conn3[elemIdx];

    // Check if element has at least one owned node
    bool hasOwnedNode = false;
    for (int i = 0; i < 4; ++i) {
        if (nodeOwnership[nodes[i]] == 1) {
            hasOwnedNode = true;
            break;
        }
    }
    if (!hasOwnedNode) return;  // Skip elements with no owned nodes
    
    // Map nodes to local DOF indices (includes both owned and ghost DOFs)
    IndexType localDofs[4];
    for (int i = 0; i < 4; ++i) {
        localDofs[i] = nodeToLocalDof[nodes[i]];
    }

    RealType elem_x[4], elem_y[4], elem_z[4];
    for (int i = 0; i < 4; ++i) {
        elem_x[i] = node_x[nodes[i]];
        elem_y[i] = node_y[nodes[i]];
        elem_z[i] = node_z[nodes[i]];
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
    // Rows: only owned nodes; Columns: owned or ghost nodes (local indices)
    for (int i = 0; i < RefElem::numNodes; ++i) {
        // Skip if this node is not owned
        if (nodeOwnership[nodes[i]] != 1) continue;
        
        IndexType localRow = localDofs[i];     // Local row index (owned DOF)
        IndexType rowStart = rowOffsets[localRow];
        IndexType rowEnd = rowOffsets[localRow + 1];
        
        for (int j = 0; j < RefElem::numNodes; ++j) {
            // Column uses local index (can reference owned or ghost DOF)
            IndexType localCol = localDofs[j];
            
            // Find column index and add atomically
            atomicAddSparseEntry(values, colIndices, rowStart, rowEnd, 
                               localCol, localK[i][j]);
        }
    }
}

// Stiffness matrix assembler class
// 
// Production-ready distributed implementation with LOCAL-LOCAL indexing:
//   - Rows: LOCAL DOF indices (0 to numLocalDofs-1) - owned DOFs only
//   - Columns: LOCAL DOF indices (0 to numLocalDofs+numGhostDofs-1) - includes ghosts
//
// Matrix structure for multi-rank:
//   - Square matrix: (numLocalDofs) x (numLocalDofs + numGhostDofs)
//   - SpMV compatible: all column indices are in local storage range
//   - Requires halo exchange after SpMV to communicate ghost updates
//
// Single-rank optimization:
//   - When numRanks==1, no ghosts exist, reduces to standard local assembly
//
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class StiffnessAssembler {
public:
    using FESpace = H1FESpace<ElementTag, RealType, KeyType, AcceleratorTag>;
    using Matrix = SparseMatrix<KeyType, RealType, AcceleratorTag>;
    using Domain = typename FESpace::Domain;
    
    void assemble(FESpace& fes, Matrix& K, const std::vector<KeyType>& nodeToLocalDof) {
        auto& domain = fes.domain();
        
        // Build sparsity pattern if not already done
        if (K.nnz() == 0) {
            buildSparsityPattern(fes, K, nodeToLocalDof);
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
        
        // Use provided node-to-DOF mapping
        cstone::DeviceVector<KeyType> d_nodeToLocalDof = nodeToLocalDof;
        
        // Get ownership map
        const auto& ownership = domain.getNodeOwnershipMap();
        
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
                d_nodeToLocalDof.data(),
                ownership.data()
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
    void buildSparsityPattern(FESpace& fes, Matrix& K, const std::vector<KeyType>& nodeToLocalDof);
};

// Build sparsity pattern for finite element matrix
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void StiffnessAssembler<ElementTag, RealType, KeyType, AcceleratorTag>::buildSparsityPattern(
    FESpace& fes, Matrix& K, const std::vector<KeyType>& nodeToLocalDof) 
{
    auto& domain = fes.domain();
    
    size_t numDofs = fes.numDofs();
    size_t numElements = domain.localElementCount();
    size_t numNodes = domain.getNodeCount();
    
    // Validate input sizes
    if (nodeToLocalDof.size() != numNodes) {
        std::cerr << "Error: nodeToLocalDof.size() " << nodeToLocalDof.size() 
                  << " != numNodes " << numNodes << std::endl;
        return;
    }
    
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
    
    // Get ownership map to check which elements are fully owned
    const auto& ownership = domain.getNodeOwnershipMap();
    thrust::host_vector<uint8_t> h_ownership(numNodes);
    thrust::copy(thrust::device_pointer_cast(ownership.data()),
                 thrust::device_pointer_cast(ownership.data() + numNodes),
                 h_ownership.begin());
    
    // Build sparsity pattern:
    // - Rows: LOCAL owned DOF indices (0 to numDofs-1)
    // - Columns: LOCAL DOF indices including ghosts (0 to numDofs+numGhostDofs-1)
    //
    // Count total local DOFs (owned + ghost)
    KeyType maxLocalDof = 0;
    for (size_t i = 0; i < nodeToLocalDof.size(); ++i) {
        if (nodeToLocalDof[i] > maxLocalDof) {
            maxLocalDof = nodeToLocalDof[i];
        }
    }
    size_t numLocalDofsWithGhosts = maxLocalDof + 1;  // Total local storage
    
    std::vector<std::vector<KeyType>> rowCols(numDofs);  // Only rows for owned DOFs
    
    // For each element, add all DOF pairs
    // Only rows corresponding to owned nodes, columns can be owned or ghost
    for (size_t elem = 0; elem < numElements; ++elem) {
        KeyType nodes[4] = {h_conn0[elem], h_conn1[elem], h_conn2[elem], h_conn3[elem]};
        
        // Validate node indices
        bool invalidNodes = false;
        for (int i = 0; i < 4; ++i) {
            if (nodes[i] >= numNodes) {
                invalidNodes = true;
                break;
            }
        }
        if (invalidNodes) continue;
        
        // Check if element has any owned nodes
        bool hasOwnedNode = false;
        for (int i = 0; i < 4; ++i) {
            if (h_ownership[nodes[i]] == 1) {
                hasOwnedNode = true;
                break;
            }
        }
        if (!hasOwnedNode) continue;
        
        // Get local DOF indices for this element (includes ghosts)
        KeyType localDofs[4];
        bool invalidDofs = false;
        
        for (int i = 0; i < 4; ++i) {
            if (nodes[i] >= nodeToLocalDof.size()) {
                invalidDofs = true;
                break;
            }
            localDofs[i] = nodeToLocalDof[nodes[i]];
        }
        if (invalidDofs) continue;
        
        // Add all pairs: for each owned node (row), add coupling to all nodes (columns)
        for (int i = 0; i < 4; ++i) {
            // Only add rows for owned nodes
            if (h_ownership[nodes[i]] != 1) continue;
            
            KeyType localRow = localDofs[i];
            if (localRow >= numDofs) continue;  // Safety check
            
            for (int j = 0; j < 4; ++j) {
                // Column can be owned or ghost - use local index
                KeyType localCol = localDofs[j];
                rowCols[localRow].push_back(localCol);
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
    
    // Allocate matrix (rectangular: rows=owned DOFs, cols=owned+ghost DOFs)
    // Ghost DOF values will be updated via halo exchange before each SpMV
    K.allocate(numDofs, numLocalDofsWithGhosts, totalNnz);
    
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
