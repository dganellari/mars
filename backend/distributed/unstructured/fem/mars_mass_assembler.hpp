#pragma once

#include "mars_sparse_matrix.hpp"
#include "mars_h1_fe_space.hpp"
#include "mars_reference_element.hpp"
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <unordered_set>

namespace mars
{
namespace fem
{

// CUDA kernel for mass matrix assembly
// Computes M_ij = ∫ φ_i φ_j dV
template<typename ElementTag, typename RealType, typename IndexType>
__global__ void assembleMassKernel(const RealType* node_x,
                                   const RealType* node_y,
                                   const RealType* node_z,
                                   const IndexType* conn0,
                                   const IndexType* conn1,
                                   const IndexType* conn2,
                                   const IndexType* conn3,
                                   IndexType numElements,
                                   const IndexType* rowOffsets,
                                   const IndexType* colIndices,
                                   RealType* values,
                                   const IndexType* nodeToDof)
{
    using RefElem = ReferenceElement<ElementTag, RealType>;

    IndexType elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Get element nodes from tuple connectivity (already local IDs)
    IndexType nodes[4];
    nodes[0] = conn0[elemIdx];
    nodes[1] = conn1[elemIdx];
    nodes[2] = conn2[elemIdx];
    nodes[3] = conn3[elemIdx];
    
    // Check if ALL nodes are owned and map to DOF indices
    bool allOwned = true;
    IndexType dofs[4];
    for (int i = 0; i < 4; ++i) {
        dofs[i] = nodeToDof[nodes[i]];
        if (dofs[i] == static_cast<IndexType>(-1)) {
            allOwned = false;
            break;
        }
    }
    
    if (!allOwned) return;  // Skip elements with ghost nodes

    RealType elem_x[4], elem_y[4], elem_z[4];
    for (int i = 0; i < 4; ++i)
    {
        elem_x[i] = node_x[dofs[i]];
        elem_y[i] = node_y[dofs[i]];
        elem_z[i] = node_z[dofs[i]];
    }

    // Compute Jacobian determinant
    RealType J[3][3];
    RefElem::computeJacobian(elem_x, elem_y, elem_z, J);
    RealType detJ = RefElem::computeJacobianDeterminant(J);

    if (detJ <= 0.0) return;

    // Local mass matrix
    RealType localM[RefElem::numNodes][RefElem::numNodes];

    for (int i = 0; i < RefElem::numNodes; ++i)
    {
        for (int j = 0; j < RefElem::numNodes; ++j)
        {
            localM[i][j] = 0.0;
        }
    }

    // Quadrature integration
    for (int q = 0; q < RefElem::numQuadPoints; ++q)
    {
        auto qpt        = RefElem::getQuadraturePoint(q);
        RealType xi     = qpt.xi;
        RealType eta    = qpt.eta;
        RealType zeta   = qpt.zeta;
        RealType weight = qpt.weight;

        for (int i = 0; i < RefElem::numNodes; ++i)
        {
            RealType phi_i = RefElem::evaluateBasis(i, xi, eta, zeta);

            for (int j = 0; j < RefElem::numNodes; ++j)
            {
                RealType phi_j = RefElem::evaluateBasis(j, xi, eta, zeta);

                localM[i][j] += weight * detJ * phi_i * phi_j;
            }
        }
    }

    // Assemble into global matrix
    for (int i = 0; i < RefElem::numNodes; ++i)
    {
        IndexType globalRow = dofs[i];
        IndexType rowStart  = rowOffsets[globalRow];
        IndexType rowEnd    = rowOffsets[globalRow + 1];

        for (int j = 0; j < RefElem::numNodes; ++j)
        {
            IndexType globalCol = dofs[j];
            atomicAddSparseEntry(values, colIndices, rowStart, rowEnd, globalCol, localM[i][j]);
        }
    }
}

// CUDA kernel for RHS vector assembly
// Computes b_i = ∫ f φ_i dV where f is source term
template<typename ElementTag, typename RealType, typename IndexType, typename SourceFunc>
__global__ void assembleRHSKernel(const RealType* node_x,
                                  const RealType* node_y,
                                  const RealType* node_z,
                                  const IndexType* conn0,
                                  const IndexType* conn1,
                                  const IndexType* conn2,
                                  const IndexType* conn3,
                                  IndexType numElements,
                                  SourceFunc sourceTerm,
                                  RealType* rhs,
                                  const IndexType* nodeToDof)
{
    using RefElem = ReferenceElement<ElementTag, RealType>;

    IndexType elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    // Get element nodes from tuple connectivity (already local IDs)
    IndexType nodes[4];
    nodes[0] = conn0[elemIdx];
    nodes[1] = conn1[elemIdx];
    nodes[2] = conn2[elemIdx];
    nodes[3] = conn3[elemIdx];
    
    // Check if ALL nodes are owned and map to DOF indices
    bool allOwned = true;
    IndexType dofs[4];
    for (int i = 0; i < 4; ++i) {
        dofs[i] = nodeToDof[nodes[i]];
        if (dofs[i] == static_cast<IndexType>(-1)) {
            allOwned = false;
            break;
        }
    }
    
    if (!allOwned) return;  // Skip elements with ghost nodes

    RealType elem_x[4], elem_y[4], elem_z[4];
    for (int i = 0; i < 4; ++i)
    {
        elem_x[i] = node_x[nodes[i]];
        elem_y[i] = node_y[nodes[i]];
        elem_z[i] = node_z[nodes[i]];
    }

    // Compute Jacobian
    RealType J[3][3];
    RefElem::computeJacobian(elem_x, elem_y, elem_z, J);
    RealType detJ = RefElem::computeJacobianDeterminant(J);

    if (detJ == 0.0) return;  // Skip degenerate elements
    detJ = std::abs(detJ);    // Use absolute value to handle inverted elements

    // Local RHS vector
    RealType localB[RefElem::numNodes];
    for (int i = 0; i < RefElem::numNodes; ++i)
    {
        localB[i] = 0.0;
    }

    // Quadrature integration
    for (int q = 0; q < RefElem::numQuadPoints; ++q)
    {
        auto qpt        = RefElem::getQuadraturePoint(q);
        RealType xi     = qpt.xi;
        RealType eta    = qpt.eta;
        RealType zeta   = qpt.zeta;
        RealType weight = qpt.weight;

        // Map quadrature point to physical coordinates
        RealType x = 0.0, y = 0.0, z = 0.0;
        for (int i = 0; i < RefElem::numNodes; ++i)
        {
            RealType phi = RefElem::evaluateBasis(i, xi, eta, zeta);
            x += phi * elem_x[i];
            y += phi * elem_y[i];
            z += phi * elem_z[i];
        }

        // Evaluate source term at quadrature point
        RealType f = sourceTerm(x, y, z);

        // Integrate
        for (int i = 0; i < RefElem::numNodes; ++i)
        {
            RealType phi_i = RefElem::evaluateBasis(i, xi, eta, zeta);
            localB[i] += weight * detJ * f * phi_i;
        }
    }

    // Assemble into global vector
    for (int i = 0; i < RefElem::numNodes; ++i)
    {
        atomicAdd(&rhs[dofs[i]], localB[i]);
    }
}

// Mass matrix assembler class
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class MassAssembler
{
public:
    using FESpace = H1FESpace<ElementTag, RealType, KeyType, AcceleratorTag>;
    using Matrix  = SparseMatrix<KeyType, RealType, AcceleratorTag>;
    using Vector  = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    using Domain  = typename FESpace::Domain;

    void assemble(FESpace& fes, Matrix& M, const std::vector<KeyType>& nodeToDof)
    {
        auto& domain = fes.domain();

        // Build sparsity pattern if needed (same as stiffness)
        if (M.nnz() == 0)
        {
            // Use same sparsity pattern as stiffness matrix
            // In practice, share this with StiffnessAssembler
            size_t numDofs      = fes.numDofs();
            size_t estimatedNnz = numDofs * 15;
            M.allocate(numDofs, numDofs, estimatedNnz);
        }

        M.zero();

        // Get mesh data
        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();

        const auto& conn_tuple = domain.getElementToNodeConnectivity();
        const auto& conn0      = std::get<0>(conn_tuple);
        const auto& conn1      = std::get<1>(conn_tuple);
        const auto& conn2      = std::get<2>(conn_tuple);
        const auto& conn3      = std::get<3>(conn_tuple);

        size_t numElements = domain.localElementCount();
        
        // Use provided node-to-DOF mapping
        cstone::DeviceVector<KeyType> d_nodeToDof = nodeToDof;

        // Launch kernel
        const int blockSize = 256;
        const int gridSize  = (numElements + blockSize - 1) / blockSize;

        assembleMassKernel<ElementTag, RealType, KeyType>
            <<<gridSize, blockSize>>>(d_x.data(), d_y.data(), d_z.data(), conn0.data(), conn1.data(), conn2.data(),
                                      conn3.data(), numElements, M.rowOffsetsPtr(), M.colIndicesPtr(), M.valuesPtr(),
                                      d_nodeToDof.data());

        cudaDeviceSynchronize();
    }

    // Assemble RHS vector with source term
    template<typename SourceFunc>
    void assembleRHS(FESpace& fes, Vector& b, SourceFunc f, const std::vector<KeyType>& nodeToDof)
    {
        auto& domain = fes.domain();

        size_t numDofs = fes.numDofs();
        b.resize(numDofs);
        thrust::fill(thrust::device_pointer_cast(b.data()), thrust::device_pointer_cast(b.data() + b.size()),
                     RealType(0));

        // Get mesh data
        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();

        const auto& conn_tuple = domain.getElementToNodeConnectivity();
        const auto& conn0      = std::get<0>(conn_tuple);
        const auto& conn1      = std::get<1>(conn_tuple);
        const auto& conn2      = std::get<2>(conn_tuple);
        const auto& conn3      = std::get<3>(conn_tuple);

        size_t numElements = domain.localElementCount();
        
        // Use provided node-to-DOF mapping
        cstone::DeviceVector<KeyType> d_nodeToDof = nodeToDof;

        // Launch kernel
        const int blockSize = 256;
        const int gridSize  = (numElements + blockSize - 1) / blockSize;

        assembleRHSKernel<ElementTag, RealType, KeyType>
            <<<gridSize, blockSize>>>(d_x.data(), d_y.data(), d_z.data(), conn0.data(), conn1.data(), conn2.data(), conn3.data(),
                                      numElements, f, thrust::raw_pointer_cast(b.data()), d_nodeToDof.data());

        cudaDeviceSynchronize();
    }
};

} // namespace fem
} // namespace mars
