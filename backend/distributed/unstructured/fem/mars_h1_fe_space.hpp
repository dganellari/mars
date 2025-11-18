#pragma once

#include <set>
#include <vector>
#include "mars.hpp"
#include "mars_sparse_matrix.hpp"
#include "domain.hpp"

namespace mars
{
namespace fem
{

// Helper kernel to flatten connectivity tuple into DOF map
template<typename KeyType>
__global__ void flattenConnectivityToDofMapKernel(KeyType* d_dofMap,
                                                  const KeyType* conn0,
                                                  const KeyType* conn1,
                                                  const KeyType* conn2,
                                                  const KeyType* conn3,
                                                  size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e < numElements)
    {
        d_dofMap[e * 4 + 0] = conn0[e];
        d_dofMap[e * 4 + 1] = conn1[e];
        d_dofMap[e * 4 + 2] = conn2[e];
        d_dofMap[e * 4 + 3] = conn3[e];
    }
}

// H1 Finite Element Space (conforming, continuous)
// For linear Lagrange elements: DOFs are at nodes
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class H1FESpace
{
public:
    using Domain           = ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>;
    using DeviceVectorKey  = typename mars::VectorSelector<KeyType, AcceleratorTag>::type;
    using DeviceVectorReal = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    using RefElement       = ReferenceElement<ElementTag, RealType>;

    H1FESpace(Domain& domain, int order = 1)
        : domain_(domain)
        , order_(order)
    {

        if (order != 1) { throw std::runtime_error("Only linear (order=1) elements currently supported"); }

        // For linear elements, DOFs = nodes
        buildDofMap();
    }

    // Get total number of degrees of freedom
    size_t numDofs() const { return numDofs_; }

    // Get polynomial order
    int order() const { return order_; }

    // Get reference element
    const RefElement& referenceElement() const { return refElement_; }

    // Get underlying mesh domain
    Domain& domain() { return domain_; }
    const Domain& domain() const { return domain_; }

    // Get DOF map (element DOF -> global DOF)
    const DeviceVectorKey& dofMap() const { return d_dofMap_; }

    // Identify boundary DOFs (nodes on domain boundary)
    std::vector<KeyType> getBoundaryDofs() const;

    // Get number of DOFs per element
    static constexpr int dofsPerElement() { return RefElement::numNodes; }

    // Build DOF map (public to allow CUDA kernel call)
    void buildDofMap()
    {
        const size_t numElements = domain_.localElementCount();
        const size_t numNodes    = domain_.getNodeCount();

        numDofs_ = numNodes;

        size_t totalDofs = numElements * dofsPerElement();
        d_dofMap_.resize(totalDofs);

        // Get local ID connectivity (built lazily) - returns tuple of vectors
        const auto& conn_tuple = domain_.getElementToNodeConnectivity();

        // Flatten tuple connectivity into DOF map using CUDA kernel
        const auto& conn0 = std::get<0>(conn_tuple);
        const auto& conn1 = std::get<1>(conn_tuple);
        const auto& conn2 = std::get<2>(conn_tuple);
        const auto& conn3 = std::get<3>(conn_tuple);

        // Launch kernel to interleave connectivity
        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        flattenConnectivityToDofMapKernel<<<numBlocks, blockSize>>>(d_dofMap_.data(), conn0.data(), conn1.data(),
                                                                    conn2.data(), conn3.data(), numElements);
        cudaDeviceSynchronize();
    }

private:
    Domain& domain_;
    int order_;
    size_t numDofs_;

    DeviceVectorKey d_dofMap_; // Element DOF layout (flat: elem0_dof0, elem0_dof1, ..., elem1_dof0, ...)
    RefElement refElement_;
};

// Implementation of getBoundaryDofs (out-of-line with template declaration)
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
std::vector<KeyType> H1FESpace<ElementTag, RealType, KeyType, AcceleratorTag>::getBoundaryDofs() const
{
    const auto& haloIndices = domain_.getHaloElementIndices();

    std::vector<KeyType> h_haloIndices(haloIndices.size());
    thrust::copy(thrust::device_pointer_cast(haloIndices.data()),
                 thrust::device_pointer_cast(haloIndices.data() + haloIndices.size()), h_haloIndices.begin());

    std::set<KeyType> boundaryDofSet;

    // Get connectivity tuple
    const auto& conn_tuple = domain_.getElementToNodeConnectivity();

    // Copy each component to host
    const auto& conn0 = std::get<0>(conn_tuple);
    const auto& conn1 = std::get<1>(conn_tuple);
    const auto& conn2 = std::get<2>(conn_tuple);
    const auto& conn3 = std::get<3>(conn_tuple);

    std::vector<KeyType> h_conn0(conn0.size());
    std::vector<KeyType> h_conn1(conn1.size());
    std::vector<KeyType> h_conn2(conn2.size());
    std::vector<KeyType> h_conn3(conn3.size());

    thrust::copy(thrust::device_pointer_cast(conn0.data()), thrust::device_pointer_cast(conn0.data() + conn0.size()),
                 h_conn0.begin());
    thrust::copy(thrust::device_pointer_cast(conn1.data()), thrust::device_pointer_cast(conn1.data() + conn1.size()),
                 h_conn1.begin());
    thrust::copy(thrust::device_pointer_cast(conn2.data()), thrust::device_pointer_cast(conn2.data() + conn2.size()),
                 h_conn2.begin());
    thrust::copy(thrust::device_pointer_cast(conn3.data()), thrust::device_pointer_cast(conn3.data() + conn3.size()),
                 h_conn3.begin());

    const int nodesPerElem = dofsPerElement();
    for (auto elemIdx : h_haloIndices)
    {
        // Add all nodes of this halo element to boundary set
        boundaryDofSet.insert(h_conn0[elemIdx]);
        boundaryDofSet.insert(h_conn1[elemIdx]);
        boundaryDofSet.insert(h_conn2[elemIdx]);
        boundaryDofSet.insert(h_conn3[elemIdx]);
    }

    return std::vector<KeyType>(boundaryDofSet.begin(), boundaryDofSet.end());
}

} // namespace fem
} // namespace mars
