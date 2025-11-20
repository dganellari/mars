#pragma once

#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include "mars.hpp"
#include "mars_sparse_matrix.hpp"
#include "domain.hpp"

namespace mars
{
namespace fem
{

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
        // Count owned nodes (where ownership == 1)
        const auto& ownership = domain_.getNodeOwnershipMap();
        size_t numNodes = domain_.getNodeCount();
            
        thrust::host_vector<uint8_t> h_ownership(numNodes);
        thrust::copy(thrust::device_pointer_cast(ownership.data()),
                    thrust::device_pointer_cast(ownership.data() + numNodes),
                    h_ownership.begin());
            
        size_t ownedCount = 0;
        for (size_t i = 0; i < numNodes; ++i) {
            if (h_ownership[i] == 1) {
                ownedCount++;
            }
        }
            
        numDofs_ = ownedCount;
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

    // Get number of DOFs per element
    static constexpr int dofsPerElement() { return RefElement::numNodes; }

private:
    Domain& domain_;
    int order_;
    size_t numDofs_;
    RefElement refElement_;
};

} // namespace fem
} // namespace mars
