#pragma once

#include "domain.hpp"
#include <map>
#include <array>

namespace mars
{

// ============================================================================
// FaceTopology: Face connectivity and geometry for CVFEM
// Provides face extraction, element-face-element connectivity, and face data
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class FaceTopology
{
public:
    // Number of faces per element (element-type dependent)
    static constexpr int FacesPerElement = std::is_same_v<ElementTag, TetTag> ? 4 : 
                                           std::is_same_v<ElementTag, HexTag> ? 6 : 0;
    static constexpr int NodesPerFace = std::is_same_v<ElementTag, TetTag> ? 3 : 
                                        std::is_same_v<ElementTag, HexTag> ? 4 : 0;

    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    // Face extraction and numbering
    KeyType numFaces_ = 0;
    KeyType numBoundaryFaces_ = 0;
    KeyType numInteriorFaces_ = 0;

    // Element-to-face connectivity (each element has FacesPerElement faces)
    DeviceVector<KeyType> d_elementToFaces_; // [elementCount * FacesPerElement]

    // Face-to-element connectivity (CSR format)
    // Interior faces: 2 elements, boundary faces: 1 element
    DeviceVector<KeyType> d_faceToElementOffsets_; // [numFaces + 1]
    DeviceVector<KeyType> d_faceToElementList_;    // [numInteriorFaces*2 + numBoundaryFaces]

    // Face node connectivity (each face defined by NodesPerFace nodes)
    DeviceVector<KeyType> d_faceNodes_; // [numFaces * NodesPerFace] - SFC keys

    // Face geometric data
    DeviceVector<RealType> d_faceNormalX_;   // [numFaces]
    DeviceVector<RealType> d_faceNormalY_;   // [numFaces]
    DeviceVector<RealType> d_faceNormalZ_;   // [numFaces]
    DeviceVector<RealType> d_faceArea_;      // [numFaces]
    DeviceVector<RealType> d_faceCentroidX_; // [numFaces]
    DeviceVector<RealType> d_faceCentroidY_; // [numFaces]
    DeviceVector<RealType> d_faceCentroidZ_; // [numFaces]

    // Face ownership for MPI (which rank owns this face)
    DeviceVector<int> d_faceOwnerRank_; // [numFaces]

    // Boundary face markers
    DeviceVector<uint8_t> d_isBoundaryFace_; // [numFaces] - 1 if boundary, 0 if interior

    // Build face topology from parent domain
    FaceTopology(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

private:
    void extractFaces(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    void buildFaceConnectivity(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    void computeFaceGeometry(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    void determineFaceOwnership(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    
    // Helper for face definition
    using Face = std::array<KeyType, NodesPerFace>;
    Face makeFace(const std::vector<KeyType>& nodes) const;
};

// ============================================================================
// Implementation
// ============================================================================

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::FaceTopology(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    // Build face topology in stages
    extractFaces(domain);
    buildFaceConnectivity(domain);
    computeFaceGeometry(domain);
    determineFaceOwnership(domain);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
typename FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::Face
FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::makeFace(const std::vector<KeyType>& nodes) const
{
    Face face;
    for (int i = 0; i < NodesPerFace; ++i) {
        face[i] = nodes[i];
    }
    // Sort for canonical representation
    std::sort(face.begin(), face.end());
    return face;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::extractFaces(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // For GPU implementation, we need to extract faces from device data
        // For now, copy to host, process, and copy back
        
        const auto& conn = domain.getElementToNodeConnectivity();
        size_t elementCount = domain.getElementCount();
        
        // Copy connectivity to host for processing
        std::vector<KeyType> h_conn0, h_conn1, h_conn2, h_conn3;
        if constexpr (std::is_same_v<ElementTag, TetTag>) {
            h_conn0.resize(elementCount);
            h_conn1.resize(elementCount);
            h_conn2.resize(elementCount);
            h_conn3.resize(elementCount);
            
            cudaMemcpy(h_conn0.data(), thrust::raw_pointer_cast(std::get<0>(conn).data()), 
                      elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_conn1.data(), thrust::raw_pointer_cast(std::get<1>(conn).data()), 
                      elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_conn2.data(), thrust::raw_pointer_cast(std::get<2>(conn).data()), 
                      elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_conn3.data(), thrust::raw_pointer_cast(std::get<3>(conn).data()), 
                      elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        }
        
        // Map from sorted face nodes to face info
        struct FaceInfo {
            KeyType faceId;
            std::vector<KeyType> elements;
            std::vector<KeyType> nodes; // Original unsorted order from first element
        };
        std::map<Face, FaceInfo> faceMap;
        
        // Extract faces from each element
        for (size_t elemIdx = 0; elemIdx < elementCount; ++elemIdx) {
            std::vector<std::vector<KeyType>> elementFaces;
            
            if constexpr (std::is_same_v<ElementTag, TetTag>) {
                // Tetrahedral faces: (0,1,2), (0,1,3), (0,2,3), (1,2,3)
                elementFaces = {
                    {h_conn0[elemIdx], h_conn1[elemIdx], h_conn2[elemIdx]},
                    {h_conn0[elemIdx], h_conn1[elemIdx], h_conn3[elemIdx]},
                    {h_conn0[elemIdx], h_conn2[elemIdx], h_conn3[elemIdx]},
                    {h_conn1[elemIdx], h_conn2[elemIdx], h_conn3[elemIdx]}
                };
            } else if constexpr (std::is_same_v<ElementTag, HexTag>) {
                // TODO: Implement hex faces
                // Hex faces: 6 quadrilateral faces
            }
            
            for (const auto& faceNodes : elementFaces) {
                Face sortedFace = makeFace(faceNodes);
                
                auto it = faceMap.find(sortedFace);
                if (it == faceMap.end()) {
                    // New face
                    FaceInfo info;
                    info.faceId = faceMap.size();
                    info.elements.push_back(elemIdx);
                    info.nodes = faceNodes; // Store original order
                    faceMap[sortedFace] = info;
                } else {
                    // Existing face - add this element
                    it->second.elements.push_back(elemIdx);
                }
            }
        }
        
        // Count boundary and interior faces
        numFaces_ = faceMap.size();
        numBoundaryFaces_ = 0;
        numInteriorFaces_ = 0;
        
        for (const auto& [face, info] : faceMap) {
            if (info.elements.size() == 1) {
                numBoundaryFaces_++;
            } else {
                numInteriorFaces_++;
            }
        }
        
        std::cout << "Extracted " << numFaces_ << " faces (" 
                  << numBoundaryFaces_ << " boundary, " 
                  << numInteriorFaces_ << " interior)" << std::endl;
        
        // Allocate device memory
        d_faceNodes_.resize(numFaces_ * NodesPerFace);
        d_elementToFaces_.resize(elementCount * FacesPerElement);
        d_faceToElementOffsets_.resize(numFaces_ + 1);
        d_faceToElementList_.resize(numInteriorFaces_ * 2 + numBoundaryFaces_);
        d_isBoundaryFace_.resize(numFaces_);
        
        // TODO: Fill device vectors with extracted face data
        // For now, this is a placeholder
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::buildFaceConnectivity(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    // TODO: Implement face-to-element connectivity building
    std::cout << "FaceTopology::buildFaceConnectivity - Not yet implemented" << std::endl;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::computeFaceGeometry(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if (numFaces_ == 0) return;
    
    // Allocate geometry arrays
    d_faceNormalX_.resize(numFaces_);
    d_faceNormalY_.resize(numFaces_);
    d_faceNormalZ_.resize(numFaces_);
    d_faceArea_.resize(numFaces_);
    d_faceCentroidX_.resize(numFaces_);
    d_faceCentroidY_.resize(numFaces_);
    d_faceCentroidZ_.resize(numFaces_);
    
    // TODO: Implement face geometry computation (normals, areas, centroids)
    std::cout << "FaceTopology::computeFaceGeometry - Not yet implemented" << std::endl;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::determineFaceOwnership(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if (numFaces_ == 0) return;
    
    d_faceOwnerRank_.resize(numFaces_);
    
    // TODO: Implement face ownership determination for MPI
    // A face is owned by the rank that owns its first (or lowest-rank) adjacent element
    std::cout << "FaceTopology::determineFaceOwnership - Not yet implemented" << std::endl;
}

} // namespace mars
