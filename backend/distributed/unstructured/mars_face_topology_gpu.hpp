#pragma once

#include "domain.hpp"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/memory.h>
#include <map>
#include <array>

namespace mars
{

// Struct to encode both face and element ID
template<typename KeyType>
struct FaceElementPair {
    KeyType faceHash;  // Unique hash for face (based on sorted nodes)
    KeyType node0, node1, node2, node3;  // Sorted face nodes (max 4 for hex)
    KeyType elementId;
    
    __device__ __host__ bool operator<(const FaceElementPair& other) const {
        if (node0 != other.node0) return node0 < other.node0;
        if (node1 != other.node1) return node1 < other.node1;
        if (node2 != other.node2) return node2 < other.node2;
        return node3 < other.node3;
    }
    
    __device__ __host__ bool operator==(const FaceElementPair& other) const {
        return node0 == other.node0 && node1 == other.node1 && 
               node2 == other.node2 && node3 == other.node3;
    }
};

// GPU kernel for extracting face-element pairs from hex elements
template<typename KeyType>
__global__ void extractHexFaceElementPairsKernel(
    const KeyType* __restrict__ conn0,
    const KeyType* __restrict__ conn1,
    const KeyType* __restrict__ conn2,
    const KeyType* __restrict__ conn3,
    const KeyType* __restrict__ conn4,
    const KeyType* __restrict__ conn5,
    const KeyType* __restrict__ conn6,
    const KeyType* __restrict__ conn7,
    size_t numElements,
    FaceElementPair<KeyType>* faceElemPairs)  // Output: [numElements * 6] (6 quad faces per hex)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;
    
    KeyType n0 = conn0[elemIdx];
    KeyType n1 = conn1[elemIdx];
    KeyType n2 = conn2[elemIdx];
    KeyType n3 = conn3[elemIdx];
    KeyType n4 = conn4[elemIdx];
    KeyType n5 = conn5[elemIdx];
    KeyType n6 = conn6[elemIdx];
    KeyType n7 = conn7[elemIdx];
    
    // Helper to sort 4 nodes (quad face)
    auto sortFace = [](KeyType a, KeyType b, KeyType c, KeyType d, 
                       KeyType& s0, KeyType& s1, KeyType& s2, KeyType& s3) {
        s0 = a; s1 = b; s2 = c; s3 = d;
        // Bubble sort for 4 elements
        if (s0 > s1) { KeyType tmp = s0; s0 = s1; s1 = tmp; }
        if (s2 > s3) { KeyType tmp = s2; s2 = s3; s3 = tmp; }
        if (s0 > s2) { KeyType tmp = s0; s0 = s2; s2 = tmp; }
        if (s1 > s3) { KeyType tmp = s1; s1 = s3; s3 = tmp; }
        if (s1 > s2) { KeyType tmp = s1; s1 = s2; s2 = tmp; }
    };
    
    int baseIdx = elemIdx * 6;
    KeyType dummy = 0; // 4th node for struct (only uses 3, but hex has 4-node faces)
    
    // Hex faces (6 quadrilateral faces)
    // Face 0: bottom (0,1,2,3)
    sortFace(n0, n1, n2, n3,
             faceElemPairs[baseIdx + 0].node0,
             faceElemPairs[baseIdx + 0].node1,
             faceElemPairs[baseIdx + 0].node2,
             dummy);
    faceElemPairs[baseIdx + 0].elementId = elemIdx;
    
    // Face 1: top (4,5,6,7)
    sortFace(n4, n5, n6, n7,
             faceElemPairs[baseIdx + 1].node0,
             faceElemPairs[baseIdx + 1].node1,
             faceElemPairs[baseIdx + 1].node2,
             dummy);
    faceElemPairs[baseIdx + 1].elementId = elemIdx;
    
    // Face 2: front (0,1,5,4)
    sortFace(n0, n1, n5, n4,
             faceElemPairs[baseIdx + 2].node0,
             faceElemPairs[baseIdx + 2].node1,
             faceElemPairs[baseIdx + 2].node2,
             dummy);
    faceElemPairs[baseIdx + 2].elementId = elemIdx;
    
    // Face 3: back (2,3,7,6)
    sortFace(n2, n3, n7, n6,
             faceElemPairs[baseIdx + 3].node0,
             faceElemPairs[baseIdx + 3].node1,
             faceElemPairs[baseIdx + 3].node2,
             dummy);
    faceElemPairs[baseIdx + 3].elementId = elemIdx;
    
    // Face 4: left (0,3,7,4)
    sortFace(n0, n3, n7, n4,
             faceElemPairs[baseIdx + 4].node0,
             faceElemPairs[baseIdx + 4].node1,
             faceElemPairs[baseIdx + 4].node2,
             dummy);
    faceElemPairs[baseIdx + 4].elementId = elemIdx;
    
    // Face 5: right (1,2,6,5)
    sortFace(n1, n2, n6, n5,
             faceElemPairs[baseIdx + 5].node0,
             faceElemPairs[baseIdx + 5].node1,
             faceElemPairs[baseIdx + 5].node2,
             dummy);
    faceElemPairs[baseIdx + 5].elementId = elemIdx;
}

// GPU kernel for extracting face-element pairs from tet elements
template<typename KeyType>
__global__ void extractTetFaceElementPairsKernel(
    const KeyType* __restrict__ conn0,
    const KeyType* __restrict__ conn1,
    const KeyType* __restrict__ conn2,
    const KeyType* __restrict__ conn3,
    size_t numElements,
    FaceElementPair<KeyType>* faceElemPairs)  // Output: [numElements * 4]
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;
    
    KeyType n0 = conn0[elemIdx];
    KeyType n1 = conn1[elemIdx];
    KeyType n2 = conn2[elemIdx];
    KeyType n3 = conn3[elemIdx];
    
    // Helper to sort 3 nodes
    auto sortFace = [](KeyType a, KeyType b, KeyType c, KeyType& s0, KeyType& s1, KeyType& s2) {
        s0 = a; s1 = b; s2 = c;
        if (s0 > s1) { KeyType tmp = s0; s0 = s1; s1 = tmp; }
        if (s1 > s2) { KeyType tmp = s1; s1 = s2; s2 = tmp; }
        if (s0 > s1) { KeyType tmp = s0; s0 = s1; s1 = tmp; }
    };
    
    int baseIdx = elemIdx * 4;
    
    // Face 0: nodes (0,1,2)
    sortFace(n0, n1, n2, 
             faceElemPairs[baseIdx + 0].node0,
             faceElemPairs[baseIdx + 0].node1,
             faceElemPairs[baseIdx + 0].node2);
    faceElemPairs[baseIdx + 0].elementId = elemIdx;
    
    // Face 1: nodes (0,1,3)
    sortFace(n0, n1, n3,
             faceElemPairs[baseIdx + 1].node0,
             faceElemPairs[baseIdx + 1].node1,
             faceElemPairs[baseIdx + 1].node2);
    faceElemPairs[baseIdx + 1].elementId = elemIdx;
    
    // Face 2: nodes (0,2,3)
    sortFace(n0, n2, n3,
             faceElemPairs[baseIdx + 2].node0,
             faceElemPairs[baseIdx + 2].node1,
             faceElemPairs[baseIdx + 2].node2);
    faceElemPairs[baseIdx + 2].elementId = elemIdx;
    
    // Face 3: nodes (1,2,3)
    sortFace(n1, n2, n3,
             faceElemPairs[baseIdx + 3].node0,
             faceElemPairs[baseIdx + 3].node1,
             faceElemPairs[baseIdx + 3].node2);
    faceElemPairs[baseIdx + 3].elementId = elemIdx;
}

// Kernel to build CSR face-to-element connectivity from sorted face-element pairs
template<typename KeyType>
__global__ void buildFaceToElementCSRKernel(
    const FaceElementPair<KeyType>* __restrict__ sortedPairs,
    size_t numPairs,
    KeyType* faceNodes,           // Output: [numFaces * 3]
    KeyType* faceToElemOffsets,   // Output: [numFaces + 1]
    KeyType* faceToElemList,      // Output: [numPairs]
    uint8_t* isBoundary,          // Output: [numFaces]
    KeyType* numFacesOut)         // Output: total unique faces
{
    // This kernel processes the sorted face-element pairs
    // Consecutive pairs with same face nodes belong to same face
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numPairs) return;
    
    bool isFirstInGroup = (idx == 0) || 
        !(sortedPairs[idx] == sortedPairs[idx-1]);
    
    bool isLastInGroup = (idx == numPairs - 1) ||
        !(sortedPairs[idx] == sortedPairs[idx+1]);
    
    if (isFirstInGroup) {
        // Start of a new face
        KeyType faceId = 0;
        // Count how many faces before this one (atomic would be needed, simplified here)
        for (int i = 0; i < idx; ++i) {
            if (i == 0 || !(sortedPairs[i] == sortedPairs[i-1])) {
                faceId++;
            }
        }
        
        // Store face nodes
        faceNodes[faceId * 3 + 0] = sortedPairs[idx].node0;
        faceNodes[faceId * 3 + 1] = sortedPairs[idx].node1;
        faceNodes[faceId * 3 + 2] = sortedPairs[idx].node2;
        
        // Count elements for this face
        int numElems = 1;
        int j = idx + 1;
        while (j < numPairs && (sortedPairs[j] == sortedPairs[idx])) {
            numElems++;
            j++;
        }
        
        // Boundary if only 1 element
        isBoundary[faceId] = (numElems == 1) ? 1 : 0;
    }
}

// GPU kernel to compute face normals and areas
template<typename KeyType, typename RealType>
__global__ void computeFaceGeometryKernel(
    const KeyType* __restrict__ faceNodes,
    const KeyType* __restrict__ sfc_to_local,
    size_t nodeCount,
    const RealType* __restrict__ x,
    const RealType* __restrict__ y,
    const RealType* __restrict__ z,
    size_t numFaces,
    int nodesPerFace,
    RealType* normalX,
    RealType* normalY,
    RealType* normalZ,
    RealType* area,
    RealType* centroidX,
    RealType* centroidY,
    RealType* centroidZ)
{
    int faceIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (faceIdx >= numFaces) return;
    
    // Get face nodes (now local indices)
    KeyType local_nodes[4];
    for (int i = 0; i < nodesPerFace; ++i) {
        local_nodes[i] = faceNodes[faceIdx * nodesPerFace + i];
    }
    
    // Get SFC keys from local indices
    KeyType sfc_nodes[4];
    for (int i = 0; i < nodesPerFace; ++i) {
        sfc_nodes[i] = sfc_to_local[local_nodes[i]];
    }
    
    // Get coordinates
    RealType faceX[4], faceY[4], faceZ[4];
    RealType centroidX_val = 0, centroidY_val = 0, centroidZ_val = 0;
    for (int i = 0; i < nodesPerFace; ++i) {
        faceX[i] = x[local_nodes[i]];
        faceY[i] = y[local_nodes[i]];
        faceZ[i] = z[local_nodes[i]];
        centroidX_val += faceX[i] / nodesPerFace;
        centroidY_val += faceY[i] / nodesPerFace;
        centroidZ_val += faceZ[i] / nodesPerFace;
    }
    
    // Compute normal and area
    RealType nx = 0, ny = 0, nz = 0;
    RealType faceArea = 0;
    
    if (nodesPerFace == 3) {
        // Triangle
        RealType ax = faceX[1] - faceX[0];
        RealType ay = faceY[1] - faceY[0];
        RealType az = faceZ[1] - faceZ[0];
        RealType bx = faceX[2] - faceX[0];
        RealType by = faceY[2] - faceY[0];
        RealType bz = faceZ[2] - faceZ[0];
        
        nx = ay * bz - az * by;
        ny = az * bx - ax * bz;
        nz = ax * by - ay * bx;
        
        RealType norm = sqrt(nx*nx + ny*ny + nz*nz);
        faceArea = 0.5 * norm;
        if (norm > 1e-14) {
            nx /= norm;
            ny /= norm;
            nz /= norm;
        }
    } else if (nodesPerFace == 4) {
        // Quad - use cross product of two adjacent edges
        RealType ax = faceX[1] - faceX[0];
        RealType ay = faceY[1] - faceY[0];
        RealType az = faceZ[1] - faceZ[0];
        RealType bx = faceX[3] - faceX[0];  // Use diagonal
        RealType by = faceY[3] - faceY[0];
        RealType bz = faceZ[3] - faceZ[0];
        
        nx = ay * bz - az * by;
        ny = az * bx - ax * bz;
        nz = ax * by - ay * bx;
        
        RealType norm = sqrt(nx*nx + ny*ny + nz*nz);
        if (norm > 1e-14) {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            // Area as sum of two triangles
            RealType cx = faceX[2] - faceX[0];
            RealType cy = faceY[2] - faceY[0];
            RealType cz = faceZ[2] - faceZ[0];
            
            RealType area1 = 0.5 * sqrt(
                ((ay*cz - az*cy)*(ay*cz - az*cy) + 
                 (az*cx - ax*cz)*(az*cx - ax*cz) + 
                 (ax*cy - ay*cx)*(ax*cy - ay*cx)));
            
            RealType area2 = 0.5 * sqrt(
                ((by*cz - bz*cy)*(by*cz - bz*cy) + 
                 (bz*cx - bx*cz)*(bz*cx - bx*cz) + 
                 (bx*cy - by*cx)*(bx*cy - by*cx)));
            
            faceArea = area1 + area2;
        }
    }
    
    // Write results
    normalX[faceIdx] = nx;
    normalY[faceIdx] = ny;
    normalZ[faceIdx] = nz;
    area[faceIdx] = faceArea;
    centroidX[faceIdx] = centroidX_val;
    centroidY[faceIdx] = centroidY_val;
    centroidZ[faceIdx] = centroidZ_val;
}

// ============================================================================
// FaceTopology: Face connectivity and geometry for CVFEM (GPU production version)
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class FaceTopology
{
public:
    static constexpr int NodesPerFace = std::is_same_v<ElementTag, TetTag> ? 3 : 4;

    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    KeyType numFaces_ = 0;
    KeyType numBoundaryFaces_ = 0;
    KeyType numInteriorFaces_ = 0;

    DeviceVector<KeyType> d_faceNodes_;              // [numFaces * NodesPerFace] - sorted SFC keys
    DeviceVector<KeyType> d_faceToElementOffsets_;   // [numFaces + 1] CSR
    DeviceVector<KeyType> d_faceToElementList_;      // adjacency list
    DeviceVector<RealType> d_faceNormalX_, d_faceNormalY_, d_faceNormalZ_;
    DeviceVector<RealType> d_faceArea_;
    DeviceVector<RealType> d_faceCentroidX_, d_faceCentroidY_, d_faceCentroidZ_;
    DeviceVector<uint8_t> d_isBoundaryFace_;

    FaceTopology(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

    void extractAndDeduplicate(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    void computeFaceGeometry(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
};

// ============================================================================
// Implementation
// ============================================================================

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::FaceTopology(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    extractAndDeduplicate(domain);
    computeFaceGeometry(domain);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>::extractAndDeduplicate(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if constexpr (!std::is_same_v<AcceleratorTag, cstone::GpuTag>) {
        throw std::runtime_error("FaceTopology only implemented for GPU");
    }
    
    const auto& conn = domain.getConnectivity();
    size_t numElements = domain.getElementCount();
    
    if constexpr (std::is_same_v<ElementTag, TetTag>) {
        // Extract all tet face-element pairs (4 per element)
        thrust::device_vector<FaceElementPair<KeyType>> d_faceElemPairs(numElements * 4);
        
        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        
        extractTetFaceElementPairsKernel<<<numBlocks, blockSize>>>(
            std::get<0>(conn).data(),
            std::get<1>(conn).data(),
            std::get<2>(conn).data(),
            std::get<3>(conn).data(),
            numElements,
            thrust::raw_pointer_cast(d_faceElemPairs.data())
        );
        cudaDeviceSynchronize();
        
        // Sort face-element pairs by face nodes
        // This brings together all pairs belonging to the same face
        thrust::sort(d_faceElemPairs.begin(), d_faceElemPairs.end());
        
        // Use thrust::reduce_by_key to count unique faces and build connectivity
        // Create keys (face signatures) and values (element IDs)
        thrust::device_vector<KeyType> d_faceKeys(numElements * 4);
        thrust::device_vector<KeyType> d_elementIds(numElements * 4);
        
        // Build face signature as combined hash
        thrust::transform(d_faceElemPairs.begin(), d_faceElemPairs.end(),
                         d_faceKeys.begin(),
                         [] __device__ (const FaceElementPair<KeyType>& p) {
                             // Simple hash: use first node as primary key
                             return p.node0;
                         });
        
        thrust::transform(d_faceElemPairs.begin(), d_faceElemPairs.end(),
                         d_elementIds.begin(),
                         [] __device__ (const FaceElementPair<KeyType>& p) {
                             return p.elementId;
                         });
        
        // Count consecutive runs of same face
        // This gives us unique faces and their element counts
        thrust::device_vector<int> d_faceCounts(numElements * 4);
        thrust::device_vector<KeyType> d_uniqueFaceNodes(numElements * 4 * 3);
        
        // Manual implementation: scan through sorted pairs to build connectivity
        std::vector<FaceElementPair<KeyType>> h_pairs(numElements * 4);
        thrust::copy(d_faceElemPairs.begin(), d_faceElemPairs.end(), h_pairs.begin());
        
        std::vector<std::array<KeyType, 3>> uniqueFaceNodes;
        std::vector<std::vector<KeyType>> faceElements;
        
        for (size_t i = 0; i < h_pairs.size(); ) {
            // Start of new face
            std::array<KeyType, 3> faceNodes = {h_pairs[i].node0, h_pairs[i].node1, h_pairs[i].node2};
            std::vector<KeyType> elems;
            elems.push_back(h_pairs[i].elementId);
            
            // Collect all elements sharing this face
            size_t j = i + 1;
            while (j < h_pairs.size() && h_pairs[j] == h_pairs[i]) {
                elems.push_back(h_pairs[j].elementId);
                j++;
            }
            
            uniqueFaceNodes.push_back(faceNodes);
            faceElements.push_back(elems);
            i = j;
        }
        
        numFaces_ = uniqueFaceNodes.size();
        numBoundaryFaces_ = 0;
        numInteriorFaces_ = 0;
        
        // Count boundary and interior faces
        for (const auto& elems : faceElements) {
            if (elems.size() == 1) numBoundaryFaces_++;
            else numInteriorFaces_++;
        }
        
        // Build device arrays
        d_faceNodes_.resize(numFaces_ * 3);
        d_isBoundaryFace_.resize(numFaces_);
        d_faceToElementOffsets_.resize(numFaces_ + 1);
        
        size_t totalConnectivity = 0;
        for (const auto& elems : faceElements) {
            totalConnectivity += elems.size();
        }
        d_faceToElementList_.resize(totalConnectivity);
        
        // Fill device arrays from host data
        std::vector<KeyType> h_faceNodes(numFaces_ * 3);
        std::vector<uint8_t> h_isBoundary(numFaces_);
        std::vector<KeyType> h_offsets(numFaces_ + 1);
        std::vector<KeyType> h_elemList(totalConnectivity);
        
        size_t offset = 0;
        for (size_t faceId = 0; faceId < numFaces_; ++faceId) {
            h_faceNodes[faceId * 3 + 0] = uniqueFaceNodes[faceId][0];
            h_faceNodes[faceId * 3 + 1] = uniqueFaceNodes[faceId][1];
            h_faceNodes[faceId * 3 + 2] = uniqueFaceNodes[faceId][2];
            
            h_isBoundary[faceId] = (faceElements[faceId].size() == 1) ? 1 : 0;
            
            h_offsets[faceId] = offset;
            for (KeyType elemId : faceElements[faceId]) {
                h_elemList[offset++] = elemId;
            }
        }
        h_offsets[numFaces_] = offset;
        
        // Copy to device
        thrust::copy(h_faceNodes.begin(), h_faceNodes.end(), 
                     thrust::device_pointer_cast(d_faceNodes_.data()));
        thrust::copy(h_isBoundary.begin(), h_isBoundary.end(), 
                     thrust::device_pointer_cast(d_isBoundaryFace_.data()));
        thrust::copy(h_offsets.begin(), h_offsets.end(), 
                     thrust::device_pointer_cast(d_faceToElementOffsets_.data()));
        thrust::copy(h_elemList.begin(), h_elemList.end(), 
                     thrust::device_pointer_cast(d_faceToElementList_.data()));
    }
    else if constexpr (std::is_same_v<ElementTag, HexTag>) {
        // Extract all hex face-element pairs (6 per element)
        thrust::device_vector<FaceElementPair<KeyType>> d_faceElemPairs(numElements * 6);
        
        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        
        extractHexFaceElementPairsKernel<<<numBlocks, blockSize>>>(
            std::get<0>(conn).data(),
            std::get<1>(conn).data(),
            std::get<2>(conn).data(),
            std::get<3>(conn).data(),
            std::get<4>(conn).data(),
            std::get<5>(conn).data(),
            std::get<6>(conn).data(),
            std::get<7>(conn).data(),
            numElements,
            thrust::raw_pointer_cast(d_faceElemPairs.data())
        );
        cudaDeviceSynchronize();
        
        // Sort face-element pairs
        thrust::sort(d_faceElemPairs.begin(), d_faceElemPairs.end());
        
        // Copy to host and build connectivity
        std::vector<FaceElementPair<KeyType>> h_pairs(numElements * 6);
        thrust::copy(d_faceElemPairs.begin(), d_faceElemPairs.end(), h_pairs.begin());
        
        std::vector<std::array<KeyType, 4>> uniqueFaceNodes;  // Store all 4 nodes
        std::vector<std::vector<KeyType>> faceElements;
        
        for (size_t i = 0; i < h_pairs.size(); ) {
            // Start of new face (quad face with 4 nodes)
            std::array<KeyType, 4> faceNodes = {h_pairs[i].node0, h_pairs[i].node1, h_pairs[i].node2, h_pairs[i].node3};
            std::vector<KeyType> elems;
            elems.push_back(h_pairs[i].elementId);
            
            // Collect all elements sharing this face
            size_t j = i + 1;
            while (j < h_pairs.size() && h_pairs[j] == h_pairs[i]) {
                elems.push_back(h_pairs[j].elementId);
                j++;
            }
            
            uniqueFaceNodes.push_back(faceNodes);
            faceElements.push_back(elems);
            i = j;
        }
        
        numFaces_ = uniqueFaceNodes.size();
        numBoundaryFaces_ = 0;
        numInteriorFaces_ = 0;
        
        for (const auto& elems : faceElements) {
            if (elems.size() == 1) numBoundaryFaces_++;
            else numInteriorFaces_++;
        }
        
        // Build device arrays
        d_faceNodes_.resize(numFaces_ * 4);  // 4 nodes per quad face
        d_isBoundaryFace_.resize(numFaces_);
        d_faceToElementOffsets_.resize(numFaces_ + 1);
        
        size_t totalConnectivity = 0;
        for (const auto& elems : faceElements) {
            totalConnectivity += elems.size();
        }
        d_faceToElementList_.resize(totalConnectivity);
        
        // Fill device arrays
        std::vector<KeyType> h_faceNodes(numFaces_ * 4);
        std::vector<uint8_t> h_isBoundary(numFaces_);
        std::vector<KeyType> h_offsets(numFaces_ + 1);
        std::vector<KeyType> h_elemList(totalConnectivity);
        
        size_t offset = 0;
        for (size_t faceId = 0; faceId < numFaces_; ++faceId) {
            h_faceNodes[faceId * 4 + 0] = uniqueFaceNodes[faceId][0];
            h_faceNodes[faceId * 4 + 1] = uniqueFaceNodes[faceId][1];
            h_faceNodes[faceId * 4 + 2] = uniqueFaceNodes[faceId][2];
            h_faceNodes[faceId * 4 + 3] = uniqueFaceNodes[faceId][3];
            
            h_isBoundary[faceId] = (faceElements[faceId].size() == 1) ? 1 : 0;
            
            h_offsets[faceId] = offset;
            for (KeyType elemId : faceElements[faceId]) {
                h_elemList[offset++] = elemId;
            }
        }
        h_offsets[numFaces_] = offset;
        
        // Copy to device
        thrust::copy(h_faceNodes.begin(), h_faceNodes.end(), 
                     thrust::device_pointer_cast(d_faceNodes_.data()));
        
        // Map SFC keys to local indices
        const auto& d_sfc_map = domain.getLocalToGlobalSfcMap();
        int blockSizeMap = 256;
        int numBlocksMap = (d_faceNodes_.size() + blockSizeMap - 1) / blockSizeMap;
        mapSfcToLocalIdKernel<<<numBlocksMap, blockSizeMap>>>(
            d_faceNodes_.data(),
            d_faceNodes_.data(),  // in place
            d_sfc_map.data(),
            d_faceNodes_.size(),
            domain.getNodeCount()
        );
        
        thrust::copy(h_isBoundary.begin(), h_isBoundary.end(), 
                     thrust::device_pointer_cast(d_isBoundaryFace_.data()));
        thrust::copy(h_offsets.begin(), h_offsets.end(), 
                     thrust::device_pointer_cast(d_faceToElementOffsets_.data()));
        thrust::copy(h_elemList.begin(), h_elemList.end(), 
                     thrust::device_pointer_cast(d_faceToElementList_.data()));
    }
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
    
    // Get coordinates
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();
    const auto& d_sfc_map = domain.getLocalToGlobalSfcMap();
    
    int blockSize = 256;
    int numBlocks = (numFaces_ + blockSize - 1) / blockSize;
    
    computeFaceGeometryKernel<<<numBlocks, blockSize>>>(
        d_faceNodes_.data(),
        d_sfc_map.data(),
        domain.getNodeCount(),
        d_x.data(),
        d_y.data(),
        d_z.data(),
        numFaces_,
        NodesPerFace,
        d_faceNormalX_.data(),
        d_faceNormalY_.data(),
        d_faceNormalZ_.data(),
        d_faceArea_.data(),
        d_faceCentroidX_.data(),
        d_faceCentroidY_.data(),
        d_faceCentroidZ_.data()
    );
    cudaDeviceSynchronize();
}

} // namespace mars
