#include "domain.hpp"

namespace mars
{
// CUDA kernels with RealType template parameter instead of Real
template<typename RealType>
__global__ void transformCharacteristicSizesKernel(RealType* d_h, size_t size, RealType meshFactor, RealType minH, RealType maxH)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size)
    {
        RealType val    = d_h[idx];
        RealType result = val * meshFactor;
        d_h[idx]    = fmaxf(minH, fminf(maxH, result));
    }
}

template<typename RealType>
__global__ void fillCharacteristicSizesKernel(RealType* d_h, size_t size, RealType value)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) { d_h[idx] = value; }
}

// CUDA kernel to calculate element characteristic sizes
template<typename ElementTag, typename RealType>
__global__ void computeCharacteristicSizesKernel(const RealType* x,
                                                 const RealType* y,
                                                 const RealType* z,
                                                 const int* indices0,
                                                 const int* indices1,
                                                 const int* indices2,
                                                 const int* indices3,
                                                 int* nodeTetCount,
                                                 RealType* h,
                                                 int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get the four nodes of this tetrahedron
            int n0 = indices0[elemIdx];
            int n1 = indices1[elemIdx];
            int n2 = indices2[elemIdx];
            int n3 = indices3[elemIdx];

            // Calculate edge lengths and contribute to characteristic size
            // Edge n0-n1
            RealType dx         = x[n0] - x[n1];
            RealType dy         = y[n0] - y[n1];
            RealType dz         = z[n0] - z[n1];
            RealType edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n1], 1);

            // Edge n0-n2
            dx         = x[n0] - x[n2];
            dy         = y[n0] - y[n2];
            dz         = z[n0] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n2], 1);

            // Edge n0-n3
            dx         = x[n0] - x[n3];
            dy         = y[n0] - y[n3];
            dz         = z[n0] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n3], 1);

            // Edge n1-n2
            dx         = x[n1] - x[n2];
            dy         = y[n1] - y[n2];
            dz         = z[n1] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1);
            atomicAdd(&nodeTetCount[n2], 1);

            // Edge n1-n3
            dx         = x[n1] - x[n3];
            dy         = y[n1] - y[n3];
            dz         = z[n1] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1);
            atomicAdd(&nodeTetCount[n3], 1);

            // Edge n2-n3
            dx         = x[n2] - x[n3];
            dy         = y[n2] - y[n3];
            dz         = z[n2] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n2], 1);
            atomicAdd(&nodeTetCount[n3], 1);
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // Hexahedral element implementation
            // This would need all 8 node indices and 12 edges
            // ...
        }
    }
}

template<typename RealType>
__global__ void finalizeCharacteristicSizesKernel(RealType* h, int* nodeTetCount, int numNodes)
{
    int nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (nodeIdx < numNodes)
    {
        if (nodeTetCount[nodeIdx] > 0) { h[nodeIdx] /= nodeTetCount[nodeIdx]; }
        else
        {
            h[nodeIdx] = 0.01; // Default for isolated nodes
        }
    }
}

// Generic kernel for finding representative nodes
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void findRepresentativeNodesKernel(const int* indices0,
                                              const int* indices1,
                                              const int* indices2,
                                              const int* indices3,
                                              const KeyType* sfcCodes,
                                              int* elemToNodeMap,
                                              int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get the four nodes of this tetrahedron
            int node0 = indices0[elemIdx];
            int node1 = indices1[elemIdx];
            int node2 = indices2[elemIdx];
            int node3 = indices3[elemIdx];

            // Get SFC codes
            unsigned sfc0 = sfcCodes[node0];
            unsigned sfc1 = sfcCodes[node1];
            unsigned sfc2 = sfcCodes[node2];
            unsigned sfc3 = sfcCodes[node3];

            // Find minimum SFC code
            int repNode     = node0;
            unsigned minSfc = sfc0;

            if (sfc1 < minSfc)
            {
                minSfc  = sfc1;
                repNode = node1;
            }

            if (sfc2 < minSfc)
            {
                minSfc  = sfc2;
                repNode = node2;
            }

            if (sfc3 < minSfc)
            {
                minSfc  = sfc3;
                repNode = node3;
            }

            // Store the representative node
            elemToNodeMap[elemIdx] = repNode;
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // For hexahedra - would need additional parameters for indices4-7
            // Implementation would be similar but with 8 nodes
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            // For triangles - would use only indices0-2
            int node0 = indices0[elemIdx];
            int node1 = indices1[elemIdx];
            int node2 = indices2[elemIdx];

            unsigned sfc0 = sfcCodes[node0];
            unsigned sfc1 = sfcCodes[node1];
            unsigned sfc2 = sfcCodes[node2];

            int repNode     = node0;
            unsigned minSfc = sfc0;

            if (sfc1 < minSfc)
            {
                minSfc  = sfc1;
                repNode = node1;
            }

            if (sfc2 < minSfc)
            {
                minSfc  = sfc2;
                repNode = node2;
            }

            elemToNodeMap[elemIdx] = repNode;
        }
    }
}

// Implementation of mapElementsToNodes for GPU
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::mapElementsToNodes()
{
    // This is a stub implementation - the actual mapping happens in the sync() method
    // Just ensure the vectors are sized correctly
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        if (d_elemToNodeMap_.size() != elementCount_)
        {
            d_elemToNodeMap_.resize(elementCount_);
        }
    }
}

// Explicit instantiations for mapElementsToNodes
template void ElementDomain<TetTag, float, unsigned int, cstone::GpuTag>::mapElementsToNodes();
template void ElementDomain<TetTag, double, unsigned int, cstone::GpuTag>::mapElementsToNodes();
template void ElementDomain<TetTag, float, uint64_t, cstone::GpuTag>::mapElementsToNodes();
template void ElementDomain<TetTag, double, uint64_t, cstone::GpuTag>::mapElementsToNodes();

template<typename KeyType, typename RealType>
void generateSfcKeys(
    const RealType* x, const RealType* y, const RealType* z, KeyType* keys, size_t numKeys, const cstone::Box<RealType>& box)
{
    // Use sfcKindPointer to match cornerstone's template instantiation
    cstone::computeSfcKeysGpu(x, y, z, cstone::sfcKindPointer(keys), numKeys, box);
    cudaCheckError();
}

// Kernel to extract representative node coordinates and compute SFC keys
template<typename ElementTag, typename RealType>
__global__ void extractRepCoordinatesKernel(const RealType* x,
                                            const RealType* y,
                                            const RealType* z,
                                            const RealType* h,
                                            const int* elemToNodeMap,
                                            RealType* elemX,
                                            RealType* elemY,
                                            RealType* elemZ,
                                            RealType* elemH,
                                            int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        int repNodeIdx = elemToNodeMap[elemIdx];

        // Extract coordinates and properties
        elemX[elemIdx] = x[repNodeIdx];
        elemY[elemIdx] = y[repNodeIdx];
        elemZ[elemIdx] = z[repNodeIdx];
        elemH[elemIdx] = h[repNodeIdx];
    }
}

// Compute tetrahedron volumes
template<typename ElementTag, typename RealType>
__global__ void computeElementVolumesKernel(const RealType* x,
                                            const RealType* y,
                                            const RealType* z,
                                            const int* indices0,
                                            const int* indices1,
                                            const int* indices2,
                                            const int* indices3,
                                            RealType* volumes,
                                            int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get node indices
            int n0 = indices0[elemIdx];
            int n1 = indices1[elemIdx];
            int n2 = indices2[elemIdx];
            int n3 = indices3[elemIdx];

            // Get node coordinates
            RealType x0 = x[n0], y0 = y[n0], z0 = z[n0];
            RealType x1 = x[n1], y1 = y[n1], z1 = z[n1];
            RealType x2 = x[n2], y2 = y[n2], z2 = z[n2];
            RealType x3 = x[n3], y3 = y[n3], z3 = z[n3];

            // Compute vectors for volume calculation
            RealType v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
            RealType v2x = x2 - x0, v2y = y2 - y0, v2z = z2 - z0;
            RealType v3x = x3 - x0, v3y = y3 - y0, v3z = z3 - z0;

            // Compute volume using the scalar triple product
            volumes[elemIdx] =
                fabs(v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) + v1z * (v2x * v3y - v2y * v3x)) /
                6.0;
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // Hexahedron volume calculation would go here
            // Would need additional indices parameters
        }
    }
}

// Explicit instantiations for float
template __global__ void transformCharacteristicSizesKernel<float>(float* d_h, size_t size, float meshFactor, float minH, float maxH);
template __global__ void fillCharacteristicSizesKernel<float>(float* d_h, size_t size, float value);
template __global__ void finalizeCharacteristicSizesKernel<float>(float* h, int* nodeTetCount, int numNodes);

// Explicit instantiations for double
template __global__ void transformCharacteristicSizesKernel<double>(double* d_h, size_t size, double meshFactor, double minH, double maxH);
template __global__ void fillCharacteristicSizesKernel<double>(double* d_h, size_t size, double value);
template __global__ void finalizeCharacteristicSizesKernel<double>(double* h, int* nodeTetCount, int numNodes);

// Explicit instantiation for each element type with float
template __global__ void computeCharacteristicSizesKernel<TetTag, float>(const float* x,
                                                                  const float* y,
                                                                  const float* z,
                                                                  const int* indices0,
                                                                  const int* indices1,
                                                                  const int* indices2,
                                                                  const int* indices3,
                                                                  int* nodeTetCount,
                                                                  float* h,
                                                                  int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, float>(const float* x,
                                                             const float* y,
                                                             const float* z,
                                                             const float* h,
                                                             const int* elemToNodeMap,
                                                             float* elemX,
                                                             float* elemY,
                                                             float* elemZ,
                                                             float* elemH,
                                                             int numElements);

template __global__ void computeElementVolumesKernel<TetTag, float>(const float* x,
                                                             const float* y,
                                                             const float* z,
                                                             const int* indices0,
                                                             const int* indices1,
                                                             const int* indices2,
                                                             const int* indices3,
                                                             float* volumes,
                                                             int numElements);

// Explicit instantiation for each element type with double
template __global__ void computeCharacteristicSizesKernel<TetTag, double>(const double* x,
                                                                  const double* y,
                                                                  const double* z,
                                                                  const int* indices0,
                                                                  const int* indices1,
                                                                  const int* indices2,
                                                                  const int* indices3,
                                                                  int* nodeTetCount,
                                                                  double* h,
                                                                  int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, double>(const double* x,
                                                             const double* y,
                                                             const double* z,
                                                             const double* h,
                                                             const int* elemToNodeMap,
                                                             double* elemX,
                                                             double* elemY,
                                                             double* elemZ,
                                                             double* elemH,
                                                             int numElements);

template __global__ void computeElementVolumesKernel<TetTag, double>(const double* x,
                                                             const double* y,
                                                             const double* z,
                                                             const int* indices0,
                                                             const int* indices1,
                                                             const int* indices2,
                                                             const int* indices3,
                                                             double* volumes,
                                                             int numElements);

// For float with unsigned keys
template __global__ void findRepresentativeNodesKernel<TetTag, unsigned, float>(
    const int* indices0, const int* indices1, const int* indices2, const int* indices3,
    const unsigned* sfcCodes, int* elemToNodeMap, int numElements);

// For double with unsigned keys
template __global__ void findRepresentativeNodesKernel<TetTag, unsigned, double>(
    const int* indices0, const int* indices1, const int* indices2, const int* indices3,
    const unsigned* sfcCodes, int* elemToNodeMap, int numElements);

// For float with uint64_t keys
template __global__ void findRepresentativeNodesKernel<TetTag, uint64_t, float>(
    const int* indices0, const int* indices1, const int* indices2, const int* indices3,
    const uint64_t* sfcCodes, int* elemToNodeMap, int numElements);

// For double with uint64_t keys
template __global__ void findRepresentativeNodesKernel<TetTag, uint64_t, double>(
    const int* indices0, const int* indices1, const int* indices2, const int* indices3,
    const uint64_t* sfcCodes, int* elemToNodeMap, int numElements);                                                          

// Explicit instantiation for computeSfcKeysGpu with common combinations
template void generateSfcKeys<unsigned, float>(
    const float* x, const float* y, const float* z, unsigned* keys, size_t numKeys, const cstone::Box<float>& box);
template void generateSfcKeys<unsigned, double>(
    const double* x, const double* y, const double* z, unsigned* keys, size_t numKeys, const cstone::Box<double>& box);
template void generateSfcKeys<uint64_t, float>(
    const float* x, const float* y, const float* z, uint64_t* keys, size_t numKeys, const cstone::Box<float>& box);
template void generateSfcKeys<uint64_t, double>(
    const double* x, const double* y, const double* z, uint64_t* keys, size_t numKeys, const cstone::Box<double>& box);
} // namespace mars