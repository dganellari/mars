#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/sfc/sfc.hpp"
#include "domain.hpp"

#define cudaCheckError()                                                                 \
    {                                                                                    \
        cudaError_t e = cudaGetLastError();                                              \
        if (e != cudaSuccess) {                                                          \
            printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
            exit(EXIT_FAILURE);                                                          \
        }                                                                                \
    }

// CUDA kernel to calculate element characteristic sizes
template <typename ElementTag>
__global__ void computeCharacteristicSizesKernel(const double* x,
                                                 const double* y,
                                                 const double* z,
                                                 const int* indices0,
                                                 const int* indices1,
                                                 const int* indices2,
                                                 const int* indices3,
                                                 int* nodeTetCount,
                                                 double* h,
                                                 int numElements) {
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements) {
        if constexpr (std::is_same_v<ElementTag, TetTag>) {
            // Get the four nodes of this tetrahedron
            int n0 = indices0[elemIdx];
            int n1 = indices1[elemIdx];
            int n2 = indices2[elemIdx];
            int n3 = indices3[elemIdx];

            // Calculate edge lengths and contribute to characteristic size
            // Edge n0-n1
            double dx = x[n0] - x[n1];
            double dy = y[n0] - y[n1];
            double dz = z[n0] - z[n1];
            double edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n1], 1);

            // Edge n0-n2
            dx = x[n0] - x[n2];
            dy = y[n0] - y[n2];
            dz = z[n0] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n2], 1);

            // Edge n0-n3
            dx = x[n0] - x[n3];
            dy = y[n0] - y[n3];
            dz = z[n0] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n3], 1);

            // Edge n1-n2
            dx = x[n1] - x[n2];
            dy = y[n1] - y[n2];
            dz = z[n1] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1);
            atomicAdd(&nodeTetCount[n2], 1);

            // Edge n1-n3
            dx = x[n1] - x[n3];
            dy = y[n1] - y[n3];
            dz = z[n1] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1);
            atomicAdd(&nodeTetCount[n3], 1);

            // Edge n2-n3
            dx = x[n2] - x[n3];
            dy = y[n2] - y[n3];
            dz = z[n2] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n2], 1);
            atomicAdd(&nodeTetCount[n3], 1);
        } else if constexpr (std::is_same_v<ElementTag, HexTag>) {
            // Hexahedral element implementation
            // This would need all 8 node indices and 12 edges
            // ...
        }
    }
}

__global__ void finalizeCharacteristicSizesKernel(double* h, int* nodeTetCount, int numNodes) {
    int nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (nodeIdx < numNodes) {
        if (nodeTetCount[nodeIdx] > 0) {
            h[nodeIdx] /= nodeTetCount[nodeIdx];
        } else {
            h[nodeIdx] = 0.01;  // Default for isolated nodes
        }
    }
}

// Generic kernel for finding representative nodes
template <typename ElementTag>
__global__ void findRepresentativeNodesKernel(const int* indices0,
                                              const int* indices1,
                                              const int* indices2,
                                              const int* indices3,
                                              const unsigned* sfcCodes,
                                              int* elemToNodeMap,
                                              int numElements) {
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements) {
        if constexpr (std::is_same_v<ElementTag, TetTag>) {
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
            int repNode = node0;
            unsigned minSfc = sfc0;

            if (sfc1 < minSfc) {
                minSfc = sfc1;
                repNode = node1;
            }

            if (sfc2 < minSfc) {
                minSfc = sfc2;
                repNode = node2;
            }

            if (sfc3 < minSfc) {
                minSfc = sfc3;
                repNode = node3;
            }

            // Store the representative node
            elemToNodeMap[elemIdx] = repNode;
        } else if constexpr (std::is_same_v<ElementTag, HexTag>) {
            // For hexahedra - would need additional parameters for indices4-7
            // Implementation would be similar but with 8 nodes
        } else if constexpr (std::is_same_v<ElementTag, TriTag>) {
            // For triangles - would use only indices0-2
            int node0 = indices0[elemIdx];
            int node1 = indices1[elemIdx];
            int node2 = indices2[elemIdx];

            unsigned sfc0 = sfcCodes[node0];
            unsigned sfc1 = sfcCodes[node1];
            unsigned sfc2 = sfcCodes[node2];

            int repNode = node0;
            unsigned minSfc = sfc0;

            if (sfc1 < minSfc) {
                minSfc = sfc1;
                repNode = node1;
            }

            if (sfc2 < minSfc) {
                minSfc = sfc2;
                repNode = node2;
            }

            elemToNodeMap[elemIdx] = repNode;
        }
    }
}

// Use cornerstone's implementation for SFC key computation
template <typename KeyType>
void computeSfcKeysGpu(const double* x,
                       const double* y,
                       const double* z,
                       KeyType* keys,
                       size_t numKeys,
                       const cstone::Box<double>& box) {
    // Call the cornerstone library's implementation directly
    cstone::computeSfcKeysGpu(x, y, z, keys, numKeys, box);

    cudaCheckError();
}

// Kernel to extract representative node coordinates and compute SFC keys
template <typename ElementTag>
__global__ void extractRepCoordinatesKernel(const double* x,
                                            const double* y,
                                            const double* z,
                                            const double* h,
                                            const int* elemToNodeMap,
                                            double* elemX,
                                            double* elemY,
                                            double* elemZ,
                                            double* elemH,
                                            int numElements) {
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements) {
        int repNodeIdx = elemToNodeMap[elemIdx];

        // Extract coordinates and properties
        elemX[elemIdx] = x[repNodeIdx];
        elemY[elemIdx] = y[repNodeIdx];
        elemZ[elemIdx] = z[repNodeIdx];
        elemH[elemIdx] = h[repNodeIdx];
    }
}

// Compute tetrahedron volumes
template <typename ElementTag>
__global__ void computeElementVolumesKernel(const double* x,
                                            const double* y,
                                            const double* z,
                                            const int* indices0,
                                            const int* indices1,
                                            const int* indices2,
                                            const int* indices3,
                                            double* volumes,
                                            int numElements) {
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements) {
        if constexpr (std::is_same_v<ElementTag, TetTag>) {
            // Get node indices
            int n0 = indices0[elemIdx];
            int n1 = indices1[elemIdx];
            int n2 = indices2[elemIdx];
            int n3 = indices3[elemIdx];

            // Get node coordinates
            double x0 = x[n0], y0 = y[n0], z0 = z[n0];
            double x1 = x[n1], y1 = y[n1], z1 = z[n1];
            double x2 = x[n2], y2 = y[n2], z2 = z[n2];
            double x3 = x[n3], y3 = y[n3], z3 = z[n3];

            // Compute vectors for volume calculation
            double v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
            double v2x = x2 - x0, v2y = y2 - y0, v2z = z2 - z0;
            double v3x = x3 - x0, v3y = y3 - y0, v3z = z3 - z0;

            // Compute volume using the scalar triple product
            volumes[elemIdx] =
                fabs(v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) + v1z * (v2x * v3y - v2y * v3x)) /
                6.0;
        } else if constexpr (std::is_same_v<ElementTag, HexTag>) {
            // Hexahedron volume calculation would go here
            // Would need additional indices parameters
        }
    }
}

// Explicit instantiation for each element type
template __global__ void computeCharacteristicSizesKernel<TetTag>(const double* x,
                                                                  const double* y,
                                                                  const double* z,
                                                                  const int* indices0,
                                                                  const int* indices1,
                                                                  const int* indices2,
                                                                  const int* indices3,
                                                                  int* nodeTetCount,
                                                                  double* h,
                                                                  int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag>(const int* indices0,
                                                               const int* indices1,
                                                               const int* indices2,
                                                               const int* indices3,
                                                               const unsigned* sfcCodes,
                                                               int* elemToNodeMap,
                                                               int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag>(const double* x,
                                                             const double* y,
                                                             const double* z,
                                                             const double* h,
                                                             const int* elemToNodeMap,
                                                             double* elemX,
                                                             double* elemY,
                                                             double* elemZ,
                                                             double* elemH,
                                                             int numElements);

template __global__ void computeElementVolumesKernel<TetTag>(const double* x,
                                                             const double* y,
                                                             const double* z,
                                                             const int* indices0,
                                                             const int* indices1,
                                                             const int* indices2,
                                                             const int* indices3,
                                                             double* volumes,
                                                             int numElements);

// Instantiation for hex elements would follow
// template __global__ void computeCharacteristicSizesKernel<HexTag>...
// etc.

// Explicit instantiation for the type used in our code
template void computeSfcKeysGpu<KeyType>(const double* x,
                                         const double* y,
                                         const double* z,
                                         KeyType* keys,
                                         size_t numKeys,
                                         const cstone::Box<double>& box);
