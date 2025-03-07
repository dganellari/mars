#include "domain.hpp"

// CUDA kernel implementations
__global__ void computeCharacteristicSizesKernel(double* nodes,
                                                 int* tets,
                                                 int* nodeTetCount,
                                                 double* h,
                                                 int localTetCount) {
    int tetIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (tetIdx < localTetCount) {
        // For each tetrahedron, process all 6 edges
        for (int i = 0; i < 4; ++i) {
            int nodeIdx1 = tets[tetIdx * 4 + i];

            for (int j = i + 1; j < 4; ++j) {
                int nodeIdx2 = tets[tetIdx * 4 + j];

                // Calculate edge length
                double dx = nodes[nodeIdx1 * 3] - nodes[nodeIdx2 * 3];
                double dy = nodes[nodeIdx1 * 3 + 1] - nodes[nodeIdx2 * 3 + 1];
                double dz = nodes[nodeIdx1 * 3 + 2] - nodes[nodeIdx2 * 3 + 2];
                double edgeLength = sqrt(dx * dx + dy * dy + dz * dz);

                // Add to each node's total (using atomic to avoid race conditions)
                atomicAdd(&h[nodeIdx1], edgeLength);
                atomicAdd(&h[nodeIdx2], edgeLength);
                atomicAdd(&nodeTetCount[nodeIdx1], 1);
                atomicAdd(&nodeTetCount[nodeIdx2], 1);
            }
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

__global__ void findRepresentativeNodesKernel(int* tets, unsigned* sfcCodes, int* tetToNodeMap, int numTets) {
    int tetIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (tetIdx < numTets) {
        // Get the four nodes of this tetrahedron
        int node0 = tets[tetIdx * 4];
        int node1 = tets[tetIdx * 4 + 1];
        int node2 = tets[tetIdx * 4 + 2];
        int node3 = tets[tetIdx * 4 + 3];

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
        tetToNodeMap[tetIdx] = repNode;
    }
}

__global__ void computeTetrahedralVolumesKernel(double* x, double* y, double* z, int* connectivity, int numTets) {
    int tetIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (tetIdx < numTets) {
        // Get the four nodes of this tetrahedron
        int n0 = connectivity[tetIdx * 4];
        int n1 = connectivity[tetIdx * 4 + 1];
        int n2 = connectivity[tetIdx * 4 + 2];
        int n3 = connectivity[tetIdx * 4 + 3];

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
        double volume =
            fabs(v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) + v1z * (v2x * v3y - v2y * v3x)) / 6.0;

        // Volume calculation result could be stored or used here
    }
}