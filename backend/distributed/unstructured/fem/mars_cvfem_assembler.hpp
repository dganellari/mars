#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_graph.hpp"
#include "mars_cvfem_hex_kernel_optimized.hpp"
#include "mars_cvfem_hex_kernel_shmem.hpp"
#include "mars_cvfem_hex_kernel_tensor.hpp"
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// Kernel variant for CVFEM assembly
enum class CvfemKernelVariant {
    Original,    // Original graph kernel with linear search
    Optimized,   // Binary search, pragma unroll optimizations
    Shmem,       // Low-register kernel with on-the-fly shape derivatives and binary search
    Team,        // Warp-cooperative (32 threads/element) with shared memory + parallel SCS integration
    Tensor       // Tensor core variant using FP64 WMMA for 8×8 matrix assembly (GH200)
};

// High-level CVFEM assembler for hex elements
// Encapsulates kernel selection and launch configuration
template<typename KeyType, typename RealType>
class CvfemHexAssembler {
public:
    struct Config {
        int blockSize = 256;
        CvfemKernelVariant variant = CvfemKernelVariant::Original;
        cudaStream_t stream = 0;
    };

    // Assemble CVFEM system with graph-based sparsity and diagonal lumping
    static void assembleGraphLump(
        // Element connectivity (8 arrays for hex8)
        const KeyType* d_conn0,
        const KeyType* d_conn1,
        const KeyType* d_conn2,
        const KeyType* d_conn3,
        const KeyType* d_conn4,
        const KeyType* d_conn5,
        const KeyType* d_conn6,
        const KeyType* d_conn7,
        size_t numElements,
        // Node coordinates
        const RealType* d_x,
        const RealType* d_y,
        const RealType* d_z,
        // Field data (per node)
        const RealType* d_gamma,
        const RealType* d_phi,
        const RealType* d_beta,
        const RealType* d_grad_phi_x,
        const RealType* d_grad_phi_y,
        const RealType* d_grad_phi_z,
        // Element data (per element * 12 SCS)
        const RealType* d_mdot,
        const RealType* d_areaVec_x,
        const RealType* d_areaVec_y,
        const RealType* d_areaVec_z,
        // DOF mapping
        const int* d_node_to_dof,
        const uint8_t* d_ownership,
        // Matrix assembly (CSR format)
        CSRMatrix<RealType>* d_matrix,
        RealType* d_rhs,
        // Configuration
        const Config& config = Config{})
    {
        int blockSize = config.blockSize;
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        switch (config.variant) {
            case CvfemKernelVariant::Original:
                cvfem_hex_assembly_kernel_graph<KeyType, RealType><<<numBlocks, blockSize, 0, config.stream>>>(
                    d_conn0, d_conn1, d_conn2, d_conn3,
                    d_conn4, d_conn5, d_conn6, d_conn7,
                    numElements,
                    d_x, d_y, d_z,
                    d_gamma, d_phi, d_beta,
                    d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                    d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                    d_node_to_dof, d_ownership,
                    d_matrix, d_rhs
                );
                break;

            case CvfemKernelVariant::Optimized:
                cvfem_hex_assembly_kernel_optimized<KeyType, RealType><<<numBlocks, blockSize, 0, config.stream>>>(
                    d_conn0, d_conn1, d_conn2, d_conn3,
                    d_conn4, d_conn5, d_conn6, d_conn7,
                    numElements,
                    d_x, d_y, d_z,
                    d_gamma, d_phi, d_beta,
                    d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                    d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                    d_node_to_dof, d_ownership,
                    d_matrix, d_rhs
                );
                break;

            case CvfemKernelVariant::Shmem:
                // Thread-per-element with on-the-fly shape derivative computation
                cvfem_hex_assembly_kernel_shmem<KeyType, RealType, 256><<<numBlocks, blockSize, 0, config.stream>>>(
                    d_conn0, d_conn1, d_conn2, d_conn3,
                    d_conn4, d_conn5, d_conn6, d_conn7,
                    numElements,
                    d_x, d_y, d_z,
                    d_gamma, d_phi, d_beta,
                    d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                    d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                    d_node_to_dof, d_ownership,
                    d_matrix, d_rhs
                );
                break;

            case CvfemKernelVariant::Team: {
                // Team kernel uses 32 threads per element (one warp per element)
                int teamBlocks = (numElements * 32 + blockSize - 1) / blockSize;
                cvfem_hex_assembly_kernel_team<KeyType, RealType, 256><<<teamBlocks, blockSize, 0, config.stream>>>(
                    d_conn0, d_conn1, d_conn2, d_conn3,
                    d_conn4, d_conn5, d_conn6, d_conn7,
                    numElements,
                    d_x, d_y, d_z,
                    d_gamma, d_phi, d_beta,
                    d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                    d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                    d_node_to_dof, d_ownership,
                    d_matrix, d_rhs
                );
                break;
            }

            case CvfemKernelVariant::Tensor:
                // Tensor core variant - assembles full 8×8 local matrix
                cvfem_hex_assembly_kernel_tensor<KeyType, RealType, 256><<<numBlocks, blockSize, 0, config.stream>>>(
                    d_conn0, d_conn1, d_conn2, d_conn3,
                    d_conn4, d_conn5, d_conn6, d_conn7,
                    numElements,
                    d_x, d_y, d_z,
                    d_gamma, d_phi, d_beta,
                    d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                    d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                    d_node_to_dof, d_ownership,
                    d_matrix, d_rhs
                );
                break;
        }
    }

    // Assemble CVFEM system with full sparsity (27 NNZ/row)
    static void assembleFull(
        const KeyType* d_conn0,
        const KeyType* d_conn1,
        const KeyType* d_conn2,
        const KeyType* d_conn3,
        const KeyType* d_conn4,
        const KeyType* d_conn5,
        const KeyType* d_conn6,
        const KeyType* d_conn7,
        size_t numElements,
        const RealType* d_x,
        const RealType* d_y,
        const RealType* d_z,
        const RealType* d_gamma,
        const RealType* d_phi,
        const RealType* d_beta,
        const RealType* d_grad_phi_x,
        const RealType* d_grad_phi_y,
        const RealType* d_grad_phi_z,
        const RealType* d_mdot,
        const RealType* d_areaVec_x,
        const RealType* d_areaVec_y,
        const RealType* d_areaVec_z,
        const int* d_node_to_dof,
        const uint8_t* d_ownership,
        CSRMatrix<RealType>* d_matrix,
        RealType* d_rhs,
        const Config& config = Config{})
    {
        int blockSize = config.blockSize;
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        // Full assembly uses original kernel (no lumping needed)
        cvfem_hex_assembly_kernel<KeyType, RealType><<<numBlocks, blockSize, 0, config.stream>>>(
            d_conn0, d_conn1, d_conn2, d_conn3,
            d_conn4, d_conn5, d_conn6, d_conn7,
            numElements,
            d_x, d_y, d_z,
            d_gamma, d_phi, d_beta,
            d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
            d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
            d_node_to_dof, d_ownership,
            d_matrix, d_rhs
        );
    }

    // Get kernel variant name for logging
    static const char* variantName(CvfemKernelVariant variant) {
        switch (variant) {
            case CvfemKernelVariant::Original: return "original";
            case CvfemKernelVariant::Optimized: return "optimized";
            case CvfemKernelVariant::Shmem: return "shmem";
            case CvfemKernelVariant::Team: return "team";
            case CvfemKernelVariant::Tensor: return "tensor";
            default: return "unknown";
        }
    }
};

} // namespace fem
} // namespace mars
