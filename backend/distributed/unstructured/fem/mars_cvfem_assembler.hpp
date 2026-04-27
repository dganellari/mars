#pragma once

#include "mars_cvfem_hex_kernel.hpp"
#include "mars_cvfem_hex_kernel_graph.hpp"
#include "mars_cvfem_hex_kernel_optimized.hpp"
#include "mars_cvfem_hex_kernel_shmem.hpp"
#include "mars_cvfem_hex_kernel_tensor.hpp"
#include "mars_cvfem_hex_kernel_tensor_wmma.hpp"
#include "mars_cvfem_hex_kernel_tensor_colored.hpp"
#include "mars_cvfem_hex_kernel_tensor_aos.hpp"
#include "mars_cvfem_hex_kernel_tensor_perip.hpp"
#include "mars_cvfem_hex_kernel_smem_cache.hpp"
#include "mars_cvfem_node_data.hpp"
#include "mars_cvfem_coloring.hpp"
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// Kernel variant for CVFEM assembly
enum class CvfemKernelVariant {
    Original,    // Original graph kernel with linear search
    Optimized,   // Binary search, pragma unroll optimizations
    Shmem,       // Low-register kernel with on-the-fly shape derivatives and binary search
    Team,        // Warp-cooperative (32 threads/element) with shared memory + parallel SCS integration
    Tensor,        // Full 8x8 local matrix, scalar (current best: 4.60ms on GH200)
    TensorColored, // Tensor kernel + graph coloring → atomic-free CSR scatter
    TensorAoS,     // Tensor kernel with AoS node data (3× less L2 traffic for scattered node reads)
    TensorPerip,   // Pre-lookup pos[64], scatter per-SCS: -64 regs vs tensor (lhs[64]→pos[64])
    TensorPeripLb2,// Same + __launch_bounds__(256,2): forces 2 blocks/SM, ~7 regs forced spill
    SmemCache,     // Smem node cache: deduplicate+load ~400 unique nodes/block, 5cyc vs 30cyc L2
    WmmaTensor,    // FP64 WMMA tensor cores: 1 warp/element, diffusion = S^T x D matmul (SM80+)
    WgmmaTensor    // Hopper WGMMA m64n8k4: 2 elem/warp (75% util) + async WGMMA for diffusion (SM90+ only)
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
        const CvfemColoringData* coloring  = nullptr; // required for TensorColored
        const NodeData*          nodeData  = nullptr; // required for TensorAoS
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
                // Scalar full 8x8 local matrix — current best on GH200
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

            case CvfemKernelVariant::TensorColored: {
                // Graph-colored tensor kernel: atomic-free CSR scatter.
                // Requires config.coloring to be set (built once from connectivity).
                // Launches one kernel per color; within a color no two elements share a node.
                const CvfemColoringData* col = config.coloring;
                for (int c = 0; c < col->numColors; ++c) {
                    const int  cSize   = col->colorSize(c);
                    const int* cList   = col->colorStart(c);
                    const int  cBlocks = (cSize + blockSize - 1) / blockSize;
                    cvfem_hex_assembly_kernel_tensor_colored<KeyType, RealType, 256>
                        <<<cBlocks, blockSize, 0, config.stream>>>(
                        cList, static_cast<size_t>(cSize),
                        d_conn0, d_conn1, d_conn2, d_conn3,
                        d_conn4, d_conn5, d_conn6, d_conn7,
                        d_x, d_y, d_z,
                        d_gamma, d_phi, d_beta,
                        d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                        d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                        d_node_to_dof, d_ownership,
                        d_matrix, d_rhs
                    );
                }
                break;
            }

            case CvfemKernelVariant::TensorAoS:
                // AoS node data: packs 9 fields per node into one 72-byte struct.
                // Reduces scattered node reads from 9 × 1 sector/field to 3 sectors/node.
                cvfem_hex_assembly_kernel_tensor_aos<KeyType, RealType, 256><<<numBlocks, blockSize, 0, config.stream>>>(
                    d_conn0, d_conn1, d_conn2, d_conn3,
                    d_conn4, d_conn5, d_conn6, d_conn7,
                    numElements,
                    config.nodeData,
                    d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                    d_node_to_dof, d_ownership,
                    d_matrix, d_rhs
                );
                break;

            case CvfemKernelVariant::TensorPerip:
                // Pre-lookup all 64 CSR positions, scatter LHS per SCS.
                // lhs[64 doubles] → pos[64 ints]: saves ~64 registers.
                // Binary search moved out of 12-iteration SCS hot loop.
                cvfem_hex_assembly_kernel_tensor_perip<KeyType, RealType, 256><<<numBlocks, blockSize, 0, config.stream>>>(
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

            case CvfemKernelVariant::TensorPeripLb2:
                // Same as TensorPerip + __launch_bounds__(256,2).
                // Forces 2 blocks/SM with ~7 regs of forced spilling (vs ~44 for tensor+lb2).
                cvfem_hex_assembly_kernel_tensor_perip_lb2<KeyType, RealType, 256><<<numBlocks, blockSize, 0, config.stream>>>(
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

            case CvfemKernelVariant::SmemCache: {
                // Smem node cache: deduplicate + cooperatively load ~400 unique nodes/block.
                // Uses dynamic smem (>48 KB limit): unlock via MaxDynamicSharedMemorySize.
                // Node reads: ~30 cycle L2 → ~5 cycle smem. rowPtr/diagPtr also smem-cached.
                constexpr int    SC_HT   = 1024;
                constexpr int    SC_MN   = 512;
                constexpr size_t smemSz  = smemCacheBytes<SC_HT, SC_MN>();
                auto* scPtr = cvfem_hex_assembly_kernel_smem_cache<KeyType, RealType, 256, SC_HT, SC_MN>;
                cudaFuncSetAttribute(scPtr, cudaFuncAttributeMaxDynamicSharedMemorySize,    (int)smemSz);
                cudaFuncSetAttribute(scPtr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
                scPtr<<<numBlocks, blockSize, smemSz, config.stream>>>(
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

            case CvfemKernelVariant::WmmaTensor: {
                // FP64 WMMA: 1 warp (32 threads) per element, diffusion = S^T x D
                // Optimised for Hopper/GH200 (SM90a): hardware smem atomics, 228 KB smem.
                // Request 100% shared memory carveout so all 228 KB go to smem not L1.
                auto* wmmaKernelPtr = cvfem_hex_assembly_kernel_tensor_wmma<KeyType, RealType, 256>;
                cudaFuncSetAttribute(wmmaKernelPtr,
                    cudaFuncAttributePreferredSharedMemoryCarveout, 100);
                constexpr int WarpsPerBlock = 256 / 32;  // 8 elements per block
                const int wmmaBlocks = (numElements + WarpsPerBlock - 1) / WarpsPerBlock;
                const size_t smemSize = WarpsPerBlock * sizeof(CvfemWmmaWarpData);
                wmmaKernelPtr<<<wmmaBlocks, 256, smemSize, config.stream>>>(
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

            case CvfemKernelVariant::WgmmaTensor: {
                // Hopper WGMMA: 2 elem/warp (75% Jacobian util), async WGMMA m64n8k4 for diffusion.
                // BlockSize=128 = 1 warpgroup = 8 elements/block.
                // ~29 KB smem/block allows up to 7 blocks/SM from smem.
                auto* wgmmaKernelPtr = cvfem_hex_assembly_kernel_wgmma<KeyType, RealType, 128>;
                cudaFuncSetAttribute(wgmmaKernelPtr,
                    cudaFuncAttributePreferredSharedMemoryCarveout, 100);
                constexpr int ElemsPerBlock = 8;  // 4 warps × 2 elements/warp
                const int wgmmaBlocks = (numElements + ElemsPerBlock - 1) / ElemsPerBlock;
                const size_t smemSize = sizeof(CvfemWgmaWarpgroupData);
                wgmmaKernelPtr<<<wgmmaBlocks, 128, smemSize, config.stream>>>(
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

        // Full assembly supports kernel variants (but no diagonal lumping)
        switch (config.variant) {
            case CvfemKernelVariant::Original:
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
                cvfem_hex_assembly_kernel_shmem<KeyType, RealType><<<numBlocks, blockSize, 0, config.stream>>>(
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

            case CvfemKernelVariant::Tensor:
                cvfem_hex_assembly_kernel_tensor<KeyType, RealType><<<numBlocks, blockSize, 0, config.stream>>>(
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

            case CvfemKernelVariant::TensorColored: {
                const CvfemColoringData* col = config.coloring;
                for (int c = 0; c < col->numColors; ++c) {
                    const int  cSize   = col->colorSize(c);
                    const int* cList   = col->colorStart(c);
                    const int  cBlocks = (cSize + blockSize - 1) / blockSize;
                    cvfem_hex_assembly_kernel_tensor_colored<KeyType, RealType, 256>
                        <<<cBlocks, blockSize, 0, config.stream>>>(
                        cList, static_cast<size_t>(cSize),
                        d_conn0, d_conn1, d_conn2, d_conn3,
                        d_conn4, d_conn5, d_conn6, d_conn7,
                        d_x, d_y, d_z,
                        d_gamma, d_phi, d_beta,
                        d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                        d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z,
                        d_node_to_dof, d_ownership,
                        d_matrix, d_rhs
                    );
                }
                break;
            }

            case CvfemKernelVariant::WmmaTensor: {
                auto* wmmaKernelPtr = cvfem_hex_assembly_kernel_tensor_wmma<KeyType, RealType, 256>;
                cudaFuncSetAttribute(wmmaKernelPtr,
                    cudaFuncAttributePreferredSharedMemoryCarveout, 100);
                constexpr int WarpsPerBlock = 256 / 32;
                const int wmmaBlocks = (numElements + WarpsPerBlock - 1) / WarpsPerBlock;
                const size_t smemSize = WarpsPerBlock * sizeof(CvfemWmmaWarpData);
                wmmaKernelPtr<<<wmmaBlocks, 256, smemSize, config.stream>>>(
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

            case CvfemKernelVariant::WgmmaTensor: {
                auto* wgmmaKernelPtr = cvfem_hex_assembly_kernel_wgmma<KeyType, RealType, 128>;
                cudaFuncSetAttribute(wgmmaKernelPtr,
                    cudaFuncAttributePreferredSharedMemoryCarveout, 100);
                constexpr int ElemsPerBlock = 8;
                const int wgmmaBlocks = (numElements + ElemsPerBlock - 1) / ElemsPerBlock;
                const size_t smemSize = sizeof(CvfemWgmaWarpgroupData);
                wgmmaKernelPtr<<<wgmmaBlocks, 128, smemSize, config.stream>>>(
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

            case CvfemKernelVariant::Team: {
                int teamBlocks = (numElements * 32 + blockSize - 1) / blockSize;
                cvfem_hex_assembly_kernel_team<KeyType, RealType><<<teamBlocks, blockSize, 0, config.stream>>>(
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
    }

    // Get kernel variant name for logging
    static const char* variantName(CvfemKernelVariant variant) {
        switch (variant) {
            case CvfemKernelVariant::Original:   return "original";
            case CvfemKernelVariant::Optimized:  return "optimized";
            case CvfemKernelVariant::Shmem:      return "shmem";
            case CvfemKernelVariant::Team:       return "team";
            case CvfemKernelVariant::Tensor:         return "tensor";
            case CvfemKernelVariant::TensorColored: return "tensor_colored";
            case CvfemKernelVariant::TensorAoS:      return "tensor_aos";
            case CvfemKernelVariant::TensorPerip:    return "tensor_perip";
            case CvfemKernelVariant::TensorPeripLb2: return "tensor_perip_lb2";
            case CvfemKernelVariant::SmemCache:      return "smem_cache";
            case CvfemKernelVariant::WmmaTensor:    return "wmma_tensor";
            case CvfemKernelVariant::WgmmaTensor: return "wgmma_tensor";
            default: return "unknown";
        }
    }
};

} // namespace fem
} // namespace mars
