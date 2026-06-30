#pragma once

#include "mars_cvfem_tet_kernel_graph.hpp"
#include "mars_cvfem_tet_kernel_full.hpp"
#include "mars_cvfem_tet_kernel_full_perip.hpp"
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// CVFEM assembler for linear tetrahedra. Counterpart to CvfemHexAssembler
// in mars_cvfem_assembler.hpp. Three paths exist: graph-lumped, full, and
// full-perip (full with pre-looked-up CSR positions). The tensor/wmma/shmem
// variants the hex assembler offers do not exist for tet.
template<typename KeyType, typename RealType>
class CvfemTetAssembler {
public:
    struct Config {
        int blockSize     = 256;
        cudaStream_t stream = 0;
    };

    static void assembleGraphLump(
        // Connectivity: 4 columns (tet has 4 corner nodes)
        const KeyType* d_conn0,
        const KeyType* d_conn1,
        const KeyType* d_conn2,
        const KeyType* d_conn3,
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
        // Element data: mdot per (element, SCS), 6 SCS per tet
        const RealType* d_mdot,
        // DOF mapping
        const int* d_node_to_dof,
        const uint8_t* d_ownership,
        // Output system
        CSRMatrix<RealType>* d_matrix,
        RealType* d_rhs,
        const Config& config = Config{})
    {
        int blockSize = config.blockSize;
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        cvfem_tet_assembly_kernel_graph<KeyType, RealType>
            <<<numBlocks, blockSize, 0, config.stream>>>(
                d_conn0, d_conn1, d_conn2, d_conn3,
                numElements,
                d_x, d_y, d_z,
                d_gamma, d_phi, d_beta,
                d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                d_mdot,
                d_node_to_dof, d_ownership,
                d_matrix, d_rhs);
    }

    // Full assembly: scatters the full 4x4 element matrix into a 16-NNZ
    // sparsity pattern. No diagonal lumping. Required when the consumer
    // needs every (i,j) pair (e.g. Poisson with CG, symmetric solvers).
    static void assembleFull(
        const KeyType* d_conn0,
        const KeyType* d_conn1,
        const KeyType* d_conn2,
        const KeyType* d_conn3,
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
        const int* d_node_to_dof,
        const uint8_t* d_ownership,
        CSRMatrix<RealType>* d_matrix,
        RealType* d_rhs,
        const Config& config = Config{})
    {
        int blockSize = config.blockSize;
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        cvfem_tet_assembly_kernel_full<KeyType, RealType>
            <<<numBlocks, blockSize, 0, config.stream>>>(
                d_conn0, d_conn1, d_conn2, d_conn3,
                numElements,
                d_x, d_y, d_z,
                d_gamma, d_phi, d_beta,
                d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                d_mdot,
                d_node_to_dof, d_ownership,
                d_matrix, d_rhs);
    }

    // Full perip: same math as assembleFull, but pre-looks-up the 16 CSR
    // positions once per element and scatters via direct atomicAdd into
    // values[pos]. Avoids the per-(i,j) findColumnIndex inside the hot loop.
    // Uses __launch_bounds__(256, 4) to pack 4 blocks/SM. Recommended path
    // for production runs; assembleFull stays as the simpler reference.
    static void assembleFullPerip(
        const KeyType* d_conn0,
        const KeyType* d_conn1,
        const KeyType* d_conn2,
        const KeyType* d_conn3,
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
        const int* d_node_to_dof,
        const uint8_t* d_ownership,
        CSRMatrix<RealType>* d_matrix,
        RealType* d_rhs,
        const Config& config = Config{})
    {
        int blockSize = config.blockSize;
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        cvfem_tet_assembly_kernel_full_perip<KeyType, RealType>
            <<<numBlocks, blockSize, 0, config.stream>>>(
                d_conn0, d_conn1, d_conn2, d_conn3,
                numElements,
                d_x, d_y, d_z,
                d_gamma, d_phi, d_beta,
                d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                d_mdot,
                d_node_to_dof, d_ownership,
                d_matrix, d_rhs);
    }

    enum class Variant { GraphLump, Full, FullPerip };
    static const char* variantName(Variant v) {
        switch (v) {
            case Variant::GraphLump: return "TetGraphLump";
            case Variant::Full:      return "TetFull";
            case Variant::FullPerip: return "TetFullPerip";
        }
        return "TetUnknown";
    }
};

} // namespace fem
} // namespace mars
