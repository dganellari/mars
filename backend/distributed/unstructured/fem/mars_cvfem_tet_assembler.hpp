#pragma once

#include "mars_cvfem_tet_kernel_graph.hpp"
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// CVFEM assembler for linear tetrahedra. Counterpart to CvfemHexAssembler
// in mars_cvfem_assembler.hpp. Only the graph-sparsity path is implemented
// here; tet does not yet have the tensor/wmma/shmem variants the hex
// assembler offers.
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

    static const char* variantName() { return "TetGraph"; }
};

} // namespace fem
} // namespace mars
