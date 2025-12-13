#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_face_kernel.hpp"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <mpi.h>
#include <iomanip>

#include <thrust/pair.h>

using namespace mars;
using namespace mars::fem;

// GPU kernel to build sparsity pattern
template<typename KeyType>
__global__ void buildPatternKernel(
    const KeyType* d_faceToElementOffsets,
    const KeyType* d_faceToElementList,
    const uint8_t* d_isBoundaryFace,
    size_t numFaces,
    const KeyType* d_conn0,
    const KeyType* d_conn1,
    const KeyType* d_conn2,
    const KeyType* d_conn3,
    const KeyType* d_conn4,
    const KeyType* d_conn5,
    const KeyType* d_conn6,
    const KeyType* d_conn7,
    size_t elementCount,
    const int* d_nodeToLocalDof,
    const int* d_ownership,
    const KeyType* d_sfc_map,
    size_t nodeCount,
    const KeyType* d_faceNodes,
    int nodesPerFace,
    thrust::pair<int, int>* d_pattern,
    int* d_counter,
    int max_entries
) {
    int faceIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (faceIdx >= numFaces) return;
    if (d_isBoundaryFace[faceIdx]) return;
    
    int faceDofs[4];
    int numDofs = 0;
    for (int n = 0; n < nodesPerFace; ++n) {
        KeyType sfc = d_faceNodes[faceIdx * nodesPerFace + n];
        // binary search for local_id
        int left = 0, right = nodeCount - 1;
        int local_id = -1;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (d_sfc_map[mid] == sfc) {
                local_id = mid;
                break;
            }
            if (d_sfc_map[mid] < sfc) left = mid + 1;
            else right = mid - 1;
        }
        if (local_id >= 0 && d_ownership[local_id] != 0) {
            faceDofs[numDofs++] = d_nodeToLocalDof[local_id];
        }
    }
    // Add pairs
    for (int i = 0; i < numDofs; ++i) {
        for (int j = 0; j < numDofs; ++j) {
            int row = faceDofs[i];
            int col = faceDofs[j];
            int idx = atomicAdd(d_counter, 1);
            if (idx < max_entries) {
                d_pattern[idx] = thrust::make_pair(row, col);
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    
    // Default to pitz_daly.mesh (hex mesh from reference CVFEM), allow override from command line
    std::string meshFile = "../../../meshes/pitz_daly.mesh";
    if (argc >= 2) {
        meshFile = argv[1];
    }
    
    if (rank == 0) {
        std::cout << "=== MARS CVFEM Example ===" << std::endl;
        std::cout << "Mesh: " << meshFile << std::endl;
        std::cout << "Ranks: " << numRanks << std::endl;
    }
    
    // Create domain using MARS ElementDomain (HexTag for hexahedral meshes)
    using ElemTag = HexTag;
    using RealType = double;
    using KeyType = unsigned long;
    using Domain = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;
    
    // Load mesh - mesh is broadcast to all ranks
    Domain domain(meshFile, rank, numRanks);
    
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    
    if (rank == 0) {
        std::cout << "Mesh loaded: " << nodeCount << " nodes, " 
                  << elementCount << " elements" << std::endl;
    }
    
    // Create FE space and DOF handler (like other MARS examples)
    using DofHandlerT = mars::fem::UnstructuredDofHandler<ElemTag, RealType, KeyType, cstone::GpuTag>;
    using FESpaceT = mars::fem::H1FESpace<ElemTag, RealType, KeyType, cstone::GpuTag>;
    using SparseMatrixType = mars::fem::SparseMatrix<KeyType, RealType, cstone::GpuTag>;
    
    FESpaceT fes(domain, 1);  // Linear elements
    DofHandlerT dof_handler(domain, rank, numRanks);
    dof_handler.initialize();
    dof_handler.enumerate_dofs();
    
    size_t numLocalDofs = dof_handler.get_num_local_dofs();
    size_t numLocalDofsWithGhosts = dof_handler.get_num_local_dofs_with_ghosts();
    
    if (rank == 0) {
        std::cout << "FE space: " << numLocalDofs << " local DOFs" << std::endl;
    }
    
    // Create sparse matrix for CVFEM system
    SparseMatrixType A(numLocalDofs, numLocalDofsWithGhosts);
    
    // Create RHS vector
    cstone::DeviceVector<RealType> b(numLocalDofs, 0.0);
    
    // Get mappings needed for kernel
    const auto& nodeToLocalDof = dof_handler.get_node_to_local_dof();
    const auto& resolvedOwnership = dof_handler.get_resolved_ownership();
    
    // Upload host data to device and cast to int for kernel
    cstone::DeviceVector<KeyType> d_nodeToLocalDof(nodeToLocalDof.size());
    thrust::copy(nodeToLocalDof.begin(), nodeToLocalDof.end(), thrust::device_pointer_cast(d_nodeToLocalDof.data()));
    cstone::DeviceVector<uint8_t> d_resolvedOwnership(resolvedOwnership.size());
    thrust::copy(resolvedOwnership.begin(), resolvedOwnership.end(), thrust::device_pointer_cast(d_resolvedOwnership.data()));
    cstone::DeviceVector<int> d_nodeToLocalDof_int(nodeCount);
    cstone::DeviceVector<int> d_ownership_int(nodeCount);
    
    thrust::transform(thrust::device_pointer_cast(d_nodeToLocalDof.data()),
                      thrust::device_pointer_cast(d_nodeToLocalDof.data() + d_nodeToLocalDof.size()),
                      thrust::device_pointer_cast(d_nodeToLocalDof_int.data()),
                      [] __host__ __device__ (KeyType k) { return static_cast<int>(k); });
    thrust::transform(thrust::device_pointer_cast(d_resolvedOwnership.data()),
                      thrust::device_pointer_cast(d_resolvedOwnership.data() + d_resolvedOwnership.size()),
                      thrust::device_pointer_cast(d_ownership_int.data()),
                      [] __host__ __device__ (uint8_t k) { return static_cast<int>(k); });
    
    // Build face topology (production CVFEM uses faces, not elements)
    if (rank == 0) {
        std::cout << "Building face topology..." << std::endl;
    }
    mars::FaceTopology<ElemTag, RealType, KeyType, cstone::GpuTag> faceTopology(domain);
    
    if (rank == 0) {
        std::cout << "Face topology built:" << std::endl;
        std::cout << "  Total faces: " << faceTopology.numFaces_ << std::endl;
        std::cout << "  Boundary faces: " << faceTopology.numBoundaryFaces_ << std::endl;
        std::cout << "  Interior faces: " << faceTopology.numInteriorFaces_ << std::endl;
        
        // Verify connectivity
        std::vector<KeyType> h_offsets(std::min(size_t(10), faceTopology.numFaces_ + 1));
        thrust::copy(thrust::device_pointer_cast(faceTopology.d_faceToElementOffsets_.data()),
                     thrust::device_pointer_cast(faceTopology.d_faceToElementOffsets_.data() + h_offsets.size()),
                     h_offsets.begin());
        std::cout << "  First few face element counts: ";
        for (size_t i = 0; i < h_offsets.size() - 1; ++i) {
            std::cout << (h_offsets[i+1] - h_offsets[i]) << " ";
        }
        std::cout << std::endl;
    }
    
    // Create int versions of sparse matrix indices for kernel (kernel expects int*)
    cstone::DeviceVector<int> d_A_rowptr_int(A.numRows() + 1);
    cstone::DeviceVector<int> d_A_colidx_int(A.nnz());
    
    // Get coordinates
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();
    
    // Get SFC to local ID mapping
    const auto& d_sfc_map = domain.getLocalToGlobalSfcMap();
    
    // Build sparsity pattern on GPU
    if (rank == 0) {
        std::cout << "Building sparsity pattern..." << std::endl;
    }
    size_t max_entries = faceTopology.numFaces_ * 16;
    thrust::device_vector<thrust::pair<int, int>> d_pattern(max_entries);
    thrust::device_vector<int> d_counter(1, 0);
    
    // Get connectivity pointers
    auto& d_conn = domain.getElementToNodeConnectivity();
    auto d_conn0 = std::get<0>(d_conn).data();
    auto d_conn1 = std::get<1>(d_conn).data();
    auto d_conn2 = std::get<2>(d_conn).data();
    auto d_conn3 = std::get<3>(d_conn).data();
    auto d_conn4 = std::get<4>(d_conn).data();
    auto d_conn5 = std::get<5>(d_conn).data();
    auto d_conn6 = std::get<6>(d_conn).data();
    auto d_conn7 = std::get<7>(d_conn).data();
    
    // Launch kernel
    int blockSizePattern = 256;
    int numBlocksPattern = (faceTopology.numFaces_ + blockSizePattern - 1) / blockSizePattern;
    buildPatternKernel<KeyType><<<numBlocksPattern, blockSizePattern>>>(
        faceTopology.d_faceToElementOffsets_.data(),
        faceTopology.d_faceToElementList_.data(),
        faceTopology.d_isBoundaryFace_.data(),
        faceTopology.numFaces_,
        d_conn0, d_conn1, d_conn2, d_conn3, d_conn4, d_conn5, d_conn6, d_conn7,
        elementCount,
        d_nodeToLocalDof_int.data(),
        d_ownership_int.data(),
        d_sfc_map.data(),
        nodeCount,
        faceTopology.d_faceNodes_.data(),
        4,
        thrust::raw_pointer_cast(d_pattern.data()),
        thrust::raw_pointer_cast(d_counter.data()),
        max_entries
    );
    
    // Get num_entries
    int num_entries;
    thrust::copy(d_counter.begin(),
                 d_counter.begin() + 1,
                 &num_entries);
    
    // Sort and unique
    thrust::sort(d_pattern.begin(),
                 d_pattern.begin() + num_entries);
    auto new_end = thrust::unique(d_pattern.begin(),
                                  d_pattern.begin() + num_entries);
    num_entries = new_end - d_pattern.begin();
    
    // Copy to host for CSR building
    std::vector<thrust::pair<int, int>> h_pattern(num_entries);
    thrust::copy(d_pattern.begin(),
                 d_pattern.begin() + num_entries,
                 h_pattern.begin());
    
    // Build CSR
    std::vector<int> rowPtr(numLocalDofs + 1, 0);
    for (auto& p : h_pattern) {
        rowPtr[p.first + 1]++;
    }
    for (size_t i = 1; i <= numLocalDofs; ++i) {
        rowPtr[i] += rowPtr[i-1];
    }
    std::vector<int> colInd(num_entries);
    std::vector<int> currentPos = rowPtr;
    for (auto& p : h_pattern) {
        int row = p.first;
        int col = p.second;
        int pos = currentPos[row]++;
        colInd[pos] = col;
    }
    
    // Allocate matrix
    A.allocate(numLocalDofs, numLocalDofsWithGhosts, num_entries);
    thrust::copy(rowPtr.begin(), rowPtr.end(), thrust::device_pointer_cast(A.rowOffsets().data()));
    thrust::copy(colInd.begin(), colInd.end(), thrust::device_pointer_cast(A.colIndices().data()));
    A.zero();
    
    // Resize int vectors
    d_A_colidx_int.resize(A.nnz());
    
    // Update int versions
    thrust::transform(thrust::device_pointer_cast(A.rowOffsets().data()),
                      thrust::device_pointer_cast(A.rowOffsets().data() + A.numRows() + 1),
                      thrust::device_pointer_cast(d_A_rowptr_int.data()),
                      [] __host__ __device__ (KeyType k) { return static_cast<int>(k); });
    thrust::transform(thrust::device_pointer_cast(A.colIndices().data()),
                      thrust::device_pointer_cast(A.colIndices().data() + A.nnz()),
                      thrust::device_pointer_cast(d_A_colidx_int.data()),
                      [] __host__ __device__ (KeyType k) { return static_cast<int>(k); });
    
    if (rank == 0) {
        std::cout << "Sparsity pattern built: " << A.nnz() << " non-zeros" << std::endl;
    }
    
    // Create fields on device (node-based) - match STL field values
    cstone::DeviceVector<RealType> d_gamma(nodeCount, 0.1);        // Diffusion coefficient (STL: 0.1)
    cstone::DeviceVector<RealType> d_phi(nodeCount);               // Solution (will be filled)
    cstone::DeviceVector<RealType> d_beta(nodeCount, 1.234);       // Convection coefficient (STL: 1.234)
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount);        // Gradient x (will be filled)
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount);        // Gradient y (will be filled)
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount);        // Gradient z (will be filled)
    
    // Mass flow per face (for advection) - match STL: 0 for no advection
    cstone::DeviceVector<RealType> d_mdot(faceTopology.numFaces_, 0.0);
    
    // Fill phi and gradPhi like STL (assuming manufactured solution phi = x + y + z)
    // Coordinates already cached above
    
    // Copy coordinates to host for filling
    std::vector<RealType> h_x(nodeCount), h_y(nodeCount), h_z(nodeCount);
    thrust::copy(thrust::device_pointer_cast(d_x.data()), thrust::device_pointer_cast(d_x.data() + nodeCount), h_x.begin());
    thrust::copy(thrust::device_pointer_cast(d_y.data()), thrust::device_pointer_cast(d_y.data() + nodeCount), h_y.begin());
    thrust::copy(thrust::device_pointer_cast(d_z.data()), thrust::device_pointer_cast(d_z.data() + nodeCount), h_z.begin());
    
    // Fill phi and gradPhi - DEBUG: set to 0 to check if kernel produces zero matrix
    // std::vector<RealType> h_phi(nodeCount), h_grad_x(nodeCount, 1.0), h_grad_y(nodeCount, 1.0), h_grad_z(nodeCount, 1.0);
    // for (size_t i = 0; i < nodeCount; ++i) {
    //     h_phi[i] = h_x[i] + h_y[i] + h_z[i];
    // }
    std::vector<RealType> h_phi(nodeCount, 0.0), h_grad_x(nodeCount, 0.0), h_grad_y(nodeCount, 0.0), h_grad_z(nodeCount, 0.0);
    
    // Copy to device
    thrust::copy(h_phi.begin(), h_phi.end(), thrust::device_pointer_cast(d_phi.data()));
    thrust::copy(h_grad_x.begin(), h_grad_x.end(), thrust::device_pointer_cast(d_grad_phi_x.data()));
    thrust::copy(h_grad_y.begin(), h_grad_y.end(), thrust::device_pointer_cast(d_grad_phi_y.data()));
    thrust::copy(h_grad_z.begin(), h_grad_z.end(), thrust::device_pointer_cast(d_grad_phi_z.data()));
    
    // Int versions already created above
    
    if (rank == 0) {
        std::cout << "Fields initialized, assembling with face-based CVFEM..." << std::endl;
    }
    
    // Launch face-based CVFEM assembly kernel (production approach)
    int blockSize = 256;
    int numBlocks = (faceTopology.numFaces_ + blockSize - 1) / blockSize;

    cvfem_face_assembly_kernel<KeyType, RealType><<<numBlocks, blockSize>>>(
        // Face topology
        faceTopology.d_faceNodes_.data(),
        faceTopology.d_faceToElementOffsets_.data(),
        faceTopology.d_faceToElementList_.data(),
        faceTopology.d_isBoundaryFace_.data(),
        faceTopology.d_faceNormalX_.data(),
        faceTopology.d_faceNormalY_.data(),
        faceTopology.d_faceNormalZ_.data(),
        faceTopology.d_faceArea_.data(),
        faceTopology.numFaces_,
        4,  // nodesPerFace for hex
        // Node data
        d_sfc_map.data(),
        nodeCount,
        d_x.data(),
        d_y.data(),
        d_z.data(),
        // Field data
        thrust::raw_pointer_cast(d_gamma.data()),
        thrust::raw_pointer_cast(d_phi.data()),
        thrust::raw_pointer_cast(d_beta.data()),
        thrust::raw_pointer_cast(d_grad_phi_x.data()),
        thrust::raw_pointer_cast(d_grad_phi_y.data()),
        thrust::raw_pointer_cast(d_grad_phi_z.data()),
        thrust::raw_pointer_cast(d_mdot.data()),
        // Sparse matrix data
        thrust::raw_pointer_cast(A.values().data()),
        thrust::raw_pointer_cast(d_A_colidx_int.data()),
        thrust::raw_pointer_cast(d_A_rowptr_int.data()),
        // DOF mappings
        thrust::raw_pointer_cast(d_nodeToLocalDof_int.data()),
        thrust::raw_pointer_cast(d_ownership_int.data()),
        // RHS vector
        thrust::raw_pointer_cast(b.data()),
        false  // include advection (false to match STL with mdot=0)
    );
    
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        std::cout << "Assembly complete" << std::endl;
    }
    
    // Compute RHS norm
    RealType local_norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(b.data()),
        thrust::device_pointer_cast(b.data() + numLocalDofs),
        [] __host__ __device__ (RealType x) { return x * x; },
        0.0,
        thrust::plus<RealType>()
    );
    
    RealType global_norm2;
    MPI_Allreduce(&local_norm2, &global_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    RealType rhs_norm = std::sqrt(global_norm2);
    
    // Compute matrix norm (L2/Frobenius norm - sqrt(sum of squares))
    RealType local_matrix_norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(A.values().data()),
        thrust::device_pointer_cast(A.values().data() + A.nnz()),
        [] __host__ __device__ (RealType x) { return x * x; },
        0.0,
        thrust::plus<RealType>()
    );
    
    RealType global_matrix_norm2;
    MPI_Allreduce(&local_matrix_norm2, &global_matrix_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    RealType matrix_norm = std::sqrt(global_matrix_norm2);
    
    if (rank == 0) {
        std::cout << "\n=== Results ===" << std::endl;
        std::cout << std::scientific;
        std::cout << "MARS CVFEM: [matrix norm: " << matrix_norm << "]" << std::endl;
        std::cout.unsetf(std::ios::scientific);
    }
    
    MPI_Finalize();
    return 0;
}
