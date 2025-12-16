#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel.hpp"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <mpi.h>
#include <iomanip>
#include <chrono>

using namespace mars;
using namespace mars::fem;

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    if (argc < 2) {
        if (rank == 0) {
            std::cout << "Usage: " << argv[0] << " <mesh_file>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string meshFile = argv[1];

    using KeyType = unsigned int;
    using RealType = double;
    using ElemTag = HexTag;

    // Create domain
    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(meshFile, rank, numRanks);

    // Build halo for node ownership (this triggers SFC map creation which updates node count)
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();

    // Get counts AFTER node ownership map is built (node count is updated during SFC map creation)
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();  // Total elements (includes halos for kernel)

    // Get element range for local (owned) elements
    // After sync, elements are ordered as: [halos_before, local, halos_after]
    auto [startIdx, endIdx] = domain.localElementRange();
    size_t localElementCount = endIdx - startIdx;

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "DEBUG: startIdx=" << startIdx << " endIdx=" << endIdx << " diff=" << localElementCount << std::endl;
        std::cout << "DEBUG: elementCount (total with halos)=" << elementCount << std::endl;
        std::cout << "Mesh Size (Nodes):    " << nodeCount << std::endl;
        std::cout << "Mesh Size (Elements): " << localElementCount << std::endl;
        std::cout << std::endl;
    }

    // Count owned DOFs
    std::vector<uint8_t> h_ownership(nodeCount);
    thrust::copy(thrust::device_pointer_cast(d_nodeOwnership.data()),
                 thrust::device_pointer_cast(d_nodeOwnership.data() + nodeCount),
                 h_ownership.begin());

    int numOwnedDofs = 0;
    for (auto o : h_ownership) {
        if (o == 1 || o == 2) numOwnedDofs++;  // owned or shared
    }

    // Create DOF mapping (simple: local node ID = DOF ID for owned nodes)
    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    std::vector<int> h_nodeToDof(nodeCount);
    int dofCounter = 0;
    for (size_t i = 0; i < nodeCount; ++i) {
        if (h_ownership[i] != 0) {  // owned or shared
            h_nodeToDof[i] = dofCounter++;
        } else {
            h_nodeToDof[i] = -1;  // ghost node
        }
    }
    thrust::copy(h_nodeToDof.begin(), h_nodeToDof.end(),
                 thrust::device_pointer_cast(d_nodeToDof.data()));

    // Build sparsity pattern from LOCAL elements only (not halos)
    if (rank == 0) {
        std::cout << "Building sparsity pattern..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::set<std::pair<int,int>> pattern_set;
    const auto& d_conn = domain.getElementToNodeConnectivity();

    // Copy LOCAL elements from [startIdx, endIdx) range
    std::vector<KeyType> h_conn[8];
    for (int i = 0; i < 8; ++i) {
        h_conn[i].resize(localElementCount);
        auto& d_conni = (i==0) ? std::get<0>(d_conn) :
                        (i==1) ? std::get<1>(d_conn) :
                        (i==2) ? std::get<2>(d_conn) :
                        (i==3) ? std::get<3>(d_conn) :
                        (i==4) ? std::get<4>(d_conn) :
                        (i==5) ? std::get<5>(d_conn) :
                        (i==6) ? std::get<6>(d_conn) :
                                 std::get<7>(d_conn);
        // Copy from startIdx to endIdx (local element range)
        thrust::copy(thrust::device_pointer_cast(d_conni.data() + startIdx),
                     thrust::device_pointer_cast(d_conni.data() + endIdx),
                     h_conn[i].begin());
    }

    for (size_t e = 0; e < localElementCount; ++e) {
        KeyType nodes[8];
        for (int i = 0; i < 8; ++i) {
            nodes[i] = h_conn[i][e];
            // Bounds check
            if (nodes[i] >= nodeCount) {
                std::cerr << "Rank " << rank << ": Element " << e << " node " << i
                          << " index " << nodes[i] << " exceeds nodeCount " << nodeCount << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        for (int i = 0; i < 8; ++i) {
            int row_dof = h_nodeToDof[nodes[i]];
            if (row_dof < 0) continue;  // skip ghost

            for (int j = 0; j < 8; ++j) {
                int col_dof = h_nodeToDof[nodes[j]];
                if (col_dof < 0) continue;  // skip ghost

                pattern_set.insert({row_dof, col_dof});
            }
        }
    }

    // Build CSR structure
    std::vector<std::pair<int,int>> pattern_vec(pattern_set.begin(), pattern_set.end());
    std::sort(pattern_vec.begin(), pattern_vec.end());

    std::vector<int> rowPtr(numOwnedDofs + 1, 0);
    for (auto& p : pattern_vec) {
        rowPtr[p.first + 1]++;
    }
    for (int i = 1; i <= numOwnedDofs; ++i) {
        rowPtr[i] += rowPtr[i-1];
    }

    int nnz = pattern_vec.size();
    std::vector<int> colInd(nnz);
    for (int i = 0; i < nnz; ++i) {
        colInd[i] = pattern_vec[i].second;
    }

    if (rank == 0) {
        std::cout << "Sparsity pattern built. NNZ = " << nnz << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate matrix and RHS on device
    cstone::DeviceVector<int> d_rowPtr(rowPtr.data(), rowPtr.data() + rowPtr.size());
    cstone::DeviceVector<int> d_colInd(colInd.data(), colInd.data() + colInd.size());
    cstone::DeviceVector<RealType> d_values(nnz, 0.0);
    cstone::DeviceVector<RealType> d_rhs(numOwnedDofs, 0.0);

    // Initialize fields
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    cstone::DeviceVector<RealType> d_gamma(nodeCount, 0.1);     // diffusion
    cstone::DeviceVector<RealType> d_phi(nodeCount, 1.0);       // solution
    cstone::DeviceVector<RealType> d_beta(nodeCount, 1.0);      // upwind
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount, 0.0);

    // Element-based data (12 SCS per hex element)
    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 1.0);
    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12, 1.0);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12, 0.0);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12, 0.0);

    // Create simple CSR wrapper (matches kernel definition)
    using MatrixType = mars::fem::CSRMatrix<RealType>;

    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{
        d_rowPtr.data(),
        d_colInd.data(),
        d_values.data(),
        numOwnedDofs,
        nnz
    };
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    // Launch assembly kernel with timing
    int blockSize = 256;
    int numBlocks = (elementCount + blockSize - 1) / blockSize;

    // Warm-up run
    fem::cvfem_hex_assembly_kernel<<<numBlocks, blockSize>>>(
        std::get<0>(d_conn).data(),
        std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(),
        std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(),
        std::get<7>(d_conn).data(),
        elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(),
        d_phi.data(),
        d_beta.data(),
        d_grad_phi_x.data(),
        d_grad_phi_y.data(),
        d_grad_phi_z.data(),
        d_mdot.data(),
        d_areaVec_x.data(),
        d_areaVec_y.data(),
        d_areaVec_z.data(),
        d_nodeToDof.data(),
        d_nodeOwnership.data(),
        d_matrix,
        d_rhs.data()
    );
    cudaDeviceSynchronize();

    // Reset for timed run
    thrust::fill(thrust::device_pointer_cast(d_values.data()),
                 thrust::device_pointer_cast(d_values.data() + nnz), 0.0);
    thrust::fill(thrust::device_pointer_cast(d_rhs.data()),
                 thrust::device_pointer_cast(d_rhs.data() + numOwnedDofs), 0.0);

    // Timed run
    auto start = std::chrono::high_resolution_clock::now();

    fem::cvfem_hex_assembly_kernel<<<numBlocks, blockSize>>>(
        std::get<0>(d_conn).data(),
        std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(),
        std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(),
        std::get<7>(d_conn).data(),
        elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(),
        d_phi.data(),
        d_beta.data(),
        d_grad_phi_x.data(),
        d_grad_phi_y.data(),
        d_grad_phi_z.data(),
        d_mdot.data(),
        d_areaVec_x.data(),
        d_areaVec_y.data(),
        d_areaVec_z.data(),
        d_nodeToDof.data(),
        d_nodeOwnership.data(),
        d_matrix,
        d_rhs.data()
    );
    cudaDeviceSynchronize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    // Compute matrix norm (L2/Frobenius) - matching STK reference
    auto square_op = [] __host__ __device__ (RealType x) -> RealType {
        return x * x;
    };

    RealType local_matrix_norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(d_values.data()),
        thrust::device_pointer_cast(d_values.data() + nnz),
        square_op,
        RealType(0.0),
        thrust::plus<RealType>()
    );

    RealType local_rhs_norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(d_rhs.data()),
        thrust::device_pointer_cast(d_rhs.data() + numOwnedDofs),
        square_op,
        RealType(0.0),
        thrust::plus<RealType>()
    );

    // Gather squared norms and timing from all ranks
    RealType global_matrix_norm2, global_rhs_norm2;
    MPI_Allreduce(&local_matrix_norm2, &global_matrix_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_rhs_norm2, &global_rhs_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    RealType matrix_norm = std::sqrt(global_matrix_norm2);
    RealType rhs_norm = std::sqrt(global_rhs_norm2);

    // Calculate throughput (bytes read/written per element assembly)
    // Estimate: 8 nodes * (3 coords + 4 fields + 3 grads) * 8 bytes = ~640 bytes per element
    // Plus output: 64 matrix entries + 8 RHS entries * 8 bytes = ~576 bytes
    size_t bytes_per_element = 8 * (3 + 4 + 3) * sizeof(RealType) +
                                (64 + 8) * sizeof(RealType);
    double total_bytes_gb = (bytes_per_element * elementCount) / 1e9;
    double throughput_gbs = total_bytes_gb / (elapsed.count() / 1000.0);

    if (rank == 0) {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Assembler: MARS CVFEM Hex..................: "
                  << elapsed.count() << " milliseconds @ "
                  << throughput_gbs << " GB/s (average of 1 samples) [matrix norm: "
                  << matrix_norm << "]" << std::endl;
        std::cout << std::endl;
    }

    cudaFree(d_matrix);
    MPI_Finalize();
    return 0;
}
