#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include <thrust/device_vector.h>
#include <set>
#include <thrust/reduce.h>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <numeric>

using namespace mars;
using namespace mars::fem;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    if (argc < 2) {
        if (rank == 0) {
            std::cout << "Usage: " << argv[0] << " <mesh_file> [options]" << std::endl;
            std::cout << "  mesh_file: .mesh or .exo format mesh file" << std::endl;
            std::cout << "\nOptions:" << std::endl;
            std::cout << "  --iterations=N      Number of assembly iterations (default: 10)" << std::endl;
            std::cout << "  --block-size=N      CUDA block size (default: 256)" << std::endl;
            std::cout << "  --kernel=VARIANT    original, optimized, team (default: original)" << std::endl;
            std::cout << "  --quiet             Suppress detailed output" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string meshFile = argv[1];

    // Parse command-line options
    int numIterations = 10;
    int blockSize = 256;
    bool quiet = false;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Original;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--iterations=") == 0) {
            numIterations = std::stoi(arg.substr(13));
        } else if (arg.find("--block-size=") == 0) {
            blockSize = std::stoi(arg.substr(13));
        } else if (arg.find("--kernel=") == 0) {
            std::string v = arg.substr(9);
            if (v == "original") kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "team") kernelVariant = CvfemKernelVariant::Team;
            else {
                if (rank == 0) {
                    std::cerr << "Error: --kernel must be original, optimized, or team" << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "--quiet") {
            quiet = true;
        }
    }

    using KeyType = uint64_t;
    using RealType = double;
    using ElemTag = HexTag;
    using Assembler = CvfemHexAssembler<KeyType, RealType>;

    // Initialize performance counters and timer
    PerfCounters perf;
    GpuTimer timer;

    const char* kernelName = Assembler::variantName(kernelVariant);

    if (rank == 0 && !quiet) {
        std::cout << "Configuration:" << std::endl;
        std::cout << "  Iterations:      " << numIterations << std::endl;
        std::cout << "  Block size:      " << blockSize << std::endl;
        std::cout << "  Kernel:          " << kernelName << std::endl;
        std::cout << std::endl;
    }

    // Load mesh
    auto meshLoadStart = std::chrono::high_resolution_clock::now();
    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(meshFile, rank, numRanks, true);
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    cudaDeviceSynchronize();
    auto meshLoadEnd = std::chrono::high_resolution_clock::now();
    perf.meshLoadTime = std::chrono::duration<float, std::milli>(meshLoadEnd - meshLoadStart).count();

    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    perf.numNodes = nodeCount;
    perf.numElements = elementCount;

    if (rank == 0 && !quiet) {
        std::cout << "Mesh loaded in " << perf.meshLoadTime << " ms" << std::endl;
        std::cout << "  Total nodes:    " << nodeCount << std::endl;
        std::cout << "  Total elements: " << elementCount << std::endl;
        std::cout << std::endl;
    }

    // Build DOF mapping
    std::vector<uint8_t> h_ownership(nodeCount);
    thrust::copy(thrust::device_pointer_cast(d_nodeOwnership.data()),
                 thrust::device_pointer_cast(d_nodeOwnership.data() + nodeCount),
                 h_ownership.begin());

    int numOwnedDofs = 0;
    for (auto o : h_ownership) {
        if (o == 1 || o == 2) numOwnedDofs++;
    }

    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    std::vector<int> h_nodeToDof(nodeCount);
    int dofCounter = 0;
    for (size_t i = 0; i < nodeCount; ++i) {
        if (h_ownership[i] != 0) {
            h_nodeToDof[i] = dofCounter++;
        } else {
            h_nodeToDof[i] = -1;
        }
    }
    thrust::copy(h_nodeToDof.begin(), h_nodeToDof.end(),
                 thrust::device_pointer_cast(d_nodeToDof.data()));

    // Build reduced sparsity pattern from element connectivity
    auto sparsityStart = std::chrono::high_resolution_clock::now();
    if (rank == 0 && !quiet) {
        std::cout << "Building reduced sparsity pattern..." << std::endl;
    }

    const auto& d_conn = domain.getElementToNodeConnectivity();

    std::vector<KeyType> h_conn[8];
    for (int i = 0; i < 8; ++i) {
        h_conn[i].resize(elementCount);
        auto& d_conni = (i==0) ? std::get<0>(d_conn) :
                        (i==1) ? std::get<1>(d_conn) :
                        (i==2) ? std::get<2>(d_conn) :
                        (i==3) ? std::get<3>(d_conn) :
                        (i==4) ? std::get<4>(d_conn) :
                        (i==5) ? std::get<5>(d_conn) :
                        (i==6) ? std::get<6>(d_conn) :
                                 std::get<7>(d_conn);
        thrust::copy(thrust::device_pointer_cast(d_conni.data()),
                     thrust::device_pointer_cast(d_conni.data() + elementCount),
                     h_conn[i].begin());
    }

    // hexLRSCV: maps SCS index to left/right nodes
    static const int hexLRSCV[24] = {
        0, 1, 1, 2, 2, 3, 0, 3,
        4, 5, 5, 6, 6, 7, 4, 7,
        0, 4, 1, 5, 2, 6, 3, 7
    };

    std::vector<std::set<int>> adj(numOwnedDofs);

    for (size_t e = 0; e < elementCount; ++e) {
        KeyType nodes[8];
        for (int i = 0; i < 8; ++i) {
            nodes[i] = h_conn[i][e];
        }

        for (int scs = 0; scs < 12; ++scs) {
            int nodeL_local = hexLRSCV[scs * 2];
            int nodeR_local = hexLRSCV[scs * 2 + 1];

            KeyType nodeL = nodes[nodeL_local];
            KeyType nodeR = nodes[nodeR_local];

            int dofL = h_nodeToDof[nodeL];
            int dofR = h_nodeToDof[nodeR];

            if (dofL >= 0) {
                adj[dofL].insert(dofL);
                if (dofR >= 0) adj[dofL].insert(dofR);
            }
            if (dofR >= 0) {
                adj[dofR].insert(dofR);
                if (dofL >= 0) adj[dofR].insert(dofL);
            }
        }
    }

    std::vector<int> rowPtr(numOwnedDofs + 1);
    std::vector<int> colInd;

    rowPtr[0] = 0;
    for (int d = 0; d < numOwnedDofs; ++d) {
        for (int col : adj[d]) {
            colInd.push_back(col);
        }
        rowPtr[d + 1] = colInd.size();
    }

    int nnz = colInd.size();
    auto sparsityEnd = std::chrono::high_resolution_clock::now();
    perf.sparsityBuildTime = std::chrono::duration<float, std::milli>(sparsityEnd - sparsityStart).count();
    perf.numDofs = numOwnedDofs;
    perf.matrixNnz = nnz;

    if (rank == 0 && !quiet) {
        std::cout << "  NNZ: " << nnz << ", avg/row: "
                  << std::setprecision(1) << std::fixed
                  << static_cast<double>(nnz) / numOwnedDofs << std::endl;
    }

    // Allocate matrix and RHS on device
    cstone::DeviceVector<int> d_rowPtr(rowPtr.data(), rowPtr.data() + rowPtr.size());
    cstone::DeviceVector<int> d_colInd(colInd.data(), colInd.data() + colInd.size());
    cstone::DeviceVector<RealType> d_values(nnz, 0.0);
    cstone::DeviceVector<RealType> d_rhs(numOwnedDofs, 0.0);

    // Initialize fields
    auto fieldInitStart = std::chrono::high_resolution_clock::now();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    std::vector<RealType> h_x(nodeCount), h_y(nodeCount), h_z(nodeCount);
    thrust::copy(thrust::device_pointer_cast(d_x.data()),
                 thrust::device_pointer_cast(d_x.data() + nodeCount),
                 h_x.begin());
    thrust::copy(thrust::device_pointer_cast(d_y.data()),
                 thrust::device_pointer_cast(d_y.data() + nodeCount),
                 h_y.begin());
    thrust::copy(thrust::device_pointer_cast(d_z.data()),
                 thrust::device_pointer_cast(d_z.data() + nodeCount),
                 h_z.begin());

    std::vector<RealType> h_gamma(nodeCount);
    std::vector<RealType> h_phi(nodeCount);
    std::vector<RealType> h_beta(nodeCount);
    std::vector<RealType> h_grad_phi_x(nodeCount);
    std::vector<RealType> h_grad_phi_y(nodeCount);
    std::vector<RealType> h_grad_phi_z(nodeCount);

    for (size_t i = 0; i < nodeCount; ++i) {
        RealType x = h_x[i];
        RealType y = h_y[i];
        RealType z = h_z[i];

        h_phi[i] = std::sin(x) + 3.0 * std::cos(y) + 4.0 * std::sin(5.0 * x * y * z);
        h_grad_phi_x[i] = std::cos(x) + 20.0 * y * z * std::cos(5.0 * x * y * z);
        h_grad_phi_y[i] = 20.0 * x * z * std::cos(5.0 * x * y * z) - 3.0 * std::sin(y);
        h_grad_phi_z[i] = 20.0 * x * y * std::cos(5.0 * x * y * z);
        h_gamma[i] = 0.1;
        h_beta[i] = 1.234;
    }

    cstone::DeviceVector<RealType> d_gamma(h_gamma.data(), h_gamma.data() + nodeCount);
    cstone::DeviceVector<RealType> d_phi(h_phi.data(), h_phi.data() + nodeCount);
    cstone::DeviceVector<RealType> d_beta(h_beta.data(), h_beta.data() + nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_x(h_grad_phi_x.data(), h_grad_phi_x.data() + nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_y(h_grad_phi_y.data(), h_grad_phi_y.data() + nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_z(h_grad_phi_z.data(), h_grad_phi_z.data() + nodeCount);

    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 0.0);
    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12, 1.0);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12, 0.0);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12, 0.0);

    cudaDeviceSynchronize();
    auto fieldInitEnd = std::chrono::high_resolution_clock::now();
    perf.fieldInitTime = std::chrono::duration<float, std::milli>(fieldInitEnd - fieldInitStart).count();

    // Calculate data sizes for bandwidth calculation
    perf.coordinateBytes = 3 * nodeCount * sizeof(RealType);
    perf.inputFieldBytes = 6 * nodeCount * sizeof(RealType);
    perf.connectivityBytes = 8 * elementCount * sizeof(KeyType);
    perf.outputMatrixBytes = nnz * sizeof(RealType) + numOwnedDofs * sizeof(RealType);

    // Create CSR matrix wrapper
    using MatrixType = CSRMatrix<RealType>;

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

    // Configure assembler
    Assembler::Config config;
    config.blockSize = blockSize;
    config.variant = kernelVariant;

    int numBlocks = (elementCount + blockSize - 1) / blockSize;

    if (rank == 0 && !quiet) {
        std::cout << "\nAssembly kernel: " << numBlocks << " blocks x " << blockSize << " threads\n";
    }

    // Lambda to run assembly
    auto runAssembly = [&]() {
        Assembler::assembleGraphLump(
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
            d_rhs.data(),
            config
        );
    };

    // Lambda to reset matrix and RHS
    auto resetSystem = [&]() {
        thrust::fill(thrust::device_pointer_cast(d_values.data()),
                     thrust::device_pointer_cast(d_values.data() + nnz), 0.0);
        thrust::fill(thrust::device_pointer_cast(d_rhs.data()),
                     thrust::device_pointer_cast(d_rhs.data() + numOwnedDofs), 0.0);
    };

    // Warmup run
    timer.start();
    runAssembly();
    timer.stop();
    perf.assemblyWarmupTime = timer.elapsedMs();

    if (rank == 0 && !quiet) {
        std::cout << "  Warmup: " << perf.assemblyWarmupTime << " ms\n";
    }

    // Timed iterations
    for (int iter = 0; iter < numIterations; ++iter) {
        resetSystem();
        cudaDeviceSynchronize();

        timer.start();
        runAssembly();
        timer.stop();

        perf.recordIteration(timer.elapsedMs());
    }

    // Show per-iteration timing if verbose
    if (rank == 0 && !quiet && numIterations <= 20) {
        std::cout << "  Per-iteration times: ";
        for (int i = 0; i < numIterations; ++i) {
            std::cout << std::fixed << std::setprecision(2) << perf.iterTimes[i];
            if (i < numIterations - 1) std::cout << ", ";
        }
        std::cout << " ms\n";
    }

    // Compute matrix norm
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

    RealType global_matrix_norm2, global_rhs_norm2;
    MPI_Allreduce(&local_matrix_norm2, &global_matrix_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_rhs_norm2, &global_rhs_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    RealType matrix_norm = std::sqrt(global_matrix_norm2);
    RealType rhs_norm = std::sqrt(global_rhs_norm2);

    // Print STK-compatible single-line summary
    if (rank == 0) {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Assembler: MARS CVFEM Hex (Graph+Lump " << kernelName << ")....: "
                  << perf.avgAssemblyTimeMs() << " milliseconds @ "
                  << perf.bandwidthGBs() << " GB/s (average of " << numIterations
                  << " samples) [matrix norm: " << matrix_norm
                  << ", rhs norm: " << rhs_norm << "]" << std::endl;
    }

    // Print detailed performance summary
    if (!quiet) {
        perf.print(rank);
    }

    cudaFree(d_matrix);
    MPI_Finalize();
    return 0;
}
