#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include <thrust/device_vector.h>
#include <set>
#include <thrust/reduce.h>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <numeric>

// Performance counter structure for detailed timing
struct PerfCounters {
    // Timing (in milliseconds)
    float meshLoadTime = 0.0f;
    float sparsityBuildTime = 0.0f;
    float fieldInitTime = 0.0f;
    float assemblyWarmupTime = 0.0f;
    float assemblyTime = 0.0f;       // Total time for all iterations
    float assemblyMinTime = 1e9f;
    float assemblyMaxTime = 0.0f;
    int   assemblyIterations = 0;

    // Data sizes (in bytes)
    size_t inputFieldBytes = 0;      // Field data read per element
    size_t outputMatrixBytes = 0;    // Matrix/RHS data written
    size_t connectivityBytes = 0;    // Connectivity data read
    size_t coordinateBytes = 0;      // Coordinate data read

    // Counts
    size_t numElements = 0;
    size_t numNodes = 0;
    size_t numDofs = 0;
    size_t matrixNnz = 0;

    // Computed metrics
    double avgAssemblyTimeMs() const {
        return assemblyIterations > 0 ? assemblyTime / assemblyIterations : 0.0;
    }

    size_t totalBytesPerAssembly() const {
        return inputFieldBytes + outputMatrixBytes + connectivityBytes + coordinateBytes;
    }

    double bandwidthGBs() const {
        double timeS = avgAssemblyTimeMs() / 1000.0;
        return timeS > 0 ? (totalBytesPerAssembly() / 1e9) / timeS : 0.0;
    }

    double elementsPerSecond() const {
        double timeS = avgAssemblyTimeMs() / 1000.0;
        return timeS > 0 ? numElements / timeS : 0.0;
    }

    // FLOP estimate for CVFEM hex assembly (approximate)
    // Per element: 12 SCS * (shape functions + Jacobian + diffusion) ~ 2000 FLOPs
    double estimatedGFLOPs() const {
        const double flopsPerElement = 2000.0;
        double timeS = avgAssemblyTimeMs() / 1000.0;
        return timeS > 0 ? (numElements * flopsPerElement / 1e9) / timeS : 0.0;
    }

    void print(int rank) const {
        if (rank != 0) return;

        std::cout << "\n";
        std::cout << "================================================================================\n";
        std::cout << "                         PERFORMANCE SUMMARY\n";
        std::cout << "================================================================================\n";
        std::cout << std::fixed << std::setprecision(3);

        std::cout << "\n--- Problem Size ---\n";
        std::cout << "  Elements:              " << std::setw(12) << numElements << "\n";
        std::cout << "  Nodes:                 " << std::setw(12) << numNodes << "\n";
        std::cout << "  DOFs:                  " << std::setw(12) << numDofs << "\n";
        std::cout << "  Matrix NNZ:            " << std::setw(12) << matrixNnz << "\n";
        std::cout << "  Avg NNZ/row:           " << std::setw(12) << std::setprecision(1)
                  << (numDofs > 0 ? static_cast<double>(matrixNnz) / numDofs : 0.0) << "\n";

        std::cout << "\n--- Data Transfer (per assembly) ---\n";
        std::cout << std::setprecision(2);
        std::cout << "  Input fields:          " << std::setw(12) << inputFieldBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Coordinates:           " << std::setw(12) << coordinateBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Connectivity:          " << std::setw(12) << connectivityBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Output (matrix+RHS):   " << std::setw(12) << outputMatrixBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Total:                 " << std::setw(12) << totalBytesPerAssembly() / 1024.0 / 1024.0 << " MB\n";

        std::cout << "\n--- Timing Breakdown ---\n";
        std::cout << std::setprecision(3);
        std::cout << "  Mesh loading:          " << std::setw(12) << meshLoadTime << " ms\n";
        std::cout << "  Sparsity build:        " << std::setw(12) << sparsityBuildTime << " ms\n";
        std::cout << "  Field initialization:  " << std::setw(12) << fieldInitTime << " ms\n";
        std::cout << "  Assembly warmup:       " << std::setw(12) << assemblyWarmupTime << " ms\n";

        std::cout << "\n--- Assembly Performance (" << assemblyIterations << " iterations) ---\n";
        std::cout << "  Average time:          " << std::setw(12) << avgAssemblyTimeMs() << " ms\n";
        std::cout << "  Min time:              " << std::setw(12) << assemblyMinTime << " ms\n";
        std::cout << "  Max time:              " << std::setw(12) << assemblyMaxTime << " ms\n";
        std::cout << "  Bandwidth:             " << std::setw(12) << bandwidthGBs() << " GB/s\n";
        std::cout << "  Throughput:            " << std::setw(12) << elementsPerSecond() / 1e6 << " M elem/s\n";
        std::cout << "  Est. compute:          " << std::setw(12) << estimatedGFLOPs() << " GFLOP/s\n";

        std::cout << "================================================================================\n";
    }
};

// CUDA event-based GPU timer
class GpuTimer {
public:
    GpuTimer() {
        cudaEventCreate(&start_);
        cudaEventCreate(&stop_);
    }

    ~GpuTimer() {
        cudaEventDestroy(start_);
        cudaEventDestroy(stop_);
    }

    void start(cudaStream_t stream = 0) {
        cudaEventRecord(start_, stream);
    }

    void stop(cudaStream_t stream = 0) {
        cudaEventRecord(stop_, stream);
    }

    float elapsedMs() {
        cudaEventSynchronize(stop_);
        float ms = 0.0f;
        cudaEventElapsedTime(&ms, start_, stop_);
        return ms;
    }

private:
    cudaEvent_t start_, stop_;
};

using namespace mars;
using namespace mars::fem;

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // Set CUDA device based on local MPI rank (for multi-GPU nodes)
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) {
        int device = rank % deviceCount;
        cudaSetDevice(device);
    }

    // Parse command-line options
    std::string meshFile;
    int keyBits = 64;
    bool useExactCoords = false;
    int numIterations = 10;
    int blockSize = 256;
    bool quiet = false;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Original;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0) {
            meshFile = arg.substr(7);
        } else if (arg.find("--key-bits=") == 0) {
            keyBits = std::stoi(arg.substr(11));
            if (keyBits != 32 && keyBits != 64) {
                if (rank == 0) {
                    std::cerr << "Error: --key-bits must be 32 or 64" << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "--exact-coords") {
            useExactCoords = true;
        } else if (arg.find("--iterations=") == 0) {
            numIterations = std::stoi(arg.substr(13));
        } else if (arg.find("--block-size=") == 0) {
            blockSize = std::stoi(arg.substr(13));
        } else if (arg.find("--kernel=") == 0) {
            std::string v = arg.substr(9);
            if (v == "original") kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "tensor") kernelVariant = CvfemKernelVariant::Tensor;
            else {
                if (rank == 0) {
                    std::cerr << "Error: Unknown kernel variant: " << v << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg[0] != '-' && meshFile.empty()) {
            // Positional argument (backward compatibility)
            meshFile = arg;
        }
    }

    if (meshFile.empty()) {
        if (rank == 0) {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "\nOptions:" << std::endl;
            std::cout << "  --mesh=FILE         Mesh file (.mesh or .exo format) [REQUIRED]" << std::endl;
            std::cout << "  --kernel=VARIANT    tensor, shmem, optimized, original (default: original)" << std::endl;
            std::cout << "  --block-size=N      CUDA block size (default: 256)" << std::endl;
            std::cout << "  --iterations=N      Number of assembly iterations (default: 10)" << std::endl;
            std::cout << "  --key-bits=32|64    SFC key size (default: 64)" << std::endl;
            std::cout << "  --exact-coords      Store exact original coordinates" << std::endl;
            std::cout << "  --quiet             Suppress detailed debug output" << std::endl;
            std::cout << "\nBackward compatible positional syntax:" << std::endl;
            std::cout << "  " << argv[0] << " <mesh_file> [options]" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Initialize performance counters
    PerfCounters perf;
    GpuTimer timer;

    if (rank == 0 && !quiet) {
        std::cout << "CVFEM Full Assembly (27 NNZ/row sparsity pattern)" << std::endl;
        std::cout << "Configuration:" << std::endl;
        std::cout << "  SFC key bits:    " << keyBits << std::endl;
        std::cout << "  Coord mode:      " << (useExactCoords ? "exact" : "SFC-decoded") << std::endl;
        std::cout << "  Iterations:      " << numIterations << std::endl;
        std::cout << "  Block size:      " << blockSize << std::endl;
        std::cout << "  Kernel variant:  ";
        switch (kernelVariant) {
            case CvfemKernelVariant::Original: std::cout << "original"; break;
            case CvfemKernelVariant::Optimized: std::cout << "optimized"; break;
            case CvfemKernelVariant::Shmem: std::cout << "shmem"; break;
            case CvfemKernelVariant::Tensor: std::cout << "tensor"; break;
            default: std::cout << "unknown"; break;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

    // Use 64-bit SFC keys by default for better coordinate precision
#ifdef USE_32BIT_KEYS
    using KeyType = unsigned int;
#else
    using KeyType = uint64_t;
#endif
    using RealType = double;
    using ElemTag = HexTag;

    // Create domain with optional exact coordinate storage
    auto meshLoadStart = std::chrono::high_resolution_clock::now();
    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(meshFile, rank, numRanks, useExactCoords);

    // Build halo for node ownership (this triggers SFC map creation which updates node count)
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    cudaDeviceSynchronize();
    auto meshLoadEnd = std::chrono::high_resolution_clock::now();
    perf.meshLoadTime = std::chrono::duration<float, std::milli>(meshLoadEnd - meshLoadStart).count();

    // Get counts AFTER node ownership map is built
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();

    // Get element range for local (owned) elements
    auto [startIdx, endIdx] = domain.localElementRange();
    size_t localElementCount = endIdx - startIdx;

    // Store counts in perf counters
    perf.numNodes = nodeCount;
    perf.numElements = elementCount;

    if (rank == 0 && !quiet) {
        std::cout << "Mesh loaded in " << perf.meshLoadTime << " ms" << std::endl;
        std::cout << "  Total nodes:    " << nodeCount << std::endl;
        std::cout << "  Total elements: " << elementCount << " (local: " << localElementCount << ")" << std::endl;
        std::cout << std::endl;
    }

    // Count owned DOFs (owned + shared nodes)
    std::vector<uint8_t> h_ownership(nodeCount);
    thrust::copy(thrust::device_pointer_cast(d_nodeOwnership.data()),
                 thrust::device_pointer_cast(d_nodeOwnership.data() + nodeCount),
                 h_ownership.begin());

    int numOwnedDofs = 0;
    for (auto o : h_ownership) {
        if (o == 1 || o == 2) numOwnedDofs++;
    }

    // Debug: Print node ownership statistics
    if (!quiet) {
        int num_ghost = 0, num_owned = 0, num_shared = 0;
        for (auto o : h_ownership) {
            if (o == 0) num_ghost++;
            else if (o == 1) num_owned++;
            else if (o == 2) num_shared++;
        }
        if (rank == 0) {
            std::cout << "Node ownership: ghost=" << num_ghost
                      << " owned=" << num_owned
                      << " shared=" << num_shared << std::endl;
        }
    }

    // Create DOF mapping
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

    // Build FULL sparsity pattern (all 8x8 = 64 entries per element, ~27 NNZ/row)
    auto sparsityStart = std::chrono::high_resolution_clock::now();
    if (rank == 0 && !quiet) {
        std::cout << "Building full sparsity pattern (27 NNZ/row)..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Get element connectivity
    const auto& d_conn = domain.getElementToNodeConnectivity();

    // Copy connectivity to host
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

    // Build full adjacency (all node pairs in each element)
    std::vector<std::set<int>> adj(numOwnedDofs);

    for (size_t e = 0; e < elementCount; ++e) {
        KeyType nodes[8];
        for (int i = 0; i < 8; ++i) {
            nodes[i] = h_conn[i][e];
        }

        // Full 8x8 coupling
        for (int i = 0; i < 8; ++i) {
            int row_dof = h_nodeToDof[nodes[i]];
            if (row_dof < 0) continue;

            for (int j = 0; j < 8; ++j) {
                int col_dof = h_nodeToDof[nodes[j]];
                if (col_dof >= 0) {
                    adj[row_dof].insert(col_dof);
                }
            }
        }
    }

    // Convert to CSR format
    std::vector<int> rowPtr(numOwnedDofs + 1);
    std::vector<int> colInd;
    std::vector<int> diagPtr(numOwnedDofs);  // Diagonal pointers for fast diagonal access

    rowPtr[0] = 0;
    for (int d = 0; d < numOwnedDofs; ++d) {
        int diagFound = -1;
        for (int col : adj[d]) {
            if (col == d) {
                diagFound = colInd.size();  // Record diagonal position
            }
            colInd.push_back(col);
        }
        diagPtr[d] = diagFound;  // Store diagonal position for this row
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

    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate matrix and RHS on device
    cstone::DeviceVector<int> d_rowPtr(rowPtr.data(), rowPtr.data() + rowPtr.size());
    cstone::DeviceVector<int> d_colInd(colInd.data(), colInd.data() + colInd.size());
    cstone::DeviceVector<int> d_diagPtr(diagPtr.data(), diagPtr.data() + diagPtr.size());
    cstone::DeviceVector<RealType> d_values(nnz, 0.0);
    cstone::DeviceVector<RealType> d_rhs(numOwnedDofs, 0.0);

    // Initialize fields to match STK reference
    auto fieldInitStart = std::chrono::high_resolution_clock::now();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    // Copy coordinates to host for field initialization
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

    // Initialize fields using STK reference functions
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

    // Copy initialized fields to device
    cstone::DeviceVector<RealType> d_gamma(h_gamma.data(), h_gamma.data() + nodeCount);
    cstone::DeviceVector<RealType> d_phi(h_phi.data(), h_phi.data() + nodeCount);
    cstone::DeviceVector<RealType> d_beta(h_beta.data(), h_beta.data() + nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_x(h_grad_phi_x.data(), h_grad_phi_x.data() + nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_y(h_grad_phi_y.data(), h_grad_phi_y.data() + nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_z(h_grad_phi_z.data(), h_grad_phi_z.data() + nodeCount);

    // Element-based data (12 SCS per hex element)
    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 0.0);
    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12);

    // Precompute area vectors for all elements and SCS integration points
    precomputeAreaVectorsGpu<KeyType, RealType>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data()
    );

    cudaDeviceSynchronize();
    auto fieldInitEnd = std::chrono::high_resolution_clock::now();
    perf.fieldInitTime = std::chrono::duration<float, std::milli>(fieldInitEnd - fieldInitStart).count();

    // Calculate data sizes for bandwidth calculation
    perf.coordinateBytes = 3 * nodeCount * sizeof(RealType);
    perf.inputFieldBytes = 6 * nodeCount * sizeof(RealType);
    perf.connectivityBytes = 8 * elementCount * sizeof(KeyType);
    perf.outputMatrixBytes = nnz * sizeof(RealType) + numOwnedDofs * sizeof(RealType);

    // Create CSR wrapper with diagPtr for tensor kernel compatibility
    using MatrixType = mars::fem::CSRMatrix<RealType>;

    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{
        d_rowPtr.data(),
        d_colInd.data(),
        d_values.data(),
        d_diagPtr.data(),  // Diagonal pointers for fast diagonal access
        numOwnedDofs,
        nnz
    };
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    // Launch assembly kernel with timing
    int numBlocks = (elementCount + blockSize - 1) / blockSize;

    if (rank == 0 && !quiet) {
        std::cout << "\nAssembly kernel: " << numBlocks << " blocks x " << blockSize << " threads\n";
    }

    // Lambda to run the assembly kernel
    auto runAssembly = [&]() {
        CvfemHexAssembler<KeyType, RealType>::Config config;
        config.blockSize = blockSize;
        config.variant = kernelVariant;

        CvfemHexAssembler<KeyType, RealType>::assembleFull(
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
    perf.assemblyIterations = numIterations;
    std::vector<float> iterTimes(numIterations);

    for (int iter = 0; iter < numIterations; ++iter) {
        resetSystem();
        cudaDeviceSynchronize();

        timer.start();
        runAssembly();
        timer.stop();

        float iterTime = timer.elapsedMs();
        iterTimes[iter] = iterTime;
        perf.assemblyTime += iterTime;
        perf.assemblyMinTime = std::min(perf.assemblyMinTime, iterTime);
        perf.assemblyMaxTime = std::max(perf.assemblyMaxTime, iterTime);
    }

    // Show per-iteration timing if verbose
    if (rank == 0 && !quiet && numIterations <= 20) {
        std::cout << "  Per-iteration times: ";
        for (int i = 0; i < numIterations; ++i) {
            std::cout << std::fixed << std::setprecision(2) << iterTimes[i];
            if (i < numIterations - 1) std::cout << ", ";
        }
        std::cout << " ms\n";
    }

    // Compute matrix norm (L2/Frobenius)
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
        std::cout << "Assembler: MARS CVFEM Hex (Full 27-NNZ)....: "
                  << perf.avgAssemblyTimeMs() << " milliseconds @ "
                  << perf.bandwidthGBs() << " GB/s (average of " << numIterations
                  << " samples) [matrix norm: " << matrix_norm
                  << ", rhs norm: " << rhs_norm << "]" << std::endl;
    }

    // Print detailed performance summary
    perf.print(rank);

    cudaFree(d_matrix);
    MPI_Finalize();
    return 0;
}
