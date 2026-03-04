#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/system/cuda/execution_policy.h>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#ifndef MARS_DISABLE_NVTX
#include <nvToolsExt.h>
#define MARS_NVTX_PUSH(name) nvtxRangePush(name)
#define MARS_NVTX_POP()      nvtxRangePop()
#else
#define MARS_NVTX_PUSH(name)
#define MARS_NVTX_POP()
#endif

using namespace mars;
using namespace mars::fem;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // Parse command-line options
    std::string meshFile;
    int numIterations = 10;
    int blockSize = 256;
    int bucketSize = 64;  // Cornerstone octree bucket size (smaller = finer partitioning)
    bool quiet = false;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Original;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0) {
            meshFile = arg.substr(7);
        } else if (arg.find("--iterations=") == 0) {
            numIterations = std::stoi(arg.substr(13));
        } else if (arg.find("--block-size=") == 0) {
            blockSize = std::stoi(arg.substr(13));
        } else if (arg.find("--bucket-size=") == 0) {
            bucketSize = std::stoi(arg.substr(14));
        } else if (arg.find("--kernel=") == 0) {
            std::string v = arg.substr(9);
            if (v == "original") kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "team") kernelVariant = CvfemKernelVariant::Team;
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
            std::cout << "Usage: " << argv[0] << " [options]\n";
            std::cout << "\nOptions:\n";
            std::cout << "  --mesh=FILE         Mesh file (.mesh or .exo format) [REQUIRED]\n";
            std::cout << "  --kernel=VARIANT    tensor, shmem, optimized, team, original (default: original)\n";
            std::cout << "  --block-size=N      CUDA block size (default: 256)\n";
            std::cout << "  --bucket-size=N     Cornerstone octree bucket size (default: 64, try 32 or 16 for 16+ ranks)\n";
            std::cout << "  --iterations=N      Number of assembly iterations (default: 10)\n";
            std::cout << "  --quiet             Suppress detailed output\n";
            std::cout << "\nBackward compatible positional syntax:\n";
            std::cout << "  " << argv[0] << " <mesh_file> [options]\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType = uint64_t;
    using RealType = double;
    using ElemTag = HexTag;
    using Assembler = CvfemHexAssembler<KeyType, RealType>;

    PerfCounters perf;
    GpuTimer timer;

    const char* kernelName = Assembler::variantName(kernelVariant);

    if (rank == 0) {
        std::cout << "Command-line parameters:" << std::endl;
        std::cout << "  Mesh: " << meshFile << std::endl;
        std::cout << "  Kernel: " << kernelName << std::endl;
        std::cout << "  Block size: " << blockSize << std::endl;
        std::cout << "  Bucket size: " << bucketSize << std::endl;
        std::cout << "  Iterations: " << numIterations << std::endl;
    }

    // Load mesh with custom bucket size
    MARS_NVTX_PUSH("Mesh Loading");
    auto meshLoadStart = std::chrono::high_resolution_clock::now();
    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(meshFile, rank, numRanks, true, bucketSize);
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    cudaDeviceSynchronize();
    auto meshLoadEnd = std::chrono::high_resolution_clock::now();
    perf.meshLoadTime = std::chrono::duration<float, std::milli>(meshLoadEnd - meshLoadStart).count();
    MARS_NVTX_POP();

    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();

    // Build DOF mapping on GPU
    // All nodes get DOFs: owned → [0, numDofs), ghost → [numDofs, nodeCount)
    // Ghost DOFs are needed as CSR column indices for correct off-diagonal assembly.
    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs = buildDofMappingGpu<KeyType>(
        d_nodeOwnership.data(),
        d_nodeToDof.data(),
        nodeCount
    );
    int numTotalDofs = static_cast<int>(nodeCount);    // owned + ghost (for CSR sizing)
    perf.numDofs    = numDofs;
    perf.numNodes   = numDofs;                         // owned nodes (DOF = node for scalar problem)
    perf.numElements = domain.localElementCount();

    // Build sparsity pattern on GPU
    auto sparsityStart = std::chrono::high_resolution_clock::now();

    const auto& d_conn = domain.getElementToNodeConnectivity();

    // Allocate CSR structure with diagonal positions (sized for all DOFs: owned + ghost)
    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1);
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr(numTotalDofs);

    // Build sparsity pattern entirely on GPU (using numTotalDofs for rows to include ghost columns)
    int nnz = CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        elementCount,
        d_nodeToDof.data(),
        numTotalDofs,
        d_rowPtr.data(),
        nullptr,  // First pass to get nnz
        nullptr,
        0
    );

    // Allocate column indices and rebuild with diagonal positions
    d_colInd.resize(nnz);
    CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        elementCount,
        d_nodeToDof.data(),
        numTotalDofs,
        d_rowPtr.data(),
        d_colInd.data(),
        d_diagPtr.data(),
        0
    );

    cudaDeviceSynchronize();
    auto sparsityEnd = std::chrono::high_resolution_clock::now();
    perf.sparsityBuildTime = std::chrono::duration<float, std::milli>(sparsityEnd - sparsityStart).count();
    perf.matrixNnz = nnz;
    cstone::DeviceVector<RealType> d_values(nnz, 0.0);
    cstone::DeviceVector<RealType> d_rhs(numTotalDofs, 0.0);

    // Initialize fields on GPU (no CPU copies!)
    auto fieldInitStart = std::chrono::high_resolution_clock::now();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    cstone::DeviceVector<RealType> d_phi(nodeCount);
    cstone::DeviceVector<RealType> d_gamma(nodeCount);
    cstone::DeviceVector<RealType> d_beta(nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount);
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount);

    // GPU kernel for field initialization
    initFieldsGpu<RealType>(
        d_x.data(), d_y.data(), d_z.data(),
        d_phi.data(), d_gamma.data(), d_beta.data(),
        d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        nodeCount
    );

    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 0.0);
    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12);

    // Pre-compute area vectors from geometry (required for correct diffusion)
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

    // Data sizes for bandwidth (owned counts only, matching STK convention:
    // total_bytes x p = constant for a fixed problem, ghost overhead shows as efficiency loss)
    perf.coordinateBytes   = 3 * numDofs * sizeof(RealType);
    perf.inputFieldBytes   = 5 * numDofs * sizeof(RealType);
    perf.connectivityBytes = 8 * domain.localElementCount() * sizeof(int);
    perf.outputMatrixBytes = nnz * sizeof(RealType) + numDofs * sizeof(RealType);

    // Create CSR matrix wrapper with diagonal positions
    // Create CUDA streams for overlapping operations
    cudaStream_t assemblyStream, resetStream;
    cudaStreamCreate(&assemblyStream);
    cudaStreamCreate(&resetStream);

    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(), numTotalDofs, nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    // Configure assembler
    Assembler::Config config;
    config.blockSize = blockSize;
    config.variant = kernelVariant;
    config.stream = assemblyStream;  // Use dedicated stream

    auto runAssembly = [&]() {
        Assembler::assembleGraphLump(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            elementCount,
            d_x.data(), d_y.data(), d_z.data(),
            d_gamma.data(), d_phi.data(), d_beta.data(),
            d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
            d_mdot.data(), d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data(),
            d_nodeToDof.data(), d_nodeOwnership.data(),
            d_matrix, d_rhs.data(),
            config
        );
    };

    auto resetSystem = [&]() {
        // Use resetStream to potentially overlap with other operations
        thrust::fill(thrust::cuda::par.on(resetStream),
                     thrust::device_pointer_cast(d_values.data()),
                     thrust::device_pointer_cast(d_values.data() + nnz), 0.0);
        thrust::fill(thrust::cuda::par.on(resetStream),
                     thrust::device_pointer_cast(d_rhs.data()),
                     thrust::device_pointer_cast(d_rhs.data() + numDofs), 0.0);
    };

    // Warmup
    timer.start();
    runAssembly();
    cudaStreamSynchronize(assemblyStream);  // Only sync assembly stream
    timer.stop();
    perf.assemblyWarmupTime = timer.elapsedMs();

    // Timed iterations with load balance tracking
    std::vector<double> perIterationTimes(numIterations);
    for (int iter = 0; iter < numIterations; ++iter) {
        MARS_NVTX_PUSH("Assembly Iteration");
        // Reset system on dedicated stream (can overlap with CPU work)
        MARS_NVTX_PUSH("Reset System");
        resetSystem();
        cudaStreamSynchronize(resetStream);  // Ensure reset completes before assembly
        MARS_NVTX_POP();
        
        MARS_NVTX_PUSH("Assembly Kernel");
        timer.start();
        runAssembly();
        cudaStreamSynchronize(assemblyStream);  // Only sync assembly stream
        timer.stop();
        MARS_NVTX_POP();
        double iterTime = timer.elapsedMs();
        perf.recordIteration(iterTime);
        perIterationTimes[iter] = iterTime;
        MARS_NVTX_POP();
    }

    // Compute DD load imbalance from node counts (owned / ghost / total)
    MARS_NVTX_PUSH("Post-processing");
    size_t ghostNodeCount = nodeCount - static_cast<size_t>(numDofs);
    double sendBuf[3] = {
        static_cast<double>(numDofs),
        static_cast<double>(ghostNodeCount),
        static_cast<double>(nodeCount)
    };
    double recvMin[3], recvMax[3], recvSum[3];
    MPI_Request req[3];
    MARS_NVTX_PUSH("MPI DD Imbalance");
    MPI_Iallreduce(sendBuf, recvMin, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &req[0]);
    MPI_Iallreduce(sendBuf, recvMax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &req[1]);
    MPI_Iallreduce(sendBuf, recvSum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req[2]);
    MARS_NVTX_POP();

    // Compute matrix norm (overlapped with load balance computation)
    MARS_NVTX_PUSH("Norm Computation");
    auto sq = [] __host__ __device__ (RealType x) { return x * x; };
    RealType norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(d_values.data()),
        thrust::device_pointer_cast(d_values.data() + nnz),
        sq, RealType(0.0), thrust::plus<RealType>()
    );
    RealType rhs_norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(d_rhs.data()),
        thrust::device_pointer_cast(d_rhs.data() + numDofs),
        sq, RealType(0.0), thrust::plus<RealType>()
    );
    MARS_NVTX_POP();

    // Start non-blocking norm reductions
    MARS_NVTX_PUSH("MPI Norm Reduction");
    RealType global_norm2, global_rhs_norm2;
    MPI_Request normReq[2];
    MPI_Iallreduce(&norm2, &global_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &normReq[0]);
    MPI_Iallreduce(&rhs_norm2, &global_rhs_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &normReq[1]);

    // Wait for DD imbalance reductions
    MPI_Waitall(3, req, MPI_STATUSES_IGNORE);
    auto ddImbalance = [&](int k, int normK) -> double {
        double mean = recvSum[normK] / numRanks;
        return mean > 0.0 ? (recvMax[k] - recvMin[k]) / mean * 100.0 : 0.0;
    };
    perf.ddOwnedImbalancePct = ddImbalance(0, 0);  // owned spread / owned mean
    perf.ddGhostImbalancePct = ddImbalance(1, 2);  // ghost spread / total mean (avoids inflated % from small denominator)
    perf.ddTotalImbalancePct = ddImbalance(2, 2);  // total spread / total mean
    perf.ddOwnedMin  = recvMin[0]; perf.ddOwnedMax  = recvMax[0]; perf.ddOwnedMean  = recvSum[0] / numRanks;
    perf.ddGhostMin  = recvMin[1]; perf.ddGhostMax  = recvMax[1]; perf.ddGhostMean  = recvSum[1] / numRanks;
    perf.ddTotalMin  = recvMin[2]; perf.ddTotalMax  = recvMax[2]; perf.ddTotalMean  = recvSum[2] / numRanks;

    // Wait for norm reductions to complete
    MPI_Waitall(2, normReq, MPI_STATUSES_IGNORE);
    MARS_NVTX_POP();
    RealType matrix_norm = std::sqrt(global_norm2);
    RealType rhs_norm = std::sqrt(global_rhs_norm2);
    MARS_NVTX_POP();

    // Print summary
    if (rank == 0) {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Assembler: MARS CVFEM Hex (Graph+Lump " << kernelName << ")....: "
                  << perf.avgAssemblyTimeMs() << " milliseconds @ "
                  << perf.bandwidthGBs() << " GB/s (average of " << numIterations
                  << " samples) [matrix norm: " << matrix_norm
                  << ", rhs norm: " << rhs_norm << "]\n";

    }

    if (!quiet) {
        perf.print(rank);
    }

    // Cleanup streams
    cudaStreamDestroy(assemblyStream);
    cudaStreamDestroy(resetStream);

    cudaFree(d_matrix);
    MPI_Finalize();
    return 0;
}
