// CVFEM graph-sparsity assembly on a distributed tetrahedral mesh.
//
// Counterpart of mars_cvfem_graph.cu (which is hex-only). Same pipeline minus
// the kernel-variant zoo and the precomputed area-vector path: the tet graph
// kernel derives the SCS area vector inline from the edge vector.
#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_tet_assembler.hpp"
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

    std::string meshFile;
    int numIterations    = 10;
    int blockSize        = 256;
    int bucketSize       = 64;
    int bucketSizeFocus  = 8;
    bool quiet           = false;

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
        } else if (arg.find("--bucket-size-focus=") == 0) {
            bucketSizeFocus = std::stoi(arg.substr(20));
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg == "--help" || arg == "-h") {
            if (rank == 0) {
                std::cout << "Usage: " << argv[0] << " --mesh=<path> [options]\n"
                          << "  --mesh=<path>          Tetrahedral binary mesh dir\n"
                          << "  --iterations=N         Timed iterations (default 10)\n"
                          << "  --block-size=N         CUDA block size (default 256)\n"
                          << "  --bucket-size=N        cstone bucket size (default 64)\n"
                          << "  --bucket-size-focus=N  cstone focus-tree bucket size (default 8)\n"
                          << "  --quiet                Skip the perf breakdown printout\n";
            }
            MPI_Finalize();
            return 0;
        }
    }

    if (meshFile.empty()) {
        if (rank == 0) std::cerr << "Error: --mesh=<path> is required\n";
        MPI_Finalize();
        return 1;
    }

    using KeyType   = uint64_t;
    using RealType  = double;
    using ElemTag   = TetTag;
    using Assembler = CvfemTetAssembler<KeyType, RealType>;

    PerfCounters perf;
    GpuTimer timer;

    if (rank == 0) {
        std::cout << "Command-line parameters:\n"
                  << "  Mesh: "                << meshFile         << "\n"
                  << "  Kernel: "              << Assembler::variantName() << "\n"
                  << "  Block size: "          << blockSize        << "\n"
                  << "  Bucket size: "         << bucketSize       << "\n"
                  << "  Bucket size (focus): " << bucketSizeFocus  << "\n"
                  << "  Iterations: "          << numIterations    << "\n";
    }

    MARS_NVTX_PUSH("Mesh Loading");
    auto meshLoadStart = std::chrono::high_resolution_clock::now();
    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(
        meshFile, rank, numRanks, true, bucketSize, static_cast<unsigned>(bucketSizeFocus));
    if (rank == 0) { std::cout << "PHASE: domain constructed" << std::endl; std::cout.flush(); }

    auto haloStart = std::chrono::high_resolution_clock::now();
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    cudaDeviceSynchronize();
    auto haloEnd = std::chrono::high_resolution_clock::now();
    perf.nodeHaloTime = std::chrono::duration<float, std::milli>(haloEnd - haloStart).count();
    if (rank == 0) { std::cout << "PHASE: domain synced" << std::endl; std::cout.flush(); }
    auto meshLoadEnd = std::chrono::high_resolution_clock::now();
    perf.meshLoadTime = std::chrono::duration<float, std::milli>(meshLoadEnd - meshLoadStart).count();

    perf.readMeshTime = domain.readMeshTimeMs;
    perf.bboxTime     = domain.bboxTimeMs;
    perf.h2dTime      = domain.h2dTimeMs;
    perf.charSizeTime = domain.charSizeTimeMs;
    perf.syncTime     = domain.syncTimeMs;
    MARS_NVTX_POP();

    size_t nodeCount    = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    if (rank == 0) { std::cout << "PHASE: nodeCount=" << nodeCount
                               << " elementCount=" << elementCount << std::endl; std::cout.flush(); }

    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs = buildDofMappingGpu<KeyType>(
        d_nodeOwnership.data(), d_nodeToDof.data(), nodeCount);
    int numTotalDofs = static_cast<int>(nodeCount);
    if (rank == 0) { std::cout << "PHASE: DOF mapping done, numDofs=" << numDofs
                               << " numTotalDofs=" << numTotalDofs << std::endl; std::cout.flush(); }
    perf.numDofs     = numDofs;
    perf.numNodes    = numDofs;
    perf.numElements = domain.localElementCount();

    auto sparsityStart = std::chrono::high_resolution_clock::now();
    const auto& d_conn = domain.getElementToNodeConnectivity();
    if (rank == 0) { std::cout << "PHASE: starting sparsity build" << std::endl; std::cout.flush(); }

    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1);
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr(numTotalDofs);

    int nnz = CvfemTetSparsityBuilder<KeyType>::buildGraphSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        elementCount, d_nodeToDof.data(), numTotalDofs,
        d_rowPtr.data(), nullptr, nullptr, 0);

    d_colInd.resize(nnz);
    CvfemTetSparsityBuilder<KeyType>::buildGraphSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        elementCount, d_nodeToDof.data(), numTotalDofs,
        d_rowPtr.data(), d_colInd.data(), d_diagPtr.data(), 0);

    cudaDeviceSynchronize();
    if (rank == 0) { std::cout << "PHASE: sparsity done, nnz=" << nnz << std::endl; std::cout.flush(); }
    auto sparsityEnd = std::chrono::high_resolution_clock::now();
    perf.sparsityBuildTime = std::chrono::duration<float, std::milli>(sparsityEnd - sparsityStart).count();
    perf.matrixNnz = nnz;

    cstone::DeviceVector<RealType> d_values(nnz, 0.0);
    cstone::DeviceVector<RealType> d_rhs(numTotalDofs, 0.0);

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

    initFieldsGpu<RealType>(
        d_x.data(), d_y.data(), d_z.data(),
        d_phi.data(), d_gamma.data(), d_beta.data(),
        d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        nodeCount);

    // mdot: tet has 6 SCS per element (vs 12 for hex).
    cstone::DeviceVector<RealType> d_mdot(elementCount * 6, 0.0);

    cudaDeviceSynchronize();
    auto fieldInitEnd = std::chrono::high_resolution_clock::now();
    perf.fieldInitTime = std::chrono::duration<float, std::milli>(fieldInitEnd - fieldInitStart).count();

    // Bandwidth accounting: same physical fields as hex, but 4 corners/element.
    perf.coordinateBytes   = 3 * numDofs * sizeof(RealType);
    perf.inputFieldBytes   = 5 * numDofs * sizeof(RealType);
    perf.connectivityBytes = 4 * domain.localElementCount() * sizeof(int);
    perf.outputMatrixBytes = nnz * sizeof(RealType) + numDofs * sizeof(RealType);

    cudaStream_t assemblyStream, resetStream;
    cudaStreamCreate(&assemblyStream);
    cudaStreamCreate(&resetStream);

    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(),
                        d_diagPtr.data(), numTotalDofs, nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    typename Assembler::Config config;
    config.blockSize = blockSize;
    config.stream    = assemblyStream;

    auto runAssembly = [&]() {
        Assembler::assembleGraphLump(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            elementCount,
            d_x.data(), d_y.data(), d_z.data(),
            d_gamma.data(), d_phi.data(), d_beta.data(),
            d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
            d_mdot.data(),
            d_nodeToDof.data(), d_nodeOwnership.data(),
            d_matrix, d_rhs.data(),
            config);
    };

    auto resetSystem = [&]() {
        thrust::fill(thrust::cuda::par.on(resetStream),
                     thrust::device_pointer_cast(d_values.data()),
                     thrust::device_pointer_cast(d_values.data() + nnz), 0.0);
        thrust::fill(thrust::cuda::par.on(resetStream),
                     thrust::device_pointer_cast(d_rhs.data()),
                     thrust::device_pointer_cast(d_rhs.data() + numDofs), 0.0);
    };

    timer.start();
    runAssembly();
    cudaStreamSynchronize(assemblyStream);
    timer.stop();
    perf.assemblyWarmupTime = timer.elapsedMs();

    std::vector<double> perIterationTimes(numIterations);
    for (int iter = 0; iter < numIterations; ++iter) {
        MARS_NVTX_PUSH("Assembly Iteration");
        MARS_NVTX_PUSH("Reset System");
        resetSystem();
        cudaStreamSynchronize(resetStream);
        MARS_NVTX_POP();

        MARS_NVTX_PUSH("Assembly Kernel");
        timer.start();
        runAssembly();
        cudaStreamSynchronize(assemblyStream);
        timer.stop();
        MARS_NVTX_POP();
        double iterTime = timer.elapsedMs();
        perf.recordIteration(iterTime);
        perIterationTimes[iter] = iterTime;
        if (rank == 0) std::cout << "iter " << iter << ": " << iterTime << " ms" << std::endl;
        MARS_NVTX_POP();
    }

    // DD imbalance + matrix/RHS norms (overlap reductions like the hex driver).
    MARS_NVTX_PUSH("Post-processing");
    size_t ghostNodeCount = nodeCount - static_cast<size_t>(numDofs);
    double sendBuf[3] = {
        static_cast<double>(numDofs),
        static_cast<double>(ghostNodeCount),
        static_cast<double>(nodeCount)
    };
    double recvMin[3], recvMax[3], recvSum[3];
    MPI_Request req[3];
    MPI_Iallreduce(sendBuf, recvMin, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &req[0]);
    MPI_Iallreduce(sendBuf, recvMax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &req[1]);
    MPI_Iallreduce(sendBuf, recvSum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req[2]);

    auto sq = [] __host__ __device__ (RealType x) { return x * x; };
    RealType norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(d_values.data()),
        thrust::device_pointer_cast(d_values.data() + nnz),
        sq, RealType(0.0), thrust::plus<RealType>());
    RealType rhs_norm2 = thrust::transform_reduce(
        thrust::device_pointer_cast(d_rhs.data()),
        thrust::device_pointer_cast(d_rhs.data() + numDofs),
        sq, RealType(0.0), thrust::plus<RealType>());

    RealType global_norm2, global_rhs_norm2;
    MPI_Request normReq[2];
    MPI_Iallreduce(&norm2,     &global_norm2,     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &normReq[0]);
    MPI_Iallreduce(&rhs_norm2, &global_rhs_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &normReq[1]);

    MPI_Waitall(3, req, MPI_STATUSES_IGNORE);
    auto ddImbalance = [&](int k, int normK) -> double {
        double mean = recvSum[normK] / numRanks;
        return mean > 0.0 ? (recvMax[k] - recvMin[k]) / mean * 100.0 : 0.0;
    };
    perf.ddOwnedImbalancePct = ddImbalance(0, 0);
    perf.ddGhostImbalancePct = ddImbalance(1, 2);
    perf.ddTotalImbalancePct = ddImbalance(2, 2);
    perf.ddOwnedMin = recvMin[0]; perf.ddOwnedMax = recvMax[0]; perf.ddOwnedMean = recvSum[0] / numRanks;
    perf.ddGhostMin = recvMin[1]; perf.ddGhostMax = recvMax[1]; perf.ddGhostMean = recvSum[1] / numRanks;
    perf.ddTotalMin = recvMin[2]; perf.ddTotalMax = recvMax[2]; perf.ddTotalMean = recvSum[2] / numRanks;

    MPI_Waitall(2, normReq, MPI_STATUSES_IGNORE);
    MARS_NVTX_POP();
    RealType matrix_norm = std::sqrt(global_norm2);
    RealType rhs_norm    = std::sqrt(global_rhs_norm2);

    if (rank == 0) {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Assembler: MARS CVFEM Tet (Graph+Lump " << Assembler::variantName() << ")....: "
                  << perf.avgAssemblyTimeMs() << " milliseconds @ "
                  << perf.bandwidthGBs() << " GB/s (average of " << numIterations
                  << " samples) [matrix norm: " << matrix_norm
                  << ", rhs norm: " << rhs_norm << "]\n";
    }
    if (!quiet) perf.print(rank);

    cudaStreamDestroy(assemblyStream);
    cudaStreamDestroy(resetStream);
    cudaFree(d_matrix);
    MPI_Finalize();
    return 0;
}
