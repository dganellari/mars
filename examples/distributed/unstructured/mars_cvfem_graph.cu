#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <mpi.h>
#include <iomanip>
#include <chrono>

using namespace mars;
using namespace mars::fem;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    if (argc < 2) {
        if (rank == 0) {
            std::cout << "Usage: " << argv[0] << " <mesh_file> [options]\n";
            std::cout << "Options:\n";
            std::cout << "  --iterations=N      Number of assembly iterations (default: 10)\n";
            std::cout << "  --block-size=N      CUDA block size (default: 256)\n";
            std::cout << "  --kernel=VARIANT    original, optimized, shmem, team (default: original)\n";
            std::cout << "  --quiet             Suppress detailed output\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string meshFile = argv[1];

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
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "team") kernelVariant = CvfemKernelVariant::Team;
            else if (v == "tensor") kernelVariant = CvfemKernelVariant::Tensor;
        } else if (arg == "--quiet") {
            quiet = true;
        }
    }

    using KeyType = uint64_t;
    using RealType = double;
    using ElemTag = HexTag;
    using Assembler = CvfemHexAssembler<KeyType, RealType>;

    PerfCounters perf;
    GpuTimer timer;

    const char* kernelName = Assembler::variantName(kernelVariant);

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

    // Build DOF mapping on GPU
    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs = buildDofMappingGpu<KeyType>(
        d_nodeOwnership.data(),
        d_nodeToDof.data(),
        nodeCount
    );
    perf.numDofs = numDofs;

    // Build sparsity pattern on GPU
    auto sparsityStart = std::chrono::high_resolution_clock::now();

    const auto& d_conn = domain.getElementToNodeConnectivity();

    // Allocate CSR structure with diagonal positions
    cstone::DeviceVector<int> d_rowPtr(numDofs + 1);
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr(numDofs);

    // Build sparsity pattern entirely on GPU
    int nnz = CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        elementCount,
        d_nodeToDof.data(),
        numDofs,
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
        numDofs,
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
    cstone::DeviceVector<RealType> d_rhs(numDofs, 0.0);

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

    // Data sizes for bandwidth
    perf.coordinateBytes = 3 * nodeCount * sizeof(RealType);
    perf.inputFieldBytes = 6 * nodeCount * sizeof(RealType);
    perf.connectivityBytes = 8 * elementCount * sizeof(KeyType);
    perf.outputMatrixBytes = nnz * sizeof(RealType) + numDofs * sizeof(RealType);

    // Create CSR matrix wrapper with diagonal positions
    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(), numDofs, nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    // Configure assembler
    Assembler::Config config;
    config.blockSize = blockSize;
    config.variant = kernelVariant;

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
        thrust::fill(thrust::device_pointer_cast(d_values.data()),
                     thrust::device_pointer_cast(d_values.data() + nnz), 0.0);
        thrust::fill(thrust::device_pointer_cast(d_rhs.data()),
                     thrust::device_pointer_cast(d_rhs.data() + numDofs), 0.0);
    };

    // Warmup
    timer.start();
    runAssembly();
    timer.stop();
    perf.assemblyWarmupTime = timer.elapsedMs();

    // Timed iterations
    for (int iter = 0; iter < numIterations; ++iter) {
        resetSystem();
        cudaDeviceSynchronize();
        timer.start();
        runAssembly();
        timer.stop();
        perf.recordIteration(timer.elapsedMs());
    }

    // Compute matrix norm
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

    RealType global_norm2, global_rhs_norm2;
    MPI_Allreduce(&norm2, &global_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&rhs_norm2, &global_rhs_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    RealType matrix_norm = std::sqrt(global_norm2);
    RealType rhs_norm = std::sqrt(global_rhs_norm2);

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

    cudaFree(d_matrix);
    MPI_Finalize();
    return 0;
}
