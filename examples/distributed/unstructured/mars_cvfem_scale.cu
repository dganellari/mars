// CVFEM graph-assembly weak-scaling driver with PROCEDURAL per-rank cube mesh
// generation -- no mesh file. Each rank generates its slice of an Ncells^3 hex
// cube in memory and constructs ElementDomain via the in-memory ctor; cstone
// redistributes by SFC exactly as for a file mesh. This removes the mesh-file
// ceiling (an O(10 TB) file at trillion-DOF scale) so weak scaling can push past
// where files die. Same assembly harness as mars_cvfem_graph (graph sparsity +
// assembleGraphLump), so the timing is directly comparable.
//
//   --ncells=N   global cube is N^3 hex elements (N^3 elements, (N+1)^3 nodes)
//   --iters=K    timed assembly iterations (default 20)
//
// Cross-check: procedural --ncells=256 should match file cube256 (same nnz,
// same assembly time/node) -- proves the generator before extrapolating.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <cstdio>
#include <string>

using namespace mars;
using namespace mars::fem;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);

    using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
    using Assembler = CvfemHexAssembler<KeyType, RealType>;
    using Domain    = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

    size_t ncells = 256; int iters = 20; int blockSize = 256; bool buildOnly = false;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--ncells=",0)==0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--iters=",0)==0) iters = std::stoi(a.substr(8));
        else if (a.rfind("--blockSize=",0)==0) blockSize = std::stoi(a.substr(12));
        else if (a=="--build-only") buildOnly = true; }

    if (rank == 0) printf("Procedural cube: %zu^3 = %zu elements, %d ranks\n",
                          ncells, ncells*ncells*ncells, numRanks);

    // --- procedural per-rank mesh generation (no file) ---
    auto [genNodes, genElems, gx, gy, gz, lconn] =
        generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks);
    (void)genNodes; (void)genElems;
    typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
    typename Domain::HostConnectivityTuple h_conn{
        std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
        std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};

    double tBuild0 = MPI_Wtime();
    Domain domain(h_coords, h_conn, rank, numRanks, 64, false, 8u);
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    if (numRanks > 1) domain.getNodeHaloTopology();   // force halo-topo build (part of DD)
    double buildSec = MPI_Wtime() - tBuild0;

    // --build-only: report the domain decomposition (SFC partition + node
    // ownership + halo topology) and STOP before the CSR. The build is
    // sync-dominated (~0.55 GiB/M-elem), so it scales far past the assembly /
    // operator per-GPU ceiling -- this is the trillion-scale DD path.
    if (buildOnly) {
        unsigned long long locE = elementCount, gE = 0, minE = 0, maxE = 0, gN = 0, locN = nodeCount;
        MPI_Reduce(&locE, &gE,   1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&locE, &minE, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&locE, &maxE, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&locN, &gN,   1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        double maxBuild = 0; MPI_Reduce(&buildSec, &maxBuild, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        size_t freeB = 0, totB = 0; cudaMemGetInfo(&freeB, &totB);
        long long locFree = (long long)(freeB >> 20), minFree = 0;
        MPI_Reduce(&locFree, &minFree, 1, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("\n== DD-ONLY domain decomposition (no assembly) ==\n");
            printf("ranks=%d | global elems=%llu  global nodes=%llu | ~%.1fM elems/rank\n",
                   numRanks, gE, gN, double(gE)/double(numRanks)/1e6);
            printf("per-rank owned elems: min=%llu max=%llu\n", minE, maxE);
            printf("domain build (SFC partition + ownership + halo): %.3f s (slowest rank)\n", maxBuild);
            printf("min GPU free across ranks after build: %lld MiB\n", minFree);
        }
        MPI_Finalize();
        return 0;
    }

    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_nodeToDof.data(), nodeCount);
    int numTotalDofs = (int)nodeCount;

    const auto& d_conn = domain.getElementToNodeConnectivity();
    auto C = [&](int i)->const KeyType* {
        return (i==0?std::get<0>(d_conn):i==1?std::get<1>(d_conn):i==2?std::get<2>(d_conn):i==3?std::get<3>(d_conn):
                i==4?std::get<4>(d_conn):i==5?std::get<5>(d_conn):i==6?std::get<6>(d_conn):std::get<7>(d_conn)).data(); };

    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1), d_diagPtr(numTotalDofs), d_colInd;
    int nnz = CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
        C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount, d_nodeToDof.data(), numTotalDofs,
        d_rowPtr.data(), nullptr, nullptr, 0);
    d_colInd.resize(nnz);
    CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
        C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount, d_nodeToDof.data(), numTotalDofs,
        d_rowPtr.data(), d_colInd.data(), d_diagPtr.data(), 0);
    cudaDeviceSynchronize();

    cstone::DeviceVector<RealType> d_values(nnz, 0.0), d_rhs(numTotalDofs, 0.0);

    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX(); const auto& d_y = domain.getNodeY(); const auto& d_z = domain.getNodeZ();
    cstone::DeviceVector<RealType> d_phi(nodeCount,0.0), d_gamma(nodeCount,1.0), d_beta(nodeCount,0.0);
    cstone::DeviceVector<RealType> d_gx(nodeCount,0.0), d_gy(nodeCount,0.0), d_gz(nodeCount,0.0);
    cstone::DeviceVector<RealType> d_mdot(elementCount*12, 0.0);
    cstone::DeviceVector<RealType> d_avx(elementCount*12), d_avy(elementCount*12), d_avz(elementCount*12);
    precomputeAreaVectorsGpu<KeyType, RealType>(C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount,
        d_x.data(), d_y.data(), d_z.data(), d_avx.data(), d_avy.data(), d_avz.data());
    cudaDeviceSynchronize();

    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix; cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(),
                        numTotalDofs, nnz, numDofs};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    Assembler::Config cfg; cfg.blockSize = blockSize; cfg.variant = CvfemKernelVariant::Tensor;
    auto assemble = [&]() {
        thrust::fill(thrust::device_pointer_cast(d_values.data()),
                     thrust::device_pointer_cast(d_values.data() + nnz), 0.0);
        Assembler::assembleGraphLump(C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount,
            d_x.data(), d_y.data(), d_z.data(), d_gamma.data(), d_phi.data(), d_beta.data(),
            d_gx.data(), d_gy.data(), d_gz.data(), d_mdot.data(), d_avx.data(), d_avy.data(), d_avz.data(),
            d_nodeToDof.data(), d_nodeOwnership.data(), d_matrix, d_rhs.data(), cfg);
    };

    assemble(); cudaDeviceSynchronize();   // warmup

    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    cudaEventRecord(t0);
    for (int it = 0; it < iters; ++it) assemble();
    cudaEventRecord(t1); cudaEventSynchronize(t1);
    float ms = 0; cudaEventElapsedTime(&ms, t0, t1); ms /= iters;

    // Global totals + slowest-rank assembly time.
    unsigned long long locElems = domain.localElementCount(), gElems = 0;
    unsigned long long locDofs = (unsigned long long)numDofs, gDofs = 0;
    MPI_Reduce(&locElems, &gElems, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&locDofs,  &gDofs,  1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    float maxMs = 0; MPI_Reduce(&ms, &maxMs, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double dofMperS = (double)gDofs / (maxMs * 1e-3) / 1e6;
        printf("ranks=%d | global elems=%llu  global DOFs=%llu | local nnz=%d nnz/row=%.1f\n",
               numRanks, gElems, gDofs, nnz, (double)nnz/numDofs);
        printf("assembly: %.3f ms/iter (slowest rank) | %.1f MDOF/s\n", maxMs, dofMperS);
    }

    cudaFree(d_matrix); cudaEventDestroy(t0); cudaEventDestroy(t1);
    MPI_Finalize();
    return 0;
}
