// Matrix-free CVFEM vs matrix-based (assembled CSR + cuSPARSE SpMV), same operator.
//
// Builds the FULL CVFEM diffusion operator two ways on the same mesh:
//   matrix-based: buildFullSparsity + assembleFull (Tensor variant, mdot=0) -> CSR,
//                 then cuSPARSE SpMV per matvec.
//   matrix-free : cvfem_hex_matfree_apply_diffusion -- recompute the per-element
//                 8x8 SCS-flux stiffness each matvec, apply, scatter. No stored matrix.
// Both are the SAME operator (matrix-free reuses the assembler's flux math), so the
// gate is matfree(u) == SpMV(A,u) to machine precision. Then we time both and report
// the matrix memory the matrix-free path avoids.
//
// p=1 (CVFEM is p=1). Expectation: SpMV wins on matvec time (reading a lean row beats
// recomputing it); matrix-free's only edge is memory. The numbers decide the discussion.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_matfree.hpp"

#include <cusparse.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace mars;
using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return 1; } } while(0)

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);

    using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
    using Assembler = CvfemHexAssembler<KeyType, RealType>;

    std::string mesh; int iters = 50;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--mesh=", 0) == 0) mesh = a.substr(7);
        else if (a.rfind("--iters=", 0) == 0) iters = std::stoi(a.substr(8)); }
    if (mesh.empty()) { if (rank==0) printf("need --mesh=<dir> [--iters=N]\n"); MPI_Finalize(); return 1; }

    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(mesh, rank, numRanks, true, 64, 8u);
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();

    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_nodeToDof.data(), nodeCount);
    int numTotalDofs = (int)nodeCount;

    const auto& d_conn = domain.getElementToNodeConnectivity();
    auto C = [&](int i)->const KeyType* {
        return (i==0?std::get<0>(d_conn):i==1?std::get<1>(d_conn):i==2?std::get<2>(d_conn):i==3?std::get<3>(d_conn):
                i==4?std::get<4>(d_conn):i==5?std::get<5>(d_conn):i==6?std::get<6>(d_conn):std::get<7>(d_conn)).data(); };

    // FULL sparsity (all intra-element couplings) so the assembled matrix matches
    // the full 8x8 the matrix-free kernel applies.
    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1), d_diagPtr(numTotalDofs), d_colInd;
    int nnz = CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount, d_nodeToDof.data(), numTotalDofs,
        d_rowPtr.data(), nullptr, nullptr);
    d_colInd.resize(nnz);
    CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount, d_nodeToDof.data(), numTotalDofs,
        d_rowPtr.data(), d_colInd.data(), d_diagPtr.data());
    cudaDeviceSynchronize();

    cstone::DeviceVector<RealType> d_values(nnz, 0.0), d_rhs(numTotalDofs, 0.0);

    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX(); const auto& d_y = domain.getNodeY(); const auto& d_z = domain.getNodeZ();

    // Pure diffusion, unit diffusivity: gamma=1, mdot=0, advection fields irrelevant (0).
    cstone::DeviceVector<RealType> d_phi(nodeCount, 0.0), d_gamma(nodeCount, 1.0), d_beta(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_gx(nodeCount, 0.0), d_gy(nodeCount, 0.0), d_gz(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 0.0);
    cstone::DeviceVector<RealType> d_avx(elementCount*12), d_avy(elementCount*12), d_avz(elementCount*12);
    precomputeAreaVectorsGpu<KeyType, RealType>(C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount,
        d_x.data(), d_y.data(), d_z.data(), d_avx.data(), d_avy.data(), d_avz.data());
    cudaDeviceSynchronize();

    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix; CK(cudaMalloc(&d_matrix, sizeof(MatrixType)));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(),
                        numTotalDofs, nnz, numDofs};
    CK(cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice));

    Assembler::Config cfg; cfg.blockSize = 256; cfg.variant = CvfemKernelVariant::Tensor;
    Assembler::assembleFull(C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount,
        d_x.data(), d_y.data(), d_z.data(), d_gamma.data(), d_phi.data(), d_beta.data(),
        d_gx.data(), d_gy.data(), d_gz.data(), d_mdot.data(), d_avx.data(), d_avy.data(), d_avz.data(),
        d_nodeToDof.data(), d_nodeOwnership.data(), d_matrix, d_rhs.data(), cfg);
    cudaDeviceSynchronize();

    // Random input u (over all DOFs; single rank => all owned).
    std::vector<RealType> h_u(numTotalDofs);
    std::mt19937 rng(7); std::uniform_real_distribution<RealType> uni(-1, 1);
    for (auto& v : h_u) v = uni(rng);
    cstone::DeviceVector<RealType> d_u(numTotalDofs), d_ymb(numTotalDofs), d_ymf(numTotalDofs);
    CK(cudaMemcpy(d_u.data(), h_u.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));

    // --- matrix-based: cuSPARSE SpMV ---
    cusparseHandle_t h; cusparseCreate(&h);
    cusparseSpMatDescr_t matA; cusparseDnVecDescr_t vX, vY;
    const RealType alpha = 1.0, beta = 0.0;
    cusparseCreateCsr(&matA, numTotalDofs, numTotalDofs, nnz, d_rowPtr.data(), d_colInd.data(),
                      d_values.data(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
    cusparseCreateDnVec(&vX, numTotalDofs, d_u.data(), CUDA_R_64F);
    cusparseCreateDnVec(&vY, numTotalDofs, d_ymb.data(), CUDA_R_64F);
    size_t bufSize = 0;
    cusparseSpMV_bufferSize(h, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vX, &beta, vY,
                            CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &bufSize);
    void* buf = nullptr; CK(cudaMalloc(&buf, bufSize));
    auto spmv = [&]() {
        cusparseSpMV(h, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vX, &beta, vY,
                     CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, buf);
    };

    // --- matrix-free apply ---
    int blk = 256, grid = (int)((elementCount + blk - 1) / blk);
    auto matfree = [&]() {
        cudaMemset(d_ymf.data(), 0, numTotalDofs*sizeof(RealType));
        cvfem_hex_matfree_apply_diffusion<KeyType, RealType, 256><<<grid, blk>>>(
            C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7), elementCount,
            d_x.data(), d_y.data(), d_z.data(), d_gamma.data(),
            d_avx.data(), d_avy.data(), d_avz.data(),
            d_nodeToDof.data(), d_nodeOwnership.data(), numDofs, d_u.data(), d_ymf.data());
    };

    // Warmup + correctness gate.
    spmv(); matfree(); CK(cudaGetLastError()); cudaDeviceSynchronize();
    std::vector<RealType> ymb(numTotalDofs), ymf(numTotalDofs);
    CK(cudaMemcpy(ymb.data(), d_ymb.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(ymf.data(), d_ymf.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    RealType nrm = 0, mx = 0;
    for (int i = 0; i < numDofs; ++i) nrm = std::max(nrm, std::abs(ymb[i]));
    for (int i = 0; i < numDofs; ++i) mx = std::max(mx, std::abs(ymf[i] - ymb[i]));
    RealType parity = mx / (nrm > 0 ? nrm : 1.0);

    // Timing.
    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    auto timeIt = [&](auto fn) { cudaEventRecord(t0); for (int i=0;i<iters;++i) fn();
        cudaEventRecord(t1); cudaEventSynchronize(t1); float ms=0; cudaEventElapsedTime(&ms,t0,t1); return ms/iters; };
    float msMB = timeIt(spmv);
    float msMF = timeIt(matfree);
    CK(cudaGetLastError());

    if (rank == 0) {
        double matMiB = (double)nnz * (8.0 + 4.0) / (1<<20);     // values + colInd
        double mfMiB  = (double)elementCount * 12 * 3 * 8.0 / (1<<20);  // area vectors (kept)
        printf("CVFEM matrix-free vs matrix-based (full diffusion operator, p=1)\n");
        printf("nodes=%zu elems=%zu numDofs=%d nnz=%d nnz/row=%.1f\n",
               nodeCount, elementCount, numDofs, nnz, (double)nnz/numDofs);
        printf("parity(matfree vs SpMV) = %.2e\n", parity);
        printf("matrix-based : %.3f ms/matvec  | stored matrix %.1f MiB\n", msMB, matMiB);
        printf("matrix-free  : %.3f ms/matvec  | stored area-vec %.1f MiB (no matrix)\n", msMF, mfMiB);
        printf("matvec speedup (matrix-based / matrix-free) = %.2fx  (%s)\n",
               msMB / msMF, msMF < msMB ? "matrix-free faster" : "matrix-based faster");
    }

    cusparseDestroySpMat(matA); cusparseDestroyDnVec(vX); cusparseDestroyDnVec(vY); cusparseDestroy(h);
    cudaFree(buf); cudaFree(d_matrix); cudaEventDestroy(t0); cudaEventDestroy(t1);
    MPI_Finalize();
    return 0;
}
