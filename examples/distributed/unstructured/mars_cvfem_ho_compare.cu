// THE order comparison: the new high-order matrix-free CVFEM apply (Knaus Alg 2,
// mars_cvfem_ho_matfree.hpp) vs the assembled CSR + cuSPARSE SpMV path, on the
// SAME mesh and the SAME p=1 DOF space. It reports (1) PARITY, (2) THROUGHPUT of
// both, (3) MEMORY footprint, and (4) an analytic cross-p crossover table.
//
// At p=1 BOTH paths are "CVFEM diffusion at order 1", but from DIFFERENT
// formulations: the assembled CvfemHexAssembler uses MARS's SCS area-vector flux
// discretization, while the HO matrix-free uses the sum-factorized tensor-product
// stiffness of Knaus Alg 2. They SHOULD match (same method/order) -- but that is
// a physics claim, not an assumption. We therefore measure matfree(u) vs SpMV(A,u)
// on a shared random u and CHARACTERIZE any difference: a constant scaling factor
// (report the mean ratio) vs a structural per-row difference (report the spread).
//
// DOF sharing is the parity prerequisite: both paths use the SAME node->DOF map
// the assembled path builds (buildDofMappingGpu). HODofHandler is NOT used here --
// it renumbers DOFs and would break the shared-DOF parity. At p=1 the HO DOF map
// is exactly the node DOF map, so elemDof is the node map gathered through the
// element connectivity, reordered from MARS hex-corner order to the HO local index
// l = i*nn + j*n + k (see hexCornerToHoLocal below).
//
// Mesh: --mesh=<dir> (file) OR --ncells=N (procedural cube, no file). Single-rank
// focus; multi-rank is guarded as a TODO (the cross-rank owned/ghost SpMV reduce
// and HO halo apply are not wired here).

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"

#include <cusparse.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace mars;
using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return 1; } } while(0)

// MARS hex corner c -> HO local DOF index l = i*4 + j*2 + k, where (i,j,k) are the
// reference-cube indices (0/1) of corner c from c_hexCornerRef signs:
// i=(S0+1)/2, j=(S1+1)/2, k=(S2+1)/2. This is the ONLY place the two corner
// conventions are reconciled; the metric kernel reads d_corners[c] indexed by MARS
// corner c (it owns c_hexCornerRef), so corners stay in MARS order, only elemDof
// is permuted to the HO local layout.
__constant__ int c_hexCornerToHoLocal[8] = {0, 4, 6, 2, 1, 5, 7, 3};

// One thread per element: gather the 8 corner coords (MARS order) into
// d_corners[e*24 + c*3 + d] and the 8 node DOFs into d_elemDof[e*8 + l] (HO order).
template<typename KeyType, typename RealType>
__global__ void ho_p1_extract_kernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* __restrict__ d_x, const RealType* __restrict__ d_y, const RealType* __restrict__ d_z,
    const int* __restrict__ d_nodeToDof, size_t numElements,
    RealType* __restrict__ d_corners, int* __restrict__ d_elemDof)
{
    size_t e = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    const KeyType node[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    #pragma unroll
    for (int c = 0; c < 8; ++c) {
        const KeyType nd = node[c];
        d_corners[e * 24 + c * 3 + 0] = d_x[nd];
        d_corners[e * 24 + c * 3 + 1] = d_y[nd];
        d_corners[e * 24 + c * 3 + 2] = d_z[nd];
        d_elemDof[e * 8 + c_hexCornerToHoLocal[c]] = d_nodeToDof[nd];
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);

    if (numRanks > 1) {
        // Multi-rank needs cross-rank owned/ghost reduction for SpMV and a halo
        // apply for the HO path; neither is wired. Single-rank is the comparison.
        if (rank == 0) printf("TODO: multi-rank not supported (single-rank only). Run with -np 1.\n");
        MPI_Finalize(); return 1;
    }

    using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
    using Assembler = CvfemHexAssembler<KeyType, RealType>;
    using Domain    = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

    std::string mesh; size_t ncells = 0; int iters = 50;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--mesh=", 0) == 0) mesh = a.substr(7);
        else if (a.rfind("--ncells=", 0) == 0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--iters=", 0) == 0) iters = std::stoi(a.substr(8)); }
    if (mesh.empty() && ncells == 0) {
        if (rank == 0) printf("need --mesh=<dir> OR --ncells=<N> [--iters=N]\n");
        MPI_Finalize(); return 1;
    }

    // --- build ElementDomain: file or procedural cube ---
    Domain* domainPtr = nullptr;
    if (!mesh.empty()) {
        domainPtr = new Domain(mesh, rank, numRanks, true, 64, 8u);
    } else {
        auto [genNodes, genElems, gx, gy, gz, lconn] =
            generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks);
        (void)genNodes; (void)genElems;
        typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
        typename Domain::HostConnectivityTuple h_conn{
            std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
            std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};
        domainPtr = new Domain(h_coords, h_conn, rank, numRanks, 64, false, 8u);
    }
    Domain& domain = *domainPtr;

    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();

    // --- SHARED p=1 node DOF map (assembled path builds it; HO reuses it) ---
    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_nodeToDof.data(), nodeCount);
    int numTotalDofs = (int)nodeCount;

    const auto& d_conn = domain.getElementToNodeConnectivity();
    auto C = [&](int i)->const KeyType* {
        return (i==0?std::get<0>(d_conn):i==1?std::get<1>(d_conn):i==2?std::get<2>(d_conn):i==3?std::get<3>(d_conn):
                i==4?std::get<4>(d_conn):i==5?std::get<5>(d_conn):i==6?std::get<6>(d_conn):std::get<7>(d_conn)).data(); };

    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX(); const auto& d_y = domain.getNodeY(); const auto& d_z = domain.getNodeZ();

    // ====================== ASSEMBLED PATH (verbatim from matfree_compare) ======================
    // FULL sparsity so the assembled matrix carries all 8x8 intra-element couplings.
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

    // Pure diffusion, unit diffusivity: gamma=1, mdot=0, advection fields 0.
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

    // ====================== MATRIX-FREE HO PATH (p=1) ======================
    // Reference operators -> constant memory (once).
    HoCvfemOperators op = buildHoCvfemOperators(1);
    // The hex-corner -> HO-local map (c_hexCornerToHoLocal) is only correct if the
    // GLL nodes are ASCENDING (zeta[0]=-1 < zeta[1]=+1), so local index i maps to
    // reference sign sign(i). gllBasis sorts ascending today; assert it so any
    // future basis-ordering change trips loudly instead of silently corrupting parity.
    if (!(op.zeta[0] < op.zeta[1])) {
        if (rank == 0) printf("FATAL: GLL zeta not ascending (%.3f,%.3f); corner map invalid.\n",
                              op.zeta[0], op.zeta[1]);
        MPI_Finalize(); return 1;
    }
    CK(ho_cvfem_upload_operators(1, op.Btil.data(), op.Dtil.data(), op.D.data(),
                                 op.W.data(), op.xi.data(), op.zeta.data()));

    // Extract corners (MARS order) + elemDof (HO order) from the shared DOF map.
    cstone::DeviceVector<RealType> d_corners(elementCount * 24);
    cstone::DeviceVector<int>      d_elemDof(elementCount * 8);
    {
        int blk = 256, grid = (int)((elementCount + blk - 1) / blk);
        ho_p1_extract_kernel<KeyType, RealType><<<grid, blk>>>(
            C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),
            d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), elementCount,
            d_corners.data(), d_elemDof.data());
        CK(cudaGetLastError()); cudaDeviceSynchronize();
    }

    // Per-element metric, precomputed ONCE (apply never touches geometry).
    // PerPoint layout: pointsPerElem = 3*P*NN vec3 per element (a vec3 per (dir,l,s,r)
    // point); at P=1 that is 3*1*4 = 12 points * 3 = 36 doubles/elem. This MUST equal
    // the metric kernel's d_G stride (3*P*NN*3) for the allocation and writes to agree.
    constexpr int kHoG = 3 * 1 * (1 + 1) * (1 + 1) * 3;   // = 36 doubles/element at P=1
    cstone::DeviceVector<RealType> d_G(elementCount * kHoG);
    CK((ho_cvfem_metric_perpoint_launch<RealType, 1>(d_corners.data(), d_G.data(), elementCount)));
    cudaDeviceSynchronize();

    // ====================== random shared input u ======================
    std::vector<RealType> h_u(numTotalDofs);
    std::mt19937 rng(7); std::uniform_real_distribution<RealType> uni(-1, 1);
    for (auto& v : h_u) v = uni(rng);
    cstone::DeviceVector<RealType> d_u(numTotalDofs), d_ymb(numTotalDofs), d_ymf(numTotalDofs);
    CK(cudaMemcpy(d_u.data(), h_u.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));

    // assembled SpMV (cuSPARSE)
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

    // HO matrix-free apply (d_y pre-zeroed; scatter is additive via atomicAdd).
    auto matfree = [&]() {
        cudaMemset(d_ymf.data(), 0, numTotalDofs*sizeof(RealType));
        ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_ymf.data(), d_elemDof.data(),
                                           d_G.data(), elementCount);
    };

    // Warmup + bring both results to host.
    spmv(); matfree(); CK(cudaGetLastError()); cudaDeviceSynchronize();
    std::vector<RealType> ymb(numTotalDofs), ymf(numTotalDofs);
    CK(cudaMemcpy(ymb.data(), d_ymb.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(ymf.data(), d_ymf.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));

    // ====================== PARITY characterization (owned rows only) ======================
    // (1) relative max difference; (2) mean ratio matfree/SpMV (constant-scaling
    // probe); (3) spread of that ratio across rows (structural probe). We only ratio
    // rows where |SpMV| is well above noise so a near-zero denominator doesn't blow
    // up the statistic; the same threshold gates both the mean and the spread.
    RealType nrm = 0;
    for (int i = 0; i < numDofs; ++i) nrm = std::max(nrm, std::abs(ymb[i]));
    const RealType rowFloor = (nrm > 0 ? nrm : 1.0) * 1e-8;

    RealType mx = 0;
    for (int i = 0; i < numDofs; ++i) mx = std::max(mx, std::abs(ymf[i] - ymb[i]));
    RealType parity = mx / (nrm > 0 ? nrm : 1.0);

    double ratioSum = 0; int ratioN = 0;
    double ratioMin = 1e300, ratioMax = -1e300;
    for (int i = 0; i < numDofs; ++i) {
        if (std::abs(ymb[i]) < rowFloor) continue;
        double rt = ymf[i] / ymb[i];
        ratioSum += rt; ++ratioN;
        ratioMin = std::min(ratioMin, rt); ratioMax = std::max(ratioMax, rt);
    }
    double ratioMean = ratioN ? ratioSum / ratioN : 0.0;
    // relative spread: how far the per-row ratio strays from its mean. ~0 => the two
    // operators differ by a single constant factor (ratioMean); large => structural.
    double ratioSpread = 0;
    if (ratioN) {
        for (int i = 0; i < numDofs; ++i) {
            if (std::abs(ymb[i]) < rowFloor) continue;
            double rt = ymf[i] / ymb[i];
            ratioSpread = std::max(ratioSpread, std::abs(rt - ratioMean));
        }
        ratioSpread = (std::abs(ratioMean) > 0) ? ratioSpread / std::abs(ratioMean) : ratioSpread;
    }

    // ====================== timing ======================
    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    auto timeIt = [&](auto fn) { cudaEventRecord(t0); for (int i=0;i<iters;++i) fn();
        cudaEventRecord(t1); cudaEventSynchronize(t1); float ms=0; cudaEventElapsedTime(&ms,t0,t1); return ms/iters; };
    float msMB = timeIt(spmv);
    float msMF = timeIt(matfree);
    CK(cudaGetLastError());

    if (rank == 0) {
        // Memory = what each path's MATVEC reads every apply (steady-state CG cost),
        // not transient setup. Assembled SpMV reads values(8B)+colInd(4B) per nnz +
        // rowPtr. Matrix-free apply reads the per-element metric d_G (kHoG doubles)
        // + the elemDof map (8 ints); it does NOT read d_corners (those feed the
        // one-time metric precompute and are freeable afterward). Both amortized
        // over numDofs. Index-only auxiliaries (d_diagPtr) are excluded from both.
        double matBytes = (double)nnz * (8.0 + 4.0) + (double)(numTotalDofs + 1) * 4.0;
        double mfBytes  = (double)elementCount * (kHoG * 8.0 + 8 * 4.0);
        double matMiB = matBytes / (1<<20), mfMiB = mfBytes / (1<<20);

        double mbDofPerS = (double)numDofs / (msMB * 1e-3);   // SpMV throughput
        double mfDofPerS = (double)numDofs / (msMF * 1e-3);   // matrix-free throughput

        printf("HO CVFEM matrix-free (Knaus Alg 2, p=1) vs assembled CSR + cuSPARSE SpMV\n");
        printf("nodes=%zu elems=%zu numDofs=%d nnz=%d nnz/row=%.1f\n",
               nodeCount, elementCount, numDofs, nnz, (double)nnz/numDofs);
        printf("\n-- PARITY (matfree vs SpMV, owned rows) --\n");
        printf("  max|matfree-SpMV|/max|SpMV| = %.3e\n", parity);
        printf("  mean ratio matfree/SpMV     = %.6f   (rows used: %d)\n", ratioMean, ratioN);
        printf("  ratio spread (max|r-mean|/mean) = %.3e   [%s]\n", ratioSpread,
               ratioSpread < 1e-6 ? "CONSTANT scaling" : "STRUCTURAL difference");
        if (parity < 1e-10)
            printf("  => bit-identical operators (same discretization).\n");
        else if (ratioSpread < 1e-6)
            printf("  => same operator up to a constant factor %.6f.\n", ratioMean);
        else
            printf("  => DIFFERENT operators (per-row structural difference).\n");

        printf("\n-- THROUGHPUT --\n");
        printf("  assembled SpMV : %.4f ms/matvec | %.1f MDOF/s\n", msMB, mbDofPerS / 1e6);
        printf("  matrix-free HO : %.4f ms/matvec | %.1f MDOF/s\n", msMF, mfDofPerS / 1e6);
        printf("  speedup (SpMV time / matfree time) = %.2fx (%s)\n",
               msMB / msMF, msMF < msMB ? "matrix-free faster" : "SpMV faster");
        printf("  note: matrix-free time includes the per-matvec d_y memset its additive\n");
        printf("        scatter requires; SpMV uses beta=0 overwrite (no memset).\n");

        printf("\n-- MEMORY (operator storage) --\n");
        printf("  assembled CSR : %8.1f MiB | %6.1f bytes/DOF\n", matMiB, matBytes / numDofs);
        printf("  matrix-free   : %8.1f MiB | %6.1f bytes/elem amortized over %d DOFs (metric+elemDof, no matrix)\n",
               mfMiB, mfBytes / numDofs, numDofs);

        // ====================== cross-p crossover (analytic) ======================
        // Assembled CSR bytes/DOF: a p-order spectral hex has (2p+1)^3 stencil
        // nnz/row (the full element-graph couplings of a tensor-product element);
        // each nnz costs 8B value + 4B colInd = 12B. This is the storage that makes
        // assembly infeasible at high p.
        // Matrix-free bytes/DOF (steady-state apply, mirrors the measured number
        // above): the precomputed PerPoint metric is 3*p*(p+1)^2 vec3 = 9*p*(p+1)^2
        // doubles/element that the apply reads, plus the elemDof map (p+1)^3 ints.
        // Corners are NOT counted (transient, feed only the one-time precompute).
        // DOFs/element on an interior tensor element is p^3 (each element "owns" p^3
        // of its (p+1)^3 nodes; corners/faces shared). bytes/DOF =
        // (9*p*(p+1)^2 * 8 + (p+1)^3 * 4) / p^3.
        printf("\n-- CROSS-p MEMORY CROSSOVER (analytic, interior tensor-hex) --\n");
        printf("   p | nnz/row | assembled B/DOF | matfree B/DOF | assembled/matfree\n");
        for (int p = 1; p <= 7; ++p) {
            double stencil = std::pow(2.0 * p + 1.0, 3.0);
            double asmBpd  = stencil * 12.0;
            double dofsPerElem = std::pow((double)p, 3.0);
            double nodesPerElem = std::pow(p + 1.0, 3.0);
            double mfBpd = (9.0 * p * (p + 1.0) * (p + 1.0) * 8.0 + nodesPerElem * 4.0) / dofsPerElem;
            printf("  %2d | %7.0f | %15.1f | %13.1f | %14.1fx\n",
                   p, stencil, asmBpd, mfBpd, asmBpd / mfBpd);
        }
        printf("\n  Measured throughput anchors (GDOF/s): matrix-free HO 4.8 / 7.2 / 4.9 at p=1/2/4\n");
        printf("  (H100, prior validated run); this run's p=1 matrix-free = %.2f GDOF/s, SpMV = %.2f GDOF/s.\n",
               mfDofPerS / 1e9, mbDofPerS / 1e9);

        // ---- Element-level stencil dump (small meshes, e.g. --ncells=1 -> single
        //      8x8 element). Build BOTH operators' dense matrices by applying to the
        //      unit vectors, so the per-entry difference between MARS's SCS-flux CVFEM
        //      (assembled) and the Knaus matrix-free CVFEM is visible, not just an
        //      aggregate norm. This is what localizes the structural difference. ----
        if (numTotalDofs <= 64) {
            const int Nd = numTotalDofs;
            std::vector<double> Kmb(Nd*Nd), Kmf(Nd*Nd), col(Nd);
            const double one = 1.0;
            for (int J = 0; J < Nd; ++J) {
                cudaMemset(d_u.data(), 0, Nd*sizeof(double));
                cudaMemcpy(d_u.data()+J, &one, sizeof(double), cudaMemcpyHostToDevice);
                spmv();    cudaDeviceSynchronize();
                cudaMemcpy(col.data(), d_ymb.data(), Nd*sizeof(double), cudaMemcpyDeviceToHost);
                for (int i=0;i<Nd;++i) Kmb[i*Nd+J] = col[i];
                matfree(); cudaDeviceSynchronize();
                cudaMemcpy(col.data(), d_ymf.data(), Nd*sizeof(double), cudaMemcpyDeviceToHost);
                for (int i=0;i<Nd;++i) Kmf[i*Nd+J] = col[i];
            }
            auto pr = [&](const char* nm, const std::vector<double>& K){ printf("%s:\n", nm);
                for (int i=0;i<Nd;++i){ printf("  "); for (int j=0;j<Nd;++j) printf("%8.4f ", K[i*Nd+j]); printf("\n"); } };
            printf("\n-- ELEMENT STENCIL DUMP (%dx%d) --\n", Nd, Nd);
            pr("assembled (MARS SCS-flux CVFEM, --kernel=tensor)", Kmb);
            pr("matrix-free (Knaus CVFEM)", Kmf);
            double num=0, den=0, mx=0, mxAbs=0;
            for (int e=0;e<Nd*Nd;++e){ num+=Kmf[e]*Kmb[e]; den+=Kmf[e]*Kmf[e]; }
            const double sc = den>0 ? num/den : 0.0;
            for (int e=0;e<Nd*Nd;++e){ mxAbs=std::max(mxAbs,std::abs(Kmb[e])); mx=std::max(mx,std::abs(sc*Kmf[e]-Kmb[e])); }
            const double rel = mx/(mxAbs>0?mxAbs:1.0);
            printf("best scale (matfree->assembled)=%.6f  max|scale*mf - asm|/max|asm|=%.3e  [%s]\n",
                   sc, rel, rel<1e-9 ? "SAME stencil up to scale" : "DIFFERENT stencil");
        }
    }

    cusparseDestroySpMat(matA); cusparseDestroyDnVec(vX); cusparseDestroyDnVec(vY); cusparseDestroy(h);
    cudaFree(buf); cudaFree(d_matrix); cudaEventDestroy(t0); cudaEventDestroy(t1);
    delete domainPtr;
    MPI_Finalize();
    return 0;
}
