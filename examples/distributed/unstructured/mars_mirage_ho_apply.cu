// mirage parity gate (p=1): the generated Laplacian apply vs the hand-written
// mars::fem::ho_cvfem_apply_kernel, on the SAME mesh, elemDof, and metric d_G.
//
// The generated kernel (backend/.../mirage/generated/mirage_laplacian_apply.cuh) is
// emitted by mirage-compiler from specs/laplacian.op. It shares the hand kernel's
// idx layout, __constant__ operator banks (own names, mirage_c_*), and 4-step
// sweep, differing only in that step-2's flux is inlined from the authored
// expression. Same math, same reduction order => the two device outputs must
// match to round-off (expected max rel diff 0). This is the go/no-go the plan's
// Stage 3 calls for; the algorithm was already proven bit-identical to the host
// oracle by mirage-compiler/tests/validate_host.py.
//
// Setup (domain, node DOF map, corner/elemDof extract, metric d_G) mirrors
// mars_cvfem_ho_compare.cu -- the known-good p=1 HO harness.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/mirage/generated/mirage_laplacian_apply.cuh"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"

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
    printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); \
    MPI_Abort(MPI_COMM_WORLD, 1);} } while(0)

// MARS hex corner c -> HO local DOF index (copied from mars_cvfem_ho_compare.cu).
__constant__ int c_hexCornerToHoLocal_mirage[8] = {0, 4, 6, 2, 1, 5, 7, 3};

template<typename KeyType, typename RealType>
__global__ void mirage_p1_extract_kernel(
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
        d_elemDof[e * 8 + c_hexCornerToHoLocal_mirage[c]] = d_nodeToDof[nd];
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);
    if (numRanks > 1) { if (rank == 0) printf("single-rank parity check; run with -np 1.\n"); MPI_Finalize(); return 1; }

    using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
    using Domain = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

    size_t ncells = 0; int iters = 50;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--ncells=", 0) == 0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--iters=", 0) == 0) iters = std::stoi(a.substr(8)); }
    if (ncells == 0) { if (rank == 0) printf("need --ncells=<N> [--iters=N]\n"); MPI_Finalize(); return 1; }

    auto [genNodes, genElems, gx, gy, gz, lconn] =
        generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks);
    (void)genNodes; (void)genElems;
    typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
    typename Domain::HostConnectivityTuple h_conn{
        std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
        std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};
    Domain domain(h_coords, h_conn, rank, numRanks, 64, false, 8u);

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

    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX(); const auto& d_y = domain.getNodeY(); const auto& d_z = domain.getNodeZ();

    // Reference operators -> BOTH constant-memory banks (hand c_*, generated mirage_c_*).
    HoCvfemOperators op = buildHoCvfemOperators(1);
    if (!(op.zeta[0] < op.zeta[1])) { if (rank == 0) printf("FATAL: GLL zeta not ascending.\n"); MPI_Finalize(); return 1; }
    CK(ho_cvfem_upload_operators(1, op.Btil.data(), op.Dtil.data(), op.D.data(),
                                 op.W.data(), op.xi.data(), op.zeta.data()));
    CK(mars::mirage::laplacian_upload_operators(1, op.Btil.data(), op.Dtil.data(),
                                             op.D.data(), op.W.data()));

    cstone::DeviceVector<RealType> d_corners(elementCount * 24);
    cstone::DeviceVector<int>      d_elemDof(elementCount * 8);
    {
        int blk = 256, grid = (int)((elementCount + blk - 1) / blk);
        mirage_p1_extract_kernel<KeyType, RealType><<<grid, blk>>>(
            C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),
            d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), elementCount,
            d_corners.data(), d_elemDof.data());
        CK(cudaGetLastError()); cudaDeviceSynchronize();
    }

    // Per-element metric d_G, shared by both applies (geometry is library infra).
    constexpr int kHoG = 3 * 1 * (1 + 1) * (1 + 1) * 3;   // 36 doubles/element at P=1
    cstone::DeviceVector<RealType> d_G(elementCount * kHoG);
    CK((ho_cvfem_metric_perpoint_launch<RealType, 1>(d_corners.data(), d_G.data(), elementCount)));
    cudaDeviceSynchronize();

    std::vector<RealType> h_u(numTotalDofs);
    std::mt19937 rng(7); std::uniform_real_distribution<RealType> uni(-1, 1);
    for (auto& v : h_u) v = uni(rng);
    cstone::DeviceVector<RealType> d_u(numTotalDofs), d_ymf(numTotalDofs), d_ymir(numTotalDofs);
    CK(cudaMemcpy(d_u.data(), h_u.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));

    // Hand kernel.
    CK(cudaMemset(d_ymf.data(), 0, numTotalDofs*sizeof(RealType)));
    CK((ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_ymf.data(), d_elemDof.data(),
                                           d_G.data(), elementCount)));
    // Generated kernel (same d_elemDof, same d_G).
    CK(cudaMemset(d_ymir.data(), 0, numTotalDofs*sizeof(RealType)));
    CK((mars::mirage::laplacian_apply_launch<RealType, 1>(d_u.data(), d_ymir.data(), d_elemDof.data(),
                                                       d_G.data(), elementCount)));
    cudaDeviceSynchronize();

    std::vector<RealType> ymf(numTotalDofs), ymir(numTotalDofs);
    CK(cudaMemcpy(ymf.data(),  d_ymf.data(),  numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(ymir.data(), d_ymir.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));

    // Parity over owned rows, plus A.1 constant-nullspace on the generated kernel.
    RealType nrm = 0, mx = 0;
    for (int i = 0; i < numDofs; ++i) { nrm = std::max(nrm, std::abs(ymf[i])); mx = std::max(mx, std::abs(ymir[i] - ymf[i])); }
    RealType rel = mx / (nrm > 0 ? nrm : 1.0);

    cstone::DeviceVector<RealType> d_uc(numTotalDofs), d_yc(numTotalDofs);
    std::vector<RealType> ones(numTotalDofs, 1.0), yc(numTotalDofs);
    CK(cudaMemcpy(d_uc.data(), ones.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));
    CK(cudaMemset(d_yc.data(), 0, numTotalDofs*sizeof(RealType)));
    CK((mars::mirage::laplacian_apply_launch<RealType, 1>(d_uc.data(), d_yc.data(), d_elemDof.data(),
                                                       d_G.data(), elementCount)));
    cudaDeviceSynchronize();
    CK(cudaMemcpy(yc.data(), d_yc.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    RealType nullmax = 0; for (int i = 0; i < numDofs; ++i) nullmax = std::max(nullmax, std::abs(yc[i]));

    bool ok = (rel < 1e-12) && (nullmax < 1e-9);
    printf("mirage parity (p=1): elems=%zu dofs=%d  gen_vs_hand_rel=%.3e  nullspace=%.3e  %s\n",
           elementCount, numDofs, rel, nullmax, ok ? "PASS" : "FAIL");

    // ---- Throughput: generated Mirage kernel vs hand kernel, same work ----
    // Each timed op is memset(y)=0 + one apply (matches mars_cvfem_ho_compare's
    // matfree lambda) so the comparison is apples-to-apples. Both read the same
    // d_G/d_elemDof, so any timing gap is the generated code, not the geometry.
    auto handOp = [&]() {
        cudaMemset(d_ymf.data(), 0, numTotalDofs * sizeof(RealType));
        ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_ymf.data(), d_elemDof.data(), d_G.data(), elementCount);
    };
    auto mirageOp = [&]() {
        cudaMemset(d_ymir.data(), 0, numTotalDofs * sizeof(RealType));
        mars::mirage::laplacian_apply_launch<RealType, 1>(d_u.data(), d_ymir.data(), d_elemDof.data(), d_G.data(), elementCount);
    };
    for (int w = 0; w < 5; ++w) { handOp(); mirageOp(); }   // warmup
    CK(cudaGetLastError()); cudaDeviceSynchronize();

    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    auto bench = [&](auto&& fn) {
        cudaEventRecord(t0);
        for (int i = 0; i < iters; ++i) fn();
        cudaEventRecord(t1); cudaEventSynchronize(t1);
        float ms = 0; cudaEventElapsedTime(&ms, t0, t1);
        return (double)ms / iters;
    };
    double handMs = bench(handOp);
    double mirageMs = bench(mirageOp);
    double handTP = numDofs / (handMs * 1e-3);       // DOF/s
    double mirageTP = numDofs / (mirageMs * 1e-3);
    printf("mirage perf  (p=1): iters=%d  hand=%.3f ms (%.3e dof/s)  mirage=%.3f ms (%.3e dof/s)  mirage/hand=%.3f\n",
           iters, handMs, handTP, mirageMs, mirageTP, mirageMs / handMs);
    cudaEventDestroy(t0); cudaEventDestroy(t1);

    MPI_Finalize();
    return ok ? 0 : 1;
}
