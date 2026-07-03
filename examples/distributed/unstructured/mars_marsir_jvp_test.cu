// marsir JVP gate (p=1): the GENERATED Jacobian-action kernel J(u).du vs a
// central finite difference of the generated primal residual F, on a genuinely
// NONLINEAR operator (specs/nldemo.op: flux = g2*deriv^2 + g1*dt1^2 + g0*dt2).
//
//     J(u).du  ~=  (F(u + eps.du) - F(u - eps.du)) / (2 eps)
//
// This exercises the full two-gather JVP path (primal contractions from u for
// the partials, perturbation contractions from du) -- the device twin of the
// host proof in marsir-compiler/tests/validate_jvp.py. F and J are BOTH marsir-
// generated (marsir_nldemo_apply.cuh); the check is self-contained.
//
// nldemo is a differentiation-machinery test operator, not physics. Setup
// (domain/elemDof/metric d_G) mirrors mars_marsir_ho_apply.cu.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/marsir/generated/marsir_nldemo_apply.cuh"
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

__constant__ int c_hexCornerToHoLocal_jvp[8] = {0, 4, 6, 2, 1, 5, 7, 3};

template<typename KeyType, typename RealType>
__global__ void jvp_p1_extract_kernel(
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
        d_elemDof[e * 8 + c_hexCornerToHoLocal_jvp[c]] = d_nodeToDof[nd];
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);
    if (numRanks > 1) { if (rank == 0) printf("single-rank JVP check; run with -np 1.\n"); MPI_Finalize(); return 1; }

    using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
    using Domain = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

    size_t ncells = 0;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--ncells=", 0) == 0) ncells = std::stoull(a.substr(9)); }
    if (ncells == 0) { if (rank == 0) printf("need --ncells=<N>\n"); MPI_Finalize(); return 1; }

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

    HoCvfemOperators op = buildHoCvfemOperators(1);
    if (!(op.zeta[0] < op.zeta[1])) { if (rank == 0) printf("FATAL: GLL zeta not ascending.\n"); MPI_Finalize(); return 1; }
    // Hand banks (xi/zeta) feed the metric kernel; marsir banks feed nldemo.
    CK(ho_cvfem_upload_operators(1, op.Btil.data(), op.Dtil.data(), op.D.data(),
                                 op.W.data(), op.xi.data(), op.zeta.data()));
    CK(mars::marsir::nldemo_upload_operators(1, op.Btil.data(), op.Dtil.data(),
                                             op.D.data(), op.W.data()));

    cstone::DeviceVector<RealType> d_corners(elementCount * 24);
    cstone::DeviceVector<int>      d_elemDof(elementCount * 8);
    {
        int blk = 256, grid = (int)((elementCount + blk - 1) / blk);
        jvp_p1_extract_kernel<KeyType, RealType><<<grid, blk>>>(
            C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),
            d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), elementCount,
            d_corners.data(), d_elemDof.data());
        CK(cudaGetLastError()); cudaDeviceSynchronize();
    }

    constexpr int kHoG = 3 * 1 * (1 + 1) * (1 + 1) * 3;
    cstone::DeviceVector<RealType> d_G(elementCount * kHoG);
    CK((ho_cvfem_metric_perpoint_launch<RealType, 1>(d_corners.data(), d_G.data(), elementCount)));
    cudaDeviceSynchronize();

    // Random u, du; build u +/- eps.du on host (central difference).
    const double eps = 1e-6;
    std::vector<RealType> h_u(numTotalDofs), h_du(numTotalDofs), h_up(numTotalDofs), h_um(numTotalDofs);
    std::mt19937 rng(11); std::uniform_real_distribution<RealType> uni(-1, 1);
    for (int i = 0; i < numTotalDofs; ++i) { h_u[i] = uni(rng); h_du[i] = uni(rng); }
    for (int i = 0; i < numTotalDofs; ++i) { h_up[i] = h_u[i] + eps * h_du[i]; h_um[i] = h_u[i] - eps * h_du[i]; }

    cstone::DeviceVector<RealType> d_u(numTotalDofs), d_du(numTotalDofs), d_up(numTotalDofs), d_um(numTotalDofs);
    cstone::DeviceVector<RealType> d_Fp(numTotalDofs), d_Fm(numTotalDofs), d_J(numTotalDofs);
    CK(cudaMemcpy(d_u.data(),  h_u.data(),  numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_du.data(), h_du.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_up.data(), h_up.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_um.data(), h_um.data(), numTotalDofs*sizeof(RealType), cudaMemcpyHostToDevice));

    // F(u+eps.du), F(u-eps.du) via the generated PRIMAL kernel; J(u).du via the JVP.
    CK(cudaMemset(d_Fp.data(), 0, numTotalDofs*sizeof(RealType)));
    CK((mars::marsir::nldemo_apply_launch<RealType, 1>(d_up.data(), d_Fp.data(), d_elemDof.data(), d_G.data(), elementCount)));
    CK(cudaMemset(d_Fm.data(), 0, numTotalDofs*sizeof(RealType)));
    CK((mars::marsir::nldemo_apply_launch<RealType, 1>(d_um.data(), d_Fm.data(), d_elemDof.data(), d_G.data(), elementCount)));
    CK(cudaMemset(d_J.data(), 0, numTotalDofs*sizeof(RealType)));
    CK((mars::marsir::nldemo_jvp_apply_launch<RealType, 1>(d_u.data(), d_du.data(), d_J.data(), d_elemDof.data(), d_G.data(), elementCount)));
    cudaDeviceSynchronize();

    std::vector<RealType> Fp(numTotalDofs), Fm(numTotalDofs), J(numTotalDofs);
    CK(cudaMemcpy(Fp.data(), d_Fp.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(Fm.data(), d_Fm.data(), numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(J.data(),  d_J.data(),  numTotalDofs*sizeof(RealType), cudaMemcpyDeviceToHost));

    double mx = 0, nrm = 0;
    for (int i = 0; i < numDofs; ++i) {
        double fd = (Fp[i] - Fm[i]) / (2.0 * eps);
        mx  = std::max(mx,  std::fabs(J[i] - fd));
        nrm = std::max(nrm, std::fabs(fd));
    }
    double rel = mx / (nrm > 0 ? nrm : 1.0);
    bool ok = rel < 1e-6;
    printf("marsir JVP gate (nldemo, p=1): elems=%zu dofs=%d  J_vs_FD_rel=%.3e  %s\n",
           elementCount, numDofs, rel, ok ? "PASS" : "FAIL");

    MPI_Finalize();
    return ok ? 0 : 1;
}
