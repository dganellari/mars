// Multi-rank gate for the HIGH-ORDER matrix-free apply (Phase B step 3).
// Builds the distributed HO DOF (buildDistributed + resolveHoDofOwnership + HoHalo),
// then runs ONE distributed matvec with the same single-rank kernel and checks the
// cross-rank assembly via the proven invariant A*1 == 0:
//   u = 1 on OWNED DOF, 0 on ghost  ->  forward()  (ghosts <- owner's 1)
//   apply owned elements (elemDof from buildDistributed is owned-only)
//   reverseAdd()  (ghost contributions summed into owners)
//   max|y| over OWNED DOF must be ~0  (the unconstrained HO Laplacian sums every
//   row to zero; a wrong cross-rank assembly leaves O(1) residual at interfaces).
// Run: srun -N1 -n4 ./mars_ho_dist_apply_test --ncells=16 --p=2

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_halo.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <array>

using namespace mars;
using namespace mars::fem;

using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
using Domain  = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

// Distributed DOF inputs pulled from the domain (host side).
struct DistDof {
    size_t nodeCount = 0, nOwnedElem = 0;
    std::vector<std::array<int,8>> elemCorners;
    std::vector<long>    cornerGid;
    std::vector<int>     cornerOwner;
    std::vector<uint8_t> sharedCorner;
    std::vector<int>     elemOwner;
    std::vector<int>     peers;
    std::vector<RealType> nodeX, nodeY, nodeZ;
};

static DistDof extractDistDof(Domain& domain, int rank, int numRanks)
{
    DistDof D;
    domain.getNodeOwnershipMap();                 // trigger halo build -> post-sync count
    const size_t nodeCount = domain.getNodeCount();
    const size_t startE = domain.startIndex(), endE = domain.endIndex();
    D.nodeCount = nodeCount; D.nOwnedElem = endE - startE;

    std::vector<KeyType> h_sfc(nodeCount);
    cudaMemcpy(h_sfc.data(), thrust::raw_pointer_cast(domain.getLocalToGlobalSfcMap().data()),
               nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    std::vector<uint8_t> h_own(nodeCount);
    cudaMemcpy(h_own.data(), thrust::raw_pointer_cast(domain.getNodeOwnershipMap().data()),
               nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    D.cornerOwner.assign(nodeCount, -1);
    for (size_t i = 0; i < nodeCount; ++i) if (h_own[i]) D.cornerOwner[i] = rank;
    D.sharedCorner.assign(nodeCount, 0);
    if (numRanks > 1) {
        const auto& topo = domain.getNodeHaloTopology();
        D.peers = topo.peers_;
        size_t nrecv = topo.recvOffsets_.empty() ? 0 : (size_t)topo.recvOffsets_.back();
        std::vector<int> h_recv(nrecv);
        if (nrecv) cudaMemcpy(h_recv.data(), thrust::raw_pointer_cast(topo.recvNodeIds_.data()),
                              nrecv * sizeof(int), cudaMemcpyDeviceToHost);
        size_t nsend = topo.sendOffsets_.empty() ? 0 : (size_t)topo.sendOffsets_.back();
        std::vector<int> h_send(nsend);
        if (nsend) cudaMemcpy(h_send.data(), thrust::raw_pointer_cast(topo.sendNodeIds_.data()),
                              nsend * sizeof(int), cudaMemcpyDeviceToHost);
        for (size_t p = 0; p < topo.peers_.size(); ++p)
            for (int i = topo.recvOffsets_[p]; i < topo.recvOffsets_[p+1]; ++i) {
                int nd = h_recv[i];
                if (nd >= 0 && nd < (int)nodeCount) D.cornerOwner[nd] = topo.peers_[p];
            }
        for (int i = 0; i < (int)nrecv; ++i) if (h_recv[i]>=0 && h_recv[i]<(int)nodeCount) D.sharedCorner[h_recv[i]] = 1;
        for (int i = 0; i < (int)nsend; ++i) if (h_send[i]>=0 && h_send[i]<(int)nodeCount) D.sharedCorner[h_send[i]] = 1;
    }

    D.cornerGid.assign(nodeCount, 0);
    for (size_t i = 0; i < nodeCount; ++i) D.cornerGid[i] = (long)h_sfc[i];

    const auto& conn = domain.getElementToNodeConnectivity();
    auto C = [&](int c)->const KeyType* {
        return thrust::raw_pointer_cast(
            (c==0?std::get<0>(conn):c==1?std::get<1>(conn):c==2?std::get<2>(conn):c==3?std::get<3>(conn):
             c==4?std::get<4>(conn):c==5?std::get<5>(conn):c==6?std::get<6>(conn):std::get<7>(conn)).data()); };
    D.elemCorners.resize(D.nOwnedElem);
    for (int c = 0; c < 8; ++c) {
        std::vector<KeyType> col(D.nOwnedElem);
        cudaMemcpy(col.data(), C(c) + startE, D.nOwnedElem * sizeof(KeyType), cudaMemcpyDeviceToHost);
        for (size_t e = 0; e < D.nOwnedElem; ++e) D.elemCorners[e][c] = (int)col[e];
    }
    D.elemOwner.assign(D.nOwnedElem, rank);

    domain.cacheNodeCoordinates();
    D.nodeX.resize(nodeCount); D.nodeY.resize(nodeCount); D.nodeZ.resize(nodeCount);
    cudaMemcpy(D.nodeX.data(), thrust::raw_pointer_cast(domain.getNodeX().data()), nodeCount*sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(D.nodeY.data(), thrust::raw_pointer_cast(domain.getNodeY().data()), nodeCount*sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(D.nodeZ.data(), thrust::raw_pointer_cast(domain.getNodeZ().data()), nodeCount*sizeof(RealType), cudaMemcpyDeviceToHost);
    return D;
}

template<int P>
static void runDistApply(const DistDof& D, int rank, int numRanks)
{
    const int n = P + 1, N3 = n*n*n;

    HODofHandler dof;
    dof.buildDistributed(D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
    resolveHoDofOwnership(dof.dofShared, dof.dofKey, dof.dofOwner, rank, D.peers);
    HoHalo<RealType> halo;
    halo.build((int)dof.numDof, dof.dofOwner, dof.dofKey, dof.dofBoundary, rank, D.peers);

    const size_t nEl  = D.nOwnedElem;
    const long   nDof = dof.numDof;

    // owned-element corner coords [nEl*24] for the metric
    std::vector<double> h_corners(nEl * 24);
    for (size_t e = 0; e < nEl; ++e)
        for (int c = 0; c < 8; ++c) {
            int lc = D.elemCorners[e][c];
            h_corners[e*24 + c*3 + 0] = D.nodeX[lc];
            h_corners[e*24 + c*3 + 1] = D.nodeY[lc];
            h_corners[e*24 + c*3 + 2] = D.nodeZ[lc];
        }

    auto op = buildHoCvfemOperators(P);
    ho_cvfem_upload_operators(P, op.Btil.data(), op.Dtil.data(),
                              op.D.data(), op.W.data(), op.xi.data(), op.zeta.data());

    int* d_elemDof = nullptr; double *d_corners=nullptr, *d_G=nullptr, *d_u=nullptr, *d_y=nullptr;
    const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;
    // Checked allocation: report OOM as the per-GPU ceiling, not a null-deref crash.
    // (-n1 ceiling probe, or balanced runs where all ranks OOM together.)
    cudaError_t me = cudaMalloc(&d_elemDof, sizeof(int)    * nEl * N3);
    if (me==cudaSuccess) me = cudaMalloc(&d_corners, sizeof(double) * nEl * 24);
    if (me==cudaSuccess) me = cudaMalloc(&d_G,       sizeof(double) * gLen);
    if (me==cudaSuccess) me = cudaMalloc(&d_u,       sizeof(double) * nDof);
    if (me==cudaSuccess) me = cudaMalloc(&d_y,       sizeof(double) * nDof);
    if (me != cudaSuccess) {
        double needGB = (sizeof(double)*gLen + sizeof(double)*2.0*nDof
                         + sizeof(int)*(double)nEl*N3 + sizeof(double)*(double)nEl*24) / 1e9;
        if (rank == 0)
            printf("p=%d  nDof=%ld nEl=%zu  DEVICE OOM (needs ~%.1f GB, d_G=%.1f GB): %s  -- per-GPU ceiling exceeded\n",
                   P, nDof, nEl, needGB, sizeof(double)*gLen/1e9, cudaGetErrorString(me));
        cudaGetLastError();
        cudaFree(d_elemDof); cudaFree(d_corners); cudaFree(d_G); cudaFree(d_u); cudaFree(d_y);
        return;
    }
    cudaMemcpy(d_elemDof, dof.elemDof.data(), sizeof(int) * nEl * N3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_corners, h_corners.data(),   sizeof(double) * nEl * 24, cudaMemcpyHostToDevice);
    ho_cvfem_metric_perpoint_launch<double, P>(d_corners, d_G, nEl);
    cudaDeviceSynchronize();

    // A*1 distributed: u=1 on OWNED, forward fills ghosts, apply, reverseAdd.
    std::vector<double> u(nDof, 0.0);
    for (long d = 0; d < nDof; ++d) if (dof.dofOwner[d] == rank) u[d] = 1.0;
    halo.forward(u);
    cudaMemcpy(d_u, u.data(), sizeof(double) * nDof, cudaMemcpyHostToDevice);
    cudaMemset(d_y, 0, sizeof(double) * nDof);
    ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl);
    cudaError_t err = cudaDeviceSynchronize();
    std::vector<double> y(nDof);
    cudaMemcpy(y.data(), d_y, sizeof(double) * nDof, cudaMemcpyDeviceToHost);
    halo.reverseAdd(y);

    double locMax = 0; long locOwned = 0;
    for (long d = 0; d < nDof; ++d)
        if (dof.dofOwner[d] == rank) { locMax = std::max(locMax, std::abs(y[d])); ++locOwned; }
    double gMax = 0; long gOwned = 0;
    MPI_Allreduce(&locMax,   &gMax,   1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locOwned, &gOwned, 1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
        printf("p=%d  owned DOF=%ld  max|A.1| over owned (host halo) = %.3e   [%s]%s\n",
               P, gOwned, gMax, gMax < 1e-8 ? "PASS" : "FAIL",
               err != cudaSuccess ? "  (CUDA error!)" : "");

    // ---- device-halo matvec: same A.1, but forward/reverseAdd on the GPU (scaling path) ----
    halo.uploadDevice();
    std::vector<double> ud(nDof, 0.0);
    for (long d = 0; d < nDof; ++d) if (dof.dofOwner[d] == rank) ud[d] = 1.0;
    cudaMemcpy(d_u, ud.data(), sizeof(double) * nDof, cudaMemcpyHostToDevice);
    halo.forwardDevice(d_u);
    cudaMemset(d_y, 0, sizeof(double) * nDof);
    ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl);
    halo.reverseAddDevice(d_y);
    cudaDeviceSynchronize();
    std::vector<double> yd(nDof);
    cudaMemcpy(yd.data(), d_y, sizeof(double) * nDof, cudaMemcpyDeviceToHost);
    double locMaxD = 0;
    for (long d = 0; d < nDof; ++d) if (dof.dofOwner[d] == rank) locMaxD = std::max(locMaxD, std::abs(yd[d]));
    double gMaxD = 0; MPI_Allreduce(&locMaxD, &gMaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // timing: full device matvec (forward+apply+reverseAdd) vs apply-only -> MDOF/s + comm%
    const int warm = 5, iters = 50;
    for (int it = 0; it < warm; ++it) {
        cudaMemset(d_y, 0, sizeof(double)*nDof);
        halo.forwardDevice(d_u); ho_cvfem_apply_launch<double,P>(d_u,d_y,d_elemDof,d_G,nEl); halo.reverseAddDevice(d_y);
    }
    cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    for (int it = 0; it < iters; ++it) {
        cudaMemset(d_y, 0, sizeof(double)*nDof);
        halo.forwardDevice(d_u); ho_cvfem_apply_launch<double,P>(d_u,d_y,d_elemDof,d_G,nEl); halo.reverseAddDevice(d_y);
    }
    cudaDeviceSynchronize(); double tFull = (MPI_Wtime() - t0) / iters;
    MPI_Barrier(MPI_COMM_WORLD); double t1 = MPI_Wtime();
    for (int it = 0; it < iters; ++it) { cudaMemset(d_y,0,sizeof(double)*nDof); ho_cvfem_apply_launch<double,P>(d_u,d_y,d_elemDof,d_G,nEl); }
    cudaDeviceSynchronize(); double tApply = (MPI_Wtime() - t1) / iters;

    double maxFull = 0, maxApply = 0;
    MPI_Allreduce(&tFull,  &maxFull,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&tApply, &maxApply, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double mdofsFull  = (double)gOwned / maxFull  / 1e6;
    double mdofsApply = (double)gOwned / maxApply / 1e6;
    double commFrac   = (maxFull - maxApply) / maxFull * 100.0;

    if (rank == 0)
        printf("p=%d  device A.1 = %.3e [%s]   full matvec %.0f MDOF/s | apply-only %.0f MDOF/s | comm %.1f%%\n",
               P, gMaxD, gMaxD < 1e-8 ? "PASS" : "FAIL", mdofsFull, mdofsApply, commFrac);

    halo.freeDevice();
    cudaFree(d_elemDof); cudaFree(d_corners); cudaFree(d_G); cudaFree(d_u); cudaFree(d_y);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);

    size_t ncells = 16; int P = 2;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--ncells=",0)==0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--p=",0)==0) P = std::stoi(a.substr(4)); }

    auto [gn, ge, gx, gy, gz, lconn] =
        generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks);
    (void)gn; (void)ge;
    typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
    typename Domain::HostConnectivityTuple h_conn{
        std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
        std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};
    Domain domain(h_coords, h_conn, rank, numRanks, 64, false, 8u);

    DistDof D = extractDistDof(domain, rank, numRanks);

    if (rank == 0)
        printf("\n== HO distributed matrix-free apply gate (A.1 = HO-Laplacian.const = 0) ==  ncells=%zu ranks=%d\n",
               (size_t)ncells, numRanks);
    switch (P) {
        case 1: runDistApply<1>(D, rank, numRanks); break;
        case 2: runDistApply<2>(D, rank, numRanks); break;
        case 3: runDistApply<3>(D, rank, numRanks); break;
        case 4: runDistApply<4>(D, rank, numRanks); break;
        default: if (rank == 0) printf("unsupported p=%d (build adds 1..4)\n", P);
    }
    MPI_Finalize();
    return 0;
}
