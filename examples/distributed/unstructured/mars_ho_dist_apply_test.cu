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
//
// OPTIONAL GPU NUMBERING (flag-gated; host path is the default and is UNCHANGED):
//   MARS_HO_GPU_NUMBERING=1 ./mars_ho_dist_apply_test --ncells=16 --p=2
//   ./mars_ho_dist_apply_test --gpu-numbering --ncells=16 --p=2
// The GPU path (buildDistributedGpu) does the SAME numbering with thrust. The
// local DOF ids may differ from host by a permutation, so correctness is checked
// by permutation-INVARIANT quantities, not elemDof element-wise.
//
// SELF-CHECK (A/B host vs GPU on the same config; prints the invariants):
//   ./mars_ho_dist_apply_test --self-check --ncells=16 --p=2
// Asserts: numDof/nEdge/nFace equal, and the MULTISET of DofKeys equal (sort both,
// compare). Then runs the A.1 apply gate with the GPU numbering -> end-to-end check.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler_gpu.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_halo.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
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

// Permutation-independent A/B comparison of the GPU vs host numbering. The local
// DOF ids legitimately differ (host = std::map insertion order, GPU = sorted-key
// order), so we compare the INVARIANTS that cross-rank correctness depends on:
// numDof/nEdge/nFace (the unique counts) and the MULTISET of DofKeys. Returns true
// if everything matches on this rank.
static bool selfCheckNumbering(const HODofHandler& host, const HODofHandler& gpu, int rank)
{
    bool ok = true;
    if (host.numDof != gpu.numDof || host.nEdge != gpu.nEdge || host.nFace != gpu.nFace) ok = false;
    if (rank == 0)
        printf("[self-check] numDof host=%ld gpu=%ld | nEdge host=%ld gpu=%ld | nFace host=%ld gpu=%ld  [%s]\n",
               host.numDof, gpu.numDof, host.nEdge, gpu.nEdge, host.nFace, gpu.nFace,
               (host.numDof==gpu.numDof && host.nEdge==gpu.nEdge && host.nFace==gpu.nFace) ? "MATCH" : "MISMATCH");

    if (host.dofKey.size() != gpu.dofKey.size()) {
        if (rank == 0) printf("[self-check] dofKey size differs host=%zu gpu=%zu\n",
                              host.dofKey.size(), gpu.dofKey.size());
        return false;
    }

    auto packed = [](const HODofHandler::DofKey& k) {
        return std::array<long,6>{ (long)k.kind, k.g0, k.g1, k.g2, k.g3, (long)k.pos }; };
    std::vector<std::array<long,6>> hk(host.dofKey.size()), gk(gpu.dofKey.size());
    for (size_t i = 0; i < host.dofKey.size(); ++i) hk[i] = packed(host.dofKey[i]);
    for (size_t i = 0; i < gpu.dofKey.size();  ++i) gk[i] = packed(gpu.dofKey[i]);
    std::sort(hk.begin(), hk.end());
    std::sort(gk.begin(), gk.end());

    long mismatch = 0; long firstBad = -1;
    for (size_t i = 0; i < hk.size(); ++i)
        if (hk[i] != gk[i]) { ++mismatch; if (firstBad < 0) firstBad = (long)i; }
    if (mismatch) ok = false;
    if (rank == 0)
        printf("[self-check] DofKey multiset: %ld / %zu mismatched (firstBad slot=%ld)  [%s]\n",
               mismatch, hk.size(), firstBad, mismatch == 0 ? "MATCH" : "MISMATCH");

    // dofShared / dofBoundary are per-DOF flags keyed on local ids -> not directly
    // comparable. But their TOTALS are permutation-invariant -> compare counts.
    long hShared=0,gShared=0,hBnd=0,gBnd=0;
    for (auto v : host.dofShared)   hShared += v;
    for (auto v : gpu.dofShared)    gShared += v;
    for (auto v : host.dofBoundary) hBnd += v;
    for (auto v : gpu.dofBoundary)  gBnd += v;
    if (hShared != gShared || hBnd != gBnd) ok = false;
    if (rank == 0)
        printf("[self-check] shared count host=%ld gpu=%ld | boundary count host=%ld gpu=%ld  [%s]\n",
               hShared, gShared, hBnd, gBnd, (hShared==gShared && hBnd==gBnd) ? "MATCH" : "MISMATCH");

    // Owner histogram (rank -> #DOF owned) is permutation-invariant too.
    {
        std::map<int,long> ho, go;
        for (int o : host.dofOwner) ho[o]++;
        for (int o : gpu.dofOwner)  go[o]++;
        bool ownMatch = (ho == go);
        if (!ownMatch) ok = false;
        if (rank == 0) printf("[self-check] owner histogram  [%s]\n", ownMatch ? "MATCH" : "MISMATCH");
    }
    return ok;
}

enum class Numbering { Host, Gpu, SelfCheck };

template<int P>
static void runDistApply(const DistDof& D, int rank, int numRanks, Numbering mode)
{
    const int n = P + 1, N3 = n*n*n;

    HODofHandler dof;       // the handler actually used downstream (host or GPU built)
    HODofHandler dofHost;   // only populated in SelfCheck, for the A/B compare

    if (mode == Numbering::Host) {
        if (rank == 0) { printf("[build] numbering = HOST (default)\n"); fflush(stdout); }
        double tb0 = MPI_Wtime();
        dof.buildDistributed(D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        double tb1 = MPI_Wtime();
        if (rank == 0) { printf("[build] buildDistributed (host) %.3fs (numDof=%ld)\n", tb1-tb0, dof.numDof); fflush(stdout); }
    } else if (mode == Numbering::Gpu) {
        if (rank == 0) { printf("[build] numbering = GPU (MARS_HO_GPU_NUMBERING / --gpu-numbering)\n"); fflush(stdout); }
        cudaDeviceSynchronize();
        double tb0 = MPI_Wtime();
        buildDistributedGpu(dof, D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        cudaDeviceSynchronize();
        double tb1 = MPI_Wtime();
        if (rank == 0) { printf("[build] buildDistributedGpu %.3fs (numDof=%ld)\n", tb1-tb0, dof.numDof); fflush(stdout); }
    } else { // SelfCheck: build BOTH, time each, then compare invariants. Use GPU for downstream.
        if (rank == 0) { printf("[build] numbering = SELF-CHECK (host vs GPU A/B)\n"); fflush(stdout); }
        double th0 = MPI_Wtime();
        dofHost.buildDistributed(D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        double th1 = MPI_Wtime();
        cudaDeviceSynchronize();
        double tg0 = MPI_Wtime();
        buildDistributedGpu(dof, D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        cudaDeviceSynchronize();
        double tg1 = MPI_Wtime();
        double thMax=0, tgMax=0, dh=th1-th0, dg=tg1-tg0;
        MPI_Allreduce(&dh, &thMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&dg, &tgMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0)
            printf("[self-check] build time: host %.3fs  gpu %.3fs  speedup %.2fx\n",
                   thMax, tgMax, tgMax > 0 ? thMax/tgMax : 0.0);
        bool ok = selfCheckNumbering(dofHost, dof, rank);
        int allOk = ok ? 1 : 0, gAllOk = 1;
        MPI_Allreduce(&allOk, &gAllOk, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (rank == 0)
            printf("[self-check] numbering invariants (all ranks): %s -- now running A.1 gate with GPU numbering\n",
                   gAllOk ? "PASS" : "FAIL");
    }

    double tb1 = MPI_Wtime();
    resolveHoDofOwnership(dof.dofShared, dof.dofKey, dof.dofOwner, rank, D.peers);
    double tb2 = MPI_Wtime();
    if (rank == 0) { printf("[build] resolveOwnership %.3fs\n", tb2-tb1); fflush(stdout); }
    HoHalo<RealType> halo;
    halo.build((int)dof.numDof, dof.dofOwner, dof.dofKey, dof.dofBoundary, rank, D.peers);
    double tb3 = MPI_Wtime();
    if (rank == 0) { printf("[build] HoHalo %.3fs\n", tb3-tb2); fflush(stdout); }

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
    if (rank == 0) { printf("[gate] metric ready, entering host A.1\n"); fflush(stdout); }

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
    if (rank == 0) { printf("[gate] host apply+halo done, entering Allreduce (waits on slowest rank)\n"); fflush(stdout); }

    double locMax = 0; long locOwned = 0;
    for (long d = 0; d < nDof; ++d)
        if (dof.dofOwner[d] == rank) { locMax = std::max(locMax, std::abs(y[d])); ++locOwned; }
    double gMax = 0; long gOwned = 0;
    MPI_Allreduce(&locMax,   &gMax,   1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locOwned, &gOwned, 1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);

    const char* tag = (mode==Numbering::Host) ? "host-numbering"
                    : (mode==Numbering::Gpu)  ? "gpu-numbering" : "gpu-numbering(self-check)";
    if (rank == 0)
        printf("p=%d  owned DOF=%ld  max|A.1| over owned (host halo, %s) = %.3e   [%s]%s\n",
               P, gOwned, tag, gMax, gMax < 1e-8 ? "PASS" : "FAIL",
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
    bool cliGpu = false, cliSelfCheck = false;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--ncells=",0)==0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--p=",0)==0) P = std::stoi(a.substr(4));
        else if (a == "--gpu-numbering") cliGpu = true;
        else if (a == "--self-check")    cliSelfCheck = true; }

    // Default = host numbering (UNCHANGED). Opt into GPU via env or CLI; --self-check
    // runs both and A/Bs the invariants. The env var matches MARS' other flags.
    const char* envGpu = std::getenv("MARS_HO_GPU_NUMBERING");
    bool useGpu = cliGpu || (envGpu && std::atoi(envGpu) != 0);
    Numbering mode = cliSelfCheck ? Numbering::SelfCheck
                   : useGpu       ? Numbering::Gpu
                                  : Numbering::Host;

    auto [gn, ge, gx, gy, gz, lconn] =
        generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks);
    (void)gn; (void)ge;
    typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
    typename Domain::HostConnectivityTuple h_conn{
        std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
        std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};
    Domain domain(h_coords, h_conn, rank, numRanks, 64, false, 8u);

    double te0 = MPI_Wtime();
    DistDof D = extractDistDof(domain, rank, numRanks);
    double te1 = MPI_Wtime();
    if (rank == 0) { printf("[build] extractDistDof %.1fs (nodeCount=%zu nOwnedElem=%zu)\n",
                            te1 - te0, D.nodeCount, D.nOwnedElem); fflush(stdout); }

    if (rank == 0)
        printf("\n== HO distributed matrix-free apply gate (A.1 = HO-Laplacian.const = 0) ==  ncells=%zu ranks=%d  mode=%s\n",
               (size_t)ncells, numRanks,
               mode==Numbering::Host ? "host" : mode==Numbering::Gpu ? "gpu" : "self-check");
    switch (P) {
        case 1: runDistApply<1>(D, rank, numRanks, mode); break;
        case 2: runDistApply<2>(D, rank, numRanks, mode); break;
        case 3: runDistApply<3>(D, rank, numRanks, mode); break;
        case 4: runDistApply<4>(D, rank, numRanks, mode); break;
        default: if (rank == 0) printf("unsupported p=%d (build adds 1..4)\n", P);
    }
    MPI_Finalize();
    return 0;
}
