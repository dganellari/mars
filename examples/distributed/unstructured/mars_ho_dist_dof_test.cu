// Multi-rank gate for the high-order DOF numbering + halo (Phase B foundation).
// Builds a small procedural cube, runs HODofHandler::buildDistributed (owner+key)
// and HoHalo (receiver-driven), then two correctness gates that DON'T need the
// operator:
//   A. global owned-DOF count == analytic (ncells*P+1)^3  (no orphans/double-own)
//   B. forward a constant -> every ghost DOF receives it    (halo matches correctly)
// Run on 2-4 ranks:  srun -N1 -n4 ./mars_ho_dist_dof_test --ncells=16 --p=2

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_halo.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <cstdio>
#include <string>
#include <vector>
#include <array>

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
    using Domain  = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

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

    // Trigger the lazy halo build BEFORE reading the node count: ensureHalo()
    // updates nodeCount_ to the post-sync total (owned + halo nodes), which is the
    // index space the connectivity and the halo recv-list live in. Reading
    // getNodeCount() earlier returns the stale pre-sync input count and undersizes
    // every node array. (mars_cvfem_scale relies on the same ordering.)
    domain.getNodeOwnershipMap();
    const size_t nodeCount = domain.getNodeCount();
    const size_t startE = domain.startIndex(), endE = domain.endIndex();
    const size_t nOwnedElem = endE - startE;

    // --- pull global-key map + P1 ownership to host ---
    std::vector<KeyType> h_sfc(nodeCount);
    cudaMemcpy(h_sfc.data(), thrust::raw_pointer_cast(domain.getLocalToGlobalSfcMap().data()),
               nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    std::vector<uint8_t> h_own(nodeCount);
    cudaMemcpy(h_own.data(), thrust::raw_pointer_cast(domain.getNodeOwnershipMap().data()),
               nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    // --- cornerOwner: owned -> myRank; recv-ghost -> the peer it comes from ---
    std::vector<int> cornerOwner(nodeCount, -1);
    for (size_t i = 0; i < nodeCount; ++i) if (h_own[i]) cornerOwner[i] = rank;
    std::vector<uint8_t> sharedCorner(nodeCount, 0);   // P1 boundary nodes (send U recv)
    std::vector<int> peers;
    if (numRanks > 1) {
        const auto& topo = domain.getNodeHaloTopology();
        peers = topo.peers_;
        size_t nrecv = topo.recvOffsets_.empty() ? 0 : (size_t)topo.recvOffsets_.back();
        std::vector<int> h_recv(nrecv);
        if (nrecv) cudaMemcpy(h_recv.data(), thrust::raw_pointer_cast(topo.recvNodeIds_.data()),
                              nrecv * sizeof(int), cudaMemcpyDeviceToHost);
        size_t nsend = topo.sendOffsets_.empty() ? 0 : (size_t)topo.sendOffsets_.back();
        std::vector<int> h_send(nsend);
        if (nsend) cudaMemcpy(h_send.data(), thrust::raw_pointer_cast(topo.sendNodeIds_.data()),
                              nsend * sizeof(int), cudaMemcpyDeviceToHost);
        long badRecv = 0;
        for (size_t p = 0; p < topo.peers_.size(); ++p)
            for (int i = topo.recvOffsets_[p]; i < topo.recvOffsets_[p+1]; ++i) {
                int nd = h_recv[i];
                if (nd >= 0 && nd < (int)nodeCount) cornerOwner[nd] = topo.peers_[p];
                else ++badRecv;
            }
        for (int i = 0; i < (int)nrecv; ++i) if (h_recv[i]>=0 && h_recv[i]<(int)nodeCount) sharedCorner[h_recv[i]] = 1;
        for (int i = 0; i < (int)nsend; ++i) if (h_send[i]>=0 && h_send[i]<(int)nodeCount) sharedCorner[h_send[i]] = 1;
        if (rank == 0) printf("[dbg] cornerOwner done  badRecv=%ld (recvNodeId out of range)\n", badRecv);
        fflush(stdout);
    }

    // --- elemCorners (owned elements, local node indices) + cornerGid + elemOwner ---
    std::vector<long> cornerGid(nodeCount);
    for (size_t i = 0; i < nodeCount; ++i) cornerGid[i] = (long)h_sfc[i];

    const auto& conn = domain.getElementToNodeConnectivity();
    auto C = [&](int c)->const KeyType* {
        return thrust::raw_pointer_cast(
            (c==0?std::get<0>(conn):c==1?std::get<1>(conn):c==2?std::get<2>(conn):c==3?std::get<3>(conn):
             c==4?std::get<4>(conn):c==5?std::get<5>(conn):c==6?std::get<6>(conn):std::get<7>(conn)).data()); };
    std::vector<std::array<int,8>> elemCorners(nOwnedElem);
    long badConn = 0;
    for (int c = 0; c < 8; ++c) {
        std::vector<KeyType> col(nOwnedElem);
        cudaMemcpy(col.data(), C(c) + startE, nOwnedElem * sizeof(KeyType), cudaMemcpyDeviceToHost);
        for (size_t e = 0; e < nOwnedElem; ++e) {
            long lc = (long)col[e];
            if (lc < 0 || lc >= (long)nodeCount) { ++badConn; lc = 0; }
            elemCorners[e][c] = (int)lc;
        }
    }
    std::vector<int> elemOwner(nOwnedElem, rank);
    if (rank == 0) printf("[dbg] nodeCount=%zu elemCount=%zu startE=%zu endE=%zu nOwned=%zu badConn=%ld\n",
                          nodeCount, domain.getElementCount(), startE, endE, nOwnedElem, badConn);
    fflush(stdout);

    // --- build distributed HO DOF + halo ---
    HODofHandler dof;
    dof.buildDistributed(elemCorners, (long)nodeCount, P, cornerGid, cornerOwner, elemOwner, rank, sharedCorner);
    long nShared = 0; for (long d = 0; d < dof.numDof; ++d) nShared += dof.dofShared[d];
    if (rank == 0) printf("[dbg] buildDistributed done  numDof=%ld nEdge=%ld nFace=%ld shared=%ld\n",
                          dof.numDof, dof.nEdge, dof.nFace, nShared); fflush(stdout);

    // Resolve edge/face ownership (min-rank-among-holders) BEFORE building the halo.
    resolveHoDofOwnership(dof.dofShared, dof.dofKey, dof.dofOwner, rank, peers);

    HoHalo<RealType> halo;
    halo.build((int)dof.numDof, dof.dofOwner, dof.dofKey, dof.dofBoundary, rank, peers);
    if (rank == 0) printf("[dbg] HoHalo done  sendDof=%zu recvDof=%zu peers=%zu\n",
                          halo.sendDof_.size(), halo.recvDof_.size(), halo.peers_.size()); fflush(stdout);

    // --- GATE A: global owned-DOF count == analytic ---
    long localOwned = 0;
    for (long d = 0; d < dof.numDof; ++d) if (dof.dofOwner[d] == rank) ++localOwned;
    long globalOwned = 0;
    MPI_Allreduce(&localOwned, &globalOwned, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    long side = (long)(ncells * P + 1);
    long analytic = side * side * side;

    // --- GATE B: forward a constant -> every ghost DOF must receive it ---
    std::vector<RealType> v(dof.numDof, 0.0);
    for (long d = 0; d < dof.numDof; ++d) if (dof.dofOwner[d] == rank) v[d] = 1.0;
    halo.forward(v);
    long badGhost = 0;
    for (long d = 0; d < dof.numDof; ++d)
        if (dof.dofOwner[d] >= 0 && dof.dofOwner[d] != rank && v[d] != 1.0) ++badGhost;
    long globalBad = 0;
    MPI_Allreduce(&badGhost, &globalBad, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\n== HO distributed DOF + halo gate ==  p=%d ncells=%zu ranks=%d\n", P, (size_t)ncells, numRanks);
        printf("GATE A  global owned DOF = %ld   analytic (ncells*P+1)^3 = %ld   [%s]\n",
               globalOwned, analytic, globalOwned == analytic ? "PASS" : "FAIL");
        printf("GATE B  forward-constant ghost mismatches = %ld   [%s]\n",
               globalBad, globalBad == 0 ? "PASS" : "FAIL");
    }
    MPI_Finalize();
    return 0;
}
