// Multi-rank gates for mars::fem::GhostRegistry (backend/distributed/unstructured/fem).
//
// Host-only + MPI on purpose: the registry deliberately depends on <mpi.h> and std alone, so its
// correctness can be pinned down without a GPU or a mesh file. Topology is a synthetic 1-D chain --
// each rank owns kOwnedPerRank entities with global keys [rank*N, rank*N+N) and ghosts its
// neighbours' kGhostDepth nearest entities.
//
// Run with several rank counts; 1 rank exercises almost nothing:
//   mpirun -np 4 ./mars_ghost_registry_mpi_test

#include <gtest/gtest.h>

#include <vector>

#include "mars_env.hpp"
// Path-relative on purpose: mars_add_test() builds a SECOND executable from this same source that
// does not inherit addMarsMpiTest()'s include directories, so a bare include resolves in one target
// and not the other.
#include "../../unstructured/fem/mars_ghost_registry.hpp"

using mars::fem::GhostRegistry;

namespace
{

constexpr int kOwnedPerRank = 100;
constexpr int kGhostDepth   = 3;

struct Layout
{
    std::vector<long>    key;
    std::vector<int>     owner;
    std::vector<uint8_t> participates;
    std::vector<int>     peers;
    int                  nLeftGhost = 0, nRightGhost = 0;

    size_t size() const { return key.size(); }
};

// Owned entities [0, kOwnedPerRank), then left-neighbour ghosts, then right-neighbour ghosts.
// Only near-edge owned entities participate, which is the realistic case (interior never ghosted).
Layout makeLayout(int rank, int size)
{
    Layout L;
    for (int i = 0; i < kOwnedPerRank; ++i)
    {
        L.key.push_back(long(rank) * kOwnedPerRank + i);
        L.owner.push_back(rank);
        const bool nearEdge = (i < kGhostDepth) || (i >= kOwnedPerRank - kGhostDepth);
        L.participates.push_back(nearEdge ? 1 : 0);
    }
    if (rank > 0)
    {
        L.peers.push_back(rank - 1);
        for (int i = 0; i < kGhostDepth; ++i)
        {
            L.key.push_back(long(rank - 1) * kOwnedPerRank + (kOwnedPerRank - kGhostDepth + i));
            L.owner.push_back(rank - 1);
            L.participates.push_back(1);
            ++L.nLeftGhost;
        }
    }
    if (rank + 1 < size)
    {
        L.peers.push_back(rank + 1);
        for (int i = 0; i < kGhostDepth; ++i)
        {
            L.key.push_back(long(rank + 1) * kOwnedPerRank + i);
            L.owner.push_back(rank + 1);
            L.participates.push_back(1);
            ++L.nRightGhost;
        }
    }
    return L;
}

double valueOf(long key) { return 1000.0 + double(key) * 0.5; }

int myRank()
{
    int r = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    return r;
}
int numRanks()
{
    int s = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &s);
    return s;
}

} // namespace

// Send/recv counts must agree with the peer's, which is what the receiver-driven build buys: a
// mismatch here is the MPI_ERR_TRUNCATE failure mode a sender-driven build invites.
TEST(GhostRegistry, SendRecvSymmetry)
{
    const int rank = myRank(), size = numRanks();
    Layout L = makeLayout(rank, size);

    GhostRegistry<long> reg("interface");
    reg.build(L.size(), L.owner, L.key, L.participates, rank, L.peers);

    EXPECT_TRUE(reg.isConsistent()) << "requested keys were unresolvable";

    for (size_t p = 0; p < reg.peers().size(); ++p)
    {
        const int peer = reg.peers()[p];
        int mySend = reg.sendCount(int(p)), myRecv = reg.recvCount(int(p));
        int peerSendToMe = -1, peerRecvFromMe = -1;
        MPI_Sendrecv(&mySend, 1, MPI_INT, peer, 77, &peerSendToMe, 1, MPI_INT, peer, 77,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myRecv, 1, MPI_INT, peer, 78, &peerRecvFromMe, 1, MPI_INT, peer, 78,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        EXPECT_EQ(mySend, peerRecvFromMe) << "my send count != peer " << peer << "'s recv count";
        EXPECT_EQ(myRecv, peerSendToMe) << "my recv count != peer " << peer << "'s send count";
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// forward(): every ghost slot must end up holding its OWNER's value.
TEST(GhostRegistry, ForwardOwnerToGhost)
{
    const int rank = myRank(), size = numRanks();
    Layout L = makeLayout(rank, size);

    GhostRegistry<long> reg("interface");
    reg.build(L.size(), L.owner, L.key, L.participates, rank, L.peers);

    std::vector<double> v(L.size());
    for (size_t e = 0; e < L.size(); ++e)
        v[e] = (L.owner[e] == rank) ? valueOf(L.key[e]) : -12345.0;   // ghosts start as garbage

    reg.forwardHost(v);

    for (size_t e = 0; e < L.size(); ++e)
        if (L.owner[e] != rank)
            EXPECT_DOUBLE_EQ(v[e], valueOf(L.key[e])) << "ghost " << e << " missed the owner value";
    MPI_Barrier(MPI_COMM_WORLD);
}

// reverseAdd(): an owned entity collects one contribution per remote ghost referencing it.
TEST(GhostRegistry, ReverseAddGhostToOwner)
{
    const int rank = myRank(), size = numRanks();
    Layout L = makeLayout(rank, size);

    GhostRegistry<long> reg("interface");
    reg.build(L.size(), L.owner, L.key, L.participates, rank, L.peers);

    std::vector<double> v(L.size(), 1.0);
    reg.reverseAddHost(v);

    for (int i = 0; i < kOwnedPerRank; ++i)
    {
        const bool ghostedLeft  = (rank > 0) && (i < kGhostDepth);
        const bool ghostedRight = (rank + 1 < size) && (i >= kOwnedPerRank - kGhostDepth);
        const double expect = 1.0 + (ghostedLeft ? 1.0 : 0.0) + (ghostedRight ? 1.0 : 0.0);
        EXPECT_DOUBLE_EQ(v[i], expect) << "owned entity " << i << " got the wrong contribution count";
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// STK's parallel_owner_rank, answered from O(interface) storage rather than a per-entity array.
TEST(GhostRegistry, OwnerRankQuery)
{
    const int rank = myRank(), size = numRanks();
    Layout L = makeLayout(rank, size);

    GhostRegistry<long> reg("interface");
    reg.build(L.size(), L.owner, L.key, L.participates, rank, L.peers);

    for (size_t e = 0; e < L.size(); ++e)
    {
        EXPECT_EQ(reg.ownerRank(int(e)), L.owner[e]);
        EXPECT_EQ(reg.isGhost(int(e)), L.owner[e] != rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// Two registries live on one communicator must not exchange each other's payloads. This is what
// the per-registry MPI_Comm_dup buys, and what fixed-tag HoHalo cannot do. Registry B selects a
// strict subset BY GLOBAL KEY, so both sides of every shared entity agree on participation.
TEST(GhostRegistry, TwoRegistriesDoNotCrossTalk)
{
    const int rank = myRank(), size = numRanks();
    Layout A = makeLayout(rank, size);
    Layout B = makeLayout(rank, size);
    for (size_t e = 0; e < B.size(); ++e)
        if (B.key[e] % 2 != 0) B.participates[e] = 0;

    GhostRegistry<long> regA("interface");
    GhostRegistry<long> regB("interface_even_keys");
    regA.build(A.size(), A.owner, A.key, A.participates, rank, A.peers);
    regB.build(B.size(), B.owner, B.key, B.participates, rank, B.peers);
    EXPECT_TRUE(regB.isConsistent());

    std::vector<double> a(A.size()), b(B.size());
    for (size_t e = 0; e < A.size(); ++e) a[e] = (A.owner[e] == rank) ? valueOf(A.key[e]) : -1.0;
    for (size_t e = 0; e < B.size(); ++e) b[e] = (B.owner[e] == rank) ? valueOf(B.key[e]) * 7.0 : -1.0;

    // Interleaved: sharing a message space would swap these payloads.
    regA.forwardHost(a);
    regB.forwardHost(b);
    regA.reverseAddHost(a);
    regB.forwardHost(b);

    for (size_t e = 0; e < B.size(); ++e)
        if (B.owner[e] != rank && B.participates[e])
            EXPECT_DOUBLE_EQ(b[e], valueOf(B.key[e]) * 7.0) << "B payload corrupted by A";
    for (size_t e = 0; e < A.size(); ++e)
        if (A.owner[e] != rank)
            EXPECT_DOUBLE_EQ(a[e], valueOf(A.key[e])) << "A payload corrupted by B";

    if (size > 1) EXPECT_LT(regB.totalRecv(), regA.totalRecv()) << "B is not a strict subset of A";
    MPI_Barrier(MPI_COMM_WORLD);
}

// Participation must agree across ranks. If a requester asks for an entity its owner did not mark,
// the owner cannot resolve the key -- that must be reported, never silently exchanged as zero.
// The counter lives on the rank being ASKED, so the collective verdict is what a requester checks.
TEST(GhostRegistry, InconsistentParticipationIsReported)
{
    const int rank = myRank(), size = numRanks();
    if (size < 2) GTEST_SKIP() << "needs >= 2 ranks";

    Layout C = makeLayout(rank, size);
    for (int i = kOwnedPerRank - kGhostDepth; i < kOwnedPerRank; ++i) C.participates[i] = 0;

    GhostRegistry<long> regC("inconsistent");
    regC.build(C.size(), C.owner, C.key, C.participates, rank, C.peers);

    if (rank + 1 < size)
        EXPECT_GT(regC.unresolvedIncomingRequests(), 0u) << "owner did not report unresolvable keys";
    EXPECT_FALSE(regC.isConsistent()) << "collective verdict missed an inconsistent ghosting";
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
