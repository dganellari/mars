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

#include <algorithm>
#include <random>
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

// Gate: a random sequence of update(add, remove) must leave the registry in EXACTLY the state
// build() would produce for the same final ghost set. Anything less shows up later as misaligned
// exchange slots rather than as an error, which is why this compares the bookkeeping AND the
// actual forward/reverseAdd results.
//
// Topology: rank r owns keys [r*N, r*N+N) and may ghost the nearest few of each neighbour. The
// OWNER always participates, so what varies per step is purely which ghosts each rank asks for.
TEST(GhostRegistryHost, IncrementalMatchesRebuild)
{
    int rank = 0, np = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    constexpr int kOwn = 40, kCand = 6;

    for (int seed = 0; seed < 12; ++seed)
    {
        std::vector<long>    key;
        std::vector<int>     owner;
        std::vector<uint8_t> part;
        std::vector<int>     ghostSlots, peers;
        for (int i = 0; i < kOwn; ++i)
        {
            key.push_back(long(rank) * kOwn + i);
            owner.push_back(rank);
            part.push_back(1);
        }
        auto addCand = [&](int nb)
        {
            if (nb < 0 || nb >= np) return;
            for (int g = 0; g < kCand; ++g)
            {
                key.push_back(long(nb) * kOwn + g);
                owner.push_back(nb);
                part.push_back(0);
                ghostSlots.push_back(int(key.size()) - 1);
            }
        };
        addCand(rank - 1);
        addCand(rank + 1);
        for (int r = 0; r < np; ++r) peers.push_back(r);

        GhostRegistry<long> incr("incr");
        incr.build(key.size(), owner, key, part, rank, peers, MPI_COMM_WORLD);

        std::mt19937 gen(900 + seed);
        std::uniform_int_distribution<int> coin(0, 1);
        for (int step = 0; step < 6; ++step)
        {
            std::vector<GhostRegistry<long>::Delta> add;
            std::vector<long> rem;
            for (int lid : ghostSlots)
            {
                const bool want = coin(gen) == 1;
                if (want && !part[lid]) { add.push_back({key[lid], lid, owner[lid]}); part[lid] = 1; }
                else if (!want && part[lid]) { rem.push_back(key[lid]); part[lid] = 0; }
            }
            incr.update(add, rem);
        }

        GhostRegistry<long> fresh("fresh");
        fresh.build(key.size(), owner, key, part, rank, peers, MPI_COMM_WORLD);

        ASSERT_EQ(incr.peers().size(), fresh.peers().size()) << "seed " << seed;
        for (size_t p = 0; p < incr.peers().size(); ++p)
        {
            EXPECT_EQ(incr.peers()[p], fresh.peers()[p]) << "seed " << seed;
            ASSERT_EQ(incr.sendCount(int(p)), fresh.sendCount(int(p))) << "seed " << seed;
            ASSERT_EQ(incr.recvCount(int(p)), fresh.recvCount(int(p))) << "seed " << seed;
            auto a = incr.sendList(int(p)), b = fresh.sendList(int(p));
            EXPECT_TRUE(std::equal(a.first, a.second, b.first)) << "send list, seed " << seed;
            auto c = incr.recvList(int(p)), d = fresh.recvList(int(p));
            EXPECT_TRUE(std::equal(c.first, c.second, d.first)) << "recv list, seed " << seed;
        }

        // The exchange itself, not just the bookkeeping.
        std::vector<double> v1(key.size(), 0.0), v2(key.size(), 0.0);
        for (size_t i = 0; i < key.size(); ++i)
            if (owner[i] == rank) v1[i] = v2[i] = double(key[i]);
        incr.forwardHost(v1);
        fresh.forwardHost(v2);
        EXPECT_EQ(v1, v2) << "forward, seed " << seed;

        std::vector<double> w1(key.size(), 0.0), w2(key.size(), 0.0);
        for (size_t i = 0; i < key.size(); ++i)
            if (owner[i] != rank && part[i]) w1[i] = w2[i] = 1.0;
        incr.reverseAddHost(w1);
        fresh.reverseAddHost(w2);
        EXPECT_EQ(w1, w2) << "reverseAdd, seed " << seed;
    }
}

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
