// Device-vs-host equivalence gates for coarseSearch and GhostRegistry.
//
// Both headers carry two implementations of the same algorithm: a host reference that the plain
// MPI gates pin down without a GPU, and the device path that production actually uses. This file
// is the only thing that compiles the device halves at all -- they are header-only and every other
// consumer is a .cpp, so without this the kernels are never handed to nvcc.
//
// The gate is equality, not plausibility: run both paths on the SAME input and require the results
// to match element-for-element. Both sort and dedup, so that comparison is well-defined and does
// not depend on message arrival order.
//
//   mpirun -np 4 ./mars_search_ghost_device_test

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "mars_env.hpp"
#include "../../unstructured/fem/mars_coarse_search.hpp"
#include "../../unstructured/fem/mars_ghost_registry.hpp"

using mars::fem::Aabb;
using mars::fem::coarseSearch;
using mars::fem::coarseSearchHost;
using mars::fem::GhostRegistry;
using mars::fem::SearchPair;

namespace
{

int worldRank()
{
    int r = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    return r;
}
int worldSize()
{
    int n = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    return n;
}

}  // namespace

// Gate: the device coarse search returns exactly what the host reference returns. Random boxes and
// random ownership, so the routing has to be right for non-neighbouring ranks too, not just for a
// tidy slab layout.
TEST(CoarseSearchDevice, MatchesHostReference)
{
    const int rank = worldRank(), np = worldSize();
    constexpr int kSeeds = 10;
    constexpr int nDomainGlobal = 40, nRangeGlobal = 60;

    for (int seed = 0; seed < kSeeds; ++seed)
    {
        std::mt19937 gen(4321 + seed);
        std::uniform_real_distribution<double> pos(0.0, 10.0);
        std::uniform_real_distribution<double> ext(0.05, 2.5);
        std::uniform_int_distribution<int> owner(0, np - 1);

        std::vector<Aabb<double>> myDom, myRan;
        std::vector<uint64_t> myDomId, myRanId;
        auto draw = [&](Aabb<double>& b)
        {
            for (int d = 0; d < 3; ++d)
            {
                double c = pos(gen), h = ext(gen);
                b.lo[d] = c - h;
                b.hi[d] = c + h;
            }
        };
        for (int i = 0; i < nDomainGlobal; ++i)
        {
            Aabb<double> b;
            draw(b);
            if (owner(gen) == rank) { myDom.push_back(b); myDomId.push_back(i); }
        }
        for (int j = 0; j < nRangeGlobal; ++j)
        {
            Aabb<double> b;
            draw(b);
            if (owner(gen) == rank) { myRan.push_back(b); myRanId.push_back(1000 + j); }
        }

        std::vector<int> peers;
        for (int r = 0; r < np; ++r) peers.push_back(r);

        std::vector<SearchPair> hostOut;
        coarseSearchHost(myDom, myDomId, myRan, myRanId, peers, MPI_COMM_WORLD, hostOut);

        thrust::device_vector<Aabb<double>> d_dom(myDom.begin(), myDom.end());
        thrust::device_vector<Aabb<double>> d_ran(myRan.begin(), myRan.end());
        thrust::device_vector<uint64_t> d_domId(myDomId.begin(), myDomId.end());
        thrust::device_vector<uint64_t> d_ranId(myRanId.begin(), myRanId.end());
        thrust::device_vector<SearchPair> d_out;
        coarseSearch<double>(thrust::raw_pointer_cast(d_dom.data()),
                             thrust::raw_pointer_cast(d_domId.data()), myDom.size(),
                             thrust::raw_pointer_cast(d_ran.data()),
                             thrust::raw_pointer_cast(d_ranId.data()), myRan.size(), peers,
                             MPI_COMM_WORLD, d_out);

        thrust::host_vector<SearchPair> devOut = d_out;
        ASSERT_EQ(devOut.size(), hostOut.size()) << "seed " << seed;
        for (size_t k = 0; k < hostOut.size(); ++k)
            EXPECT_TRUE(devOut[k] == hostOut[k]) << "seed " << seed << " entry " << k;
    }
}

namespace
{

// 1-D chain, same shape as the host registry gates: each rank owns a block of keys and ghosts its
// neighbours' nearest few.
struct Chain
{
    std::vector<long>    key;
    std::vector<int>     owner;
    std::vector<uint8_t> participates;
    std::vector<int>     peers;
};

Chain makeChain(int rank, int np)
{
    constexpr int kOwned = 50, kGhost = 3;
    Chain c;
    for (int i = 0; i < kOwned; ++i)
    {
        c.key.push_back(long(rank) * kOwned + i);
        c.owner.push_back(rank);
        c.participates.push_back((i < kGhost || i >= kOwned - kGhost) ? 1 : 0);
    }
    if (rank > 0)
        for (int g = 0; g < kGhost; ++g)
        {
            c.key.push_back(long(rank - 1) * kOwned + (kOwned - kGhost + g));
            c.owner.push_back(rank - 1);
            c.participates.push_back(1);
        }
    if (rank + 1 < np)
        for (int g = 0; g < kGhost; ++g)
        {
            c.key.push_back(long(rank + 1) * kOwned + g);
            c.owner.push_back(rank + 1);
            c.participates.push_back(1);
        }
    for (int r = 0; r < np; ++r) c.peers.push_back(r);
    return c;
}

}  // namespace

// Gate: device forward writes the same field as forwardHost.
TEST(GhostRegistryDevice, ForwardMatchesHost)
{
    const int rank = worldRank(), np = worldSize();
    Chain c = makeChain(rank, np);

    GhostRegistry<long> reg("dev-fwd");
    reg.build(c.key.size(), c.owner, c.key, c.participates, rank, c.peers, MPI_COMM_WORLD);
    ASSERT_TRUE(reg.isConsistent());

    std::vector<double> v0(c.key.size(), 0.0);
    for (size_t i = 0; i < c.key.size(); ++i)
        if (c.owner[i] == rank) v0[i] = double(c.key[i]);

    std::vector<double> hostV = v0;
    reg.forwardHost(hostV);

    thrust::device_vector<double> d_v(v0.begin(), v0.end());
    reg.forward(thrust::raw_pointer_cast(d_v.data()));
    thrust::host_vector<double> devV = d_v;

    ASSERT_EQ(devV.size(), hostV.size());
    for (size_t i = 0; i < hostV.size(); ++i) EXPECT_DOUBLE_EQ(devV[i], hostV[i]) << "slot " << i;
}

// Gate: device reverseAdd sums the same as reverseAddHost. This is the one that would catch the
// atomicAdd being wrong -- an owned entity ghosted by several peers appears once per peer, so a
// plain store would drop contributions where the serial host loop accumulates them.
TEST(GhostRegistryDevice, ReverseAddMatchesHost)
{
    const int rank = worldRank(), np = worldSize();
    Chain c = makeChain(rank, np);

    GhostRegistry<long> reg("dev-rev");
    reg.build(c.key.size(), c.owner, c.key, c.participates, rank, c.peers, MPI_COMM_WORLD);
    ASSERT_TRUE(reg.isConsistent());

    std::vector<double> w0(c.key.size(), 0.0);
    for (size_t i = 0; i < c.key.size(); ++i)
        if (c.owner[i] != rank) w0[i] = 1.0;

    std::vector<double> hostW = w0;
    reg.reverseAddHost(hostW);

    thrust::device_vector<double> d_w(w0.begin(), w0.end());
    reg.reverseAdd(thrust::raw_pointer_cast(d_w.data()));
    thrust::host_vector<double> devW = d_w;

    ASSERT_EQ(devW.size(), hostW.size());
    for (size_t i = 0; i < hostW.size(); ++i) EXPECT_DOUBLE_EQ(devW[i], hostW[i]) << "slot " << i;
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    mars::Env env(argc, argv);
    return RUN_ALL_TESTS();
}
