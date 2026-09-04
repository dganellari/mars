// Multi-rank gates for mars::fem::coarseSearch (backend/distributed/unstructured/fem).
//
// Host-only + MPI on purpose, same reasoning as mars_ghost_registry_mpi.cpp: the hard part of a
// parallel search is the communication plan, and that can be pinned down without a GPU or a mesh
// file. The CUDA path in the same header is written to return the identical sorted result, so
// this suite is also the oracle it gets compared against once a GPU is available.
//
// Topology is a synthetic 1-D slab stack: rank r owns range boxes inside [r, r+1) and domain
// boxes that straddle its right-hand seam, so every rank has a partner on a NEIGHBOURING rank.
// A single rank exercises almost nothing -- run several counts:
//   mpirun -np 4 ./mars_coarse_search_mpi_test

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "mars_env.hpp"
// Path-relative on purpose -- see the note in mars_ghost_registry_mpi.cpp.
#include "../../unstructured/fem/mars_coarse_search.hpp"

using mars::fem::Aabb;
using mars::fem::coarseSearchHost;
using mars::fem::IdentProc;
using mars::fem::SearchPair;

namespace
{

constexpr int kBoxesPerRank = 8;

// Range box j on rank r spans [r + j/N, r + (j+1)/N) in x, the full unit square in y/z.
// Ids are globally unique so a pair identifies its owner without ambiguity.
Aabb<double> rangeBox(int rank, int j)
{
    const double w = 1.0 / kBoxesPerRank;
    Aabb<double> b;
    b.lo[0] = rank + j * w;
    b.hi[0] = rank + (j + 1) * w;
    b.lo[1] = b.lo[2] = 0.0;
    b.hi[1] = b.hi[2] = 1.0;
    return b;
}

uint64_t rangeId(int rank, int j) { return uint64_t(rank) * 1000 + uint64_t(j); }

// A domain box centred on rank r's right seam: it reaches into the last range box of rank r and
// the first range box of rank r+1. That is the case the whole thing exists for -- a query whose
// partner is not on the querying rank.
Aabb<double> seamBox(int rank)
{
    const double w = 1.0 / kBoxesPerRank;
    Aabb<double> b;
    b.lo[0] = (rank + 1) - 0.5 * w;
    b.hi[0] = (rank + 1) + 0.5 * w;
    b.lo[1] = b.lo[2] = 0.25;
    b.hi[1] = b.hi[2] = 0.75;
    return b;
}

uint64_t seamId(int rank) { return 900000 + uint64_t(rank); }

struct Fixture
{
    int rank = 0, nRanks = 1;
    std::vector<Aabb<double>> domain, range;
    std::vector<uint64_t> domainIds, rangeIds;
    // A superset is legal by contract, and passing every rank also exercises the drop-out path
    // for peers that turn out to have no traffic. Production supplies cstone's peer list.
    std::vector<int> peers;
};

Fixture makeFixture()
{
    Fixture f;
    MPI_Comm_rank(MPI_COMM_WORLD, &f.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &f.nRanks);
    for (int j = 0; j < kBoxesPerRank; ++j)
    {
        f.range.push_back(rangeBox(f.rank, j));
        f.rangeIds.push_back(rangeId(f.rank, j));
    }
    // The last rank's seam points outside the stack; keep it anyway, it is the "query with a
    // partner on only one side" case.
    f.domain.push_back(seamBox(f.rank));
    f.domainIds.push_back(seamId(f.rank));
    for (int r = 0; r < f.nRanks; ++r) f.peers.push_back(r);
    return f;
}

bool hasPair(const std::vector<SearchPair>& v, uint64_t dId, int dProc, uint64_t rId, int rProc)
{
    for (const auto& p : v)
        if (p.domain.id == dId && p.domain.proc == dProc && p.range.id == rId && p.range.proc == rProc)
            return true;
    return false;
}

}  // namespace

// Gate 1: the seam query finds its partner on THIS rank (last range box).
TEST(CoarseSearch, FindsPartnerOnOwnRank)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearchHost(f.domain, f.domainIds, f.range, f.rangeIds, f.peers, MPI_COMM_WORLD, out);
    EXPECT_TRUE(hasPair(out, seamId(f.rank), f.rank, rangeId(f.rank, kBoxesPerRank - 1), f.rank));
}

// Gate 2: and its partner on the NEXT rank. This is the one that fails if the routing or the
// reply leg is wrong, and it is the whole point of the exercise.
TEST(CoarseSearch, FindsPartnerOnNeighbourRank)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearchHost(f.domain, f.domainIds, f.range, f.rangeIds, f.peers, MPI_COMM_WORLD, out);
    if (f.rank + 1 < f.nRanks)
        EXPECT_TRUE(hasPair(out, seamId(f.rank), f.rank, rangeId(f.rank + 1, 0), f.rank + 1));
    else
        SUCCEED() << "last rank has no right neighbour";
}

// Gate 3: a caller only gets pairs it has a stake in. With enforceSymmetry (STK's default) that
// means owning EITHER side, not just the domain side.
TEST(CoarseSearch, ResultIsOwnedByCaller)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearchHost(f.domain, f.domainIds, f.range, f.rangeIds, f.peers, MPI_COMM_WORLD, out);
    for (const auto& p : out) EXPECT_TRUE(p.domain.proc == f.rank || p.range.proc == f.rank);
}

// Gate 4: no spurious pairs. My seam box touches the two range boxes either side of the seam, and
// under symmetry I also keep the pair rank-1's seam makes with MY first range box.
TEST(CoarseSearch, NoSpuriousPairs)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearchHost(f.domain, f.domainIds, f.range, f.rangeIds, f.peers, MPI_COMM_WORLD, out);
    const size_t fromMySeam   = (f.rank + 1 < f.nRanks) ? 2u : 1u;  // own last box, plus neighbour's first
    const size_t fromPrevSeam = (f.rank > 0) ? 1u : 0u;             // symmetry: I own that range box
    const size_t expected     = fromMySeam + fromPrevSeam;
    EXPECT_EQ(out.size(), expected);
}

// Gate 5: sorted and deduplicated, so the answer does not depend on message arrival order.
TEST(CoarseSearch, SortedAndUnique)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearchHost(f.domain, f.domainIds, f.range, f.rangeIds, f.peers, MPI_COMM_WORLD, out);
    EXPECT_TRUE(std::is_sorted(out.begin(), out.end()));
    EXPECT_EQ(std::adjacent_find(out.begin(), out.end()), out.end());
}

// Gate 6: a rank contributing nothing must not hang the collective, and must not corrupt the
// answer for everyone else. This is the failure mode that cost a night on MARS_OFS_DBG.
TEST(CoarseSearch, EmptyRankDoesNotHang)
{
    Fixture f = makeFixture();
    std::vector<Aabb<double>> emptyBoxes;
    std::vector<uint64_t> emptyIds;
    std::vector<SearchPair> out;
    const bool mute = (f.nRanks > 1) && (f.rank == f.nRanks - 1);
    coarseSearchHost(mute ? emptyBoxes : f.domain, mute ? emptyIds : f.domainIds,
                     mute ? emptyBoxes : f.range, mute ? emptyIds : f.rangeIds, f.peers,
                     MPI_COMM_WORLD, out);
    if (mute) EXPECT_TRUE(out.empty());
    else SUCCEED();
}

// Gate 7: randomized brute-force oracle. The coarse-search result is DEFINED, not chosen -- it is
// the set of (domain, range) pairs whose AABBs intersect -- so an O(N^2) serial pass over the
// global boxes is the ground truth, and it is the same set stk::search::coarse_search returns
// under the conventions listed at the bottom of this file. This catches routing bugs the hand-built
// slab fixture cannot: boxes spanning three or more ranks, overlaps between NON-neighbouring ranks,
// and ranks holding no domain or no range boxes at all.
TEST(CoarseSearch, MatchesBruteForceOracle)
{
    int rank = 0, np = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    constexpr int kSeeds        = 25;
    constexpr int nDomainGlobal = 40;
    constexpr int nRangeGlobal  = 60;

    for (int seed = 0; seed < kSeeds; ++seed)
    {
        // Every rank draws from the SAME stream and keeps only what it owns, so the oracle can be
        // rebuilt locally without communicating the global input.
        std::mt19937 gen(1234 + seed);
        std::uniform_real_distribution<double> pos(0.0, 10.0);
        std::uniform_real_distribution<double> ext(0.05, 2.5);
        std::uniform_int_distribution<int> owner(0, np - 1);

        std::vector<Aabb<double>> gDom(nDomainGlobal), gRan(nRangeGlobal);
        std::vector<int> gDomOwner(nDomainGlobal), gRanOwner(nRangeGlobal);
        auto draw = [&](Aabb<double>& b)
        {
            for (int d = 0; d < 3; ++d)
            {
                double c = pos(gen), h = ext(gen);
                b.lo[d] = c - h;
                b.hi[d] = c + h;
            }
        };
        for (int i = 0; i < nDomainGlobal; ++i) { draw(gDom[i]); gDomOwner[i] = owner(gen); }
        for (int j = 0; j < nRangeGlobal; ++j) { draw(gRan[j]); gRanOwner[j] = owner(gen); }

        std::vector<Aabb<double>> myDom, myRan;
        std::vector<uint64_t> myDomId, myRanId;
        for (int i = 0; i < nDomainGlobal; ++i)
            if (gDomOwner[i] == rank) { myDom.push_back(gDom[i]); myDomId.push_back(i); }
        for (int j = 0; j < nRangeGlobal; ++j)
            if (gRanOwner[j] == rank) { myRan.push_back(gRan[j]); myRanId.push_back(1000 + j); }

        std::vector<int> peers;
        for (int r = 0; r < np; ++r) peers.push_back(r);

        std::vector<SearchPair> got;
        coarseSearchHost(myDom, myDomId, myRan, myRanId, peers, MPI_COMM_WORLD, got);

        // enforceSymmetry defaults to true, as in STK: I keep a pair if I own EITHER side.
        std::vector<SearchPair> want;
        for (int i = 0; i < nDomainGlobal; ++i)
            for (int j = 0; j < nRangeGlobal; ++j)
            {
                if (gDomOwner[i] != rank && gRanOwner[j] != rank) continue;
                if (gDom[i].overlaps(gRan[j]))
                    want.push_back(SearchPair{IdentProc{uint64_t(i), gDomOwner[i]},
                                              IdentProc{uint64_t(1000 + j), gRanOwner[j]}});
            }
        std::sort(want.begin(), want.end());
        want.erase(std::unique(want.begin(), want.end()), want.end());

        ASSERT_EQ(got.size(), want.size()) << "seed " << seed;
        EXPECT_TRUE(std::equal(got.begin(), got.end(), want.begin())) << "seed " << seed;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    mars::Env env(argc, argv);
    return RUN_ALL_TESTS();
}

// STK parity, checked against Trilinos develop rather than from memory:
//   1. OVERLAP PREDICATE -- identical. stk_search/BoundingBox.hpp tests disjointness as
//      `amax[d] < bmin[d] || bmax[d] < amin[d]`, exactly what Aabb::overlaps computes. Note this
//      is CLOSED: at amax == bmin the strict `<` is false, so touching boxes DO intersect. That is
//      what a conformal seam needs, where opposing faces share an edge exactly.
//   2. SYMMETRY -- matched. stk_search/CoarseSearch.hpp takes `enforceSearchResultSymmetry`,
//      DEFAULT TRUE, and CommonSearchUtil.hpp's communicate_vector then sends each pair to the
//      range box's processor as well as the domain box's. Our `enforceSymmetry` mirrors it, same
//      default. The DG-IP assembly needs it: master rows carry opposing-element columns and vice
//      versa, so both sides must know about the pair.
//   3. SORTING -- we are stricter on purpose. STK's `sortSearchResults` defaults to FALSE; we
//      always sort and dedup, so the result cannot depend on message arrival order. That is what
//      lets the device path be diffed against the host one element-for-element.
//   4. Self-rank pairs are included on both sides.
//
// An earlier version of this comment claimed the routing was done by a `communicateRangeBoxInfo`
// parameter. There is no such parameter; it is `enforceSearchResultSymmetry`. That note came from
// the design doc, not the source -- which is why this block now cites files.
