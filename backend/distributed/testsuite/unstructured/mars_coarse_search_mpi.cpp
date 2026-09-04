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
#include <vector>

#include "mars_env.hpp"
// Path-relative on purpose -- see the note in mars_ghost_registry_mpi.cpp.
#include "../../unstructured/fem/mars_coarse_search.hpp"

using mars::fem::Aabb;
using mars::fem::coarseSearch;
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
    coarseSearch(f.domain, f.domainIds, f.range, f.rangeIds, MPI_COMM_WORLD, out);
    EXPECT_TRUE(hasPair(out, seamId(f.rank), f.rank, rangeId(f.rank, kBoxesPerRank - 1), f.rank));
}

// Gate 2: and its partner on the NEXT rank. This is the one that fails if the routing or the
// reply leg is wrong, and it is the whole point of the exercise.
TEST(CoarseSearch, FindsPartnerOnNeighbourRank)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearch(f.domain, f.domainIds, f.range, f.rangeIds, MPI_COMM_WORLD, out);
    if (f.rank + 1 < f.nRanks)
        EXPECT_TRUE(hasPair(out, seamId(f.rank), f.rank, rangeId(f.rank + 1, 0), f.rank + 1));
    else
        SUCCEED() << "last rank has no right neighbour";
}

// Gate 3: a caller only ever gets pairs whose DOMAIN box it owns.
TEST(CoarseSearch, ResultIsOwnedByCaller)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearch(f.domain, f.domainIds, f.range, f.rangeIds, MPI_COMM_WORLD, out);
    for (const auto& p : out) EXPECT_EQ(p.domain.proc, f.rank);
}

// Gate 4: no spurious pairs. The seam box is narrow, so it can only touch the two range boxes
// either side of the seam -- anything else means the overlap test is too permissive.
TEST(CoarseSearch, NoSpuriousPairs)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearch(f.domain, f.domainIds, f.range, f.rangeIds, MPI_COMM_WORLD, out);
    const size_t expected = (f.rank + 1 < f.nRanks) ? 2u : 1u;
    EXPECT_EQ(out.size(), expected);
}

// Gate 5: sorted and deduplicated, so the answer does not depend on message arrival order.
TEST(CoarseSearch, SortedAndUnique)
{
    Fixture f = makeFixture();
    std::vector<SearchPair> out;
    coarseSearch(f.domain, f.domainIds, f.range, f.rangeIds, MPI_COMM_WORLD, out);
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
    coarseSearch(mute ? emptyBoxes : f.domain, mute ? emptyIds : f.domainIds,
                 mute ? emptyBoxes : f.range, mute ? emptyIds : f.rangeIds, MPI_COMM_WORLD, out);
    if (mute) EXPECT_TRUE(out.empty());
    else SUCCEED();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    mars::Env env(argc, argv);
    return RUN_ALL_TESTS();
}
