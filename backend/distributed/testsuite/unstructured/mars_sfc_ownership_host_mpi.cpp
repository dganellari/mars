// Multi-rank gates for SFC node ownership and element-star completion (mars_sfc_ownership.hpp).
//
// Host-only + MPI: runs the same owner function, per-element destinations and sparse exchange as the GPU sync,
// on cornerstone's CPU domain with a structured hex cube laid out like ElementDomain lays out elements (element
// position = corner with the smallest node key, h = local edge length). The property checked is the one assembly
// needs: every rank holds all elements that touch the nodes it owns, and every node has exactly one owner.
//
// Run with several rank counts; 1 rank exercises almost nothing:
//   mpirun -np 4 ./mars_sfc_ownership_host_mpi_test

#include <gtest/gtest.h>

#include <map>
#include <tuple>
#include <vector>

#include "cstone/domain/domain.hpp"
#include "mars_env.hpp"
// Path-relative on purpose, see mars_ghost_registry_host_mpi.cpp
#include "../../unstructured/mars_sfc_ownership.hpp"

using namespace mars;

namespace
{

template<class KeyType, class RealType>
struct StarResult
{
    long ownedNodes      = 0; // nodes owned by this rank
    long incompleteStars = 0; // owned nodes missing an incident element
    std::vector<KeyType> ownedKeys;
};

// n^3 hexes on an anisotropic box; each rank checks the element star of every node it owns
template<class KeyType, class RealType>
StarResult<KeyType, RealType> runCube(int n, RealType hScale, bool complete)
{
    int rank = 0, numRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    const RealType L[3] = {2.0, 1.0, 1.5};
    const RealType d[3] = {L[0] / n, L[1] / n, L[2] / n};
    // padded like ElementDomain's bounding box
    cstone::Box<RealType> meshBox(-0.05 * L[0], 1.05 * L[0], -0.05 * L[1], 1.05 * L[1], -0.05 * L[2], 1.05 * L[2]);

    auto nodeKey = [&](int i, int j, int k)
    { return KeyType(cstone::sfc3D<cstone::SfcKind<KeyType>>(i * d[0], j * d[1], k * d[2], meshBox)); };

    std::map<KeyType, std::array<int, 3>> nodeIndex;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            for (int k = 0; k <= n; ++k)
                nodeIndex[nodeKey(i, j, k)] = {i, j, k};
    EXPECT_EQ(nodeIndex.size(), size_t(n + 1) * (n + 1) * (n + 1)) << "node keys not unique";

    std::vector<RealType> x, y, z, h;
    std::array<std::vector<KeyType>, 8> conn;
    for (int e = rank; e < n * n * n; e += numRanks)
    {
        int i = e / (n * n), j = (e / n) % n, k = e % n;
        KeyType minKey = std::numeric_limits<KeyType>::max();
        int rep        = 0;
        for (int c = 0; c < 8; ++c)
        {
            KeyType key = nodeKey(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
            conn[c].push_back(key);
            if (key < minKey) { minKey = key, rep = c; }
        }
        x.push_back((i + (rep & 1)) * d[0]);
        y.push_back((j + ((rep >> 1) & 1)) * d[1]);
        z.push_back((k + ((rep >> 2) & 1)) * d[2]);
        h.push_back(hScale * std::min({d[0], d[1], d[2]}));
    }

    cstone::Domain<KeyType, RealType> domain(rank, numRanks, 64, 8, 0.5, meshBox);
    std::vector<KeyType> keys(x.size());
    std::vector<RealType> s1, s2, s3;
    std::vector<KeyType> k1, k2, k3, order;
    auto props   = std::tie(conn[0], conn[1], conn[2], conn[3], conn[4], conn[5], conn[6], conn[7]);
    auto scratch = std::tie(s1, s2, s3, k1, k2, k3, order);
    domain.sync(keys, x, y, z, h, props, scratch);

    std::vector<KeyType> bounds(numRanks + 1);
    for (int r = 0; r <= numRanks; ++r)
    {
        bounds[r] = domain.assignment()[r];
    }
    SfcNodeOwner<KeyType, RealType> owner{meshBox, domain.box(), bounds.data(), numRanks};

    if (complete)
    {
        auto ranges = flattenHaloSendRanges(domain.outgoingHaloIndices(), numRanks);
        HaloSendRanges sent{ranges.rankOffsets.data(), ranges.start.data(), ranges.end.data()};

        std::vector<std::pair<int, KeyType>> pushes;
        for (size_t e = domain.startIndex(); e < domain.endIndex(); ++e)
        {
            KeyType corners[8];
            for (int c = 0; c < 8; ++c)
            {
                corners[c] = conn[c][e];
            }
            int dests[8];
            int numDests = missingStarRanks<8>(corners, cstone::LocalIndex(e), owner, sent, rank, dests);
            for (int i = 0; i < numDests; ++i)
            {
                pushes.emplace_back(dests[i], keys[e]);
            }
        }
        std::sort(pushes.begin(), pushes.end());

        std::vector<int> dests, counts;
        std::vector<KeyType> sendKeys;
        for (auto [dest, key] : pushes)
        {
            if (dests.empty() || dests.back() != dest)
            {
                dests.push_back(dest);
                counts.push_back(0);
            }
            counts.back()++;
            sendKeys.push_back(key);
        }

        int any = !pushes.empty();
        MPI_Allreduce(MPI_IN_PLACE, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any)
        {
            std::vector<KeyType> haloKeys;
            sparseExchange<KeyType>(dests, counts, sendKeys.data(), haloKeys, 0x5354, MPI_COMM_WORLD);
            domain.addHalos(haloKeys, keys, x, y, z, h, props, scratch);
        }
    }
    domain.exchangeHalos(props, s1, s2);

    // element star size of every node present on this rank
    std::map<KeyType, int> star;
    for (size_t e = 0; e < x.size(); ++e)
    {
        for (int c = 0; c < 8; ++c)
        {
            star[conn[c][e]]++;
        }
    }

    StarResult<KeyType, RealType> result;
    for (auto [key, count] : star)
    {
        if (owner(key) != rank) { continue; }
        auto [i, j, k]    = nodeIndex.at(key);
        auto side         = [n](int a) { return (a > 0) + (a < n); };
        int expectedCount = side(i) * side(j) * side(k);
        result.ownedNodes++;
        result.ownedKeys.push_back(key);
        result.incompleteStars += count != expectedCount;
        EXPECT_LE(count, expectedCount) << "element counted twice at node " << key;
    }
    return result;
}

template<class KeyType>
std::vector<KeyType> gatherAll(const std::vector<KeyType>& local)
{
    int numRanks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int count = int(local.size());
    std::vector<int> counts(numRanks), displs(numRanks + 1, 0);
    MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int r = 0; r < numRanks; ++r)
    {
        displs[r + 1] = displs[r] + counts[r];
    }
    std::vector<KeyType> all(displs.back());
    MPI_Allgatherv(local.data(), count, ::MpiType<KeyType>{}, all.data(), counts.data(), displs.data(),
                   ::MpiType<KeyType>{}, MPI_COMM_WORLD);
    return all;
}

template<class KeyType, class RealType>
void checkCompleteStars(int n, RealType hScale)
{
    auto result = runCube<KeyType, RealType>(n, hScale, true);
    EXPECT_EQ(result.incompleteStars, 0);

    // every node owned by exactly one rank, and that rank holds it
    auto all = gatherAll(result.ownedKeys);
    std::sort(all.begin(), all.end());
    EXPECT_EQ(std::adjacent_find(all.begin(), all.end()), all.end()) << "node owned twice";
    EXPECT_EQ(all.size(), size_t(n + 1) * (n + 1) * (n + 1)) << "node without owner";
}

} // namespace

TEST(SfcOwnership, completeStarsDefaultReach)
{
    checkCompleteStars<uint64_t, double>(12, 1.0);
    checkCompleteStars<unsigned, float>(12, 1.0);
}

// A reach far below the element size leaves most seam stars to the completion
TEST(SfcOwnership, completeStarsShortReach)
{
    checkCompleteStars<uint64_t, double>(12, 0.05);
    checkCompleteStars<unsigned, float>(12, 0.05);
}

// Negative control: the distance search alone does not complete the stars at short reach, so the gates above test
// the completion and not the search
TEST(SfcOwnership, searchAloneLeavesStarsIncomplete)
{
    int numRanks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    if (numRanks == 1) { GTEST_SKIP() << "no seams on one rank"; }

    auto result     = runCube<uint64_t, double>(12, 0.05, false);
    long incomplete = result.incompleteStars;
    MPI_Allreduce(MPI_IN_PLACE, &incomplete, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_GT(incomplete, 0);
}

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
