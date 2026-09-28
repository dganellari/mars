// Locality of the cornerstone decomposition on a hex mesh that is one cell thick.
//
// Each hex element is one cornerstone particle, placed at its corner with the smallest SFC key (as MARS does) or at
// its centroid. After Domain::sync, every mesh node is given to the rank whose SFC range holds the node position.
// That rank should hold elements next to the node; an FEM code then needs little extra communication to complete the
// element star of each node it owns. Per run, rank 0 prints:
//   max owner distance: largest distance, in cells, from a node to the nearest element assigned to its owner
//   far nodes:          nodes whose owner holds no assigned element within 2 cells
//   missing star:       (node, element) incidences where the node's owner holds the element neither as assigned
//                       element nor as halo; these must be sent to the owner before it can assemble its rows
//   halos:              halo particles over all ranks after sync
//
// Only public cornerstone API is used. Build against a cornerstone checkout and run with several rank counts:
//   mpicxx -std=c++20 -O2 -DNDEBUG -I <cornerstone>/include thin_mesh_locality.cpp -o thin_mesh_locality
//   mpirun -np 16 ./thin_mesh_locality 320 32 1 10 1 0.06 corner

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "cstone/domain/domain.hpp"

using KeyType  = uint64_t;
using RealType = double;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    const int n[3]         = {argc > 1 ? std::atoi(argv[1]) : 320, argc > 2 ? std::atoi(argv[2]) : 32,
                              argc > 3 ? std::atoi(argv[3]) : 1};
    const RealType L[3]    = {argc > 4 ? std::atof(argv[4]) : 10.0, argc > 5 ? std::atof(argv[5]) : 1.0,
                              argc > 6 ? std::atof(argv[6]) : 0.06};
    const bool centroid    = argc > 7 && std::string(argv[7]) == "centroid";
    const RealType d[3]    = {L[0] / n[0], L[1] / n[1], L[2] / n[2]};
    const long numElements = long(n[0]) * n[1] * n[2];

    // Smallest-key corner, with keys from a box padded by 5% like the MARS mesh box
    cstone::Box<RealType> meshBox(-0.05 * L[0], 1.05 * L[0], -0.05 * L[1], 1.05 * L[1], -0.05 * L[2], 1.05 * L[2]);
    auto nodeId = [&](int i, int j, int k) { return (KeyType(i) * (n[1] + 1) + j) * (n[2] + 1) + k; };

    // Contiguous element blocks per rank, as a mesh reader would give them
    std::vector<RealType> x, y, z, h;
    std::array<std::vector<KeyType>, 8> corners;
    for (long e = rank * numElements / numRanks; e < (rank + 1) * numElements / numRanks; ++e)
    {
        int i = int(e / (long(n[1]) * n[2])), j = int((e / n[2]) % n[1]), k = int(e % n[2]);
        KeyType minKey = ~KeyType(0);
        int rep        = 0;
        for (int c = 0; c < 8; ++c)
        {
            int ci = i + (c & 1), cj = j + ((c >> 1) & 1), ck = k + ((c >> 2) & 1);
            corners[c].push_back(nodeId(ci, cj, ck));
            KeyType key = cstone::sfc3D<cstone::SfcKind<KeyType>>(ci * d[0], cj * d[1], ck * d[2], meshBox);
            if (key < minKey) { minKey = key, rep = c; }
        }
        RealType off[3] = {RealType(rep & 1), RealType((rep >> 1) & 1), RealType((rep >> 2) & 1)};
        if (centroid) { off[0] = off[1] = off[2] = 0.5; }
        x.push_back((i + off[0]) * d[0]);
        y.push_back((j + off[1]) * d[1]);
        z.push_back((k + off[2]) * d[2]);
        h.push_back((d[0] + d[1] + d[2]) / 3);
    }

    cstone::Domain<KeyType, RealType> domain(cstone::execution::Cpu{}, rank, numRanks, 64, 8, 0.5, MPI_COMM_WORLD,
                                             meshBox);
    std::vector<KeyType> keys(x.size());
    std::vector<RealType> s1, s2, s3;
    std::vector<KeyType> k1, k2, k3, order;
    auto props   = std::tie(corners[0], corners[1], corners[2], corners[3], corners[4], corners[5], corners[6],
                            corners[7]);
    auto scratch = std::tie(s1, s2, s3, k1, k2, k3, order);
    domain.sync(keys, x, y, z, h, props, scratch);
    domain.exchangeHalos(props, s1, s2);

    // Owner of a point: the rank whose SFC range holds the key of the point, clamped into the domain box
    std::vector<KeyType> starts(numRanks);
    KeyType myStart = domain.assignmentStart();
    MPI_Allgather(&myStart, 1, MPI_UINT64_T, starts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
    const auto& box = domain.box();
    auto owner      = [&](int i, int j, int k)
    {
        RealType p[3] = {std::clamp(i * d[0], box.xmin(), box.xmax()), std::clamp(j * d[1], box.ymin(), box.ymax()),
                         std::clamp(k * d[2], box.zmin(), box.zmax())};
        KeyType key   = cstone::sfc3D<cstone::SfcKind<KeyType>>(p[0], p[1], p[2], box);
        return int(std::upper_bound(starts.begin(), starts.end(), key) - starts.begin()) - 1;
    };

    // Cells assigned to this rank, and how many local elements (assigned or halo) touch each node
    std::vector<char> myCell(numElements, 0);
    for (size_t e = domain.startIndex(); e < domain.endIndex(); ++e)
    {
        KeyType id = corners[0][e];
        int k = int(id % (n[2] + 1)), j = int((id / (n[2] + 1)) % (n[1] + 1)), i = int(id / ((n[1] + 1) * (n[2] + 1)));
        myCell[(long(i) * n[1] + j) * n[2] + k] = 1;
    }
    std::vector<int> starCount((n[0] + 1L) * (n[1] + 1) * (n[2] + 1), 0);
    for (size_t e = 0; e < x.size(); ++e)
        for (int c = 0; c < 8; ++c)
            starCount[corners[c][e]]++;

    long farNodes = 0, missingStar = 0, maxDistance = 0;
    const int maxSearch = 64;
    for (int i = 0; i <= n[0]; ++i)
        for (int j = 0; j <= n[1]; ++j)
            for (int k = 0; k <= n[2]; ++k)
            {
                if (owner(i, j, k) != rank) { continue; }
                auto side = [](int a, int m) { return (a > 0) + (a < m); };
                missingStar += side(i, n[0]) * side(j, n[1]) * side(k, n[2]) - starCount[nodeId(i, j, k)];

                // Chebyshev distance in cells from the node to the nearest assigned cell; 0 = the node touches it
                int distance = maxSearch;
                for (int r = 0; r < maxSearch && distance == maxSearch; ++r)
                    for (int ci = i - 1 - r; ci <= i + r && distance == maxSearch; ++ci)
                        for (int cj = j - 1 - r; cj <= j + r && distance == maxSearch; ++cj)
                            for (int ck = k - 1 - r; ck <= k + r && distance == maxSearch; ++ck)
                            {
                                if (ci < 0 || cj < 0 || ck < 0 || ci >= n[0] || cj >= n[1] || ck >= n[2]) { continue; }
                                if (myCell[(long(ci) * n[1] + cj) * n[2] + ck]) { distance = r; }
                            }
                maxDistance = std::max(maxDistance, long(distance));
                farNodes += distance > 2;
            }

    long halos = long(domain.nParticlesWithHalos() - domain.nParticles());
    long sums[3] = {farNodes, missingStar, halos};
    MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maxDistance, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::printf("%dx%dx%d cells, box %gx%gx%g, %s, %d ranks: max owner distance %s%ld cells, far nodes %ld of %ld, "
                    "missing star %ld of %ld incidences, halos %ld\n",
                    n[0], n[1], n[2], L[0], L[1], L[2], centroid ? "centroid" : "corner", numRanks,
                    maxDistance >= maxSearch ? ">=" : "", maxDistance, sums[0], long(starCount.size()), sums[1],
                    8 * numElements, sums[2]);
        std::printf("  box after sync: [%g, %g] x [%g, %g] x [%g, %g]\n", box.xmin(), box.xmax(), box.ymin(), box.ymax(),
                    box.zmin(), box.zmax());
    }
    MPI_Finalize();
    return 0;
}
