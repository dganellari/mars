#include "cstone/domain/domain.hpp"
#include "backend/distributed/unstructured/mars_sfc_ownership.hpp"

#include <array>
#include <iostream>
#include <map>
#include <vector>

namespace
{
template<class KeyType>
void check_layout(int rank, int ranks)
{
    cstone::Box<double> box(-.05, 1.05);
    constexpr int n             = 4;
    constexpr int pattern[5][4] = {{1, 2, 4, 7}, {0, 1, 2, 4}, {3, 1, 2, 7}, {5, 1, 4, 7}, {6, 2, 4, 7}};
    std::vector<double> x, y, z, h;
    std::array<std::vector<KeyType>, 4> conn;
    std::map<KeyType, int> expected_stars;
    int element = 0;
    for (int iz = 0; iz < n; ++iz)
        for (int iy = 0; iy < n; ++iy)
            for (int ix = 0; ix < n; ++ix)
                for (const auto& cell : pattern)
                {
                    std::array<KeyType, 4> keys;
                    std::array<std::array<double, 3>, 4> points;
                    for (int c = 0; c < 4; ++c)
                    {
                        int corner = cell[c] ^ ((ix + iy + iz) & 1);
                        points[c]  = {double(ix + (corner & 1)) / n, double(iy + ((corner >> 1) & 1)) / n,
                                      double(iz + ((corner >> 2) & 1)) / n};
                        keys[c] =
                            cstone::sfc3D<cstone::SfcKind<KeyType>>(points[c][0], points[c][1], points[c][2], box);
                        ++expected_stars[keys[c]];
                    }
                    if (element++ % ranks != rank) continue;
                    int representative = std::min_element(keys.begin(), keys.end()) - keys.begin();
                    x.push_back(points[representative][0]);
                    y.push_back(points[representative][1]);
                    z.push_back(points[representative][2]);
                    double shortest = 1;
                    for (int a = 0; a < 4; ++a)
                        for (int b = a + 1; b < 4; ++b)
                        {
                            double distance = 0;
                            for (int d = 0; d < 3; ++d)
                                distance += (points[a][d] - points[b][d]) * (points[a][d] - points[b][d]);
                            shortest = std::min(shortest, std::sqrt(distance));
                        }
                    h.push_back(shortest);
                    for (int c = 0; c < 4; ++c)
                        conn[c].push_back(keys[c]);
                }

    cstone::Domain<KeyType, double> domain(cstone::execution::cpu, rank, ranks, 8, 1, .5, MPI_COMM_WORLD, box);
    std::vector<KeyType> keys(x.size()), k1, k2, k3, k4;
    std::vector<double> s1, s2, s3;
    auto properties = std::tie(conn[0], conn[1], conn[2], conn[3]);
    auto scratch    = std::tie(s1, s2, s3, k1, k2, k3, k4);
    domain.sync(keys, x, y, z, h, properties, scratch);

    int invalid = 0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        bool owned    = i >= domain.startIndex() && i < domain.endIndex();
        int suppliers = owned;
        for (const auto& range : domain.incomingHaloIndices())
            suppliers += i >= size_t(range.start()) && i < size_t(range.end());
        invalid += suppliers != 1;
        // A missing receive must not pass because old buffer contents look like valid connectivity.
        if (!owned)
            for (auto& column : conn)
                column[i] = ~KeyType(0);
    }
    domain.exchangeHalos(properties, s1, s2);

    std::vector<KeyType> bounds(ranks + 1);
    for (int r = 0; r <= ranks; ++r)
        bounds[r] = domain.assignment()[r];
    mars::SfcNodeOwner<KeyType, double> owner{box, domain.box(), bounds.data(), ranks};
    std::map<KeyType, int> stars;
    for (const auto& column : conn)
        for (KeyType key : column)
        {
            invalid += !expected_stars.count(key);
            ++stars[key];
        }
    for (auto [key, count] : expected_stars)
        if (owner(key) == rank) invalid += stars[key] != count;
    if (invalid) std::cerr << "rank " << rank << ": invalid halo entries/stars=" << invalid << '\n';
    MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (invalid) throw std::runtime_error("Cornerstone left unfilled halo slots or incomplete node stars");
    if (rank == 0)
        std::cout << "Tet domain layout " << sizeof(KeyType) * 8 << "-bit keys: halo coverage and node stars PASS\n";
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    char host[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name(host, &length);
    std::cout << "CPU rank=" << rank << " host=" << std::string(host, length) << '\n';
    try
    {
        check_layout<uint32_t>(rank, ranks);
        check_layout<uint64_t>(rank, ranks);
    }
    catch (const std::exception& error)
    {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
