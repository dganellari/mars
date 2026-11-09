#pragma once

#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace amr_test
{
using KeyType  = uint64_t;
using RealType = double;
using Point    = std::array<RealType, 3>;
using Edge     = std::array<KeyType, 2>;
using Face     = std::array<Point, 3>;

struct Tet
{
    std::array<KeyType, 4> keys{};
    std::array<Point, 4> points{};
    int rank = 0;
};

struct Node
{
    KeyType key = 0;
    Point point{};
    int owner = -1;
};

void require(bool condition, const char* message);
RealType determinant(const std::array<Point, 4>& points);
RealType check_mesh(const std::vector<Tet>& global);
RealType affine(const Point& point);
Edge longest_edge(const Tet& tet);
Face tet_face(const Tet& tet, int opposite);
int boundary_plane(const Face& face);
void select_device(int rank);
void run_gpu_checks(int rank, int ranks, bool staged);

template<class T>
std::vector<T> gather(const std::vector<T>& local)
{
    static_assert(std::is_trivially_copyable_v<T>);
    int ranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    int bytes = int(local.size() * sizeof(T));
    std::vector<int> counts(ranks), offsets(ranks + 1, 0);
    MPI_Allgather(&bytes, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int r = 0; r < ranks; ++r)
        offsets[r + 1] = offsets[r] + counts[r];
    require(offsets.back() % sizeof(T) == 0, "Invalid oracle gather size");
    std::vector<T> result(offsets.back() / sizeof(T));
    MPI_Allgatherv(local.data(), bytes, MPI_BYTE, result.data(), counts.data(), offsets.data(), MPI_BYTE,
                   MPI_COMM_WORLD);
    return result;
}

} // namespace amr_test
