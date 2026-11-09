#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_TET_RED_HD __host__ __device__
#else
#define MARS_TET_RED_HD
#endif

namespace mars::amr
{
// Nodes 0..3 are the parent vertices; 4..9 are midpoints 01,02,03,12,13,23.
// The octahedron diagonal is 02--13. Every child keeps the parent's orientation.
template<class KeyType>
MARS_TET_RED_HD void red_child_connectivity(int child, const KeyType* nodes, KeyType* output)
{
    constexpr int corners[8][4] = {{0, 4, 5, 6}, {4, 1, 7, 8}, {5, 7, 2, 9}, {6, 8, 9, 3},
                                   {4, 5, 6, 8}, {4, 5, 8, 7}, {5, 6, 8, 9}, {5, 8, 7, 9}};
    for (int c = 0; c < 4; ++c) output[c] = nodes[corners[child][c]];
}
} // namespace mars::amr

#undef MARS_TET_RED_HD
