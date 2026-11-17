#pragma once
#include "mars_amr_lepp_boundary.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_RED_BOUNDARY_HD __host__ __device__
#else
#define MARS_RED_BOUNDARY_HD
#endif

namespace mars::amr
{
// The refiner numbers new nodes after the old nodes, in sorted edge order.
template<class KeyType>
struct RedBoundaryMidpoint
{
    const uint64_t* edges;
    size_t count;
    KeyType first;
    MARS_RED_BOUNDARY_HD KeyType operator()(KeyType a, KeyType b) const
    {
        const uint64_t lo = a < b ? a : b, hi = a < b ? b : a;
        if (hi > UINT32_MAX) return ~KeyType(0);
        const uint64_t key = (hi << 32) | lo;
        size_t l = 0, r = count;
        while (l < r)
        {
            size_t m = l + (r - l) / 2;
            if (edges[m] < key)
                l = m + 1;
            else
                r = m;
        }
        return l < count && edges[l] == key ? first + l : ~KeyType(0);
    }
};

// Preserve labels by parent-face ancestry, using the refiner's actual midpoint indices.
template<class Execution, class KeyType>
typename Execution::template Vector<LeppFaceRecord<KeyType>> save_red_boundary(const typename Execution::template Vector<LeppTaggedFace>& labels,
                       std::array<const KeyType*, 4> conn,
                       const uint8_t* marks,
                       RedBoundaryMidpoint<KeyType> midpoint,
                       MPI_Comm comm)
{
    const size_t n = labels.size();
    typename Execution::template Vector<uint64_t> counts(n + 1), offsets(n + 1), invalid(n);
    auto count = Execution::data(counts), offset = Execution::data(offsets), bad = Execution::data(invalid);
    auto faces = Execution::data(labels);
    Execution::for_each(n + 1, [=] MARS_RED_BOUNDARY_HD(size_t i)
                        { count[i] = i < n ? (marks[faces[i].element] ? 4 : 1) : 0; });
    Execution::scan(count, n + 1, offset);
    typename Execution::template Vector<LeppFaceRecord<KeyType>> saved(Execution::read(offset + n));
    auto out = Execution::data(saved);
    Execution::for_each(n,
                        [=] MARS_RED_BOUNDARY_HD(size_t i)
                        {
                            KeyType v[3];
                            for (int c = 0, j = 0; c < 4; ++c)
                                if (c != faces[i].face) v[j++] = conn[c][faces[i].element];
                            bad[i] = 0;
                            if (count[i] == 1)
                            {
                                out[offset[i]] = {LeppFaceKey<KeyType>(v[0], v[1], v[2]), faces[i].tag};
                                return;
                            }
                            KeyType m[3] = {midpoint(v[0], v[1]), midpoint(v[1], v[2]), midpoint(v[2], v[0])};
                            for (int j = 0; j < 3; ++j)
                                bad[i] |= m[j] == ~KeyType(0);
                            if (bad[i]) return;
                            for (int j = 0; j < 3; ++j)
                                out[offset[i] + j] = {LeppFaceKey<KeyType>(v[j], m[j], m[(j + 2) % 3]), faces[i].tag};
                            out[offset[i] + 3] = {LeppFaceKey<KeyType>(m[0], m[1], m[2]), faces[i].tag};
                        });
    require_collectively(Execution::sum(bad, n) == 0, "Red boundary face has no emitted edge midpoint", comm);
    return saved;
}
} // namespace mars::amr
#undef MARS_RED_BOUNDARY_HD
