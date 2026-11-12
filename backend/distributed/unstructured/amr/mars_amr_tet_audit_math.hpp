#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_AUDIT_HD __host__ __device__
#else
#define MARS_AUDIT_HD
#endif

namespace mars::amr
{
MARS_AUDIT_HD inline uint64_t audit_mix(uint64_t x)
{
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

template<class KeyType>
struct TetAuditRecord
{
    KeyType keys[4]{};
    double volume         = 0;
    double roundoff_scale = 0;
    uint64_t invalid      = 0;
};

struct TetAuditTotals
{
    uint64_t count = 0, hash_a = 0, hash_b = 0, invalid = 0, inverted = 0;
    double volume = 0, absolute_volume = 0, min_abs_volume = std::numeric_limits<double>::max();
};

MARS_AUDIT_HD inline double audit_tet_volume(const double (&p)[4][3])
{
    double a[3], b[3], c[3];
    for (int d = 0; d < 3; ++d)
    {
        a[d] = p[1][d] - p[0][d];
        b[d] = p[2][d] - p[0][d];
        c[d] = p[3][d] - p[0][d];
    }
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])) /
           6.0;
}

template<class KeyType>
MARS_AUDIT_HD TetAuditRecord<KeyType> audit_record(const KeyType (&keys)[4], const double (&p)[4][3])
{
    TetAuditRecord<KeyType> r;
    r.volume      = audit_tet_volume(p);
    r.invalid     = !std::isfinite(r.volume) || r.volume == 0;
    double extent = 0, coordinate = 0;
    for (int i = 0; i < 4; ++i)
    {
        r.keys[i] = keys[i];
        for (int d = 0; d < 3; ++d)
        {
            r.invalid |= !std::isfinite(p[i][d]);
            coordinate = ::fabs(p[i][d]) > coordinate ? ::fabs(p[i][d]) : coordinate;
            for (int j = 0; j < i; ++j)
                extent = ::fabs(p[i][d] - p[j][d]) > extent ? ::fabs(p[i][d] - p[j][d]) : extent;
        }
    }
    r.roundoff_scale = (coordinate > extent ? coordinate : extent) * extent * extent;
    for (int i = 1; i < 4; ++i)
        for (int j = i; j > 0 && r.keys[j] < r.keys[j - 1]; --j)
        {
            KeyType t     = r.keys[j];
            r.keys[j]     = r.keys[j - 1];
            r.keys[j - 1] = t;
        }
    for (int i = 1; i < 4; ++i)
        r.invalid |= r.keys[i] == r.keys[i - 1];
    return r;
}

struct TetAuditFold
{
    MARS_AUDIT_HD TetAuditTotals operator()(TetAuditTotals a, TetAuditTotals b) const
    {
        a.count += b.count;
        a.hash_a += b.hash_a;
        a.hash_b += b.hash_b;
        a.invalid += b.invalid;
        a.inverted += b.inverted;
        a.volume += b.volume;
        a.absolute_volume += b.absolute_volume;
        a.min_abs_volume = a.min_abs_volume < b.min_abs_volume ? a.min_abs_volume : b.min_abs_volume;
        return a;
    }
};

template<class KeyType>
MARS_AUDIT_HD bool audit_red_children_valid(const TetAuditRecord<KeyType>& parent,
                                            const TetAuditRecord<KeyType>* children,
                                            size_t size,
                                            size_t first,
                                            bool marked)
{
    const size_t count = marked ? 8 : 1;
    if (parent.invalid || first > size || count > size - first) return false;
    const double expected = parent.volume / double(count);
    for (size_t j = 0; j < count; ++j)
    {
        const auto& child = children[first + j];
        const double tolerance =
            128 * std::numeric_limits<double>::epsilon() * (parent.roundoff_scale + child.roundoff_scale);
        if (child.invalid || parent.volume * child.volume <= 0 || ::fabs(child.volume - expected) > tolerance)
            return false;
        for (size_t k = 0; k < j; ++k)
        {
            bool duplicate = true;
            for (int vertex = 0; vertex < 4; ++vertex)
                duplicate &= child.keys[vertex] == children[first + k].keys[vertex];
            if (duplicate) return false;
        }
    }
    return true;
}

template<class KeyType>
struct TetAuditValue
{
    MARS_AUDIT_HD TetAuditTotals operator()(const TetAuditRecord<KeyType>& r) const
    {
        TetAuditTotals t;
        t.count    = 1;
        t.invalid  = r.invalid;
        uint64_t a = UINT64_C(0x243f6a8885a308d3), b = UINT64_C(0x13198a2e03707344);
        for (int i = 0; i < 4; ++i)
        {
            a = audit_mix(a ^ uint64_t(r.keys[i]));
            b = audit_mix(b + uint64_t(r.keys[i]) + uint64_t(i));
        }
        t.hash_a   = a;
        t.hash_b   = b;
        t.inverted = r.volume < 0;
        if (std::isfinite(r.volume))
        {
            t.volume          = r.volume;
            t.absolute_volume = ::fabs(r.volume);
            t.min_abs_volume  = ::fabs(r.volume);
        }
        return t;
    }
};

inline bool audit_same_cells(const TetAuditTotals& a, const TetAuditTotals& b)
{
    return a.count == b.count && a.hash_a == b.hash_a && a.hash_b == b.hash_b && !a.invalid && !b.invalid;
}

// Count/identity comparisons are exact. This bound covers summation of positive cell volumes in double.
inline double audit_volume_tolerance(const TetAuditTotals& a, const TetAuditTotals& b)
{
    double scale = a.absolute_volume > b.absolute_volume ? a.absolute_volume : b.absolute_volume;
    return 16 * std::numeric_limits<double>::epsilon() * double(a.count + b.count + 1) * scale;
}
} // namespace mars::amr

#undef MARS_AUDIT_HD
