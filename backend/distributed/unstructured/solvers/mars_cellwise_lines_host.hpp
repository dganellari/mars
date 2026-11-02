#pragma once
// The element lines of the node-line smoother (mars_cellwise_lines.hpp), built on the host from
// one rank's elements: no CUDA or MPI dependency, so tests/cellwise_lines_host_check.cpp runs it
// on the host (scripts/test_node_lines.py).

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <utility>
#include <vector>

namespace mars {
namespace cellwise {

// ---- element lines (host) --------------------------------------------------------------

// One rank's local elements as lines, longest first: line l is entries start[l] .. start[l+1)
// of elem and code. code bits 0-1: the element's thin axis; bit 2: the line runs from node P
// down to 0 along it; bits 3-5: the map from the line's transverse node position, given in its
// first element's frame, to this element's (bit 3 swaps the two transverse axes, bits 4 and 5
// reverse the first and the second).
// The lateral neighbours of entry m, by its transverse axes k1 < k2 (the two other than its thin
// axis). Face f (0, 1: across k1 at its low / high side; 2, 3: across k2): face_elem[m * 4 + f]
// the local element there (-1: none on this rank), face_code its frame: bits 0-1 its axis along
// the line, bits 2-3 its axis across the face, bit 4 its side of that axis, bit 5 its
// coordinate along the line runs the other way, bit 6 so does its coordinate along the face's
// other axis. Edge g along the line (bit 0: the high side of k1, bit 1 of k2): the other local
// elements of the edge are entries edge_start[m * 4 + g] .. edge_start[m * 4 + g + 1) of
// edge_elem and edge_code: bits 0-1 their axis along the edge, bits 2 and 3 the edge's corner
// bits on their two other axes (ascending), bit 4 their coordinate along the edge runs the
// other way.
struct ElementLines {
    std::vector<int> start{0}, elem;
    std::vector<unsigned char> code;
    std::vector<int> face_elem, edge_start{0}, edge_elem;
    std::vector<unsigned char> face_code, edge_code;
    long long lines() const { return (long long)start.size() - 1; }
};

namespace line_detail {

MARS_CW_HD inline void transverse(int code, int P, int& i, int& j)
{
    if (code & 1) {
        const int t = i;
        i = j;
        j = t;
    }
    if (code & 2) i = P - i;
    if (code & 4) j = P - j;
}

// The two reference axes other than a, ascending.
MARS_CW_HD inline void other_axes(int a, int& k1, int& k2)
{
    k1 = a == 0 ? 1 : 0;
    k2 = a == 2 ? 1 : 2;
}

// The code c with transverse(c) = transverse(outer) after transverse(inner).
inline int compose(int outer, int inner)
{
    for (int c = 0; c < 8; ++c) {
        bool same = true;
        for (int i = 0; i <= 1 && same; ++i)
            for (int j = 0; j <= 1 && same; ++j) {
                int a = i, b = j, x = i, y = j;
                transverse(inner, 1, a, b);
                transverse(outer, 1, a, b);
                transverse(c, 1, x, y);
                same = a == x && b == y;
            }
        if (same) return c;
    }
    return -1;
}

}  // namespace line_detail

// corners[(e * 8 + c) * 3 + a]: coordinate a of corner c (corner_bits order) of local element
// e; vid[e * 8 + c]: an id of that corner's vertex, equal ids for equal vertices. As
// element_lines of the reference: elements in order, a line started at every element with no
// chain neighbour on one side, closed loops opened at their lowest element, then sorted by
// length (stable).
inline ElementLines build_element_lines(const std::vector<double>& corners, const std::vector<long long>& vid,
                                        long long E)
{
    using line_detail::other_axes;
    std::vector<int> axis(E);
    for (long long e = 0; e < E; ++e) {
        double len[3];
        for (int k = 0; k < 3; ++k) {
            double sum = 0.0;
            for (int c = 0; c < 8; ++c) {
                int b[3];
                corner_bits(c, b[0], b[1], b[2]);
                if (b[k]) continue;
                b[k] = 1;
                const int d = corner_from_bits(b[0], b[1], b[2]);
                double s2 = 0.0;
                for (int a = 0; a < 3; ++a) {
                    const double t = corners[(e * 8 + d) * 3 + a] - corners[(e * 8 + c) * 3 + a];
                    s2 += t * t;
                }
                sum += std::sqrt(s2);
            }
            len[k] = sum / 4.0;
        }
        int a = 0;
        if (len[1] < len[a]) a = 1;
        if (len[2] < len[a]) a = 2;
        axis[e] = a;
    }
    // The faces across each element's thin axis, sorted by their vertices: a face listed twice
    // chains its two elements.
    struct Face {
        std::array<long long, 4> v;
        long long e;
        int side;
    };
    std::vector<Face> faces;
    faces.reserve(2 * E);
    for (long long e = 0; e < E; ++e)
        for (int side = 0; side < 2; ++side) {
            Face f{{}, e, side};
            int m = 0;
            for (int c = 0; c < 8; ++c) {
                int b[3];
                corner_bits(c, b[0], b[1], b[2]);
                if (b[axis[e]] == side) f.v[m++] = vid[e * 8 + c];
            }
            std::sort(f.v.begin(), f.v.end());
            faces.push_back(f);
        }
    std::sort(faces.begin(), faces.end(), [](const Face& x, const Face& y) {
        return x.v != y.v ? x.v < y.v : (x.e != y.e ? x.e < y.e : x.side < y.side);
    });
    std::vector<long long> nxt(2 * E, -1);   // nxt[e * 2 + side]: neighbour * 2 + its side
    for (size_t i = 0; i + 1 < faces.size(); ++i)
        if (faces[i].v == faces[i + 1].v) {
            nxt[faces[i].e * 2 + faces[i].side] = faces[i + 1].e * 2 + faces[i + 1].side;
            nxt[faces[i + 1].e * 2 + faces[i + 1].side] = faces[i].e * 2 + faces[i].side;
            ++i;
        }
    // Walk: a line entered through side `entry` (0: the axis' low face, 1: its high face).
    std::vector<char> seen(E, 0);
    std::vector<std::vector<std::pair<long long, int>>> lines;
    auto walk = [&](long long e, int entry) {
        std::vector<std::pair<long long, int>> line;
        while (!seen[e]) {
            seen[e] = 1;
            line.emplace_back(e, entry);
            const long long n = nxt[e * 2 + (1 - entry)];
            if (n < 0) break;
            e = n / 2;
            entry = (int)(n % 2);
        }
        lines.push_back(std::move(line));
    };
    for (long long e = 0; e < E; ++e)
        for (int side = 0; side < 2; ++side)
            if (!seen[e] && nxt[e * 2 + side] < 0) walk(e, side);
    for (long long e = 0; e < E; ++e)
        if (!seen[e]) walk(e, 0);
    std::stable_sort(lines.begin(), lines.end(), [](const auto& x, const auto& y) { return x.size() > y.size(); });

    ElementLines out;
    out.elem.reserve(E);
    out.code.reserve(E);
    for (const auto& line : lines) {
        int cum = 0;
        for (size_t m = 0; m < line.size(); ++m) {
            const long long e = line[m].first;
            const int entry = line[m].second, a = axis[e];
            if (m > 0) {
                // The previous element's exit face and this one's entry face are the same 4
                // vertices; their transverse corner positions give the map between the frames.
                const long long p = line[m - 1].first;
                const int pa = axis[p], pexit = 1 - line[m - 1].second;
                int pk1, pk2, k1, k2;
                other_axes(pa, pk1, pk2);
                other_axes(a, k1, k2);
                auto corner_vid = [&](long long el, int ax, int bit, int t1, int t2, int i, int j) {
                    int b[3];
                    b[ax] = bit;
                    b[t1] = i;
                    b[t2] = j;
                    return vid[el * 8 + corner_from_bits(b[0], b[1], b[2])];
                };
                int rel = -1;
                for (int c = 0; c < 8 && rel < 0; ++c) {
                    bool ok = true;
                    for (int i = 0; i <= 1 && ok; ++i)
                        for (int j = 0; j <= 1 && ok; ++j) {
                            int x = i, y = j;
                            line_detail::transverse(c, 1, x, y);
                            ok = corner_vid(p, pa, pexit, pk1, pk2, i, j) == corner_vid(e, a, entry, k1, k2, x, y);
                        }
                    if (ok) rel = c;
                }
                cum = line_detail::compose(rel, cum);
            }
            out.elem.push_back((int)e);
            out.code.push_back((unsigned char)(a | (entry << 2) | (cum << 3)));
        }
        out.start.push_back((int)out.elem.size());
    }

    // The lateral neighbours of every entry, matched by vertex ids.
    auto bits_of = [&](long long el, long long v, int (&b)[3]) {
        for (int c = 0; c < 8; ++c)
            if (vid[el * 8 + c] == v) {
                corner_bits(c, b[0], b[1], b[2]);
                return;
            }
    };
    auto corner_vid = [&](long long el, int (&b)[3]) { return vid[el * 8 + corner_from_bits(b[0], b[1], b[2])]; };
    struct AnyFace {
        std::array<long long, 4> v;
        long long e;
    };
    std::vector<AnyFace> all;
    all.reserve(6 * E);
    for (long long e = 0; e < E; ++e)
        for (int k = 0; k < 3; ++k)
            for (int side = 0; side < 2; ++side) {
                AnyFace f{{}, e};
                int m = 0;
                for (int c = 0; c < 8; ++c) {
                    int b[3];
                    corner_bits(c, b[0], b[1], b[2]);
                    if (b[k] == side) f.v[m++] = vid[e * 8 + c];
                }
                std::sort(f.v.begin(), f.v.end());
                all.push_back(f);
            }
    std::sort(all.begin(), all.end(),
              [](const AnyFace& x, const AnyFace& y) { return x.v != y.v ? x.v < y.v : x.e < y.e; });
    std::vector<std::pair<long long, long long>> vel;   // (vertex, element), sorted
    vel.reserve(8 * E);
    for (long long e = 0; e < E; ++e)
        for (int c = 0; c < 8; ++c) vel.emplace_back(vid[e * 8 + c], e);
    std::sort(vel.begin(), vel.end());
    auto elements_of = [&](long long v) {
        std::vector<long long> r;
        for (auto it = std::lower_bound(vel.begin(), vel.end(), std::make_pair(v, (long long)-1));
             it != vel.end() && it->first == v; ++it)
            r.push_back(it->second);
        return r;
    };
    const size_t M = out.elem.size();
    out.face_elem.assign(4 * M, -1);
    out.face_code.assign(4 * M, 0);
    for (size_t m = 0; m < M; ++m) {
        const long long e = out.elem[m];
        const int a = out.code[m] & 3;
        int k1, k2;
        other_axes(a, k1, k2);
        for (int f = 0; f < 4; ++f) {
            const int k = f < 2 ? k1 : k2, t = f < 2 ? k2 : k1, side = f & 1;
            int b00[3], b10[3], b01[3];
            b00[a] = 0, b00[t] = 0, b00[k] = side;
            b10[a] = 1, b10[t] = 0, b10[k] = side;
            b01[a] = 0, b01[t] = 1, b01[k] = side;
            const long long v00 = corner_vid(e, b00), v10 = corner_vid(e, b10), v01 = corner_vid(e, b01);
            AnyFace key{{}, -1};
            int q = 0;
            for (int c = 0; c < 8; ++c) {
                int b[3];
                corner_bits(c, b[0], b[1], b[2]);
                if (b[k] == side) key.v[q++] = vid[e * 8 + c];
            }
            std::sort(key.v.begin(), key.v.end());
            long long nb = -1;
            for (auto it = std::lower_bound(all.begin(), all.end(), key,
                                            [](const AnyFace& x, const AnyFace& y) { return x.v < y.v; });
                 it != all.end() && it->v == key.v; ++it)
                if (it->e != e) nb = it->e;
            if (nb < 0) continue;
            int c00[3], c10[3], c01[3];
            bits_of(nb, v00, c00);
            bits_of(nb, v10, c10);
            bits_of(nb, v01, c01);
            int an = 0, bn = 0;
            for (int d = 0; d < 3; ++d) {
                if (c00[d] != c10[d]) an = d;
                if (c00[d] != c01[d]) bn = d;
            }
            const int kn = 3 - an - bn;
            out.face_elem[m * 4 + f] = (int)nb;
            out.face_code[m * 4 + f] = (unsigned char)(an | kn << 2 | c00[kn] << 4 | c00[an] << 5 | c00[bn] << 6);
        }
        for (int g = 0; g < 4; ++g) {
            int b0[3], b1[3];
            b0[a] = 0, b0[k1] = g & 1, b0[k2] = g >> 1;
            b1[a] = 1, b1[k1] = g & 1, b1[k2] = g >> 1;
            const long long v0 = corner_vid(e, b0), v1 = corner_vid(e, b1);
            const std::vector<long long> s0 = elements_of(v0), s1 = elements_of(v1);
            std::vector<long long> star;
            std::set_intersection(s0.begin(), s0.end(), s1.begin(), s1.end(), std::back_inserter(star));
            for (long long nb : star) {
                if (nb == e) continue;
                int c0[3], c1[3];
                bits_of(nb, v0, c0);
                bits_of(nb, v1, c1);
                int an = 0;
                for (int d = 0; d < 3; ++d)
                    if (c0[d] != c1[d]) an = d;
                int o1, o2;
                other_axes(an, o1, o2);
                out.edge_elem.push_back((int)nb);
                out.edge_code.push_back((unsigned char)(an | c0[o1] << 2 | c0[o2] << 3 | c0[an] << 4));
            }
            out.edge_start.push_back((int)out.edge_elem.size());
        }
    }
    return out;
}

}  // namespace cellwise
}  // namespace mars
