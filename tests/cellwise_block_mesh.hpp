#pragma once
// A test mesh for the block-structured DSS, shared by its host check
// (tests/cellwise_block_dss_host_check.cpp) and the CPU gate of the MARSIR-generated
// kernels (marsir-mlir/tools/emu_block_dss.cpp): an M^3 lattice of blocks of B^3
// elements, each block in a random proper rotation (so the skeleton sees every face
// frame and reversed edges), random values on every copy, and the assembled sum of every
// copy's node, found by the node's global coordinates. The skeleton tables are built by
// brute force in the format of build_topology (face frame codes from the corner keys,
// reversed bit on edges).
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <random>
#include <vector>

#include "backend/distributed/unstructured/solvers/mars_cellwise_block_dss.hpp"

namespace block_test {

namespace cw = mars::cellwise;

// A proper rotation of the block frame: global axis perm[d] = local axis d, reversed
// when flip[d].
struct Rotation {
    int perm[3], flip[3];
    void apply(const int (&u)[3], int extent, int (&g)[3]) const
    {
        for (int d = 0; d < 3; ++d) g[perm[d]] = flip[d] ? extent - u[d] : u[d];
    }
};

inline std::vector<Rotation> proper_rotations()
{
    std::vector<Rotation> r;
    int p[3] = {0, 1, 2};
    do {
        const int inversions = (p[0] > p[1]) + (p[0] > p[2]) + (p[1] > p[2]);
        for (int f = 0; f < 8; ++f) {
            const int flips = (f & 1) + (f >> 1 & 1) + (f >> 2 & 1);
            if ((inversions + flips) % 2) continue;
            r.push_back({{p[0], p[1], p[2]}, {f & 1, f >> 1 & 1, f >> 2 & 1}});
        }
    } while (std::next_permutation(p, p + 3));
    return r;
}

template <int P, int B>
struct BlockMesh {
    using S = cw::BlockShape<P, B>;
    int M = 0;
    long long blocks = 0, N = 0;
    std::vector<int> face_nbr, edge_off{0}, edge_ent, vert_off{0}, vert_ent;
    std::vector<unsigned char> face_code;
    std::vector<double> in, sum;   // sum[t]: the assembled value of copy t's node

    cw::BlockSkeleton skeleton() const
    {
        return {face_nbr.data(), edge_off.data(), edge_ent.data(), vert_off.data(), vert_ent.data(), face_code.data(),
                blocks,          (long long)edge_off.size() - 1, (long long)vert_off.size() - 1};
    }

    // Largest |y - sum| over all copies; infinite when any value is not finite.
    double error(const std::vector<double>& y) const
    {
        double err = 0.0;
        for (long long t = 0; t < N; ++t) {
            if (!std::isfinite(y[t]) || !std::isfinite(sum[t])) return INFINITY;
            err = std::fmax(err, std::fabs(y[t] - sum[t]));
        }
        return err;
    }

    static BlockMesh build(int M, unsigned seed)
    {
        constexpr int Pm = S::Pm;
        BlockMesh m;
        m.M = M;
        m.blocks = (long long)M * M * M;
        m.N = m.blocks * S::B3 * S::n3;
        const std::vector<Rotation> rots = proper_rotations();
        std::mt19937 rng(seed);
        std::vector<Rotation> rot(m.blocks);
        for (auto& r : rot) r = rots[rng() % rots.size()];
        auto origin = [&](long long blk, int (&o)[3]) {
            o[0] = (int)(blk / ((long long)M * M));
            o[1] = (int)((blk / M) % M);
            o[2] = (int)(blk % M);
        };
        // Block corner c's key is its global lattice vertex.
        auto corner_key = [&](long long blk, int c) {
            int b[3], g[3], o[3];
            cw::corner_bits(c, b[0], b[1], b[2]);
            rot[blk].apply(b, 1, g);
            origin(blk, o);
            return ((unsigned long long)(o[0] + g[0]) * (M + 1) + o[1] + g[1]) * (M + 1) + o[2] + g[2];
        };
        m.face_nbr.assign(m.blocks * 6, -1);
        m.face_code.assign(m.blocks * 6, 1 << 6);
        std::map<std::array<unsigned long long, 4>, std::vector<int>> faces;
        std::map<std::pair<unsigned long long, unsigned long long>, std::vector<int>> edges;
        std::map<unsigned long long, std::vector<int>> verts;
        for (long long blk = 0; blk < m.blocks; ++blk) {
            for (int f = 0; f < 6; ++f) {
                std::array<unsigned long long, 4> k;
                for (int s = 0; s < 4; ++s) k[s] = corner_key(blk, cw::face_corner(f, s));
                std::sort(k.begin(), k.end());
                faces[k].push_back((int)(blk * 6 + f));
            }
            for (int e = 0; e < 12; ++e) {
                const unsigned long long a = corner_key(blk, cw::edge_corner(e, 0)),
                                         b = corner_key(blk, cw::edge_corner(e, 1));
                edges[{std::min(a, b), std::max(a, b)}].push_back((int)(blk * 12 + e) | (a > b ? (int)0x80000000u : 0));
            }
            for (int c = 0; c < 8; ++c) verts[corner_key(blk, c)].push_back((int)(blk * 8 + c));
        }
        for (auto& [k, v] : faces) {
            if (v.size() != 2) continue;
            unsigned long long ka[4], kb[4];
            for (int s = 0; s < 4; ++s) {
                ka[s] = corner_key(v[0] / 6, cw::face_corner(v[0] % 6, s));
                kb[s] = corner_key(v[1] / 6, cw::face_corner(v[1] % 6, s));
            }
            const int ca = cw::face_frame_code(ka), cb = cw::face_frame_code(kb);
            m.face_nbr[v[0]] = v[1];
            m.face_nbr[v[1]] = v[0];
            m.face_code[v[0]] = (unsigned char)(ca | cb << 3);
            m.face_code[v[1]] = (unsigned char)(cb | ca << 3);
        }
        for (auto& [k, v] : edges) {
            m.edge_ent.insert(m.edge_ent.end(), v.begin(), v.end());
            m.edge_off.push_back((int)m.edge_ent.size());
        }
        for (auto& [k, v] : verts) {
            m.vert_ent.insert(m.vert_ent.end(), v.begin(), v.end());
            m.vert_off.push_back((int)m.vert_ent.size());
        }

        std::uniform_real_distribution<double> val(-1.0, 1.0);
        m.in.resize(m.N);
        for (auto& x : m.in) x = val(rng);
        auto node_of = [&](long long t) {
            const long long e = t / S::n3, blk = e / S::B3;
            const int l = (int)(t % S::n3), le = (int)(e % S::B3);
            const int u[3] = {(le / (B * B)) * P + l / (S::n * S::n), ((le / B) % B) * P + (l / S::n) % S::n,
                              (le % B) * P + l % S::n};
            int g[3], o[3];
            rot[blk].apply(u, Pm, g);
            origin(blk, o);
            return std::array<long long, 3>{(long long)o[0] * Pm + g[0], (long long)o[1] * Pm + g[1],
                                            (long long)o[2] * Pm + g[2]};
        };
        std::map<std::array<long long, 3>, double> total;
        for (long long t = 0; t < m.N; ++t) total[node_of(t)] += m.in[t];
        m.sum.resize(m.N);
        for (long long t = 0; t < m.N; ++t) m.sum[t] = total[node_of(t)];
        return m;
    }
};

}  // namespace block_test
