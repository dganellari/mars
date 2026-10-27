// Host check of the block-structured DSS bodies (backend/distributed/unstructured/
// solvers/mars_cellwise_block_dss.hpp): an M^3 lattice of blocks, each block in a random
// proper rotation (so the skeleton sees every face frame and reversed edges), random
// values on every copy. The passes and the skeleton, run thread by thread in launch
// order, must give every copy the sum over all copies of its node, found here by the
// node's global coordinates. The skeleton tables are built by brute force in the format
// of build_topology (face frame codes from the corner keys, reversed bit on edges).
//
// Build and run from the MARS root:
//   c++ -std=c++20 -O2 -I. tests/cellwise_block_dss_host_check.cpp -o /tmp/block_dss_check && /tmp/block_dss_check
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <vector>

#include "backend/distributed/unstructured/solvers/mars_cellwise_block_dss.hpp"

namespace cw = mars::cellwise;

namespace {

// A proper rotation of the block frame: global axis perm[d] = local axis d, reversed
// when flip[d].
struct Rotation {
    int perm[3], flip[3];
    void apply(const int (&u)[3], int extent, int (&g)[3]) const
    {
        for (int d = 0; d < 3; ++d) g[perm[d]] = flip[d] ? extent - u[d] : u[d];
    }
};

std::vector<Rotation> proper_rotations()
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
int check(int M, unsigned seed)
{
    using S = cw::BlockShape<P, B>;
    constexpr int Pm = S::Pm;
    const long long blocks = (long long)M * M * M, N = blocks * S::B3 * S::n3;
    const std::vector<Rotation> rots = proper_rotations();
    std::mt19937 rng(seed);
    std::vector<Rotation> rot(blocks);
    for (auto& r : rot) r = rots[rng() % rots.size()];
    auto origin = [&](long long blk, int (&o)[3]) {
        o[0] = (int)(blk / ((long long)M * M));
        o[1] = (int)((blk / M) % M);
        o[2] = (int)(blk % M);
    };

    // The skeleton: block corner c's key is its global lattice vertex.
    auto corner_key = [&](long long blk, int c) {
        int b[3], g[3], o[3];
        cw::corner_bits(c, b[0], b[1], b[2]);
        rot[blk].apply(b, 1, g);
        origin(blk, o);
        return ((unsigned long long)(o[0] + g[0]) * (M + 1) + o[1] + g[1]) * (M + 1) + o[2] + g[2];
    };
    std::vector<int> face_nbr(blocks * 6, -1);
    std::vector<unsigned char> face_code(blocks * 6, 1 << 6);
    std::map<std::array<unsigned long long, 4>, std::vector<int>> faces;
    std::map<std::pair<unsigned long long, unsigned long long>, std::vector<int>> edges;
    std::map<unsigned long long, std::vector<int>> verts;
    for (long long blk = 0; blk < blocks; ++blk) {
        for (int f = 0; f < 6; ++f) {
            std::array<unsigned long long, 4> k;
            for (int s = 0; s < 4; ++s) k[s] = corner_key(blk, cw::face_corner(f, s));
            std::sort(k.begin(), k.end());
            faces[k].push_back((int)(blk * 6 + f));
        }
        for (int e = 0; e < 12; ++e) {
            const unsigned long long a = corner_key(blk, cw::edge_corner(e, 0)), b = corner_key(blk, cw::edge_corner(e, 1));
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
        face_nbr[v[0]] = v[1];
        face_nbr[v[1]] = v[0];
        face_code[v[0]] = (unsigned char)(ca | cb << 3);
        face_code[v[1]] = (unsigned char)(cb | ca << 3);
    }
    std::vector<int> edge_off{0}, edge_ent, vert_off{0}, vert_ent;
    for (auto& [k, v] : edges) {
        edge_ent.insert(edge_ent.end(), v.begin(), v.end());
        edge_off.push_back((int)edge_ent.size());
    }
    for (auto& [k, v] : verts) {
        vert_ent.insert(vert_ent.end(), v.begin(), v.end());
        vert_off.push_back((int)vert_ent.size());
    }
    const cw::BlockSkeleton sk{face_nbr.data(), edge_off.data(), edge_ent.data(), vert_off.data(), vert_ent.data(),
                               face_code.data(),  blocks,          (long long)edges.size(), (long long)verts.size()};

    // Values, and the assembled sum of every node by its global coordinates.
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    std::vector<double> in(N), out(N, std::nan(""));
    for (auto& x : in) x = val(rng);
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
    for (long long t = 0; t < N; ++t) total[node_of(t)] += in[t];

    for (long long t = 0; t < N; ++t) cw::block_pass_z<P, B>(in.data(), out.data(), t);
    if constexpr (B > 1) {
        const long long pairs = blocks * cw::block_pairs_per_block<P, B>;
        for (long long t = 0; t < pairs; ++t) cw::block_pass_xy<P, B, 0>(out.data(), t);
        for (long long t = 0; t < pairs; ++t) cw::block_pass_xy<P, B, 1>(out.data(), t);
    }
    const long long threads = cw::block_skeleton_threads<P, B>(sk);
    for (long long t = 0; t < threads; ++t) cw::block_skeleton_node<P, B>(out.data(), sk, t);

    double err = 0.0;
    for (long long t = 0; t < N; ++t) err = std::fmax(err, std::fabs(out[t] - total[node_of(t)]));
    if (std::isnan(err)) err = INFINITY;
    const bool ok = err <= 1e-13;
    printf("block DSS p=%d B=%d, %d^3 rotated blocks (%lld copies, %zu nodes): max|err| = %.2e  %s\n", P, B, M, N,
           total.size(), err, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

}  // namespace

int main()
{
    int fails = 0;
    fails += check<1, 1>(3, 1);
    fails += check<1, 2>(3, 2);
    fails += check<1, 4>(2, 3);
    fails += check<2, 1>(3, 4);
    fails += check<2, 3>(3, 5);
    fails += check<3, 2>(3, 6);
    fails += check<3, 4>(2, 7);
    fails += check<7, 2>(2, 8);
    return fails;
}
