// CPU gate for the pass-generated DSS kernel (dss_gpu_pipeline.py --emulate):
// runs the SAME kernel body the GPU gets with one thread per lane, 32 per
// element, and compares every written copy BIT FOR BIT with the formula of the
// hand-written MARS kernel (for_each_copy_sum in
// backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp), written
// out here on the host with MARS's own reference-hex functions.
//
// The tables are random but valid as indices: face pairs, frame codes, Dirichlet
// bits, reversed edges, stars of 1 to 9 entries, and copies on ghost elements.
// They need not describe a real mesh: the gate is about indexing and summation
// order, which a mesh would exercise only partly (a cube has no 5-edge stars,
// for instance). Real meshes are the GPU check against the hand kernel.
//
// Two passes: every element in a shuffled order, then a subset, whose other
// copies must stay untouched (NaN).
//
// Build + run: python3 test/dss_gpu_pipeline.py [p] --emulate [--emu-tsan]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <random>
#include <thread>
#include <vector>

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"
#include "emu_warp_runtime.h"

#ifndef MIR_DSS_P
#define MIR_DSS_P 7
#endif
#define MIR_CAT2(a, b) a##b
#define MIR_CAT(a, b) MIR_CAT2(a, b)
#define MIR_DSS_FN MIR_CAT(_mlir_ciface_dss_p, MIR_DSS_P)

using D64 = Desc<double>;
using D32 = Desc<int32_t>;
using D8 = Desc<int8_t>;
extern "C" void MIR_DSS_FN(D64*, D64*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, D64*);

namespace {

namespace cw = mars::cellwise;
constexpr int P = MIR_DSS_P, n = P + 1, n3 = n * n * n, F = (P - 1) * (P - 1), Ed = P - 1;

template <typename T>
Desc<T> desc(std::vector<T>& v)
{
    return {v.data(), v.data(), 0, {(int64_t)v.size()}, {1}};
}

struct Mesh {
    long long E, G;   // local and ghost elements
    std::vector<double> values, ghost;
    std::vector<int32_t> face_nbr, edge_of, edge_off, edge_ent, vert_of, vert_off, vert_ent;
    std::vector<int8_t> face_code;

    double value(long long i) const { return i < (long long)values.size() ? values[i] : ghost[i - values.size()]; }
};

Mesh random_mesh(long long E, long long G, unsigned seed)
{
    std::mt19937 rng(seed);
    auto uni = [&](long long lo, long long hi) { return std::uniform_int_distribution<long long>(lo, hi)(rng); };
    Mesh m{E, G};
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    m.values.resize(E * n3);
    m.ghost.resize(G * n3);
    for (auto& v : m.values) v = val(rng);
    for (auto& v : m.ghost) v = val(rng);
    const long long copies = E + G;
    for (long long i = 0; i < E * 6; ++i) {
        m.face_nbr.push_back(uni(0, 4) == 0 ? -1 : (int32_t)uni(0, copies * 6 - 1));
        m.face_code.push_back((int8_t)uni(0, 127));
    }
    auto star = [&](long long entities, long long slots, std::vector<int32_t>& off, std::vector<int32_t>& ent,
                    bool reversal) {
        off.push_back(0);
        for (long long g = 0; g < entities; ++g) {
            const long long size = uni(1, 9);
            for (long long k = 0; k < size; ++k) {
                int32_t c = (int32_t)uni(0, copies * slots - 1);
                if (reversal && uni(0, 1)) c = (int32_t)((uint32_t)c | 0x80000000u);
                ent.push_back(c);
            }
            off.push_back((int32_t)ent.size());
        }
    };
    const long long edges = 3 * E + 1, verts = E + 1;
    star(edges, 12, m.edge_off, m.edge_ent, true);
    star(verts, 8, m.vert_off, m.vert_ent, false);
    for (long long i = 0; i < E * 12; ++i) {
        int32_t g = (int32_t)uni(0, edges - 1);
        if (uni(0, 1)) g = (int32_t)((uint32_t)g | 0x80000000u);
        m.edge_of.push_back(g);
    }
    for (long long i = 0; i < E * 8; ++i) m.vert_of.push_back((int32_t)uni(0, verts - 1));
    return m;
}

double star_sum(const Mesh& m, const std::vector<int32_t>& ent, int begin, int end, auto index)
{
    double s = 0.0;
    for (int k = begin; k < end; ++k) s += m.value(index(ent[k]));
    return s;
}

// for_each_copy_sum<P>, one element at a time.
void reference(const Mesh& m, const std::vector<int32_t>& list, std::vector<double>& out)
{
    int face_l[6 * (F > 0 ? F : 1)], face_map[64 * (F > 0 ? F : 1)], edge_l[12 * n], corner_l[8];
    if constexpr (P > 1) {
        for (int i = 0; i < 64 * F; ++i) {
            const int q = i % F, code = i / F;
            int I, J, i2, j2;
            cw::face_to_canonical<P>(code & 7, 1 + q / (P - 1), 1 + q % (P - 1), I, J);
            cw::canonical_to_face<P>(code >> 3, I, J, i2, j2);
            face_map[i] = (i2 - 1) * (P - 1) + j2 - 1;
        }
        for (int i = 0; i < 6 * F; ++i) face_l[i] = cw::face_node<P>(i / F, 1 + (i % F) / (P - 1), 1 + (i % F) % (P - 1));
    }
    for (int i = 0; i < 12 * n; ++i) edge_l[i] = cw::edge_node<P>(i / n, i % n);
    for (int c = 0; c < 8; ++c) corner_l[c] = cw::corner_node<P>(c);

    std::vector<double> stage(n3);
    for (int32_t e : list) {
        const long long base = (long long)e * n3;
        for (int l = 0; l < n3; ++l) stage[l] = m.values[base + l];
        for (int p = 0; p < 6 * F; ++p) {
            const int f = p / F, nbr = m.face_nbr[e * 6 + f], code = (uint8_t)m.face_code[e * 6 + f];
            if (nbr < 0) continue;
            const int q = face_map[(code & 63) * F + p - f * F];
            stage[face_l[p]] += m.value((long long)(nbr / 6) * n3 + face_l[(nbr % 6) * F + q]);
        }
        for (int p = 0; p < 12 * Ed; ++p) {
            const int k = p / Ed, t = 1 + p - k * Ed, own = m.edge_of[e * 12 + k], g = own & 0x7fffffff;
            const int tc = own < 0 ? P - t : t;
            stage[edge_l[k * n + t]] = star_sum(m, m.edge_ent, m.edge_off[g], m.edge_off[g + 1], [&](int32_t c) {
                const int mm = c & 0x7fffffff;
                return (long long)(mm / 12) * n3 + edge_l[(mm % 12) * n + (c < 0 ? P - tc : tc)];
            });
        }
        for (int c = 0; c < 8; ++c) {
            const int g = m.vert_of[e * 8 + c];
            stage[corner_l[c]] = star_sum(m, m.vert_ent, m.vert_off[g], m.vert_off[g + 1],
                                          [&](int32_t v) { return (long long)(v / 8) * n3 + corner_l[v % 8]; });
        }
        for (int l = 0; l < n3; ++l) out[base + l] = stage[l];
    }
}

// The emulated kernel on `list`: grid = list size, 32 lanes per block, one
// block at a time (the workgroup buffer is one process-wide global here).
void emulated(Mesh& m, std::vector<int32_t>& list, std::vector<double>& out)
{
    D64 dv = desc(m.values), dg = desc(m.ghost), dout = desc(out);
    D32 dl = desc(list), dfn = desc(m.face_nbr), deo = desc(m.edge_of), deoff = desc(m.edge_off),
        deent = desc(m.edge_ent), dvo = desc(m.vert_of), dvoff = desc(m.vert_off), dvent = desc(m.vert_ent);
    D8 dfc = desc(m.face_code);
    Warp w;
    warp = &w;
    for (long long b = 0; b < (long long)list.size(); ++b) {
        std::vector<std::thread> lanes;
        for (int L = 0; L < kLanes; ++L)
            lanes.emplace_back([&, L, b] {
                tlLane = L;
                tlBlock = b;
                MIR_DSS_FN(&dv, &dg, &dl, &dfn, &dfc, &deo, &deoff, &deent, &dvo, &dvoff, &dvent, &dout);
            });
        for (auto& t : lanes) t.join();
    }
}

// Copies that differ in their bits (NaN equals NaN: unwritten on both sides).
long long mismatches(const std::vector<double>& a, const std::vector<double>& b)
{
    long long bad = 0;
    for (size_t i = 0; i < a.size(); ++i)
        if (memcmp(&a[i], &b[i], sizeof(double)) != 0) ++bad;
    return bad;
}

}  // namespace

int main(int argc, char** argv)
{
    const long long E = argc > 1 ? atoll(argv[1]) : 12, G = argc > 2 ? atoll(argv[2]) : 5;
    Mesh m = random_mesh(E, G, 2026);
    std::vector<int32_t> all(E);
    for (long long e = 0; e < E; ++e) all[e] = (int32_t)e;
    std::shuffle(all.begin(), all.end(), std::mt19937(7));
    std::vector<int32_t> some(all.begin(), all.begin() + E / 3);

    int failed = 0;
    for (auto* list : {&all, &some}) {
        std::vector<double> ref(E * n3, std::nan("")), got(E * n3, std::nan(""));
        reference(m, *list, ref);
        emulated(m, *list, got);
        const long long bad = mismatches(ref, got);
        printf("  p=%d, %zu of %lld elements listed, %lld ghosts: %lld of %lld copies differ\n", P, list->size(), E, G,
               bad, E * n3);
        failed |= bad != 0;
    }
    printf("DSS WARP EMULATION (p=%d, bit-exact vs the hand kernel's formula): %s\n", P, failed ? "FAIL" : "PASS");
    return failed;
}
