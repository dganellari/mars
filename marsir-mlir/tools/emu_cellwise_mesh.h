#pragma once
// Shared by the CPU gates of the generated cell-wise kernels (emu_dss.cpp,
// emu_pmg.cpp): random tables of the unstructured cell-wise topology that are valid
// as indices (face pairs, frame codes, Dirichlet bits, reversed edges, stars of 1 to
// 9 entries, copies on ghost elements; they need not describe a real mesh), and the
// formula of the hand-written MARS DSS (for_each_copy_sum in
// backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp) written out on
// the host with MARS's own reference-hex functions. Real: the value type (float for
// the single-precision kernels; the random draws are the same, rounded).
#include <cstdint>
#include <random>
#include <vector>

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"

template <int P, typename Real = double>
struct CellwiseMesh {
    static constexpr int n = P + 1, n3 = n * n * n, F = (P - 1) * (P - 1), Ed = P - 1;
    long long E = 0, G = 0;   // local and ghost elements
    std::vector<Real> values, ghost;
    std::vector<Real> minus;   // when set, local copy i reads values[i] - minus[i]
    std::vector<int32_t> face_nbr, edge_of, edge_off, edge_ent, vert_of, vert_off, vert_ent, counting;
    std::vector<int8_t> face_code, edge_dir, vert_dir;
    std::vector<Real> diag, a;

    Real value(long long i) const
    {
        if (i < (long long)values.size()) return minus.empty() ? values[i] : values[i] - minus[i];
        return ghost[i - values.size()];
    }

    static int bit(int l)
    {
        namespace cw = mars::cellwise;
        return cw::entity_bit(cw::local_entity<P>(l / (n * n), (l / n) % n, l % n));
    }

    static CellwiseMesh random(long long E, long long G, unsigned seed)
    {
        namespace cw = mars::cellwise;
        std::mt19937 rng(seed);
        auto uni = [&](long long lo, long long hi) { return std::uniform_int_distribution<long long>(lo, hi)(rng); };
        CellwiseMesh m;
        m.E = E;
        m.G = G;
        std::uniform_real_distribution<double> val(-1.0, 1.0);
        m.values.resize(E * n3);
        m.ghost.resize(G * n3);
        for (auto& v : m.values) v = (Real)val(rng);
        for (auto& v : m.ghost) v = (Real)val(rng);
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
        for (long long e = 0; e < E; ++e)
            m.counting.push_back((int32_t)((uint32_t)uni(0, 0xfffffff) | (1u << cw::kInteriorBit)));
        for (long long g = 0; g < edges; ++g) m.edge_dir.push_back((int8_t)(uni(0, 3) == 0 ? uni(1, 2) : 0));
        for (long long g = 0; g < verts; ++g) m.vert_dir.push_back((int8_t)(uni(0, 3) == 0 ? 1 : 0));
        std::uniform_real_distribution<double> pos(1.0, 2.0);
        for (long long i = 0; i < E * n3; ++i) m.diag.push_back((Real)pos(rng));
        for (long long i = 0; i < E * n3; ++i) m.a.push_back((Real)val(rng));
        return m;
    }

    // The DSS sum of every copy of element e, in the hand kernel's order.
    void stage(int32_t e, std::vector<Real>& out) const
    {
        namespace cw = mars::cellwise;
        static const Tables t;
        out.resize(n3);
        const long long base = (long long)e * n3;
        for (int l = 0; l < n3; ++l) out[l] = value(base + l);
        for (int p = 0; p < 6 * F; ++p) {
            const int f = p / F, nbr = face_nbr[e * 6 + f], code = (uint8_t)face_code[e * 6 + f];
            if (nbr < 0) continue;
            const int q = t.face_map[(code & 63) * F + p - f * F];
            out[t.face_l[p]] += value((long long)(nbr / 6) * n3 + t.face_l[(nbr % 6) * F + q]);
        }
        for (int p = 0; p < 12 * Ed; ++p) {
            const int k = p / Ed, tt = 1 + p - k * Ed, own = edge_of[e * 12 + k], g = own & 0x7fffffff;
            const int tc = own < 0 ? P - tt : tt;
            Real s = 0;
            for (int i = edge_off[g]; i < edge_off[g + 1]; ++i) {
                const int c = edge_ent[i], mm = c & 0x7fffffff;
                s += value((long long)(mm / 12) * n3 + t.edge_l[(mm % 12) * n + (c < 0 ? P - tc : tc)]);
            }
            out[t.edge_l[k * n + tt]] = s;
        }
        for (int c = 0; c < 8; ++c) {
            const int g = vert_of[e * 8 + c];
            Real s = 0;
            for (int i = vert_off[g]; i < vert_off[g + 1]; ++i) {
                const int v = vert_ent[i];
                s += value((long long)(v / 8) * n3 + t.corner_l[v % 8]);
            }
            out[t.corner_l[c]] = s;
        }
    }

    // Bit b set when entity b of element e is a Dirichlet entity (the hand kernel's ballot).
    unsigned dirichlet_mask(int32_t e) const
    {
        namespace cw = mars::cellwise;
        unsigned mask = 0;
        for (int f = 0; f < 6; ++f) mask |= (((uint8_t)face_code[e * 6 + f] >> 6) != 0 ? 1u : 0u) << f;
        for (int k = 0; k < 12; ++k)
            mask |= (edge_dir[edge_of[e * 12 + k] & 0x7fffffff] != 0 ? 1u : 0u) << (cw::kEdgeBit + k);
        for (int c = 0; c < 8; ++c) mask |= (vert_dir[vert_of[e * 8 + c]] != 0 ? 1u : 0u) << (cw::kVertexBit + c);
        return mask;
    }

private:
    struct Tables {   // build_hex_tables<P>
        int face_l[6 * (F > 0 ? F : 1)], face_map[64 * (F > 0 ? F : 1)], edge_l[12 * n], corner_l[8];
        Tables()
        {
            namespace cw = mars::cellwise;
            if constexpr (P > 1) {
                for (int i = 0; i < 64 * F; ++i) {
                    const int q = i % F, code = i / F;
                    int I, J, i2, j2;
                    cw::face_to_canonical<P>(code & 7, 1 + q / (P - 1), 1 + q % (P - 1), I, J);
                    cw::canonical_to_face<P>(code >> 3, I, J, i2, j2);
                    face_map[i] = (i2 - 1) * (P - 1) + j2 - 1;
                }
                for (int i = 0; i < 6 * F; ++i)
                    face_l[i] = cw::face_node<P>(i / F, 1 + (i % F) / (P - 1), 1 + (i % F) % (P - 1));
            }
            for (int i = 0; i < 12 * n; ++i) edge_l[i] = cw::edge_node<P>(i / n, i % n);
            for (int c = 0; c < 8; ++c) corner_l[c] = cw::corner_node<P>(c);
        }
    };
};

// The fixed-order reduction of the generated kernels over term(i), i in [lo, hi):
// lane L sums i = lo + L, lo + L + 32, ...; then an xor tree over the 32 lanes.
template <typename Term, typename Real = double>
Real lane_sum(long long lo, long long hi, Term term)
{
    Real acc[32] = {};
    for (long long i = lo; i < hi; ++i) acc[(i - lo) % 32] += term(i);
    for (int off = 16; off >= 1; off /= 2) {
        Real next[32];
        for (int L = 0; L < 32; ++L) next[L] = acc[L] + acc[L ^ off];
        for (int L = 0; L < 32; ++L) acc[L] = next[L];
    }
    return acc[0];
}
