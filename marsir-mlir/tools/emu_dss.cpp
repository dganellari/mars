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
// Two lists: every element in a shuffled order, then a subset, whose other
// copies must stay untouched (NaN). Each with one block per element and with 3
// blocks that stride over the list.
//
// -DMIR_DSS_PRECOND: the preconditioner kernel (mlir_ir.emit_dss_precond, Jacobi
// with <a, z> and <z, z>) instead, with random Dirichlet entities and counting
// bits: z must match unstructured_precondition_kernel's formula bit for bit, and
// each element's partials the kernel's reduction order (lane L sums the copies
// L, L + 32, ...; then an xor-shuffle tree over the 32 lanes).
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
#include "emu_cellwise_mesh.h"
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
#ifdef MIR_DSS_PRECOND
#define MIR_PRECOND_FN MIR_CAT(MIR_CAT(_mlir_ciface_precond_p, MIR_DSS_P), _az_zz)
extern "C" void MIR_PRECOND_FN(D64*, D64*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, D32*, D8*, D8*, D64*,
                               D64*, D64*, D64*);
#else
extern "C" void MIR_DSS_FN(D64*, D64*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, D64*);
#endif

namespace {

namespace cw = mars::cellwise;
constexpr int P = MIR_DSS_P, n = P + 1, n3 = n * n * n, F = (P - 1) * (P - 1), Ed = P - 1;

template <typename T>
Desc<T> desc(std::vector<T>& v)
{
    return {v.data(), v.data(), 0, {(int64_t)v.size()}, {1}};
}

using Mesh = CellwiseMesh<P>;

// for_each_copy_sum<P>, one element at a time.
void reference(const Mesh& m, const std::vector<int32_t>& list, std::vector<double>& out, std::vector<double>& partials)
{
    std::vector<long long> position(m.E, -1);   // place of each element in the list
    for (size_t w = 0; w < list.size(); ++w) position[list[w]] = (long long)w;
    std::vector<double> stage(n3);
    for (int32_t e : list) {
        const long long base = (long long)e * n3;
        m.stage(e, stage);
#ifdef MIR_DSS_PRECOND
        // unstructured_precondition_kernel's visit, in the generated kernel's
        // reduction order.
        const unsigned dmask = m.dirichlet_mask(e);
        double acc[32][2] = {};
        for (int l = 0; l < n3; ++l) {
            const int bit = Mesh::bit(l);
            const bool dir = (dmask >> bit) & 1u, counted = ((uint32_t)m.counting[e] >> bit) & 1u;
            const double v = dir ? 0.0 : stage[l] / m.diag[base + l];
            const double w = counted ? 1.0 : 0.0;
            out[base + l] = v;
            acc[l % 32][0] += m.a[base + l] * v * w;
            acc[l % 32][1] += v * v * w;
        }
        for (int off = 16; off >= 1; off /= 2) {
            double next[32][2];
            for (int L = 0; L < 32; ++L)
                for (int j = 0; j < 2; ++j) next[L][j] = acc[L][j] + acc[L ^ off][j];
            memcpy(acc, next, sizeof(acc));
        }
        for (int j = 0; j < 2; ++j) partials[(size_t)position[e] * 2 + j] = acc[0][j];
#else
        for (int l = 0; l < n3; ++l) out[base + l] = stage[l];
#endif
    }
}

// The emulated kernel on `list` with `grid` blocks of 32 lanes, one block at a
// time (the workgroup buffers are process-wide globals here).
void emulated(Mesh& m, std::vector<int32_t>& list, std::vector<double>& out, std::vector<double>& partials,
              long long grid)
{
    D64 dv = desc(m.values), dg = desc(m.ghost), dout = desc(out), dpart = desc(partials);
    D32 dcnt = desc(m.counting);
    D8 dedir = desc(m.edge_dir), dvdir = desc(m.vert_dir);
    D64 ddiag = desc(m.diag), da = desc(m.a);
    D32 dl = desc(list), dfn = desc(m.face_nbr), deo = desc(m.edge_of), deoff = desc(m.edge_off),
        deent = desc(m.edge_ent), dvo = desc(m.vert_of), dvoff = desc(m.vert_off), dvent = desc(m.vert_ent);
    D8 dfc = desc(m.face_code);
    Warp w;
    warp = &w;
    for (long long b = 0; b < grid; ++b) {
        std::vector<std::thread> lanes;
        for (int L = 0; L < kLanes; ++L)
            lanes.emplace_back([&, L, b] {
                tlLane = L;
                tlBlock = b;
                tlGrid = grid;
#ifdef MIR_DSS_PRECOND
                MIR_PRECOND_FN(&dv, &dg, &dl, &dfn, &dfc, &deo, &deoff, &deent, &dvo, &dvoff, &dvent, &dcnt, &dedir,
                               &dvdir, &ddiag, &da, &dout, &dpart);
#else
                (void)dpart, (void)dcnt, (void)dedir, (void)dvdir, (void)ddiag, (void)da;
                MIR_DSS_FN(&dv, &dg, &dl, &dfn, &dfc, &deo, &deoff, &deent, &dvo, &dvoff, &dvent, &dout);
#endif
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
    Mesh m = Mesh::random(E, G, 2026);
    std::vector<int32_t> all(E);
    for (long long e = 0; e < E; ++e) all[e] = (int32_t)e;
    std::shuffle(all.begin(), all.end(), std::mt19937(7));
    std::vector<int32_t> some(all.begin(), all.begin() + E / 3);

    int failed = 0;
    for (auto* list : {&all, &some})
        for (long long grid : {(long long)list->size(), 3LL}) {
            std::vector<double> ref(E * n3, std::nan("")), got(E * n3, std::nan(""));
            std::vector<double> ref_part(list->size() * 2, std::nan("")), got_part(list->size() * 2, std::nan(""));
            reference(m, *list, ref, ref_part);
            emulated(m, *list, got, got_part, grid);
            const long long bad = mismatches(ref, got);
#ifdef MIR_DSS_PRECOND
            const long long bad_part = mismatches(ref_part, got_part);
            printf("  p=%d, %zu of %lld elements listed, %lld ghosts, %lld blocks: %lld of %lld copies, %lld of %zu "
                   "partials differ\n", P, list->size(), E, G, grid, bad, E * n3, bad_part, ref_part.size());
            failed |= bad != 0 || bad_part != 0;
#else
            printf("  p=%d, %zu of %lld elements listed, %lld ghosts, %lld blocks: %lld of %lld copies differ\n", P,
                   list->size(), E, G, grid, bad, E * n3);
            failed |= bad != 0;
#endif
        }
#ifdef MIR_DSS_PRECOND
    const char* what = "PRECONDITIONER";
#else
    const char* what = "DSS";
#endif
    printf("%s WARP EMULATION (p=%d, bit-exact vs the hand kernel's formula): %s\n", what, P, failed ? "FAIL" : "PASS");
    return failed;
}
