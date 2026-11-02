// CPU gate for the pass-generated node-line smoother kernels (dss_gpu_pipeline.py p --lines
// --emulate): the same kernel bodies the GPU gets, one thread per lane, compared BIT FOR BIT
// with host formulas:
//   lines_residual0 / lines_residual  the DSS of b, or of b - y (for_each_copy_sum's sum), on
//                                     random tables with ghosts and Dirichlet entities
//                                     (emu_cellwise_mesh.h);
//   lines_update_first0 / first / next  z_t = DSS(z), zero on Dirichlet copies, then the
//                                     Chebyshev update of d and x (update_kernel's formula);
//   lines_solve                       the banded forward and backward solves of solve_kernel
//                                     on a lane-interleaved layout of blocks of varying
//                                     length over several warps, padded with zero-length
//                                     blocks to the kernel's 128 per thread block
//                                     (backend/distributed/unstructured/solvers/
//                                     mars_cellwise_lines.hpp).
// Built with -ffp-contract=off: the generated code fuses only its explicit math.fma (the solve),
// and the reference uses std::fma there.
//
// Build + run: python3 test/dss_gpu_pipeline.py 3 --lines --emulate [--emu-tsan]
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <thread>
#include <vector>

#include "emu_cellwise_mesh.h"
#include "emu_warp_runtime.h"

#ifndef MIR_DSS_P
#define MIR_DSS_P 3
#endif
constexpr int P = MIR_DSS_P;

#define CAT_(a, b) a##b
#define CAT(a, b) CAT_(a, b)
#define KP(name) CAT(name##_p, MIR_DSS_P)
#define CI(name) CAT(_mlir_ciface_, KP(name))
#define TH(name) CAT(KP(name), __threads)

using DR = Desc<double>;
using D64 = Desc<int64_t>;
using D32 = Desc<int32_t>;
using D8 = Desc<int8_t>;
using Residual0 = void(DR*, DR*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, DR*);
using Residual = void(DR*, DR*, DR*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, DR*);
using Update = void(DR*, DR*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, D32*, D8*, D8*, DR*, DR*, double,
                    double);
using Solve = void(DR*, D64*, D32*, D64*, DR*, DR*);
extern "C" {
Residual0 CI(lines_residual0);
Residual CI(lines_residual);
Update CI(lines_update_first0), CI(lines_update_first), CI(lines_update_next);
Solve CI(lines_solve);
int64_t TH(lines_residual0)(), TH(lines_update_first0)(), TH(lines_solve)();
}

namespace {

template <typename T>
Desc<T> desc(std::vector<T>& v)
{
    return {v.data(), v.data(), 0, {(int64_t)v.size()}, {1}};
}

template <typename F>
void emulate(int64_t threads, long long grid, F&& kernel)
{
    Warp w;
    w.bar.threads = (int)threads;
    warp = &w;
    for (long long b = 0; b < grid; ++b) {
        std::vector<std::thread> lanes;
        for (int L = 0; L < threads; ++L)
            lanes.emplace_back([&, L, b] {
                tlLane = L;
                tlBlock = b;
                tlGrid = grid;
                kernel();
            });
        for (auto& t : lanes) t.join();
    }
}

template <typename T>
long long differ(const std::vector<T>& a, const std::vector<T>& b)
{
    long long bad = 0;
    for (size_t i = 0; i < a.size(); ++i)
        if (memcmp(&a[i], &b[i], sizeof(T)) != 0) ++bad;
    return bad;
}

int failed = 0;
void report(const char* what, long long bad, size_t total)
{
    printf("  %-20s p=%d  %lld of %zu values differ\n", what, P, bad, total);
    failed |= bad != 0;
}

std::vector<double> random_vector(long long size, unsigned seed, double lo = -1.0, double hi = 1.0)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> val(lo, hi);
    std::vector<double> v(size);
    for (auto& x : v) x = val(rng);
    return v;
}

using Mesh = CellwiseMesh<P, double>;
constexpr int n3 = Mesh::n3;

struct Tables {
    D32 fn, eo, eoff, eent, vo, voff, vent, cnt;
    D8 fc, edir, vdir;
    explicit Tables(Mesh& m)
        : fn(desc(m.face_nbr)), eo(desc(m.edge_of)), eoff(desc(m.edge_off)), eent(desc(m.edge_ent)),
          vo(desc(m.vert_of)), voff(desc(m.vert_off)), vent(desc(m.vert_ent)), cnt(desc(m.counting)),
          fc(desc(m.face_code)), edir(desc(m.edge_dir)), vdir(desc(m.vert_dir))
    {
    }
};

void residual()
{
    Mesh m = Mesh::random(10, 4, 31 + P);
    std::vector<int32_t> list = {6, 0, 9, 2, 4};
    const std::vector<double> y = random_vector(m.E * n3, 41);
    for (const bool minus : {false, true}) {
        std::vector<double> s(m.E * n3, std::nan("")), ref = s, stage;
        m.minus = minus ? y : std::vector<double>();
        for (int32_t e : list) {
            m.stage(e, stage);
            for (int l = 0; l < n3; ++l) ref[(long long)e * n3 + l] = stage[l];
        }
        std::vector<double> yy = y;
        DR db = desc(m.values), dg = desc(m.ghost), dy = desc(yy), ds = desc(s);
        D32 dl = desc(list);
        Tables t(m);
        emulate(TH(lines_residual0)(), (long long)list.size(), [&] {
            if (minus)
                CI(lines_residual)(&db, &dg, &dy, &dl, &t.fn, &t.fc, &t.eo, &t.eoff, &t.eent, &t.vo, &t.voff, &t.vent,
                                   &ds);
            else
                CI(lines_residual0)(&db, &dg, &dl, &t.fn, &t.fc, &t.eo, &t.eoff, &t.eent, &t.vo, &t.voff, &t.vent, &ds);
        });
        report(minus ? "lines_residual" : "lines_residual0", differ(s, ref), s.size());
    }
}

void update()
{
    Mesh m = Mesh::random(10, 4, 57 + P);
    std::vector<int32_t> list = {3, 0, 9, 5, 7, 1};
    const std::vector<double> d0 = random_vector(m.E * n3, 12), x0 = random_vector(m.E * n3, 13);
    const double c1 = 0.37, c2 = 1.9;
    const char* names[3] = {"lines_update_first0", "lines_update_first", "lines_update_next"};
    Update* kernels[3] = {CI(lines_update_first0), CI(lines_update_first), CI(lines_update_next)};
    for (int kind = 0; kind < 3; ++kind) {
        std::vector<double> d = d0, x = x0, dref = d0, xref = x0, stage;
        for (int32_t e : list) {
            m.stage(e, stage);
            const unsigned dmask = m.dirichlet_mask(e);
            for (int l = 0; l < n3; ++l) {
                const long long t = (long long)e * n3 + l;
                const double zt = (dmask >> Mesh::bit(l)) & 1u ? 0.0 : stage[l];
                const double dn = kind == 2 ? c1 * d0[t] + c2 * zt : c2 * zt;
                dref[t] = dn;
                xref[t] = kind == 0 ? dn : x0[t] + dn;
            }
        }
        DR dz = desc(m.values), dg = desc(m.ghost), dd = desc(d), dx = desc(x);
        D32 dl = desc(list);
        Tables t(m);
        emulate(TH(lines_update_first0)(), (long long)list.size(), [&] {
            kernels[kind](&dz, &dg, &dl, &t.fn, &t.fc, &t.eo, &t.eoff, &t.eent, &t.vo, &t.voff, &t.vent, &t.cnt, &t.edir,
                          &t.vdir, &dd, &dx, c1, c2);
        });
        report(names[kind], differ(d, dref) + differ(x, xref), d.size() + x.size());
    }
}

// solve_kernel's arithmetic on the host, block by block, each multiply-add one fma as nvcc
// contracts it (this file is built with -ffp-contract=off).
void solve_reference(const std::vector<double>& L, const std::vector<int64_t>& copy, const std::vector<int32_t>& len,
                     const std::vector<int64_t>& warp_off, const std::vector<double>& s, std::vector<double>& z)
{
    for (size_t q = 0; q < len.size(); ++q) {
        const int lane = (int)(q & 31), m = len[q];
        const int64_t base = warp_off[q >> 5];
        double ring[P];
        for (double& r : ring) r = 0.0;
        for (int r = 0; r < m; ++r) {
            const double* Lr = L.data() + (base + r) * (P + 1) * 32 + lane;
            const int64_t t = copy[(base + r) * 32 + lane];
            double acc = s[t];
            for (int dd = 1; dd <= P; ++dd) acc = std::fma(-Lr[dd * 32], ring[dd - 1], acc);
            const double v = acc / Lr[0];
            z[t] = v;
            for (int dd = P - 1; dd > 0; --dd) ring[dd] = ring[dd - 1];
            ring[0] = v;
        }
        for (double& r : ring) r = 0.0;
        for (int r = m - 1; r >= 0; --r) {
            const double* Lr = L.data() + (base + r) * (P + 1) * 32 + lane;
            const int64_t t = copy[(base + r) * 32 + lane];
            const double v = (z[t] - ring[0]) / Lr[0];
            z[t] = v;
            for (int dd = 1; dd <= P; ++dd) ring[dd - 1] = std::fma(Lr[dd * 32], v, dd < P ? ring[dd] : 0.0);
        }
    }
}

void solve()
{
    std::mt19937 rng(97 + P);
    const int blocks = 3 * 32 + 11;   // three full warps and a partial one
    std::vector<int32_t> len(blocks);
    for (auto& m : len) m = std::uniform_int_distribution<int>(1, 4 * P + 9)(rng);
    std::sort(len.begin(), len.end(), std::greater<int>());
    // The kernel's contract: a multiple of its elements per block (128), padded with
    // zero-length blocks, as NodeLines lays them out.
    len.resize((blocks + 127) / 128 * 128, 0);
    const int warps = (int)len.size() / 32;
    std::vector<int64_t> warp_off(warps + 1, 0);
    for (int w = 0; w < warps; ++w) warp_off[w + 1] = warp_off[w] + len[w * 32];
    const int64_t slots = warp_off[warps];
    long long rows = 0;
    for (int m : len) rows += m;
    // Every row its own copy, scattered over a larger vector (the counting copies).
    std::vector<int64_t> pool(3 * rows);
    for (size_t i = 0; i < pool.size(); ++i) pool[i] = (int64_t)i;
    std::shuffle(pool.begin(), pool.end(), rng);
    std::vector<int64_t> copy(slots * 32, -1);
    std::vector<double> L = random_vector(slots * 32 * (P + 1), 5, -0.4, 0.4);
    long long next = 0;
    for (int q = 0; q < blocks; ++q)
        for (int r = 0; r < len[q]; ++r) {
            const int64_t slot = warp_off[q >> 5] + r;
            copy[slot * 32 + (q & 31)] = pool[next++];
            L[slot * (P + 1) * 32 + (q & 31)] = 1.0 + std::uniform_real_distribution<double>(0.0, 1.0)(rng);
        }
    std::vector<double> s = random_vector(3 * rows, 7), z = random_vector(3 * rows, 8), ref = z;
    solve_reference(L, copy, len, warp_off, s, ref);
    DR dL = desc(L), ds = desc(s), dz = desc(z);
    D64 dc = desc(copy), dw = desc(warp_off);
    D32 dlen = desc(len);
    const int64_t threads = TH(lines_solve)();
    emulate(threads, (long long)len.size() / threads, [&] { CI(lines_solve)(&dL, &dc, &dlen, &dw, &ds, &dz); });
    report("lines_solve", differ(z, ref), z.size());
}

}  // namespace

int main()
{
    residual();
    update();
    solve();
    printf("  block threads: residual %lld, update %lld, solve %lld\n", (long long)TH(lines_residual0)(),
           (long long)TH(lines_update_first0)(), (long long)TH(lines_solve)());
    printf("NODE-LINE KERNELS WARP EMULATION (bit-exact vs host formulas, p = %d): %s\n", P, failed ? "FAIL" : "PASS");
    return failed;
}
