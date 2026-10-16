// CPU gate for the pass-generated p-multigrid kernels (dss_gpu_pipeline.py --pmg
// --emulate): the same kernel bodies the GPU gets, one thread per lane, compared BIT
// FOR BIT with host formulas:
//   cheb_*_p{7,3,1}  unstructured_chebyshev_kernel's formula (DSS of b - y, Jacobi,
//                    d and x updates), on random tables with ghosts and Dirichlet
//                    entities (emu_cellwise_mesh.h);
//   masked_copy_p*   r = b, zero on Dirichlet copies;
//   interp_*         M (x) M (x) M per element, in the generated kernels' sum-factorized
//                    order (the hand transfer kernels sum in another order);
//   pmg_dots_*       per-element partials in the generated reduction order.
// Built with -ffp-contract=off: the generated code does not fuse multiply-adds.
//
// Build + run: python3 test/dss_gpu_pipeline.py 7 --pmg --emulate [--emu-tsan]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <thread>
#include <vector>

#include "emu_cellwise_mesh.h"
#include "emu_warp_runtime.h"

using D64 = Desc<double>;
using D32 = Desc<int32_t>;
using D8 = Desc<int8_t>;
using Cheb = void(D64*, D64*, D64*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, D32*, D8*, D8*, D64*, D64*,
                  D64*, double, double);
using Masked = void(D32*, D8*, D32*, D32*, D8*, D8*, D64*, D64*);
using Interp = void(D32*, D64*, D64*, D64*);
extern "C" {
Cheb _mlir_ciface_cheb_first0_p7, _mlir_ciface_cheb_first_p7, _mlir_ciface_cheb_next_p7;
Cheb _mlir_ciface_cheb_first0_p3, _mlir_ciface_cheb_first_p3, _mlir_ciface_cheb_next_p3;
Cheb _mlir_ciface_cheb_first0_p1, _mlir_ciface_cheb_first_p1, _mlir_ciface_cheb_next_p1;
Masked _mlir_ciface_masked_copy_p7, _mlir_ciface_masked_copy_p3;
Interp _mlir_ciface_interp_4x8, _mlir_ciface_interp_2x4, _mlir_ciface_interp_8x4, _mlir_ciface_interp_8x4_add,
    _mlir_ciface_interp_4x2, _mlir_ciface_interp_4x2_add;
void _mlir_ciface_pmg_dots_p7_az_zz(D32*, D32*, D64*, D64*, D64*, D64*);
}

namespace {

template <typename T>
Desc<T> desc(std::vector<T>& v)
{
    return {v.data(), v.data(), 0, {(int64_t)v.size()}, {1}};
}

template <typename F>
void emulate(long long grid, F&& kernel)
{
    Warp w;
    warp = &w;
    for (long long b = 0; b < grid; ++b) {
        std::vector<std::thread> lanes;
        for (int L = 0; L < kLanes; ++L)
            lanes.emplace_back([&, L, b] {
                tlLane = L;
                tlBlock = b;
                tlGrid = grid;
                kernel();
            });
        for (auto& t : lanes) t.join();
    }
}

long long differ(const std::vector<double>& a, const std::vector<double>& b)
{
    long long bad = 0;
    for (size_t i = 0; i < a.size(); ++i)
        if (memcmp(&a[i], &b[i], sizeof(double)) != 0) ++bad;
    return bad;
}

int failed = 0;
void report(const char* what, int p, long long bad, size_t total)
{
    printf("  %-14s p=%d  %lld of %zu values differ\n", what, p, bad, total);
    failed |= bad != 0;
}

std::vector<double> random_vector(long long size, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    std::vector<double> v(size);
    for (auto& x : v) x = val(rng);
    return v;
}

template <int P>
void chebyshev(Cheb* first0, Cheb* first, Cheb* next)
{
    using Mesh = CellwiseMesh<P>;
    constexpr int n3 = Mesh::n3;
    Mesh m = Mesh::random(10, 4, 2026 + P);
    std::vector<int32_t> list = {3, 0, 9, 5, 7, 1};
    std::vector<double> y = random_vector(m.E * n3, 11), d0 = random_vector(m.E * n3, 12),
                        x0 = random_vector(m.E * n3, 13);
    const double c1 = 0.37, c2 = 1.9;
    const char* names[3] = {"cheb_first0", "cheb_first", "cheb_next"};
    Cheb* kernels[3] = {first0, first, next};
    for (int kind = 0; kind < 3; ++kind) {
        std::vector<double> d = d0, x = x0, dref = d0, xref = x0, stage;
        m.minus = kind == 0 ? std::vector<double>() : y;
        for (int32_t e : list) {
            m.stage(e, stage);
            const unsigned dmask = m.dirichlet_mask(e);
            for (int l = 0; l < n3; ++l) {
                const long long t = (long long)e * n3 + l;
                const double z = (dmask >> Mesh::bit(l)) & 1u ? 0.0 : stage[l] / m.diag[t];
                const double dn = kind == 2 ? c1 * d0[t] + c2 * z : c2 * z;
                dref[t] = dn;
                xref[t] = kind == 0 ? dn : x0[t] + dn;
            }
        }
        D64 db = desc(m.values), dg = desc(m.ghost), dy = desc(y), ddiag = desc(m.diag), dd = desc(d), dx = desc(x);
        D32 dl = desc(list), fn = desc(m.face_nbr), eo = desc(m.edge_of), eoff = desc(m.edge_off),
            eent = desc(m.edge_ent), vo = desc(m.vert_of), voff = desc(m.vert_off), vent = desc(m.vert_ent),
            cnt = desc(m.counting);
        D8 fc = desc(m.face_code), edir = desc(m.edge_dir), vdir = desc(m.vert_dir);
        emulate((long long)list.size(), [&] {
            kernels[kind](&db, &dg, &dy, &dl, &fn, &fc, &eo, &eoff, &eent, &vo, &voff, &vent, &cnt, &edir, &vdir, &ddiag,
                          &dd, &dx, c1, c2);
        });
        report(names[kind], P, differ(d, dref) + differ(x, xref), d.size() + x.size());
    }
}

template <int P>
void masked_copy(Masked* kernel)
{
    using Mesh = CellwiseMesh<P>;
    constexpr int n3 = Mesh::n3;
    Mesh m = Mesh::random(10, 0, 77 + P);
    std::vector<int32_t> list = {8, 2, 4};
    std::vector<double> r(m.E * n3, std::nan("")), ref = r;
    for (int32_t e : list) {
        const unsigned dmask = m.dirichlet_mask(e);
        for (int l = 0; l < n3; ++l) {
            const long long t = (long long)e * n3 + l;
            ref[t] = (dmask >> Mesh::bit(l)) & 1u ? 0.0 : m.values[t];
        }
    }
    D32 dl = desc(list), eo = desc(m.edge_of), vo = desc(m.vert_of);
    D8 fc = desc(m.face_code), edir = desc(m.edge_dir), vdir = desc(m.vert_dir);
    D64 db = desc(m.values), dr = desc(r);
    emulate((long long)list.size(), [&] { kernel(&dl, &fc, &eo, &vo, &edir, &vdir, &db, &dr); });
    report("masked_copy", P, differ(r, ref), r.size());
}

void interpolate(const char* what, Interp* kernel, int R, int C, bool add)
{
    const long long E = 5;
    std::vector<int32_t> list = {4, 1, 3};
    std::vector<double> M = random_vector(R * C, 21), in = random_vector(E * C * C * C, 22),
                        out = random_vector(E * R * R * R, 23), ref = out;
    for (int32_t e : list) {
        const double* u = &in[(size_t)e * C * C * C];
        std::vector<double> s1(C * C * R), s2(C * R * R);
        for (int q = 0; q < C * C * R; ++q) {
            const int jk = q / R, c = q % R;
            double s = 0.0;
            for (int mm = 0; mm < C; ++mm) s = s + M[c * C + mm] * u[jk * C + mm];
            s1[q] = s;
        }
        for (int q = 0; q < C * R * R; ++q) {
            const int j = q / (R * R), b = (q / R) % R, c = q % R;
            double s = 0.0;
            for (int k = 0; k < C; ++k) s = s + M[b * C + k] * s1[(j * C + k) * R + c];
            s2[q] = s;
        }
        for (int q = 0; q < R * R * R; ++q) {
            const int a = q / (R * R), bc = q % (R * R);
            double s = 0.0;
            for (int j = 0; j < C; ++j) s = s + M[a * C + j] * s2[j * R * R + bc];
            const size_t t = (size_t)e * R * R * R + q;
            ref[t] = add ? out[t] + s : s;
        }
    }
    D32 dl = desc(list);
    D64 dm = desc(M), din = desc(in), dout = desc(out);
    emulate((long long)list.size(), [&] { kernel(&dl, &dm, &din, &dout); });
    printf("  %-14s %dx%d  %lld of %zu values differ\n", what, R, C, differ(out, ref), out.size());
    failed |= differ(out, ref) != 0;
}

void dots()
{
    using Mesh = CellwiseMesh<7>;
    constexpr int n3 = Mesh::n3;
    Mesh m = Mesh::random(6, 0, 99);
    std::vector<int32_t> list = {5, 0, 2};
    std::vector<double> z = random_vector(m.E * n3, 31), sc(7, 0.0), part(list.size() * 2, std::nan("")), ref(part.size());
    for (size_t w = 0; w < list.size(); ++w) {
        const long long e = list[w], base = e * n3;
        auto weight = [&](long long l) { return ((uint32_t)m.counting[e] >> Mesh::bit(l)) & 1u ? 1.0 : 0.0; };
        ref[2 * w] = lane_sum(0, n3, [&](long long l) { return m.a[base + l] * z[base + l] * weight(l); });
        ref[2 * w + 1] = lane_sum(0, n3, [&](long long l) { return z[base + l] * z[base + l] * weight(l); });
    }
    D32 dl = desc(list), cnt = desc(m.counting);
    D64 dz = desc(z), da = desc(m.a), dsc = desc(sc), dp = desc(part);
    emulate((long long)list.size(), [&] { _mlir_ciface_pmg_dots_p7_az_zz(&dl, &cnt, &dz, &da, &dsc, &dp); });
    report("pmg_dots", 7, differ(part, ref), part.size());
}

}  // namespace

int main()
{
    chebyshev<7>(_mlir_ciface_cheb_first0_p7, _mlir_ciface_cheb_first_p7, _mlir_ciface_cheb_next_p7);
    chebyshev<3>(_mlir_ciface_cheb_first0_p3, _mlir_ciface_cheb_first_p3, _mlir_ciface_cheb_next_p3);
    chebyshev<1>(_mlir_ciface_cheb_first0_p1, _mlir_ciface_cheb_first_p1, _mlir_ciface_cheb_next_p1);
    masked_copy<7>(_mlir_ciface_masked_copy_p7);
    masked_copy<3>(_mlir_ciface_masked_copy_p3);
    interpolate("restrict 7->3", _mlir_ciface_interp_4x8, 4, 8, false);
    interpolate("restrict 3->1", _mlir_ciface_interp_2x4, 2, 4, false);
    interpolate("prolong 3->7", _mlir_ciface_interp_8x4, 8, 4, false);
    interpolate("prolong 3->7 +", _mlir_ciface_interp_8x4_add, 8, 4, true);
    interpolate("prolong 1->3", _mlir_ciface_interp_4x2, 4, 2, false);
    interpolate("prolong 1->3 +", _mlir_ciface_interp_4x2_add, 4, 2, true);
    dots();
    printf("P-MULTIGRID KERNELS WARP EMULATION (bit-exact vs host formulas): %s\n", failed ? "FAIL" : "PASS");
    return failed;
}
