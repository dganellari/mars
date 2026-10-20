// CPU gate for the pass-generated p-multigrid kernels (dss_gpu_pipeline.py --pmg
// --emulate): the same kernel bodies the GPU gets, one thread per lane, compared BIT
// FOR BIT with host formulas:
//   cheb_*_p{7,3,1}  unstructured_chebyshev_kernel's formula (DSS of b - y, Jacobi,
//                    d and x updates), on random tables with ghosts and Dirichlet
//                    entities (emu_cellwise_mesh.h);
//   masked_copy_p*   r = b, zero on Dirichlet copies;
//   interp_*         M (x) M (x) M per element, in the generated kernels' sum-factorized
//                    order (the hand transfer kernels sum in another order);
//   pmg_dots_*       per-element partials in the generated reduction order;
//   -DMIR_EMU_F32 only: to_f32/to_f64_p* (the conversions at p = 7 and 1),
//                    interp_8x4_counted (zero off the counting copies) and dss_p7.
// Built with -ffp-contract=off: the generated code does not fuse multiply-adds.
// -DMIR_EMU_F32: the single-precision kernels (names ending in _f32), against the same
// formulas in float.
//
// Build + run: python3 test/dss_gpu_pipeline.py 7 --pmg [--f32] --emulate [--emu-tsan]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <thread>
#include <vector>

#include "emu_cellwise_mesh.h"
#include "emu_warp_runtime.h"

#ifdef MIR_EMU_F32
using Real = float;
#define K(name) name##_f32
#else
using Real = double;
#define K(name) name
#endif
#define CI(name) K(_mlir_ciface_##name)
#define CAT_(a, b) a##b
#define CAT(a, b) CAT_(a, b)
#define TH(name) CAT(K(name), __threads)

using DR = Desc<Real>;
using D64 = Desc<double>;
using D32 = Desc<int32_t>;
using D8 = Desc<int8_t>;
using Cheb = void(DR*, DR*, DR*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, D32*, D8*, D8*, DR*, DR*,
                  DR*, Real, Real);
using Masked = void(D32*, D8*, D32*, D32*, D8*, D8*, DR*, DR*);
using Interp = void(D32*, DR*, DR*, DR*);
extern "C" {
Cheb CI(cheb_first0_p7), CI(cheb_first_p7), CI(cheb_next_p7);
Cheb CI(cheb_first0_p3), CI(cheb_first_p3), CI(cheb_next_p3);
Cheb CI(cheb_first0_p1), CI(cheb_first_p1), CI(cheb_next_p1);
Masked CI(masked_copy_p7), CI(masked_copy_p3);
Interp CI(interp_4x8), CI(interp_2x4), CI(interp_8x4), CI(interp_8x4_add), CI(interp_4x2), CI(interp_4x2_add);
int64_t TH(cheb_first0_p7)(), TH(cheb_first0_p3)(), TH(cheb_first0_p1)();
int64_t TH(masked_copy_p7)(), TH(masked_copy_p3)(), TH(interp_4x8)();
#ifdef MIR_EMU_F32
using ToF32 = void(D32*, D64*, DR*);
using ToF64 = void(D32*, DR*, D64*);
ToF32 _mlir_ciface_to_f32_p7, _mlir_ciface_to_f32_p1;
ToF64 _mlir_ciface_to_f64_p7, _mlir_ciface_to_f64_p1;
int64_t to_f32_p7__threads(), to_f32_p1__threads();
void _mlir_ciface_interp_8x4_counted_f32(D32*, D32*, DR*, DR*, DR*);
void _mlir_ciface_dss_p7_f32(DR*, DR*, D32*, D32*, D8*, D32*, D32*, D32*, D32*, D32*, D32*, DR*);
int64_t interp_8x4_counted_f32__threads(), dss_p7_f32__threads();
#else
void _mlir_ciface_pmg_dots_p7_az_zz(D32*, D32*, D64*, D64*, D64*, D64*);
int64_t pmg_dots_p7_az_zz__threads();
#endif
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
void report(const char* what, int p, long long bad, size_t total)
{
    printf("  %-14s p=%d  %lld of %zu values differ\n", what, p, bad, total);
    failed |= bad != 0;
}

template <typename T = Real>
std::vector<T> random_vector(long long size, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    std::vector<T> v(size);
    for (auto& x : v) x = (T)val(rng);
    return v;
}

template <int P>
void chebyshev(int64_t threads, Cheb* first0, Cheb* first, Cheb* next)
{
    using Mesh = CellwiseMesh<P, Real>;
    constexpr int n3 = Mesh::n3;
    Mesh m = Mesh::random(10, 4, 2026 + P);
    std::vector<int32_t> list = {3, 0, 9, 5, 7, 1};
    std::vector<Real> y = random_vector(m.E * n3, 11), d0 = random_vector(m.E * n3, 12),
                      x0 = random_vector(m.E * n3, 13);
    const Real c1 = (Real)0.37, c2 = (Real)1.9;
    const char* names[3] = {"cheb_first0", "cheb_first", "cheb_next"};
    Cheb* kernels[3] = {first0, first, next};
    for (int kind = 0; kind < 3; ++kind) {
        std::vector<Real> d = d0, x = x0, dref = d0, xref = x0, stage;
        m.minus = kind == 0 ? std::vector<Real>() : y;
        for (int32_t e : list) {
            m.stage(e, stage);
            const unsigned dmask = m.dirichlet_mask(e);
            for (int l = 0; l < n3; ++l) {
                const long long t = (long long)e * n3 + l;
                const Real z = (dmask >> Mesh::bit(l)) & 1u ? (Real)0 : stage[l] / m.diag[t];
                const Real dn = kind == 2 ? c1 * d0[t] + c2 * z : c2 * z;
                dref[t] = dn;
                xref[t] = kind == 0 ? dn : x0[t] + dn;
            }
        }
        DR db = desc(m.values), dg = desc(m.ghost), dy = desc(y), ddiag = desc(m.diag), dd = desc(d), dx = desc(x);
        D32 dl = desc(list), fn = desc(m.face_nbr), eo = desc(m.edge_of), eoff = desc(m.edge_off),
            eent = desc(m.edge_ent), vo = desc(m.vert_of), voff = desc(m.vert_off), vent = desc(m.vert_ent),
            cnt = desc(m.counting);
        D8 fc = desc(m.face_code), edir = desc(m.edge_dir), vdir = desc(m.vert_dir);
        emulate(threads, (long long)list.size(), [&] {
            kernels[kind](&db, &dg, &dy, &dl, &fn, &fc, &eo, &eoff, &eent, &vo, &voff, &vent, &cnt, &edir, &vdir, &ddiag,
                          &dd, &dx, c1, c2);
        });
        report(names[kind], P, differ(d, dref) + differ(x, xref), d.size() + x.size());
    }
}

template <int P>
void masked_copy(int64_t threads, Masked* kernel)
{
    using Mesh = CellwiseMesh<P, Real>;
    constexpr int n3 = Mesh::n3;
    Mesh m = Mesh::random(10, 0, 77 + P);
    std::vector<int32_t> list = {8, 2, 4};
    std::vector<Real> r(m.E * n3, (Real)std::nan("")), ref = r;
    for (int32_t e : list) {
        const unsigned dmask = m.dirichlet_mask(e);
        for (int l = 0; l < n3; ++l) {
            const long long t = (long long)e * n3 + l;
            ref[t] = (dmask >> Mesh::bit(l)) & 1u ? (Real)0 : m.values[t];
        }
    }
    D32 dl = desc(list), eo = desc(m.edge_of), vo = desc(m.vert_of);
    D8 fc = desc(m.face_code), edir = desc(m.edge_dir), vdir = desc(m.vert_dir);
    DR db = desc(m.values), dr = desc(r);
    emulate(threads, (long long)list.size(), [&] { kernel(&dl, &fc, &eo, &vo, &edir, &vdir, &db, &dr); });
    report("masked_copy", P, differ(r, ref), r.size());
}

// counting: as mir.interpolate's counting operand (R = 8 only; null: plain).
void interpolate(const char* what, Interp* kernel, int R, int C, bool add,
                 void (*counted)(D32*, D32*, DR*, DR*, DR*) = nullptr, int64_t counted_threads = 0)
{
    const long long E = 5;
    std::vector<int32_t> list = {4, 1, 3};
    std::vector<int32_t> counting = {0x0f0f0f0f, 0x7ffffff1, 0x12345678, (int32_t)0x80000000u, 0x7fffffff};
    std::vector<Real> M = random_vector(R * C, 21), in = random_vector(E * C * C * C, 22),
                      out = random_vector(E * R * R * R, 23), ref = out;
    for (int32_t e : list) {
        const Real* u = &in[(size_t)e * C * C * C];
        std::vector<Real> s1(C * C * R), s2(C * R * R);
        for (int q = 0; q < C * C * R; ++q) {
            const int jk = q / R, c = q % R;
            Real s = 0;
            for (int mm = 0; mm < C; ++mm) s = s + M[c * C + mm] * u[jk * C + mm];
            s1[q] = s;
        }
        for (int q = 0; q < C * R * R; ++q) {
            const int j = q / (R * R), b = (q / R) % R, c = q % R;
            Real s = 0;
            for (int k = 0; k < C; ++k) s = s + M[b * C + k] * s1[(j * C + k) * R + c];
            s2[q] = s;
        }
        for (int q = 0; q < R * R * R; ++q) {
            const int a = q / (R * R), bc = q % (R * R);
            Real s = 0;
            for (int j = 0; j < C; ++j) s = s + M[a * C + j] * s2[j * R * R + bc];
            const size_t t = (size_t)e * R * R * R + q;
            ref[t] = add ? out[t] + s : s;
            if (counted && !(((uint32_t)counting[e] >> CellwiseMesh<7>::bit(q)) & 1u)) ref[t] = 0;
        }
    }
    D32 dl = desc(list), dc = desc(counting);
    DR dm = desc(M), din = desc(in), dout = desc(out);
    if (counted)
        emulate(counted_threads, (long long)list.size(), [&] { counted(&dl, &dc, &dm, &din, &dout); });
    else
        emulate(TH(interp_4x8)(), (long long)list.size(), [&] { kernel(&dl, &dm, &din, &dout); });
    printf("  %-14s %dx%d  %lld of %zu values differ\n", what, R, C, differ(out, ref), out.size());
    failed |= differ(out, ref) != 0;
}

#ifdef MIR_EMU_F32
void dss()
{
    using Mesh = CellwiseMesh<7, float>;
    constexpr int n3 = Mesh::n3;
    Mesh m = Mesh::random(10, 4, 2027);
    std::vector<int32_t> list = {3, 0, 9, 5, 7, 1};
    std::vector<float> out(m.E * n3, std::nanf("")), ref = out, stage;
    for (int32_t e : list) {
        m.stage(e, stage);
        for (int l = 0; l < n3; ++l) ref[(size_t)e * n3 + l] = stage[l];
    }
    DR dv = desc(m.values), dg = desc(m.ghost), dout = desc(out);
    D32 dl = desc(list), fn = desc(m.face_nbr), eo = desc(m.edge_of), eoff = desc(m.edge_off),
        eent = desc(m.edge_ent), vo = desc(m.vert_of), voff = desc(m.vert_off), vent = desc(m.vert_ent);
    D8 fc = desc(m.face_code);
    emulate(dss_p7_f32__threads(), (long long)list.size(),
            [&] { _mlir_ciface_dss_p7_f32(&dv, &dg, &dl, &fn, &fc, &eo, &eoff, &eent, &vo, &voff, &vent, &dout); });
    report("dss", 7, differ(out, ref), out.size());
}

template <int P>
void convert(int64_t threads, ToF32* to32, ToF64* to64)
{
    const long long E = 5, n3 = (P + 1) * (P + 1) * (P + 1);
    std::vector<int32_t> list = {4, 0, 2};
    std::vector<double> src = random_vector<double>(E * n3, 41), back = random_vector<double>(E * n3, 42),
                        back_ref = back;
    std::vector<float> low = random_vector<float>(E * n3, 43), low_ref = low;
    for (int32_t e : list)
        for (long long l = 0; l < n3; ++l) {
            low_ref[e * n3 + l] = (float)src[e * n3 + l];
            back_ref[e * n3 + l] = (double)low_ref[e * n3 + l];
        }
    D32 dl = desc(list);
    D64 ds = desc(src), db = desc(back);
    DR dlow = desc(low);
    emulate(threads, (long long)list.size(), [&] { to32(&dl, &ds, &dlow); });
    emulate(threads, (long long)list.size(), [&] { to64(&dl, &dlow, &db); });
    report("to_f32", P, differ(low, low_ref), low.size());
    report("to_f64", P, differ(back, back_ref), back.size());
}
#else
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
    emulate(pmg_dots_p7_az_zz__threads(), (long long)list.size(), [&] { _mlir_ciface_pmg_dots_p7_az_zz(&dl, &cnt, &dz, &da, &dsc, &dp); });
    report("pmg_dots", 7, differ(part, ref), part.size());
}
#endif

}  // namespace

int main()
{
    chebyshev<7>(TH(cheb_first0_p7)(), CI(cheb_first0_p7), CI(cheb_first_p7), CI(cheb_next_p7));
    chebyshev<3>(TH(cheb_first0_p3)(), CI(cheb_first0_p3), CI(cheb_first_p3), CI(cheb_next_p3));
    chebyshev<1>(TH(cheb_first0_p1)(), CI(cheb_first0_p1), CI(cheb_first_p1), CI(cheb_next_p1));
    masked_copy<7>(TH(masked_copy_p7)(), CI(masked_copy_p7));
    masked_copy<3>(TH(masked_copy_p3)(), CI(masked_copy_p3));
    interpolate("restrict 7->3", CI(interp_4x8), 4, 8, false);
    interpolate("restrict 3->1", CI(interp_2x4), 2, 4, false);
    interpolate("prolong 3->7", CI(interp_8x4), 8, 4, false);
    interpolate("prolong 3->7 +", CI(interp_8x4_add), 8, 4, true);
    interpolate("prolong 1->3", CI(interp_4x2), 4, 2, false);
    interpolate("prolong 1->3 +", CI(interp_4x2_add), 4, 2, true);
#ifdef MIR_EMU_F32
    interpolate("prolong counted", nullptr, 8, 4, false, _mlir_ciface_interp_8x4_counted_f32,
                interp_8x4_counted_f32__threads());
    dss();
    convert<7>(to_f32_p7__threads(), _mlir_ciface_to_f32_p7, _mlir_ciface_to_f64_p7);
    convert<1>(to_f32_p1__threads(), _mlir_ciface_to_f32_p1, _mlir_ciface_to_f64_p1);
#else
    dots();
#endif
    printf("  block threads: cheb p7 %lld, p3 %lld, p1 %lld\n", (long long)TH(cheb_first0_p7)(),
           (long long)TH(cheb_first0_p3)(), (long long)TH(cheb_first0_p1)());
    printf("P-MULTIGRID KERNELS WARP EMULATION (bit-exact vs host formulas, %s): %s\n",
           sizeof(Real) == 4 ? "f32" : "f64", failed ? "FAIL" : "PASS");
    return failed;
}
