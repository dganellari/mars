// CPU gate for the pass-generated BiCGStab kernels (dss_gpu_pipeline.py --bicgstab
// --emulate): runs the same kernel bodies the GPU gets with one thread per lane, 32
// per block, and compares BIT FOR BIT with host formulas in the same operation order:
//   vector updates   p = r + beta (p - omega v), s = r - alpha v, and x/r with the
//                    per-element partials of <r, r> and <rh, r> (lane L sums the
//                    copies L, L + 32, ...; then an xor-shuffle tree over the lanes);
//   fold_nv2         per-element partials -> block partials (contiguous ranges);
//   finish_*, totals_nv*, step_*  the Krylov scalar steps of mars_cellwise_krylov.hpp.
// Built with -ffp-contract=off: the generated code does not fuse multiply-adds.
//
// Build + run: python3 test/dss_gpu_pipeline.py [p] --bicgstab --emulate [--emu-tsan]
#include <cmath>
#include <cstdint>
#include <cstdio>
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

using D64 = Desc<double>;
using D32 = Desc<int32_t>;
extern "C" {
void MIR_CAT(_mlir_ciface_bicg_update_p_p, MIR_DSS_P)(D32*, D64*, D64*, D64*, D64*);
void MIR_CAT(_mlir_ciface_bicg_update_s_p, MIR_DSS_P)(D32*, D64*, D64*, D64*, D64*);
void MIR_CAT(_mlir_ciface_bicg_update_xr_p, MIR_DSS_P)(D32*, D32*, D64*, D64*, D64*, D64*, D64*, D64*, D64*, D64*);
void _mlir_ciface_fold_nv2(D64*, D64*, int64_t);
void _mlir_ciface_finish_start(D64*, D64*, int64_t);
void _mlir_ciface_finish_alpha(D64*, D64*, int64_t);
void _mlir_ciface_finish_omega(D64*, D64*, int64_t);
void _mlir_ciface_finish_rho(D64*, D64*, int64_t);
void _mlir_ciface_totals_nv2(D64*, D64*, int64_t);
void _mlir_ciface_step_rho(D64*);
}

namespace {

namespace cw = mars::cellwise;
constexpr int P = MIR_DSS_P, n = P + 1, n3 = n * n * n;
enum { kRho, kAlpha, kOmega, kBeta, kRR, kT0, kT1, kScalars };

template <typename T>
Desc<T> desc(std::vector<T>& v)
{
    return {v.data(), v.data(), 0, {(int64_t)v.size()}, {1}};
}

// `grid` blocks of 32 lanes, one block at a time.
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

// The fixed-order reduction of the generated kernels over values term(i), i in
// [lo, hi): lane L sums i = lo + L, lo + L + 32, ...; then the xor tree.
template <typename Term>
double lane_sum(long long lo, long long hi, Term term)
{
    double acc[32] = {};
    for (long long i = lo; i < hi; ++i) acc[(i - lo) % 32] += term(i);
    for (int off = 16; off >= 1; off /= 2) {
        double next[32];
        for (int L = 0; L < 32; ++L) next[L] = acc[L] + acc[L ^ off];
        memcpy(acc, next, sizeof(acc));
    }
    return acc[0];
}

long long differ(const std::vector<double>& a, const std::vector<double>& b)
{
    long long bad = 0;
    for (size_t i = 0; i < a.size(); ++i)
        if (memcmp(&a[i], &b[i], sizeof(double)) != 0) ++bad;
    return bad;
}

int failed = 0;
void report(const char* what, long long bad, size_t total)
{
    printf("  %-22s %lld of %zu values differ\n", what, bad, total);
    failed |= bad != 0;
}

}  // namespace

int main()
{
    const long long E = 9;
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    auto vec = [&](long long size) {
        std::vector<double> v(size);
        for (auto& x : v) x = val(rng);
        return v;
    };
    std::vector<int32_t> list = {4, 0, 7, 2, 8, 5};   // not every element, not in order
    std::vector<int32_t> counting(E);
    for (auto& c : counting) c = (int32_t)((uint32_t)std::uniform_int_distribution<long long>(0, 0xfffffff)(rng) |
                                           (1u << cw::kInteriorBit));
    std::vector<double> sc = vec(kScalars);
    D32 dl = desc(list), dcnt = desc(counting);
    D64 dsc = desc(sc);

    {   // p = r + beta (p - omega v), s = r - alpha v
        std::vector<double> r = vec(E * n3), v = vec(E * n3), p = vec(E * n3), s(E * n3, std::nan(""));
        std::vector<double> pref = p, sref = s;
        for (int32_t e : list)
            for (int l = 0; l < n3; ++l) {
                const long long t = (long long)e * n3 + l;
                pref[t] = r[t] + sc[kBeta] * (p[t] - sc[kOmega] * v[t]);
                sref[t] = r[t] - sc[kAlpha] * v[t];
            }
        D64 dr = desc(r), dv = desc(v), dp = desc(p), ds = desc(s);
        emulate(32, (long long)list.size(), [&] { MIR_CAT(_mlir_ciface_bicg_update_p_p, MIR_DSS_P)(&dl, &dr, &dv, &dp, &dsc); });
        emulate(32, (long long)list.size(), [&] { MIR_CAT(_mlir_ciface_bicg_update_s_p, MIR_DSS_P)(&dl, &dr, &dv, &ds, &dsc); });
        report("update_p", differ(p, pref), p.size());
        report("update_s", differ(s, sref), s.size());
    }
    {   // x += alpha p + omega s, r = s - omega t, partials of <r, r>, <rh, r>
        std::vector<double> x = vec(E * n3), p = vec(E * n3), s = vec(E * n3), t = vec(E * n3), rh = vec(E * n3);
        std::vector<double> r(E * n3, std::nan("")), part(list.size() * 2, std::nan(""));
        std::vector<double> xref = x, rref = r, pref(part.size());
        for (size_t w = 0; w < list.size(); ++w) {
            const long long e = list[w], base = e * n3;
            auto counted = [&](long long l) {
                const int bit = cw::entity_bit(cw::local_entity<P>(l / (n * n), (l / n) % n, l % n));
                return (((uint32_t)counting[e] >> bit) & 1u) ? 1.0 : 0.0;
            };
            for (int l = 0; l < n3; ++l) {
                xref[base + l] = x[base + l] + (sc[kAlpha] * p[base + l] + sc[kOmega] * s[base + l]);
                rref[base + l] = s[base + l] - sc[kOmega] * t[base + l];
            }
            pref[2 * w] = lane_sum(0, n3, [&](long long l) { return rref[base + l] * rref[base + l] * counted(l); });
            pref[2 * w + 1] = lane_sum(0, n3, [&](long long l) { return rh[base + l] * rref[base + l] * counted(l); });
        }
        D64 dx = desc(x), dp = desc(p), ds = desc(s), dt = desc(t), drh = desc(rh), dr = desc(r), dpart = desc(part);
        emulate(32, (long long)list.size(), [&] {
            MIR_CAT(_mlir_ciface_bicg_update_xr_p, MIR_DSS_P)(&dl, &dcnt, &dx, &dp, &ds, &dt, &drh, &dr, &dsc, &dpart);
        });
        report("update_xr x", differ(x, xref), x.size());
        report("update_xr r", differ(r, rref), r.size());
        report("update_xr partials", differ(part, pref), part.size());
    }
    {   // fold: 101 partial pairs into 4 block pairs
        const long long count = 101, blocks = 4, per = (count + blocks - 1) / blocks;
        std::vector<double> src = vec(count * 2), dst(blocks * 2, std::nan("")), ref(blocks * 2);
        for (long long b = 0; b < blocks; ++b)
            for (int j = 0; j < 2; ++j)
                ref[2 * b + j] = lane_sum(b * per, std::min(count, b * per + per),
                                          [&](long long i) { return src[2 * i + j]; });
        D64 dsrc = desc(src), ddst = desc(dst);
        emulate(32, blocks, [&] { _mlir_ciface_fold_nv2(&dsrc, &ddst, count); });
        report("fold_nv2", differ(dst, ref), dst.size());
    }
    {   // finish / totals / step on 37 block partials
        const long long count = 37;
        std::vector<double> src = vec(count * 2);
        D64 dsrc = desc(src);
        const double v0 = lane_sum(0, count, [&](long long i) { return src[2 * i]; });
        const double v1 = lane_sum(0, count, [&](long long i) { return src[2 * i + 1]; });
        std::vector<double> src1(count);
        for (long long i = 0; i < count; ++i) src1[i] = src[2 * i];
        D64 dsrc1 = desc(src1);
        const double u0 = lane_sum(0, count, [&](long long i) { return src1[i]; });
        auto run = [&](const char* what, auto kernel, auto expect) {
            std::vector<double> s0 = sc, ref = sc;
            D64 ds0 = desc(s0);
            expect(ref);
            emulate(32, 1, [&] { kernel(&ds0); });
            report(what, differ(s0, ref), s0.size());
        };
        run("finish_start", [&](D64* d) { _mlir_ciface_finish_start(&dsrc1, d, count); }, [&](std::vector<double>& s) {
            s[kRR] = u0; s[kRho] = u0; s[kAlpha] = 1.0; s[kOmega] = 1.0; s[kBeta] = 0.0;
        });
        run("finish_alpha", [&](D64* d) { _mlir_ciface_finish_alpha(&dsrc1, d, count); },
            [&](std::vector<double>& s) { s[kAlpha] = s[kRho] / u0; });
        run("finish_omega", [&](D64* d) { _mlir_ciface_finish_omega(&dsrc, d, count); },
            [&](std::vector<double>& s) { s[kOmega] = v0 / v1; });
        run("finish_rho", [&](D64* d) { _mlir_ciface_finish_rho(&dsrc, d, count); }, [&](std::vector<double>& s) {
            const double beta = (v1 / s[kRho]) * (s[kAlpha] / s[kOmega]);
            s[kRR] = v0; s[kBeta] = beta; s[kRho] = v1;
        });
        run("totals_nv2", [&](D64* d) { _mlir_ciface_totals_nv2(&dsrc, d, count); },
            [&](std::vector<double>& s) { s[kT0] = v0; s[kT1] = v1; });
        run("step_rho", [&](D64* d) { _mlir_ciface_step_rho(d); }, [&](std::vector<double>& s) {
            const double beta = (s[kT1] / s[kRho]) * (s[kAlpha] / s[kOmega]);
            s[kRR] = s[kT0]; s[kBeta] = beta; s[kRho] = s[kT1];
        });
    }
    printf("BICGSTAB KERNELS WARP EMULATION (p=%d, bit-exact vs host formulas): %s\n", P, failed ? "FAIL" : "PASS");
    return failed;
}
