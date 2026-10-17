// CPU gate for the pass-generated warp kernel: runs the SAME kernel body the GPU
// gets (hl_gpu_pipeline.py --emulate, which applies --mir-emulate-warp instead of
// --mir-gpu-wrap) with one thread per lane, 32 per element, and checks EVERY
// element against the Knaus oracle.
//
// The warp semantics come from emu_warp_runtime.h.
//
// -DMIR_EMU_REAL=float (hl_gpu_pipeline.py --f32): the kernel in single precision; the
// inputs are rounded to float first and the oracle runs in double on them, so the
// difference is the kernel's own rounding, checked relative to the result's size.
//
// Build + run: python3 test/hl_gpu_pipeline.py [p] --emulate [--emu-e=E]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

#include "emu_warp_runtime.h"
#include "knaus_oracle.h"

#ifndef MIR_EMU_REAL
#define MIR_EMU_REAL double
#endif
using Real = MIR_EMU_REAL;

// Argument order is the operator's: (u, Btil, Dtil, W, D, G) then the out-param.
extern "C" void _mlir_ciface_laplacian_apply(Desc<Real, 4>*, Desc<Real, 2>*, Desc<Real, 2>*,
                                             Desc<Real, 2>*, Desc<Real, 2>*, Desc<Real, 6>*,
                                             Desc<Real, 4>*);
extern "C" int64_t laplacian_apply__threads();   // the block size the kernel was built for

int main(int argc, char** argv)
{
    const long long E = argc > 1 ? atoll(argv[1]) : 3;
    const int p = argc > 2 ? atoi(argv[2]) : 7;
    const int n = p + 1, nn = n * n, n3 = nn * n;
    const long long gElem = 3LL * p * nn * 3;

    // Btil/Dtil n x n with a zero last row, as the GPU harness passes them.
    std::vector<double> hBt((size_t)nn, 0.0), hDt((size_t)nn, 0.0), hDm(nn), hW(nn);
    srand(42);
    auto rnd = [] { return 2.0 * rand() / RAND_MAX - 1.0; };
    for (int i = 0; i < p * n; ++i) hBt[i] = rnd();
    for (int i = 0; i < p * n; ++i) hDt[i] = rnd();
    for (auto& x : hDm) x = rnd();
    for (auto& x : hW) x = rnd();

    std::vector<double> U((size_t)E * n3), G((size_t)E * gElem);
    for (long long e = 0; e < E; ++e) {
        for (int i = 0; i < n3; ++i) U[(size_t)e * n3 + i] = hashVal(e, i, 0u);
        for (long long i = 0; i < gElem; ++i)
            G[(size_t)(e * gElem + i)] = hashVal(e, (size_t)i, 0x9e37u);
    }
    // The kernel's inputs in its precision; the oracle sees the same (rounded) values.
    auto narrow = [](std::vector<double>& v) {
        std::vector<Real> r(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            r[i] = (Real)v[i];
            v[i] = (double)r[i];
        }
        return r;
    };
    std::vector<Real> rU = narrow(U), rG = narrow(G), rBt = narrow(hBt), rDt = narrow(hDt), rDm = narrow(hDm),
                      rW = narrow(hW);
    // NaN, not zero: an entry of Y the kernel never writes must not pass.
    std::vector<Real> Y((size_t)E * n3, (Real)std::nan(""));

    Desc<Real, 4> dU{rU.data(), rU.data(), 0, {E, n, n, n}, {n3, nn, n, 1}};
    Desc<Real, 4> dY{Y.data(), Y.data(), 0, {E, n, n, n}, {n3, nn, n, 1}};
    Desc<Real, 2> dBt{rBt.data(), rBt.data(), 0, {n, n}, {n, 1}};
    Desc<Real, 2> dDt{rDt.data(), rDt.data(), 0, {n, n}, {n, 1}};
    Desc<Real, 2> dW{rW.data(), rW.data(), 0, {n, n}, {n, 1}};
    Desc<Real, 2> dDm{rDm.data(), rDm.data(), 0, {n, n}, {n, 1}};
    Desc<Real, 6> dG{rG.data(), rG.data(), 0, {E, 3, p, 3, n, n},   // [dir][face][component][row][col]
               {gElem, (int64_t)p * 3 * nn, 3LL * nn, (int64_t)nn, n, 1}};

    // One block at a time: the kernel's workgroup buffers are process-wide
    // globals in the emulated build, so two blocks in flight would share them.
    Warp w;
    w.bar.threads = (int)laplacian_apply__threads();
    warp = &w;
    for (long long e = 0; e < E; ++e) {   // one block per element is enough for any shape
        std::vector<std::thread> lanes;
        for (int L = 0; L < w.bar.threads; ++L)
            lanes.emplace_back([&, L, e] {
                tlLane = L;
                tlBlock = e;
                tlGrid = E;
                _mlir_ciface_laplacian_apply(&dU, &dBt, &dDt, &dW, &dDm, &dG, &dY);
            });
        for (auto& t : lanes) t.join();
#ifdef MIR_EMU_TRACE
        if (e == 0)
            trace_block(w.bar.threads);
        else
            for (int L = 0; L < w.bar.threads; ++L) g_trace[L].clear();
#endif
    }
#ifdef MIR_EMU_TRACE
    trace_report();
#endif

    double err = 0, scale = 0;
    std::vector<double> ref(n3);
    for (long long e = 0; e < E; ++e) {
        oracle(p, &U[(size_t)e * n3], hBt.data(), hDt.data(), hDm.data(), hW.data(),
               &G[(size_t)(e * gElem)], ref.data());
        double eo = 0, mo = 0, mr = 0;
        for (int i = 0; i < n3; ++i) {
            const double o = Y[(size_t)e * n3 + i];
            eo = std::isnan(o) ? INFINITY : fmax(eo, fabs(o - ref[i]));
            if (std::isinf(eo)) break;
            mo = fmax(mo, fabs(o));
            mr = fmax(mr, fabs(ref[i]));
        }
        printf("  e=%-3lld max|out|=%.3e  max|ref|=%.3e  max|out-ref|=%.3e\n",
               e, mo, mr, eo);
        err = fmax(err, eo);
        scale = fmax(scale, mr);
    }
    // double: the absolute bound of the FP64 gates; float: relative to the result.
    const bool single = sizeof(Real) == sizeof(float);
    const bool ok = single ? err < 1e-5 * scale : err < 1e-11;
    printf("WARP EMULATION (p=%d, E=%lld, %d threads per block, %s): max|err| = %.3e%s  %s\n", p, E, w.bar.threads,
           single ? "f32" : "f64", err, single ? " (relative bound 1e-5)" : "", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
