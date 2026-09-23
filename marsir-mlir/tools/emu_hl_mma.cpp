// CPU gate for the pass-generated warp kernel: runs the SAME kernel body the GPU
// gets (hl_gpu_pipeline.py --emulate, which applies --mir-emulate-warp instead of
// --mir-gpu-wrap) with one thread per lane, 32 per element, and checks EVERY
// element against the Knaus oracle.
//
// The runtime below gives the kernel the warp semantics it relies on:
//   shuffle   every lane posts its value, all wait, each reads its source lane;
//   mma.sync  every lane posts its A/B/C fragment entries, all wait, each
//             computes its own two D entries from the assembled 8x4 @ 4x8 + 8x8;
//   barrier   all 32 lanes wait.
// Private memory is private (each thread's stack), shared memory is shared, and
// the lanes really run concurrently, so a missing barrier is a real data race.
//
// But every exchange here is a real rendezvous, which in C++ also ORDERS memory,
// so a race across an mma or a shuffle stays invisible. On the GPU neither
// shfl.sync nor mma.sync orders memory -- only a barrier does. The --emu-tsan
// build (-DMIR_EMU_TSAN, -fsanitize=thread) hides the exchanges from
// ThreadSanitizer, which then reports every pair of accesses from different
// lanes that no gpu.barrier orders: a missing barrier, under the GPU's rules.
//
// Build + run: python3 test/hl_gpu_pipeline.py [p] --emulate [--emu-e=E]
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <vector>

#include "knaus_oracle.h"

template <int R> struct Desc {
    double* allocated;
    double* aligned;
    int64_t offset;
    int64_t sizes[R];
    int64_t strides[R];
};

// Argument order is the operator's: (u, Btil, Dtil, W, D, G) then the out-param.
extern "C" void _mlir_ciface_laplacian_apply(Desc<4>*, Desc<2>*, Desc<2>*,
                                             Desc<2>*, Desc<2>*, Desc<6>*,
                                             Desc<4>*);

namespace {
constexpr int kLanes = 32;

// A mutex + condition variable, not std::barrier: libc++ keeps the barrier's
// atomics in the (uninstrumented) system dylib, so ThreadSanitizer would not see
// it order anything and would report every barrier-protected access as a race.
struct Barrier {
    std::mutex m;
    std::condition_variable cv;
    int count = 0;
    unsigned gen = 0;
    void arrive_and_wait()
    {
        std::unique_lock<std::mutex> lk(m);
        const unsigned g = gen;
        if (++count == kLanes) {
            count = 0;
            ++gen;
            cv.notify_all();
        } else {
            cv.wait(lk, [&] { return gen != g; });
        }
    }
};

struct Warp {
    Barrier bar;
    double xchg[kLanes];
    double A[8][4], B[4][8], C[8][8];
};
Warp* warp = nullptr;
thread_local int64_t tlLane = 0, tlBlock = 0;
thread_local double tlHi = 0;
}  // namespace

#ifdef MIR_EMU_TSAN
extern "C" {
void AnnotateIgnoreSyncBegin(const char*, int);
void AnnotateIgnoreSyncEnd(const char*, int);
void AnnotateIgnoreReadsBegin(const char*, int);
void AnnotateIgnoreReadsEnd(const char*, int);
void AnnotateIgnoreWritesBegin(const char*, int);
void AnnotateIgnoreWritesEnd(const char*, int);
}
// An exchange TSan must not see: neither its synchronization (it would order the
// kernel's memory) nor its own buffers.
struct HiddenExchange {
    HiddenExchange()
    {
        AnnotateIgnoreSyncBegin(__FILE__, __LINE__);
        AnnotateIgnoreReadsBegin(__FILE__, __LINE__);
        AnnotateIgnoreWritesBegin(__FILE__, __LINE__);
    }
    ~HiddenExchange()
    {
        AnnotateIgnoreWritesEnd(__FILE__, __LINE__);
        AnnotateIgnoreReadsEnd(__FILE__, __LINE__);
        AnnotateIgnoreSyncEnd(__FILE__, __LINE__);
    }
};
#else
struct HiddenExchange {};
#endif

extern "C" {
int64_t mir_emu_lane() { return tlLane; }
int64_t mir_emu_block() { return tlBlock; }
void mir_emu_barrier() { warp->bar.arrive_and_wait(); }

double mir_emu_shfl(double v, int32_t src)
{
    HiddenExchange hidden;
    warp->xchg[tlLane] = v;
    warp->bar.arrive_and_wait();
    const double r = warp->xchg[src & (kLanes - 1)];
    warp->bar.arrive_and_wait();   // nobody overwrites xchg before all have read
    return r;
}

// m8n8k4 f64 fragments, lane L, i = L/4, k = L%4:
//   A[i][k] = a,  B[k][i] = b,  C[i][2k] = c0,  C[i][2k+1] = c1.
double mir_emu_mma(double a, double b, double c0, double c1)
{
    HiddenExchange hidden;
    const int i = (int)tlLane / 4, k = (int)tlLane % 4;
    warp->A[i][k] = a;
    warp->B[k][i] = b;
    warp->C[i][2 * k] = c0;
    warp->C[i][2 * k + 1] = c1;
    warp->bar.arrive_and_wait();
    double d0 = warp->C[i][2 * k], d1 = warp->C[i][2 * k + 1];
    for (int q = 0; q < 4; ++q) {
        d0 += warp->A[i][q] * warp->B[q][2 * k];
        d1 += warp->A[i][q] * warp->B[q][2 * k + 1];
    }
    warp->bar.arrive_and_wait();
    tlHi = d1;
    return d0;
}
double mir_emu_mma_hi() { return tlHi; }
}  // extern "C"

int main(int argc, char** argv)
{
    const long long E = argc > 1 ? atoll(argv[1]) : 3;
    const int p = argc > 2 ? atoi(argv[2]) : 7;
    const int n = p + 1, nn = n * n, n3 = nn * n;
    const long long gElem = 3LL * p * nn * 3;

    std::vector<double> hBt((size_t)p * n), hDt((size_t)p * n), hDm(nn), hW(nn);
    srand(42);
    auto rnd = [] { return 2.0 * rand() / RAND_MAX - 1.0; };
    for (auto& x : hBt) x = rnd();
    for (auto& x : hDt) x = rnd();
    for (auto& x : hDm) x = rnd();
    for (auto& x : hW) x = rnd();

    std::vector<double> U((size_t)E * n3), G((size_t)E * gElem);
    for (long long e = 0; e < E; ++e) {
        for (int i = 0; i < n3; ++i) U[(size_t)e * n3 + i] = hashVal(e, i, 0u);
        for (long long i = 0; i < gElem; ++i)
            G[(size_t)(e * gElem + i)] = hashVal(e, (size_t)i, 0x9e37u);
    }
    // NaN, not zero: an entry of Y the kernel never writes must not pass.
    std::vector<double> Y((size_t)E * n3, std::nan(""));

    Desc<4> dU{U.data(), U.data(), 0, {E, n, n, n}, {n3, nn, n, 1}};
    Desc<4> dY{Y.data(), Y.data(), 0, {E, n, n, n}, {n3, nn, n, 1}};
    Desc<2> dBt{hBt.data(), hBt.data(), 0, {p, n}, {n, 1}};
    Desc<2> dDt{hDt.data(), hDt.data(), 0, {p, n}, {n, 1}};
    Desc<2> dW{hW.data(), hW.data(), 0, {n, n}, {n, 1}};
    Desc<2> dDm{hDm.data(), hDm.data(), 0, {n, n}, {n, 1}};
    Desc<6> dG{G.data(), G.data(), 0, {E, 3, p, n, n, 3},
               {gElem, (int64_t)p * nn * 3, (int64_t)nn * 3, (int64_t)n * 3, 3, 1}};

    // One block at a time: the kernel's workgroup buffers are process-wide
    // globals in the emulated build, so two blocks in flight would share them.
    Warp w;
    warp = &w;
    for (long long e = 0; e < E; ++e) {
        std::vector<std::thread> lanes;
        for (int L = 0; L < kLanes; ++L)
            lanes.emplace_back([&, L, e] {
                tlLane = L;
                tlBlock = e;
                _mlir_ciface_laplacian_apply(&dU, &dBt, &dDt, &dW, &dDm, &dG, &dY);
            });
        for (auto& t : lanes) t.join();
    }

    double err = 0;
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
    }
    printf("WARP EMULATION (p=%d, E=%lld, 32 lanes/elem): max|err| = %.3e  %s\n",
           p, E, err, err < 1e-11 ? "PASS" : "FAIL");
    return err < 1e-11 ? 0 : 1;
}
