#pragma once
// Host runtime of the warp emulation (--mir-emulate-warp): the kernel body runs
// with one thread per lane, and these calls give it the warp semantics it relies on:
//   shuffle   every lane posts its value, all wait, each reads its source lane;
//   mma.sync  every lane posts its A/B/C fragment entries, all wait, each
//             computes its own two D entries from the assembled 8x4 @ 4x8 + 8x8;
//   barrier   all 32 lanes wait.
// Private memory is private (each thread's stack), shared memory is shared, and
// the lanes really run concurrently, so a missing barrier is a real data race.
//
// But every exchange here is a real rendezvous, which in C++ also ORDERS memory,
// so a race across an mma or a shuffle stays invisible. On the GPU neither
// shfl.sync nor mma.sync orders memory -- only a barrier does. The TSan build
// (-DMIR_EMU_TSAN, -fsanitize=thread) hides the exchanges from ThreadSanitizer,
// which then reports every pair of accesses from different lanes that no
// gpu.barrier orders: a missing barrier, under the GPU's rules.
//
// Include it in exactly one translation unit of a harness. The harness sets
// `warp`, and tlLane / tlBlock / tlGrid on every lane thread, before calling the kernel.
#include <condition_variable>
#include <cstdint>
#include <mutex>

template <typename T, int R = 1> struct Desc {
    T* allocated;
    T* aligned;
    int64_t offset;
    int64_t sizes[R];
    int64_t strides[R];
};

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
    int32_t xchg_i32[kLanes];
    double A[8][4], B[4][8], C[8][8];
};
Warp* warp = nullptr;
thread_local int64_t tlLane = 0, tlBlock = 0, tlGrid = 1;
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
int64_t mir_emu_grid() { return tlGrid; }
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

int32_t mir_emu_shfl_i32(int32_t v, int32_t src)
{
    HiddenExchange hidden;
    warp->xchg_i32[tlLane] = v;
    warp->bar.arrive_and_wait();
    const int32_t r = warp->xchg_i32[src & (kLanes - 1)];
    warp->bar.arrive_and_wait();
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
