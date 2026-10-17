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
// Include it in exactly one translation unit of a harness. The harness sets `warp`
// (with bar.threads = the kernel's block size, <kernel>__threads()), and tlLane (the
// thread of the block) / tlBlock / tlGrid on every thread, before calling the kernel.
// Blocks of several warps: barriers span the block, shuffles stay in the warp; the
// mma emulation assumes one warp.
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
    int count = 0, threads = kLanes;   // threads of the block
    unsigned gen = 0;
    void arrive_and_wait()
    {
        std::unique_lock<std::mutex> lk(m);
        const unsigned g = gen;
        if (++count == threads) {
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
    double xchg[1024];       // by thread of the block; a shuffle reads within its warp
    int32_t xchg_i32[1024];
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
    const double r = warp->xchg[(tlLane & ~(int64_t)(kLanes - 1)) | (src & (kLanes - 1))];
    warp->bar.arrive_and_wait();   // nobody overwrites xchg before all have read
    return r;
}

int32_t mir_emu_shfl_i32(int32_t v, int32_t src)
{
    HiddenExchange hidden;
    warp->xchg_i32[tlLane] = v;
    warp->bar.arrive_and_wait();
    const int32_t r = warp->xchg_i32[(tlLane & ~(int64_t)(kLanes - 1)) | (src & (kLanes - 1))];
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

#ifdef MIR_EMU_TRACE
// --mir-trace-shared: every lane logs its shared accesses in order; after a block,
// trace_block() groups the k-th execution of each access by the lanes of a warp
// into one warp instruction and counts its wavefronts: 32 banks of 4 bytes; an
// 8-byte access is served per half warp, a 16-byte one per quarter; lanes reading
// the same word share it. ideal is the count without bank conflicts.
#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace {
struct SharedAccess {
    int64_t id, byte, bytes, store;
};
struct TraceCount {
    int64_t instructions = 0, wavefronts = 0, ideal = 0, bytes = 0, store = 0;
};
std::vector<SharedAccess> g_trace[1024];
std::map<int64_t, TraceCount> g_counts;

void trace_instruction(const std::vector<std::pair<int, int64_t>>& lanes, int64_t bytes, TraceCount& c)
{
    const int span = bytes <= 4 ? 32 : 128 / (int)bytes;   // lanes served per pass
    const int64_t words = bytes <= 4 ? 1 : bytes / 4;
    std::map<int, std::map<int, std::set<int64_t>>> pass;   // pass -> bank -> words
    for (auto [lane, byte] : lanes)
        for (int64_t w = 0; w < words; ++w) {
            const int64_t word = byte / 4 + w;
            pass[lane % 32 / span][(int)(word % 32)].insert(word);
        }
    for (auto& [p, banks] : pass) {
        size_t worst = 0;
        for (auto& [bank, set] : banks) worst = std::max(worst, set.size());
        c.wavefronts += (int64_t)worst;
        c.ideal += 1;
    }
    c.instructions += 1;
}
}  // namespace

extern "C" void __mir_shared_trace(int64_t id, int64_t offset, int64_t bytes, int64_t store)
{
    g_trace[tlLane].push_back({id, offset * bytes, bytes, store});
}

void trace_block(int threads)
{
    for (int w = 0; w * kLanes < threads; ++w) {
        std::map<std::pair<int64_t, int64_t>, std::vector<std::pair<int, int64_t>>> groups;
        std::map<int64_t, std::pair<int64_t, int64_t>> kind;   // id -> (bytes, store)
        for (int l = w * kLanes; l < std::min(threads, (w + 1) * kLanes); ++l) {
            std::map<int64_t, int64_t> seen;
            for (const SharedAccess& a : g_trace[l]) {
                groups[{a.id, seen[a.id]++}].push_back({l, a.byte});
                kind[a.id] = {a.bytes, a.store};
            }
        }
        for (auto& [key, lanes] : groups) {
            TraceCount& c = g_counts[key.first];
            c.bytes = kind[key.first].first;
            c.store = kind[key.first].second;
            trace_instruction(lanes, c.bytes, c);
        }
    }
    for (int l = 0; l < threads; ++l) g_trace[l].clear();
}

// Per block: the accesses with conflicts, worst first, then the totals.
void trace_report()
{
    std::vector<std::pair<int64_t, TraceCount>> rows(g_counts.begin(), g_counts.end());
    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) {
        return a.second.wavefronts - a.second.ideal > b.second.wavefronts - b.second.ideal;
    });
    TraceCount load, store;
    printf("  shared accesses of one block (trace id: warp instructions, wavefronts, conflict-free)\n");
    for (auto& [id, c] : rows) {
        TraceCount& t = c.store ? store : load;
        t.instructions += c.instructions;
        t.wavefronts += c.wavefronts;
        t.ideal += c.ideal;
        if (c.wavefronts > c.ideal)
            printf("    id %4lld %s %lldB: %6lld instr  %7lld wavefronts  %7lld ideal\n", (long long)id,
                   c.store ? "st" : "ld", (long long)c.bytes, (long long)c.instructions, (long long)c.wavefronts,
                   (long long)c.ideal);
    }
    printf("  loads:  %lld instr, %lld wavefronts, %lld conflict-free\n", (long long)load.instructions,
           (long long)load.wavefronts, (long long)load.ideal);
    printf("  stores: %lld instr, %lld wavefronts, %lld conflict-free\n", (long long)store.instructions,
           (long long)store.wavefronts, (long long)store.ideal);
}
#endif
