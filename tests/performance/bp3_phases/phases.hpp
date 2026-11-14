#pragma once
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <mpi.h>

// Optional completed-work intervals. Disabled runs add no waits or file operations.
namespace bp3_phase
{
enum Id
{
    Solve,
    CoarseEntry,
    Coarse,
    Hypre,
    Operator,
    Dss,
    Smoother,
    Cycle,
    Count
};
inline const char* names[] = {
    "solve",        "coarse_entry_wait", "coarse_complete", "hypre_complete", "operator_complete",
    "outer_dss_complete", "smoother_complete", "cycle_complete"};
struct Stat
{
    long long inclusive = 0, exclusive = 0, calls = 0;
};
struct Frame
{
    Id id;
    long long start, children;
};
inline std::array<Stat, Count> stats{};
inline std::array<Frame, 64> frames{};
inline int depth   = 0;
inline bool active = false, detail = false;
inline std::string prefix;
inline void (*sync_device)() = nullptr;
inline long long now()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch())
        .count();
}
[[noreturn]] inline void fail(const char* why)
{
    std::fprintf(stderr, "BP3 phase error: %s\n", why);
    MPI_Abort(MPI_COMM_WORLD, 93);
    std::abort();
}
inline void enter(Id id)
{
    if (depth == int(frames.size())) fail("stack overflow");
    frames[depth++] = {id, now(), 0};
}
inline void leave(Id id)
{
    const auto stop = now();
    if (!depth || frames[depth - 1].id != id) fail("unbalanced intervals");
    const auto frame   = frames[--depth];
    const auto elapsed = stop - frame.start;
    if (elapsed < frame.children) fail("invalid child accounting");
    auto& s = stats[id];
    s.inclusive += elapsed;
    s.exclusive += elapsed - frame.children;
    ++s.calls;
    if (depth) frames[depth - 1].children += elapsed;
}
struct Scope
{
    Id id;
    bool enabled, drain;
    explicit Scope(Id value, bool detailed = false, bool completed = false)
        : id(value)
        , enabled(active && (!detailed || detail))
        , drain(completed)
    {
        if (enabled)
        {
            if (drain) sync_device();
            enter(id);
        }
    }
    ~Scope()
    {
        if (enabled)
        {
            if (drain) sync_device();
            leave(id);
        }
    }
    Scope(const Scope&)            = delete;
    Scope& operator=(const Scope&) = delete;
};
inline void begin(void (*sync)())
{
    const char* path = std::getenv("BP3_PHASE_PREFIX");
    if (!path) return;
    if (!*path || active || depth) fail("invalid capture start");
    prefix      = path;
    sync_device = sync;
    detail      = std::getenv("BP3_PHASE_DETAIL") && !std::strcmp(std::getenv("BP3_PHASE_DETAIL"), "1");
    stats       = {};
    sync_device();
    active = true;
    enter(Solve);
}
inline void end()
{
    if (!active) return;
    sync_device();
    leave(Solve);
    if (depth) fail("incomplete capture");
    active = false;
}
inline void write()
{
    if (prefix.empty()) return;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string name = prefix + "-rank" + std::to_string(rank) + ".json";
    FILE* f                = std::fopen(name.c_str(), "w");
    if (!f) fail("cannot open timing receipt");
    std::fprintf(f, "{\"version\":2,\"rank\":%d,\"detail\":%s,\"regions\":{", rank, detail ? "true" : "false");
    for (int i = 0; i < Count; ++i)
        std::fprintf(f, "%s\"%s\":{\"inclusive_ns\":%lld,\"exclusive_ns\":%lld,\"calls\":%lld}", i ? "," : "", names[i],
                     stats[i].inclusive, stats[i].exclusive, stats[i].calls);
    std::fprintf(f, "}}\n");
    const bool bad = std::ferror(f);
    if (std::fclose(f) || bad) fail("cannot write timing receipt");
}
inline void entry_wait()
{
    if (active)
    {
        Scope scope(CoarseEntry);
        sync_device();
    }
}
} // namespace bp3_phase
