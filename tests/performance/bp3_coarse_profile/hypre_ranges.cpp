#include <HYPRE_utilities.h>
#include <nvtx3/nvToolsExt.h>

#include <array>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Hypre 2.33's profiling implementation adds a stream wait on each pop. This
// diagnostic replaces its disabled hooks with annotations only, without waits.
namespace {
constexpr int max_depth = 1024;
constexpr int max_levels = 64;
struct Frame { int level; bool cycle; bool emitted; };
thread_local std::array<Frame, max_depth> frames{};
thread_local int depth = 0, level = -1;
thread_local bool in_cycle = false;
std::atomic<long long> pushes{0}, pops{0}, cycles{0}, errors{0};
std::array<std::atomic<long long>, max_levels> visits{};

nvtxDomainHandle_t domain()
{
    static auto handle = nvtxDomainCreateA("bp3.coarse");
    return handle;
}

void push(const char* name)
{
    ++pushes;
    if (depth >= max_depth) { ++errors; ++depth; return; }
    Frame previous{level, in_cycle, false};
    if (!std::strcmp(name, "AMGCycle")) { in_cycle = true; level = -1; ++cycles; }
    char label[512];
    int next_level = -1;
    if (in_cycle && std::sscanf(name, "AMG Level-%d", &next_level) == 1) {
        if (next_level < 0 || next_level >= max_levels) ++errors;
        else { level = next_level; ++visits[level]; }
        std::snprintf(label, sizeof(label), "bp3.amg/L%02d", level);
    } else if (!std::strcmp(name, "AMGCycle")) {
        std::snprintf(label, sizeof(label), "bp3.amg/cycle");
    } else {
        std::snprintf(label, sizeof(label), "bp3.amg/L%02d/%s", level, name);
    }
    if (in_cycle) {
        nvtxEventAttributes_t attr{};
        attr.version = NVTX_VERSION;
        attr.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
        attr.messageType = NVTX_MESSAGE_TYPE_ASCII;
        attr.message.ascii = label;
        nvtxDomainRangePushEx(domain(), &attr);
        previous.emitted = true;
    }
    frames[depth++] = previous;
}

__attribute__((destructor)) void audit()
{
    const char* path = std::getenv("BP3_RANGE_AUDIT");
    if (!path) return;
    FILE* file = std::fopen(path, "w");
    if (!file) { std::perror(path); return; }
    std::fprintf(file, "{\"pushes\":%lld,\"pops\":%lld,\"cycles\":%lld,\"errors\":%lld,\"levels\":[",
                 pushes.load(), pops.load(), cycles.load(), errors.load());
    for (int i = 0; i < max_levels; ++i)
        std::fprintf(file, "%s%lld", i ? "," : "", visits[i].load());
    std::fprintf(file, "]}\n");
    std::fclose(file);
}
}

extern "C" void hypre_GpuProfilingPushRange(const char* name) { push(name); }
extern "C" void hypre_GpuProfilingPushRangeColor(const char* name, HYPRE_Int) { push(name); }
extern "C" void hypre_GpuProfilingPopRange()
{
    ++pops;
    if (depth <= 0) { ++errors; return; }
    if (--depth >= max_depth) return;
    const auto previous = frames[depth];
    if (previous.emitted) nvtxDomainRangePop(domain());
    level = previous.level;
    in_cycle = previous.cycle;
}
