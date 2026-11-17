#include <mpi.h>
#include <impl/Kokkos_Profiling_C_Interface.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <unistd.h>
#include <cstdint>
#include <atomic>
#ifdef MARS_DIAG_CUDA
#include <cuda_runtime_api.h>
#include <dlfcn.h>
#endif

namespace {
std::atomic<bool> stopped{false};
using Clock = std::chrono::steady_clock;
double seconds() { return std::chrono::duration<double>(Clock::now().time_since_epoch()).count(); }
bool enabled(const char* name) { const char* v = std::getenv(name); return v && std::string(v) == "1"; }
std::string quote(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else out << c;
    }
    return out.str() + '"';
}
struct Event {
    std::string kind, name;
    double start, duration;
    long long count;
    int peer, tag, error;
    double value;
    std::uintptr_t pointer;
};
struct State {
    std::mutex mutex;
    std::vector<Event> events;
    std::map<uint64_t, std::tuple<std::string, std::string, double>> active;
    uint64_t next = 0;
    size_t dropped = 0;
    int rank = -1;
    bool tools = false, tools_finalized = false;
};
State& state() { static State value; return value; }
int rank() {
    auto& s = state();
    if (s.rank < 0) {
        int initialized = 0, finalized = 0;
        PMPI_Initialized(&initialized);
        if (initialized) PMPI_Finalized(&finalized);
        if (initialized && !finalized) PMPI_Comm_rank(MPI_COMM_WORLD, &s.rank);
    }
    return s.rank;
}
void event(std::string kind, std::string name, double start, double duration = 0,
           long long count = 0, int peer = -1, int tag = -1, int error = 0,
           double value = 0, const void* pointer = nullptr) {
    if (stopped.load(std::memory_order_relaxed)) return;
    auto& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.events.size() < 500000)
        s.events.push_back({kind, name, start, duration, count, peer, tag, error, value,
                            reinterpret_cast<std::uintptr_t>(pointer)});
    else ++s.dropped;
}
void flush(const char* status) {
    const char* prefix = std::getenv("MARS_DIAG_PREFIX");
    if (!prefix) return;
    int process_rank = rank();
    auto& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    std::ofstream out(std::string(prefix) + ".rank" + std::to_string(process_rank) + ".json");
    char hostname[256] = {};
    gethostname(hostname, sizeof(hostname) - 1);
    out << std::setprecision(17) << "{\"status\":" << quote(status)
        << ",\"rank\":" << process_rank << ",\"hostname\":" << quote(hostname)
        << ",\"pid\":" << getpid() << ",\"tools_loaded\":" << (s.tools ? "true" : "false")
        << ",\"tools_finalized\":" << (s.tools_finalized ? "true" : "false")
        << ",\"global_fencing\":" << (enabled("MARS_DIAG_FENCE") ? "true" : "false")
        << ",\"waitall_retains_statuses\":" << (enabled("MARS_DIAG_RETAIN_STATUSES") ? "true" : "false")
        << ",\"cuda_trace_requested\":" << (enabled("MARS_DIAG_CUDA") ? "true" : "false")
        << ",\"dropped\":" << s.dropped << ",\"events\":[";
    bool comma = false;
    for (const auto& e : s.events) {
        if (comma) out << ',';
        comma = true;
        out << "{\"kind\":" << quote(e.kind) << ",\"name\":" << quote(e.name)
            << ",\"start\":" << e.start << ",\"seconds\":" << e.duration
            << ",\"count\":" << e.count << ",\"peer\":" << e.peer
            << ",\"tag\":" << e.tag << ",\"error\":" << e.error
            << ",\"value\":" << e.value << ",\"pointer\":" << quote(std::to_string(e.pointer)) << '}';
    }
    out << "]}\n";
    if (!out) std::fprintf(stderr, "MARS_DIAGNOSTIC_WRITE_FAILED rank=%d\n", process_rank);
}
void mpi_error(const char* operation, int error, long long count, int peer, int tag) {
    char buffer[MPI_MAX_ERROR_STRING] = {};
    int length = 0, error_class = 0;
    PMPI_Error_class(error, &error_class);
    PMPI_Error_string(error, buffer, &length);
    std::string message(buffer, length);
    std::fprintf(stderr, "MARS_MPI_ERROR {\"rank\":%d,\"operation\":%s,\"code\":%d,\"class\":%d,\"count\":%lld,\"peer\":%d,\"tag\":%d,\"message\":%s}\n",
                 rank(), quote(operation).c_str(), error, error_class, count, peer, tag, quote(message).c_str());
    std::fflush(stderr);
    event("mpi_error", std::string(operation) + ": " + message, seconds(), 0, count, peer, tag, error);
    flush("mpi_error");
}
double before(const char* operation, long long count = 0, int peer = -1, int tag = -1) {
    if (enabled("MARS_DIAG_MPI_TRACE")) {
        std::fprintf(stderr, "MARS_MPI_BEGIN rank=%d operation=%s count=%lld peer=%d tag=%d\n", rank(), operation, count, peer, tag);
        std::fflush(stderr);
    }
    return seconds();
}
int after(const char* operation, int error, double start, long long count = 0, int peer = -1, int tag = -1,
          const void* pointer = nullptr) {
    event("mpi", operation, start, seconds() - start, count, peer, tag, error, 0, pointer);
    if (error != MPI_SUCCESS) mpi_error(operation, error, count, peer, tag);
    return error;
}
void begin(const char* kind, const char* name, uint64_t* id) {
    auto& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    *id = ++s.next;
    s.active[*id] = {kind, name, seconds()};
}
void end(uint64_t id) {
    double stop = seconds();
    std::tuple<std::string, std::string, double> item;
    {
        auto& s = state();
        std::lock_guard<std::mutex> lock(s.mutex);
        auto found = s.active.find(id);
        if (found == s.active.end()) return;
        item = found->second;
        s.active.erase(found);
    }
    event(std::get<0>(item), std::get<1>(item), std::get<2>(item), stop - std::get<2>(item));
}
thread_local std::vector<uint64_t> copy_stack, region_stack;
}

extern "C" {
void kokkosp_init_library(const int, const uint64_t, const uint32_t, Kokkos_Profiling_KokkosPDeviceInfo*) {
    rank();
    state().tools = true;
    event("tool", "init", seconds());
}
void kokkosp_request_tool_settings(const uint32_t, Kokkos_Tools_ToolSettings* settings) {
    settings->requires_global_fencing = enabled("MARS_DIAG_FENCE");
}
void kokkosp_finalize_library() {
    state().tools_finalized = true;
    event("tool", "finalize", seconds());
    flush("kokkos_finalized");
}
#define KERNEL_CALLBACKS(kind) \
void kokkosp_begin_##kind(const char* name, const uint32_t, uint64_t* id) { begin(#kind, name, id); } \
void kokkosp_end_##kind(uint64_t id) { end(id); }
KERNEL_CALLBACKS(parallel_for)
KERNEL_CALLBACKS(parallel_reduce)
KERNEL_CALLBACKS(parallel_scan)
KERNEL_CALLBACKS(fence)
#undef KERNEL_CALLBACKS
void kokkosp_begin_deep_copy(Kokkos_Profiling_SpaceHandle dst, const char* dst_name, const void*,
                            Kokkos_Profiling_SpaceHandle src, const char* src_name, const void*, uint64_t bytes) {
    uint64_t id;
    std::string label = std::string(src.name) + ":" + src_name + " -> " + dst.name + ":" + dst_name;
    begin("deep_copy", label.c_str(), &id);
    copy_stack.push_back(id);
    event("copy_bytes", label, seconds(), 0, bytes);
}
void kokkosp_end_deep_copy() {
    if (!copy_stack.empty()) { end(copy_stack.back()); copy_stack.pop_back(); }
}
void kokkosp_allocate_data(Kokkos_Profiling_SpaceHandle space, const char* name, const void*, uint64_t bytes) {
    event("allocate", std::string(space.name) + ":" + name, seconds(), 0, bytes);
}
void kokkosp_deallocate_data(Kokkos_Profiling_SpaceHandle space, const char* name, const void*, uint64_t bytes) {
    event("deallocate", std::string(space.name) + ":" + name, seconds(), 0, bytes);
}
void kokkosp_push_profile_region(const char* name) {
    uint64_t id;
    begin("region", name, &id);
    region_stack.push_back(id);
}
void kokkosp_pop_profile_region() {
    if (!region_stack.empty()) { end(region_stack.back()); region_stack.pop_back(); }
}
int MPI_Isend(const void* buf, int count, MPI_Datatype type, int peer, int tag, MPI_Comm comm, MPI_Request* request) {
    double t = before("MPI_Isend", count, peer, tag);
    return after("MPI_Isend", PMPI_Isend(buf, count, type, peer, tag, comm, request), t, count, peer, tag, buf);
}
int MPI_Irecv(void* buf, int count, MPI_Datatype type, int peer, int tag, MPI_Comm comm, MPI_Request* request) {
    double t = before("MPI_Irecv", count, peer, tag);
    return after("MPI_Irecv", PMPI_Irecv(buf, count, type, peer, tag, comm, request), t, count, peer, tag, buf);
}
double MPI_Wtime() {
    double start = seconds();
    double value = PMPI_Wtime();
    event("timer", "MPI_Wtime", start, seconds() - start, 0, -1, -1, 0, value);
    return value;
}
int MPI_Waitall(int count, MPI_Request* requests, MPI_Status* statuses) {
    // Status replacement is opt-in: it can change the MPI library's wait path.
    std::vector<MPI_Status> retained;
    if (count > 0 && statuses == MPI_STATUSES_IGNORE && enabled("MARS_DIAG_RETAIN_STATUSES")) {
        retained.resize(count);
        statuses = retained.data();
    }
    double t = before("MPI_Waitall", count);
    int error = PMPI_Waitall(count, requests, statuses);
    int error_class = MPI_SUCCESS;
    if (error != MPI_SUCCESS) PMPI_Error_class(error, &error_class);
    if (error_class == MPI_ERR_IN_STATUS && statuses != MPI_STATUSES_IGNORE) {
        for (int i = 0; i < count; ++i) {
            if (statuses[i].MPI_ERROR != MPI_SUCCESS && statuses[i].MPI_ERROR != MPI_ERR_PENDING) {
                std::string operation = "MPI_Waitall request " + std::to_string(i);
                mpi_error(operation.c_str(), statuses[i].MPI_ERROR, 1, -1, -1);
            }
        }
    }
    return after("MPI_Waitall", error, t, count);
}
int MPI_Allreduce(const void* send, void* recv, int count, MPI_Datatype type, MPI_Op op, MPI_Comm comm) {
    double t = before("MPI_Allreduce", count);
    return after("MPI_Allreduce", PMPI_Allreduce(send, recv, count, type, op, comm), t, count);
}
int MPI_Barrier(MPI_Comm comm) {
    double t = before("MPI_Barrier");
    return after("MPI_Barrier", PMPI_Barrier(comm), t);
}
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm* other) {
    double t = before("MPI_Comm_dup");
    return after("MPI_Comm_dup", PMPI_Comm_dup(comm, other), t);
}
int MPI_Comm_free(MPI_Comm* comm) {
    double t = before("MPI_Comm_free");
    return after("MPI_Comm_free", PMPI_Comm_free(comm), t);
}
int MPI_Dist_graph_create(MPI_Comm comm, int n, const int sources[], const int degrees[],
                          const int destinations[], const int weights[], MPI_Info info, int reorder, MPI_Comm* graph) {
    double t = before("MPI_Dist_graph_create", n);
    return after("MPI_Dist_graph_create", PMPI_Dist_graph_create(comm, n, sources, degrees, destinations,
                 weights, info, reorder, graph), t, n);
}
int MPI_Abort(MPI_Comm comm, int error) {
    mpi_error("MPI_Abort", error, 0, -1, -1);
    flush("aborting");
    return PMPI_Abort(comm, error);
}
int MPI_Finalize() {
    rank();
    int error = PMPI_Finalize();
    event("mpi", "MPI_Finalize", seconds(), 0, 0, -1, -1, error);
    flush(error == MPI_SUCCESS ? "complete" : "finalize_failed");
    stopped.store(true, std::memory_order_relaxed);
    return error;
}
}

#ifdef MARS_DIAG_CUDA
// Forward the original calls without adding CUDA operations or synchronization.
#define CUDA_CALL(name, arguments, bytes, pointer) \
    static const auto call = reinterpret_cast<decltype(&name)>(dlsym(RTLD_NEXT, #name)); \
    if (!call) { std::fprintf(stderr, "Missing CUDA symbol: %s\n", #name); std::abort(); } \
    static const bool trace = enabled("MARS_DIAG_CUDA"); \
    const double start = trace ? seconds() : 0; \
    const auto error = call arguments; \
    if (trace) event("cuda", #name, start, seconds() - start, bytes, -1, -1, error, 0, pointer); \
    return error

extern "C" {
cudaError_t cudaMalloc(void** pointer, size_t bytes) {
    CUDA_CALL(cudaMalloc, (pointer, bytes), bytes, error == cudaSuccess ? *pointer : nullptr);
}
cudaError_t cudaFree(void* pointer) { CUDA_CALL(cudaFree, (pointer), 0, pointer); }
cudaError_t cudaMallocHost(void** pointer, size_t bytes) {
    CUDA_CALL(cudaMallocHost, (pointer, bytes), bytes, error == cudaSuccess ? *pointer : nullptr);
}
cudaError_t cudaHostAlloc(void** pointer, size_t bytes, unsigned int flags) {
    CUDA_CALL(cudaHostAlloc, (pointer, bytes, flags), bytes, error == cudaSuccess ? *pointer : nullptr);
}
cudaError_t cudaFreeHost(void* pointer) { CUDA_CALL(cudaFreeHost, (pointer), 0, pointer); }
cudaError_t cudaMallocAsync(void** pointer, size_t bytes, cudaStream_t stream) {
    CUDA_CALL(cudaMallocAsync, (pointer, bytes, stream), bytes, error == cudaSuccess ? *pointer : nullptr);
}
cudaError_t cudaFreeAsync(void* pointer, cudaStream_t stream) {
    CUDA_CALL(cudaFreeAsync, (pointer, stream), 0, pointer);
}
cudaError_t cudaDeviceSynchronize() { CUDA_CALL(cudaDeviceSynchronize, (), 0, nullptr); }
cudaError_t cudaStreamSynchronize(cudaStream_t stream) { CUDA_CALL(cudaStreamSynchronize, (stream), 0, stream); }
cudaError_t cudaEventSynchronize(cudaEvent_t value) { CUDA_CALL(cudaEventSynchronize, (value), 0, value); }
cudaError_t cudaMemcpy(void* dst, const void* src, size_t bytes, cudaMemcpyKind kind) {
    CUDA_CALL(cudaMemcpy, (dst, src, bytes, kind), bytes, dst);
}
cudaError_t cudaMemcpyAsync(void* dst, const void* src, size_t bytes, cudaMemcpyKind kind, cudaStream_t stream) {
    CUDA_CALL(cudaMemcpyAsync, (dst, src, bytes, kind, stream), bytes, dst);
}
cudaError_t cudaMemset(void* dst, int value, size_t bytes) { CUDA_CALL(cudaMemset, (dst, value, bytes), bytes, dst); }
cudaError_t cudaMemsetAsync(void* dst, int value, size_t bytes, cudaStream_t stream) {
    CUDA_CALL(cudaMemsetAsync, (dst, value, bytes, stream), bytes, dst);
}
cudaError_t cudaLaunchKernel(const void* function, dim3 grid, dim3 block, void** arguments, size_t shared, cudaStream_t stream) {
    CUDA_CALL(cudaLaunchKernel, (function, grid, block, arguments, shared, stream), 0, function);
}
}
#undef CUDA_CALL
#endif
