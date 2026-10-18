// Opt-in metadata trace. Message buffers are never inspected or changed.
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

extern "C" {
#include "setup-protocol-exchange.h"
}

#if HYPRE_RELEASE_NUMBER != 23300 || defined(HYPRE_SEQUENTIAL)
#error This diagnostic requires the captured parallel Hypre 2.33 headers.
#endif

namespace {
struct Trace {
    FILE* file = nullptr;
    FILE* partition_file = nullptr;
    bool setup = false, active = false, failed = false;
    bool partition_failed = false;
    unsigned long exchange = 0, next = 0, lines = 0;
    unsigned long lookup = 0, partition_lines = 0;
    MPI_Comm comm = MPI_COMM_NULL;
    std::map<MPI_Request*, unsigned long> requests;
    std::string wait;
};
thread_local Trace trace;

template<class Function> Function next_function(const char* name) {
    // Forward to the captured library, never to a second Hypre found on the search path.
    static void* library = [] {
        const char* path = std::getenv("MARS_SETUP_PROTOCOL_HYPRE");
        return path ? dlopen(path, RTLD_NOW | RTLD_NOLOAD) : nullptr;
    }();
    void* address = library ? dlsym(library, name) : nullptr;
    if (!address) {
        std::fprintf(stderr, "Setup protocol: required symbol unavailable: %s\n", name);
        std::_Exit(125);
    }
    return reinterpret_cast<Function>(address);
}

void emit(const std::string& value) {
    if (!trace.file || trace.failed) return;
    if (++trace.lines > 100000) {
        std::fputs("{\"event\":\"truncated\"}\n", trace.file);
        std::fflush(trace.file);
        trace.failed = true;
        return;
    }
    if (std::fprintf(trace.file, "%s\n", value.c_str()) < 0 || std::fflush(trace.file)) trace.failed = true;
}

void emit_partition(const std::string& value) {
    if (!trace.partition_file || trace.partition_failed) return;
    if (++trace.partition_lines > 100000) {
        std::fputs("{\"event\":\"truncated\"}\n", trace.partition_file);
        std::fflush(trace.partition_file);
        trace.partition_failed = true;
        return;
    }
    if (std::fprintf(trace.partition_file, "%s\n", value.c_str()) < 0 || std::fflush(trace.partition_file))
        trace.partition_failed = true;
}

std::string number(long long value) { return std::to_string(value); }
const char* flag(bool value) { return value ? "true" : "false"; }
std::string event(const char* name) {
    return "{\"event\":\"" + std::string(name) + "\",\"exchange\":" + number(trace.exchange);
}

void begin_setup() {
    if (trace.setup) { emit(event("nested_setup") + "}"); return; }
    const char* directory = std::getenv("MARS_SETUP_PROTOCOL_DIR");
    const char* rank = std::getenv("SLURM_PROCID");
    if (!directory || !rank || rank[0] < '0' || rank[0] > '3' || rank[1]) return;
    if (!trace.file) {
        std::string path = std::string(directory) + "/protocol-" + rank + ".jsonl";
        trace.file = std::fopen(path.c_str(), "wx");
    }
    if (!trace.partition_file) {
        std::string path = std::string(directory) + "/partition-" + rank + ".jsonl";
        trace.partition_file = std::fopen(path.c_str(), "wx");
    }
    trace.setup = true;
    emit(event("setup_entry") + "}");
    emit_partition("{\"event\":\"setup_entry\"}");
}

unsigned long message(const char* name, HYPRE_Int peer, HYPRE_Int tag,
                      HYPRE_Int count, MPI_Comm comm, MPI_Request* request = nullptr) {
    if (!trace.active) return 0;
    const auto id = ++trace.next;
    if (request) trace.requests[request] = id;
    emit(event(name) + ",\"id\":" + number(id) + ",\"peer\":" + number(peer) +
         ",\"tag\":" + number(tag) + ",\"zero\":" + (count == 0 ? "true" : "false") +
         ",\"same_comm\":" + (comm == trace.comm ? "true" : "false") + "}");
    return id;
}

void returned(unsigned long id, HYPRE_Int code) {
    if (id) emit(event("call_return") + ",\"id\":" + number(id) +
                 ",\"ok\":" + (code == MPI_SUCCESS ? "true" : "false") + "}");
}

std::string waiting(const char* name, HYPRE_Int count, MPI_Request* requests) {
    if (!trace.active) return {};
    if (count < 0 || count > 100000) { emit(event("invalid_request_count") + "}"); return {}; }
    std::string ids = "[";
    for (HYPRE_Int i = 0; i < count; ++i) {
        if (i) ids += ',';
        auto found = trace.requests.find(requests + i);
        ids += number(found == trace.requests.end() ? 0 : found->second);
    }
    ids += ']';
    std::string value = event(name) + ",\"requests\":" + ids + "}";
    if (value != trace.wait) { emit(value); trace.wait = value; }
    return ids;
}

void tested(const std::string& ids, HYPRE_Int code, bool complete) {
    if (!ids.empty() && (code != MPI_SUCCESS || complete)) {
        emit(event("tested") + ",\"requests\":" + ids +
             ",\"ok\":" + (code == MPI_SUCCESS ? "true" : "false") +
             ",\"complete\":" + (complete ? "true" : "false") + "}");
        trace.wait.clear();
    }
}
}

extern "C" HYPRE_Int HYPRE_ParCSRGMRESSetup(HYPRE_Solver solver, HYPRE_ParCSRMatrix matrix,
                                           HYPRE_ParVector rhs, HYPRE_ParVector x) {
    static auto real = next_function<decltype(&HYPRE_ParCSRGMRESSetup)>("HYPRE_ParCSRGMRESSetup");
    begin_setup();
    const auto code = real(solver, matrix, rhs, x);
    emit_partition("{\"event\":\"setup_return\"}");
    emit(event("setup_return") + "}"); trace.setup = false;
    return code;
}

extern "C" HYPRE_Int HYPRE_ParCSRFlexGMRESSetup(HYPRE_Solver solver, HYPRE_ParCSRMatrix matrix,
                                               HYPRE_ParVector rhs, HYPRE_ParVector x) {
    static auto real = next_function<decltype(&HYPRE_ParCSRFlexGMRESSetup)>("HYPRE_ParCSRFlexGMRESSetup");
    begin_setup();
    const auto code = real(solver, matrix, rhs, x);
    emit_partition("{\"event\":\"setup_return\"}");
    emit(event("setup_return") + "}"); trace.setup = false;
    return code;
}

extern "C" HYPRE_Int hypre_GetAssumedPartitionProcFromRow(MPI_Comm comm, HYPRE_BigInt row,
    HYPRE_BigInt first, HYPRE_BigInt count, HYPRE_Int* owner) {
    using Function = HYPRE_Int(*)(MPI_Comm,HYPRE_BigInt,HYPRE_BigInt,HYPRE_BigInt,HYPRE_Int*);
    static auto real = next_function<Function>("hypre_GetAssumedPartitionProcFromRow");
    if (!trace.setup) return real(comm,row,first,count,owner);
    int size = 0, comparison = MPI_UNEQUAL;
    const bool size_ok = PMPI_Comm_size(comm, &size) == MPI_SUCCESS && size > 0;
    const bool world = PMPI_Comm_compare(comm, MPI_COMM_WORLD, &comparison) == MPI_SUCCESS && comparison == MPI_IDENT;
    // Unsigned subtraction avoids overflow even for malformed signed indices.
    const bool in_range = count > 0 && row >= first &&
        static_cast<unsigned long long>(row) - static_cast<unsigned long long>(first) <
        static_cast<unsigned long long>(count);
    const auto id = ++trace.lookup;
    const std::string key = ",\"id\":" + number(id);
    emit_partition("{\"event\":\"lookup\"" + key + ",\"input_in_range\":" + flag(in_range) +
                   ",\"size_query_ok\":" + flag(size_ok) + ",\"comm_world\":" + flag(world) + "}");
    const auto code = real(comm,row,first,count,owner);
    emit_partition("{\"event\":\"lookup_return\"" + key + ",\"owner_in_range\":" +
                   flag(size_ok && owner && *owner >= 0 && *owner < size) +
                   ",\"hypre_error_nonzero\":" + flag(code != 0) + "}");
    return code;
}

extern "C" HYPRE_Int hypre_DataExchangeList(HYPRE_Int contacts, HYPRE_Int* peers, void* buffer,
    HYPRE_Int* starts, HYPRE_Int contact_size, HYPRE_Int response_size, hypre_DataExchangeResponse* response,
    HYPRE_Int maximum, HYPRE_Int round, MPI_Comm comm, void** received, HYPRE_Int** offsets) {
    static auto real = next_function<decltype(&hypre_DataExchangeList)>("hypre_DataExchangeList");
    if (trace.active) {
        emit(event("nested_exchange") + "}");
        return real(contacts,peers,buffer,starts,contact_size,response_size,response,maximum,round,comm,received,offsets);
    }
    const bool recording = trace.setup;
    if (recording) {
        ++trace.exchange; trace.active = true; trace.comm = comm;
        trace.requests.clear(); trace.wait.clear();
        int same = MPI_UNEQUAL, rank = -1, size = -1;
        bool ok = PMPI_Comm_compare(comm, MPI_COMM_WORLD, &same) == MPI_SUCCESS;
        ok = PMPI_Comm_rank(comm, &rank) == MPI_SUCCESS && ok;
        ok = PMPI_Comm_size(comm, &size) == MPI_SUCCESS && ok;
        emit(event("exchange_entry") + ",\"round\":" + number(round) + ",\"rank\":" + number(rank) +
             ",\"ranks\":" + number(size) + ",\"comm_world\":" +
             (ok && same == MPI_IDENT ? "true" : "false") + "}");
    }
    const auto code = real(contacts,peers,buffer,starts,contact_size,response_size,response,maximum,round,comm,received,offsets);
    if (recording) { emit(event("exchange_return") + "}"); trace.active = false; }
    return code;
}

#define ASYNC_MESSAGE(NAME, EVENT) \
extern "C" HYPRE_Int NAME(void* b,HYPRE_Int n,MPI_Datatype t,HYPRE_Int p,HYPRE_Int tag,MPI_Comm c,MPI_Request* r) { \
    static auto real = next_function<HYPRE_Int(*)(void*,HYPRE_Int,MPI_Datatype,HYPRE_Int,HYPRE_Int,MPI_Comm,MPI_Request*)>(#NAME); \
    auto id = message(EVENT,p,tag,n,c,r); auto code = real(b,n,t,p,tag,c,r); returned(id,code); return code; \
}
ASYNC_MESSAGE(hypre_MPI_Isend, "isend")
ASYNC_MESSAGE(hypre_MPI_Irecv, "irecv")

extern "C" HYPRE_Int hypre_MPI_Send(void* b,HYPRE_Int n,MPI_Datatype t,HYPRE_Int p,HYPRE_Int tag,MPI_Comm c) {
    static auto real = next_function<HYPRE_Int(*)(void*,HYPRE_Int,MPI_Datatype,HYPRE_Int,HYPRE_Int,MPI_Comm)>("hypre_MPI_Send");
    auto id = message("send",p,tag,n,c); auto code = real(b,n,t,p,tag,c); returned(id,code); return code;
}

extern "C" HYPRE_Int hypre_MPI_Recv(void* b,HYPRE_Int n,MPI_Datatype t,HYPRE_Int p,HYPRE_Int tag,MPI_Comm c,MPI_Status* s) {
    static auto real = next_function<HYPRE_Int(*)(void*,HYPRE_Int,MPI_Datatype,HYPRE_Int,HYPRE_Int,MPI_Comm,MPI_Status*)>("hypre_MPI_Recv");
    auto id = message("recv",p,tag,n,c); auto code = real(b,n,t,p,tag,c,s); returned(id,code); return code;
}

extern "C" HYPRE_Int hypre_MPI_Testall(HYPRE_Int n,MPI_Request* r,HYPRE_Int* f,MPI_Status* s) {
    static auto real = next_function<HYPRE_Int(*)(HYPRE_Int,MPI_Request*,HYPRE_Int*,MPI_Status*)>("hypre_MPI_Testall");
    auto ids = waiting("testall",n,r); auto code = real(n,r,f,s); tested(ids,code,code == MPI_SUCCESS && *f); return code;
}

extern "C" HYPRE_Int hypre_MPI_Test(MPI_Request* r,HYPRE_Int* f,MPI_Status* s) {
    static auto real = next_function<HYPRE_Int(*)(MPI_Request*,HYPRE_Int*,MPI_Status*)>("hypre_MPI_Test");
    auto ids = waiting("test",1,r); auto code = real(r,f,s); tested(ids,code,code == MPI_SUCCESS && *f); return code;
}

extern "C" HYPRE_Int hypre_MPI_Waitall(HYPRE_Int n,MPI_Request* r,MPI_Status* s) {
    static auto real = next_function<HYPRE_Int(*)(HYPRE_Int,MPI_Request*,MPI_Status*)>("hypre_MPI_Waitall");
    auto ids = waiting("waitall",n,r); auto code = real(n,r,s); tested(ids,code,code == MPI_SUCCESS); return code;
}

extern "C" HYPRE_Int hypre_MPI_Wait(MPI_Request* r,MPI_Status* s) {
    static auto real = next_function<HYPRE_Int(*)(MPI_Request*,MPI_Status*)>("hypre_MPI_Wait");
    auto ids = waiting("wait",1,r); auto code = real(r,s); tested(ids,code,code == MPI_SUCCESS); return code;
}
