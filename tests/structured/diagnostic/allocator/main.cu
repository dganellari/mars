#include <Kokkos_Core.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
int rank = -1;
void mpi_check(int error, const char* operation) {
    if (error == MPI_SUCCESS) return;
    char message[MPI_MAX_ERROR_STRING]{};
    int length = 0;
    MPI_Error_string(error, message, &length);
    std::ostringstream out;
    out << "MARS_ALLOCATOR_ERROR {\"rank\":" << rank << ",\"operation\":" << std::quoted(operation)
        << ",\"message\":" << std::quoted(std::string(message, length)) << "}\n";
    std::cerr << out.str() << std::flush;
    MPI_Abort(MPI_COMM_WORLD, error);
    std::abort();
}
void require(bool valid, const char* message) {
    if (valid) return;
    std::fprintf(stderr, "allocator probe rank=%d: %s\n", rank, message);
    MPI_Abort(MPI_COMM_WORLD, 2);
    std::abort();
}
void cuda_check(cudaError_t error) { require(error == cudaSuccess, cudaGetErrorString(error)); }
constexpr unsigned long sentinel = ~0UL;
unsigned long value(int sender, int iteration, int index) {
    return 1000000UL * sender + 10000UL * iteration + index;
}

struct Buffer {
    std::string mode;
    void* allocation = nullptr;
    unsigned long* data = nullptr;
    Kokkos::View<unsigned long*, Kokkos::CudaSpace> view;
    Buffer(const std::string& kind, int words) : mode(kind) {
        const size_t bytes = words * sizeof(unsigned long);
        if (mode == "view") {
            view = decltype(view)("mpi_allocator_probe", words);
            data = view.data();
        } else if (mode == "host") {
            allocation = std::malloc(bytes);
            data = static_cast<unsigned long*>(allocation);
            require(data != nullptr, "host allocation failed");
        } else {
            // Match the 128-byte header offset of the pinned managed Kokkos View.
            if (mode == "async") cuda_check(cudaMallocAsync(&allocation, bytes + 128, nullptr));
            else cuda_check(cudaMalloc(&allocation, bytes + 128));
            data = reinterpret_cast<unsigned long*>(static_cast<char*>(allocation) + 128);
        }
    }
    ~Buffer() {
        if (mode == "host") std::free(allocation);
        else if (mode == "async") {
            cuda_check(cudaFreeAsync(allocation, nullptr));
            cuda_check(cudaDeviceSynchronize());
        } else if (mode == "cuda") cuda_check(cudaFree(allocation));
    }
    void put(const std::vector<unsigned long>& values) {
        if (mode == "host") std::copy(values.begin(), values.end(), data);
        else cuda_check(cudaMemcpy(data, values.data(), values.size() * sizeof(unsigned long), cudaMemcpyHostToDevice));
    }
    std::vector<unsigned long> get(int words) const {
        std::vector<unsigned long> values(words);
        if (mode == "host") std::copy(data, data + words, values.begin());
        else cuda_check(cudaMemcpy(values.data(), data, words * sizeof(unsigned long), cudaMemcpyDeviceToHost));
        return values;
    }
    void receipt(const char* label, const std::string& context, int cycle, int words) const {
        int memory_type = 0, device = -1, pool_query = -1;
        CUmemoryPool pool = nullptr;
        if (mode != "host") {
            cudaPointerAttributes attributes{};
            cuda_check(cudaPointerGetAttributes(&attributes, data));
            memory_type = attributes.type;
            device = attributes.device;
            pool_query = cuPointerGetAttribute(&pool, CU_POINTER_ATTRIBUTE_MEMPOOL_HANDLE,
                                               reinterpret_cast<CUdeviceptr>(data));
            require(memory_type == cudaMemoryTypeDevice, "buffer is not CUDA device memory");
        }
        std::ostringstream out;
        out << "MARS_ALLOCATOR_BUFFER {" << context << ",\"rank\":" << rank
            << ",\"cycle\":" << cycle << ",\"label\":" << std::quoted(label)
            << ",\"words\":" << words << ",\"pointer\":" << std::quoted(pointer(data))
            << ",\"memory_type\":" << memory_type << ",\"device\":" << device
            << ",\"pool_query\":" << pool_query << ",\"pooled\":" << (pool ? "true" : "false") << "}\n";
        std::cout << out.str() << std::flush;
    }
    static std::string pointer(const void* data) { std::ostringstream out; out << data; return out.str(); }
};
} // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mpi_check(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN), "world error handler");
    require(argc == 5, "usage: allocator_probe host|cuda|async|view words count offset");
    const std::string mode = argv[1];
    const int words = std::atoi(argv[2]), count = std::atoi(argv[3]), offset = std::atoi(argv[4]);
    int ranks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    require((mode == "host" || mode == "cuda" || mode == "async" || mode == "view") &&
            (ranks == 2 || ranks == 4) && words > 0 && words <= 65536 && count > 0 && offset >= 0 &&
            static_cast<long long>(offset) + static_cast<long long>(ranks) * count <= words,
            "invalid mode or exchange bounds");
    MPI_Comm local = MPI_COMM_NULL;
    mpi_check(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local), "local communicator");
    int local_rank = 0, devices = 0;
    MPI_Comm_rank(local, &local_rank);
    cuda_check(cudaGetDeviceCount(&devices));
    require(devices > 0, "no visible CUDA device");
    const int device = devices == 1 ? 0 : local_rank % devices;
    cuda_check(cudaSetDevice(device));
    if (mode == "view") Kokkos::initialize(Kokkos::InitializationSettings().set_device_id(device));
    char hostname[MPI_MAX_PROCESSOR_NAME]{}, bus[64]{};
    int hostname_length = 0;
    MPI_Get_processor_name(hostname, &hostname_length);
    cuda_check(cudaDeviceGetPCIBusId(bus, sizeof(bus), device));
    std::ostringstream context;
    context << "\"mode\":" << std::quoted(mode) << ",\"hostname\":" << std::quoted(std::string(hostname, hostname_length))
            << ",\"device_id\":" << std::quoted(bus) << ",\"ranks\":" << ranks
            << ",\"count\":" << count << ",\"offset\":" << offset
            << ",\"kokkos_async\":"
#ifdef KOKKOS_ENABLE_IMPL_CUDA_MALLOC_ASYNC
            << "true";
#else
            << "false";
#endif
    std::vector<int> peers(ranks), weights(ranks, count);
    for (int p = 0; p < ranks; ++p) peers[p] = p;
    MPI_Comm graph = MPI_COMM_NULL;
    mpi_check(MPI_Dist_graph_create(MPI_COMM_WORLD, 1, &rank, &ranks, peers.data(), weights.data(),
                                    MPI_INFO_NULL, 0, &graph), "graph communicator");
    mpi_check(MPI_Comm_set_errhandler(graph, MPI_ERRORS_RETURN), "graph error handler");
    for (int cycle = 0; cycle < 3; ++cycle) {
        Buffer send(mode, words), recv(mode, words);
        cuda_check(cudaDeviceSynchronize());
        send.receipt("send", context.str(), cycle, words);
        recv.receipt("recv", context.str(), cycle, words);
        for (int iteration = 0; iteration < 3; ++iteration) {
            const int step = cycle * 3 + iteration;
            std::vector<unsigned long> values(words);
            for (int i = 0; i < words; ++i) values[i] = value(rank, step, i);
            send.put(values);
            recv.put(std::vector<unsigned long>(words, sentinel));
            // Exclude producer ordering as a cause before either MPI request is posted.
            if (mode == "view") Kokkos::fence();
            cuda_check(cudaDeviceSynchronize());
            mpi_check(MPI_Barrier(graph), "barrier");
            std::vector<MPI_Request> requests(2 * ranks, MPI_REQUEST_NULL);
            for (int p = 0; p < ranks; ++p)
                mpi_check(MPI_Irecv(recv.data + offset + p * count, count, MPI_UNSIGNED_LONG, p, 0, graph,
                                   &requests[p]), "Irecv");
            for (int p = 0; p < ranks; ++p)
                mpi_check(MPI_Isend(send.data + offset + p * count, count, MPI_UNSIGNED_LONG, p, 0, graph,
                                   &requests[ranks + p]), "Isend");
            std::vector<MPI_Status> statuses(2 * ranks);
            mpi_check(MPI_Waitall(2 * ranks, requests.data(), statuses.data()), "Waitall");
            cuda_check(cudaDeviceSynchronize());
            auto received = recv.get(words);
            int errors = 0, total_errors = 0;
            for (int i = 0; i < words; ++i) {
                const bool active = i >= offset && i < offset + ranks * count;
                const auto expected = active ? value((i - offset) / count, step, offset + rank * count + (i - offset) % count) : sentinel;
                errors += received[i] != expected;
            }
            mpi_check(MPI_Allreduce(&errors, &total_errors, 1, MPI_INT, MPI_SUM, graph), "check values");
            require(total_errors == 0, "received values or untouched padding differ");
        }
    }
    std::cout << "MARS_ALLOCATOR_PASS {" << context.str() << ",\"rank\":" << rank
              << ",\"words\":" << words << ",\"exchanges\":9,\"errors\":0}\n" << std::flush;
    mpi_check(MPI_Comm_free(&graph), "free graph");
    if (mode == "view") Kokkos::finalize();
    mpi_check(MPI_Comm_free(&local), "free local");
    MPI_Finalize();
}
