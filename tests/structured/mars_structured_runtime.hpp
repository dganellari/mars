#pragma once

#include <mars.hpp>
#include "mars_structured_mpi_check.hpp"
#include "mars_structured_view_check.hpp"

#include <array>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#endif
#ifdef KOKKOS_ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

namespace mars::test
{

namespace detail {
struct TransportFill {
    ViewVectorType<int> outgoing;
    int sender;
    MARS_INLINE_FUNCTION void operator()(int i) const { outgoing(i) = 100 * sender + i; }
};
struct TransportCheck {
    ViewVectorType<int> incoming;
    int source;
    MARS_INLINE_FUNCTION void operator()(int i, int& failures) const {
        failures += incoming(i) != 100 * source + i;
    }
};
} // namespace detail

// The tests establish device assignment and MPI transport before exercising mesh APIs.
class StructuredRuntime
{
public:
    StructuredRuntime(int& argc, char**& argv) : mpi_(argc, argv, false)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &ranks_);
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank_, MPI_INFO_NULL, &local_);
        MPI_Comm_rank(local_, &local_rank_);
        MPI_Comm_size(local_, &local_ranks_);

        int devices = 0;
        int device_query_ok = 1;
#ifdef KOKKOS_ENABLE_CUDA
        if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::Cuda>)
        {
            backend_ = "CUDA";
            device_query_ok = cudaGetDeviceCount(&devices) == cudaSuccess && devices > 0;
        }
#endif
#ifdef KOKKOS_ENABLE_HIP
        if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::HIP>)
        {
            backend_ = "HIP";
            device_query_ok = hipGetDeviceCount(&devices) == hipSuccess && devices > 0;
        }
#endif
        require(device_query_ok, "GPU backend has no usable visible device");
        const char* expected = std::getenv("MARS_STRUCTURED_EXPECT_BACKEND");
        require(!expected || backend_ == expected, "Kokkos execution space differs from the requested backend");

        bool explicit_device = std::getenv("KOKKOS_DEVICE_ID") || std::getenv("KOKKOS_MAP_DEVICE_ID_BY");
        for (int i = 1; i < argc; ++i)
            explicit_device = explicit_device || std::string(argv[i]).rfind("--kokkos-device-id", 0) == 0 ||
                              std::string(argv[i]).rfind("--kokkos-map-device-id-by", 0) == 0;
        // Respect per-rank visibility and explicit Kokkos choices. MPI's shared
        // communicator also works with launchers whose local-rank variable is unknown to Kokkos.
        bool temporary_device = devices > 0 && !explicit_device;
        bool device_environment_ok = true;
        if (temporary_device)
        {
            std::string choice = std::to_string(devices == 1 ? 0 : local_rank_ % devices);
            device_environment_ok = setenv("KOKKOS_DEVICE_ID", choice.c_str(), 1) == 0;
        }
        require(device_environment_ok, "Cannot set the test GPU device");
        env_ = std::make_unique<mars::Env>(argc, argv, MPI_COMM_WORLD);
        { mars::Env nested(argc, argv, MPI_COMM_WORLD); }
        require(Kokkos::is_initialized(), "Nested MARS environment finalized the caller's Kokkos runtime");
        if (temporary_device) unsetenv("KOKKOS_DEVICE_ID");

        std::array<char, 64> device_id{};
        std::strcpy(device_id.data(), "cpu");
        std::string device_name = "cpu";
        int runtime_version = 0, driver_version = 0;
#ifdef KOKKOS_ENABLE_CUDA
        if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::Cuda>)
        {
            require(cudaGetDevice(&device_) == cudaSuccess, "Cannot read selected CUDA device");
            require(cudaDeviceGetPCIBusId(device_id.data(), int(device_id.size()), device_) == cudaSuccess,
                    "Cannot identify selected CUDA device");
            cudaDeviceProp properties{};
            require(cudaGetDeviceProperties(&properties, device_) == cudaSuccess,
                    "Cannot read selected CUDA device properties");
            device_name = properties.name;
            require(cudaRuntimeGetVersion(&runtime_version) == cudaSuccess &&
                    cudaDriverGetVersion(&driver_version) == cudaSuccess,
                    "Cannot read CUDA runtime and driver versions");
        }
#endif
#ifdef KOKKOS_ENABLE_HIP
        if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::HIP>)
        {
            require(hipGetDevice(&device_) == hipSuccess, "Cannot read selected HIP device");
            require(hipDeviceGetPCIBusId(device_id.data(), int(device_id.size()), device_) == hipSuccess,
                    "Cannot identify selected HIP device");
            hipDeviceProp_t properties{};
            require(hipGetDeviceProperties(&properties, device_) == hipSuccess,
                    "Cannot read selected HIP device properties");
            device_name = properties.name;
            require(hipRuntimeGetVersion(&runtime_version) == hipSuccess &&
                    hipDriverGetVersion(&driver_version) == hipSuccess,
                    "Cannot read HIP runtime and driver versions");
        }
#endif
        std::vector<char> local_devices(device_id.size() * local_ranks_);
        MPI_Allgather(device_id.data(), int(device_id.size()), MPI_CHAR, local_devices.data(),
                      int(device_id.size()), MPI_CHAR, local_);
        bool unique = true;
        if (devices > 0)
            for (int peer = 0; peer < local_ranks_; ++peer)
                if (peer != local_rank_ && std::strcmp(device_id.data(), local_devices.data() + peer * device_id.size()) == 0)
                    unique = false;
        require(unique, "Two ranks on one host selected the same GPU; assign one GPU per rank");

        check_transport();
        char hostname[MPI_MAX_PROCESSOR_NAME]{};
        int hostname_size = 0;
        MPI_Get_processor_name(hostname, &hostname_size);
        char mpi_version[MPI_MAX_LIBRARY_VERSION_STRING]{};
        int mpi_version_size = 0;
        MPI_Get_library_version(mpi_version, &mpi_version_size);
        std::string mpi_library(mpi_version, mpi_version_size);
        for (char& character : mpi_library)
            if (static_cast<unsigned char>(character) < 32) character = ' ';
        std::ostringstream receipt;
        receipt << "MARS_STRUCTURED_RECEIPT {\"backend\":" << std::quoted(backend_)
                << ",\"rank\":" << rank_ << ",\"ranks\":" << ranks_
                << ",\"local_rank\":" << local_rank_ << ",\"local_ranks\":" << local_ranks_
                << ",\"device\":" << device_ << ",\"device_id\":" << std::quoted(std::string(device_id.data()))
                << ",\"device_name\":" << std::quoted(device_name)
                << ",\"runtime_version\":" << runtime_version << ",\"driver_version\":" << driver_version
                << ",\"hostname\":" << std::quoted(std::string(hostname, hostname_size))
                << ",\"mpi_library\":" << std::quoted(mpi_library)
                << ",\"kokkos_version\":" << KOKKOS_VERSION
                << ",\"concurrency\":" << Kokkos::DefaultExecutionSpace().concurrency()
#ifdef MARS_NO_RDMA
                << ",\"transport\":\"staged\""
#else
                << ",\"transport\":\"direct\""
#endif
                << ",\"mpi_device_buffer\":" << (device_mpi_ ? "true" : "false") << "}";
        // Ordered complete lines keep launcher output machine-readable.
        for (int writer = 0; writer < ranks_; ++writer)
        {
            if (writer == rank_) std::cout << receipt.str() << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    ~StructuredRuntime()
    {
        env_.reset();
        MPI_Comm_free(&local_);
    }

private:
    void require(bool valid, const char* message) const
    {
        int value = valid ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (!value)
        {
            if (!valid) std::cerr << "Structured preflight rank " << rank_ << ": " << message << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
            std::abort();
        }
    }

    void check_transport()
    {
        static_assert(std::is_same_v<typename ViewVectorType<int>::memory_space,
                                     typename Kokkos::DefaultExecutionSpace::memory_space>,
                      "Structured views must use the selected Kokkos execution space's memory");
        constexpr int size = 32;
        ViewVectorType<int> outgoing("preflight_send", size), incoming("preflight_recv", size);
        const int source = (rank_ + ranks_ - 1) % ranks_;
        const int destination = (rank_ + 1) % ranks_;
        Kokkos::parallel_for("preflight_fill", size, detail::TransportFill{outgoing, rank_});
        Kokkos::fence();
        int status;
#ifdef MARS_NO_RDMA
        auto host_outgoing = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), outgoing);
        auto host_incoming = Kokkos::create_mirror_view(Kokkos::HostSpace(), incoming);
        status = MPI_Sendrecv(host_outgoing.data(), size, MPI_INT, destination, 719,
                              host_incoming.data(), size, MPI_INT, source, 719, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        Kokkos::deep_copy(incoming, host_incoming);
#else
        status = MPI_Sendrecv(outgoing.data(), size, MPI_INT, destination, 719,
                              incoming.data(), size, MPI_INT, source, 719, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        device_mpi_ = device_ >= 0;
#endif
        require(status == MPI_SUCCESS, "MPI ring exchange failed");
        int errors = 0;
        Kokkos::parallel_reduce("preflight_check", size, detail::TransportCheck{incoming, source}, errors);
        require(errors == 0, "Kokkos kernel read incorrect values after MPI transfer");
        require(check_mpi_view_collectives(MPI_COMM_WORLD) == 0,
                "Structured MPI view collective returned incorrect values");
        require(check_structured_view_lifetimes(MPI_COMM_WORLD) == 0,
                "Structured public view wrapper lost its data or allocation");
    }

    marsenv::mpi_guard mpi_;
    std::unique_ptr<mars::Env> env_;
    MPI_Comm local_ = MPI_COMM_NULL;
    std::string backend_ = Kokkos::DefaultExecutionSpace::name();
    int rank_ = 0, ranks_ = 1, local_rank_ = 0, local_ranks_ = 1, device_ = -1;
    bool device_mpi_ = false;
};

} // namespace mars::test
