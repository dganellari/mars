/* Copyright (c) 2016, Eidgenössische Technische Hochschule Zürich and
Forschungszentrum Jülich GmbH.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#pragma once
#include <algorithm>
#include <assert.h>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mars_mpi_error.hpp>
#include <mpi.h>
#include <type_traits>
#include <vector>

#include "mars_gathered_vector.hpp"

namespace mars
{
namespace mpi
{

// prototypes
int rank(MPI_Comm);
int size(MPI_Comm);
void barrier(MPI_Comm);

#define MPI_OR_THROW(fn, ...)        \
    while (int r_ = fn(__VA_ARGS__)) \
    throw mpi_error(r_, #fn)

// Type traits for automatically setting MPI_Datatype information for C++ types.
template <typename T>
struct mpi_traits
{
    constexpr static size_t count()
    {
        return sizeof(T);
    }
    constexpr static MPI_Datatype mpi_type()
    {
        return MPI_CHAR;
    }
    constexpr static bool is_mpi_native_type()
    {
        return false;
    }
};

#define MAKE_TRAITS(T, M)                                            \
    template <>                                                      \
    struct mpi_traits<T>                                             \
    {                                                                \
        constexpr static size_t count() { return 1; }                \
        /* constexpr */ static MPI_Datatype mpi_type() { return M; } \
        constexpr static bool is_mpi_native_type() { return true; }  \
    };

MAKE_TRAITS(float, MPI_FLOAT)
MAKE_TRAITS(double, MPI_DOUBLE)
MAKE_TRAITS(char, MPI_CHAR)
MAKE_TRAITS(int, MPI_INT)
MAKE_TRAITS(unsigned, MPI_UNSIGNED)
MAKE_TRAITS(long, MPI_LONG)
MAKE_TRAITS(unsigned long, MPI_UNSIGNED_LONG)
MAKE_TRAITS(long long, MPI_LONG_LONG)
MAKE_TRAITS(unsigned long long, MPI_UNSIGNED_LONG_LONG)

static_assert(std::is_same<std::size_t, unsigned long>::value ||
                  std::is_same<std::size_t, unsigned long long>::value,
              "size_t is not the same as unsigned long or unsigned long long");

inline int SizeToInt(unsigned int u)
{
    if (u > static_cast<unsigned int>((std::numeric_limits<int>::max)()))
    {
        /*  throw std::overflow_error(
            "unsinged int value cannot be stored in a variable of type int."); */
        std::cout << "limit: " << std::numeric_limits<unsigned char>::max() << "limit 2: " << std::numeric_limits<unsigned short>::max() << std::endl;
    }

    return static_cast<int>(u);
}
// Gather individual values of type T from each rank into a std::vector on
// the root rank.
// T must be trivially copyable.
template <typename T>
std::vector<T> gather(T value, int root, MPI_Comm comm)
{
    using traits = mpi_traits<T>;
    auto buffer_size = (rank(comm) == root) ? size(comm) : 0;
    std::vector<T> buffer(buffer_size);

    MPI_OR_THROW(MPI_Gather,
                 &value, traits::count(), traits::mpi_type(),        // send buffer
                 buffer.data(), traits::count(), traits::mpi_type(), // receive buffer
                 root, comm);

    return buffer;
}

// Scatter sfc values of type T from root into each rank
template <typename T>
ViewVectorType<T> scatter(const ViewVectorType<T> global,
                          const ViewVectorType<T> local, MPI_Comm comm)
{
    using traits = mpi_traits<T>;
    Integer chunk_size = local.extent(0);

    MPI_OR_THROW(MPI_Scatter,
                 global.data(), chunk_size, traits::mpi_type(), // send buffer
                 local.data(), chunk_size, traits::mpi_type(),  // receive buffer
                 0, comm);

    return local;
}

template <typename T>
void scatterv(const ViewVectorType<T> global,
              const ViewVectorType<T> local, const std::vector<int> &counts, MPI_Comm comm)
{
    using traits = mpi_traits<T>;
    Integer chunk_size = local.extent(0);

    auto displs = make_scan_index(counts);

    MPI_OR_THROW(MPI_Scatterv,
                 global.data(), counts.data(), displs.data(), traits::mpi_type(), // send buffer
                 local.data(), chunk_size, traits::mpi_type(),                    // receive buffer
                 0, comm);
}

template <typename T>
void i_send(const ViewVectorType<T> local, const Integer offset, const Integer count,
            const Integer remote_proc, int tag, MPI_Comm comm, MPI_Request *request)
{
    using traits = mpi_traits<T>;

    MPI_OR_THROW(MPI_Isend,
                 local.data() + offset, count, traits::mpi_type(), // send buffer
                 remote_proc, tag, comm, request);
}

template <typename T>
void i_receive(const ViewVectorType<T> local, const Integer offset, const Integer count,
               const Integer source_proc, int tag, MPI_Comm comm, MPI_Request *request)
{
    using traits = mpi_traits<T>;

    MPI_OR_THROW(MPI_Irecv,
                 local.data() + offset, count, traits::mpi_type(), // send buffer
                 source_proc, tag, comm, request);
}

template <typename T>
void i_send_v(const std::vector<T> &local, const Integer offset, const Integer count,
              const Integer remote_proc, int tag, MPI_Comm comm, MPI_Request *request)
{
    using traits = mpi_traits<T>;

    MPI_OR_THROW(MPI_Isend,
                 local.data() + offset, count, traits::mpi_type(), // send buffer
                 remote_proc, tag, comm, request);
}

template <typename T>
void i_receive_v(std::vector<T> &local, const Integer offset, const Integer count,
                 const Integer source_proc, int tag, MPI_Comm comm, MPI_Request *request)
{
    using traits = mpi_traits<T>;

    MPI_OR_THROW(MPI_Irecv,
                 local.data() + offset, count, traits::mpi_type(), // send buffer
                 source_proc, tag, comm, request);
}

inline void wait_all_send_recv(const int count, std::vector<MPI_Request> &array_req)
{
    MPI_OR_THROW(MPI_Waitall,
                 count, array_req.data(), MPI_STATUSES_IGNORE);
}

template <typename T>
void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count,
                     MPI_Comm comm)
{
    Integer proc_count = std::count_if(send_count.begin(), send_count.end(),
                                       [](T count) { return count > 0; });

    Integer proc_count_r = std::count_if(receive_count.begin(), receive_count.end(),
                                       [](T count) { return count > 0; });

    printf("proc count proc: %li - %li - %li - %li\n", proc_count, proc_count_r, receive_count.size(), send_count.size());
    std::vector<MPI_Request> send_req(proc_count);
    std::vector<MPI_Request> receive_req(proc_count_r);

    int recv_proc = 0;
    for (int i = 0; i < receive_count.size(); ++i)
    {
        if (receive_count[i] > 0)
        {
            i_receive_v(receive_count, i, 1, i, 0, comm, &receive_req[recv_proc++]);
        }
    }

    int send_proc = 0;
    for (int i = 0; i < send_count.size(); ++i)
    {
        if (send_count[i] > 0)
        {
            i_send_v(send_count, i, 1, i, 0, comm, &send_req[send_proc++]);
        }
    }

    if (proc_count > 0)
    {
        wait_all_send_recv(proc_count, send_req);
    }
    if(proc_count_r > 0)
    {
        wait_all_send_recv(proc_count_r, receive_req);
    }
}

template <typename T>
void i_send_recv_view(const ViewVectorType<T> &dest, const Integer* dest_displ,
                      const ViewVectorType<T> &src, const Integer* src_displ,
                      MPI_Comm comm)
{
    auto nranks = size(comm);

    int proc_count= 0;
    int proc_count_r = 0;
    for (int i = 0; i < nranks; ++i)
    {
        Integer count_r = dest_displ[i + 1] - dest_displ[i];
        Integer count = src_displ[i + 1] - src_displ[i];
        if (count > 0)
        {
            ++proc_count;
        }
        if(count_r > 0)
        {
            ++proc_count_r;
        }
    }

    std::vector<MPI_Request> send_req(proc_count);
    std::vector<MPI_Request> receive_req(proc_count_r);


    int recv_proc = 0;
    for (int i = 0; i < nranks; ++i)
    {
        Integer count = dest_displ[i + 1] - dest_displ[i];
        if (count > 0)
        {
            i_receive(dest, dest_displ[i], count, i, 0, comm, &receive_req[recv_proc++]);
        }
    }

    int send_proc = 0;
    for (int i = 0; i < nranks; ++i)
    {
        Integer count = src_displ[i + 1] - src_displ[i];
        if (count > 0)
        {
            i_send(src, src_displ[i], count, i, 0, comm, &send_req[send_proc++]);
        }
    }

    if (proc_count > 0)
    {
        wait_all_send_recv(proc_count, send_req);
    }
    if(proc_count_r > 0)
    {
        wait_all_send_recv(proc_count_r, receive_req);
    }
}

template <typename T>
void broadcast(const ViewVectorType<T> global, MPI_Comm comm)
{
    static_assert(std::is_trivially_copyable<T>::value,
                  "broadcast can only be performed on trivally copyable types");

    using traits = mpi_traits<T>;
    Integer count = global.extent(0);

    MPI_OR_THROW(MPI_Bcast,
                 global.data(), count, traits::mpi_type(), 0, comm);
}

/* 
Gather individual values of type T from each rank
into a view on  the every rank.
T must be trivially copyable */
template <typename T>
void gather_all_view(T value, const ViewVectorType<T>& buffer, MPI_Comm comm)
{
    using traits = mpi_traits<T>;

    MPI_OR_THROW(MPI_Allgather,
                 &value, traits::count(), traits::mpi_type(),        // send buffer
                 buffer.data(), traits::count(), traits::mpi_type(), // receive buffer
                 comm);
}

// Gather individual values of type T from each rank into a std::vector on
// the every rank.
// T must be trivially copyable
template <typename T>
std::vector<T> gather_all(T value, MPI_Comm comm)
{
    using traits = mpi_traits<T>;
    std::vector<T> buffer(size(comm));

    MPI_OR_THROW(MPI_Allgather,
                 &value, traits::count(), traits::mpi_type(),        // send buffer
                 buffer.data(), traits::count(), traits::mpi_type(), // receive buffer
                 comm);

    return buffer;
}

// Specialize gather for std::string.
inline std::vector<std::string> gather(std::string str, int root, MPI_Comm comm)
{
    using traits = mpi_traits<char>;

    auto counts = gather_all(int(str.size()), comm);
    auto displs = make_scan_index(counts);

    std::vector<char> buffer(displs.back());

    // const_cast required for MPI implementations that don't use const* in
    // their interfaces.
    std::string::value_type *ptr = const_cast<std::string::value_type *>(str.data());
    MPI_OR_THROW(MPI_Gatherv,
                 ptr, counts[rank(comm)], traits::mpi_type(),                     // send
                 buffer.data(), counts.data(), displs.data(), traits::mpi_type(), // receive
                 root, comm);

    // Unpack the raw string data into a vector of strings.
    std::vector<std::string> result;
    auto nranks = size(comm);
    result.reserve(nranks);
    for (auto i = 0; i < nranks; ++i)
    {
        result.push_back(std::string(buffer.data() + displs[i], counts[i]));
    }
    return result;
}

template <typename T>
std::vector<T> gather_all(const std::vector<T> &values, MPI_Comm comm)
{

    using traits = mpi_traits<T>;
    auto counts = gather_all(int(values.size()), comm);
    for (auto &c : counts)
    {
        c *= traits::count();
    }
    auto displs = make_scan_index(counts);

    std::vector<T> buffer(displs.back() / traits::count());
    MPI_OR_THROW(MPI_Allgatherv,
                 // const_cast required for MPI implementations that don't use const* in their interfaces
                 const_cast<T *>(values.data()), counts[rank(comm)], traits::mpi_type(), // send buffer
                 buffer.data(), counts.data(), displs.data(), traits::mpi_type(),        // receive buffer
                 comm);

    return buffer;
}

/// Gather all of a distributed vector
/// Retains the meta data (i.e. vector partition)
template <typename T>
gathered_vector<T> gather_all_with_partition(const std::vector<T> &values, MPI_Comm comm)
{
    using gathered_type = gathered_vector<T>;
    using count_type = typename gathered_vector<T>::count_type;
    using traits = mpi_traits<T>;

    // We have to use int for the count and displs vectors instead
    // of count_type because these are used as arguments to MPI_Allgatherv
    // which expects int arguments.
    auto counts = gather_all(int(values.size()), comm);
    for (auto &c : counts)
    {
        c *= traits::count();
    }
    auto displs = make_scan_index(counts);

    std::vector<T> buffer(displs.back() / traits::count());

    MPI_OR_THROW(MPI_Allgatherv,
                 // const_cast required for MPI implementations that don't use const* in their interfaces
                 const_cast<T *>(values.data()), counts[rank(comm)], traits::mpi_type(), // send buffer
                 buffer.data(), counts.data(), displs.data(), traits::mpi_type(),        // receive buffer
                 comm);

    for (auto &d : displs)
    {
        d /= traits::count();
    }

    return gathered_type(
        std::move(buffer),
        std::vector<count_type>(displs.begin(), displs.end()));
}

template <typename T>
ViewObject<T> reduce(const ViewObject<T> value, MPI_Op op, MPI_Comm comm)
{
    using traits = mpi_traits<T>;

    ViewObject<T> result("result");
    MPI_Allreduce(value.data(), result.data(), 1, traits::mpi_type(), op, comm);

    return result;
}

template <typename T>
T reduce(T value, MPI_Op op, int root, MPI_Comm comm)
{
    using traits = mpi_traits<T>;
    static_assert(traits::is_mpi_native_type(),
                  "can only perform reductions on MPI native types");

    T result;

    MPI_OR_THROW(MPI_Reduce,
                 &value, &result, 1, traits::mpi_type(), op, root, comm);

    return result;
}

template <typename T>
T reduce(T value, MPI_Op op, MPI_Comm comm)
{
    using traits = mpi_traits<T>;
    static_assert(traits::is_mpi_native_type(),
                  "can only perform reductions on MPI native types");

    T result;

    MPI_Allreduce(&value, &result, 1, traits::mpi_type(), op, comm);

    return result;
}

template <typename T>
std::pair<T, T> minmax(T value)
{
    return {reduce<T>(value, MPI_MIN), reduce<T>(value, MPI_MAX)};
}

template <typename T>
std::pair<T, T> minmax(T value, int root)
{
    return {reduce<T>(value, MPI_MIN, root), reduce<T>(value, MPI_MAX, root)};
}

template <typename T>
T broadcast(T value, int root, MPI_Comm comm)
{
    static_assert(std::is_trivially_copyable<T>::value,
                  "broadcast can only be performed on trivally copyable types");

    using traits = mpi_traits<T>;

    MPI_OR_THROW(MPI_Bcast,
                 &value, traits::count(), traits::mpi_type(), root, comm);

    return value;
}

template <typename T>
T broadcast(int root, MPI_Comm comm)
{
    static_assert(std::is_trivially_copyable<T>::value,
                  "broadcast can only be performed on trivally copyable types");

    using traits = mpi_traits<T>;
    T value;

    MPI_OR_THROW(MPI_Bcast,
                 &value, traits::count(), traits::mpi_type(), root, comm);

    return value;
}

} // namespace mpi
} // namespace mars
