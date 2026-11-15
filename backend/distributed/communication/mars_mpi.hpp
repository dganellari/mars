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
#include <assert.h>
#include <mpi.h>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mars_mpi_error.hpp>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "mars_gathered_vector.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "mars_utils_kokkos.hpp"
#endif

namespace mars {
    namespace mpi {

        // prototypes
        int rank(MPI_Comm);
        int size(MPI_Comm);
        void barrier(MPI_Comm);

#define MPI_OR_THROW(fn, ...) \
    while (int r_ = fn(__VA_ARGS__)) throw mpi_error(r_, #fn)

        // Type traits for automatically setting MPI_Datatype information for C++ types.
        template <typename T>
        struct mpi_traits {
            constexpr static size_t count() { return sizeof(T); }
            constexpr static MPI_Datatype mpi_type() { return MPI_CHAR; }
            constexpr static bool is_mpi_native_type() { return false; }
        };

#define MAKE_TRAITS(T, M)                                            \
    template <>                                                      \
    struct mpi_traits<T> {                                           \
        constexpr static size_t count() { return 1; }                \
        /* constexpr */ static MPI_Datatype mpi_type() { return M; } \
        constexpr static bool is_mpi_native_type() { return true; }  \
    };

        MAKE_TRAITS(float, MPI_FLOAT)
        MAKE_TRAITS(double, MPI_DOUBLE)
        MAKE_TRAITS(char, MPI_CHAR)
        MAKE_TRAITS(short, MPI_SHORT)
        MAKE_TRAITS(unsigned short, MPI_UNSIGNED_SHORT)
        MAKE_TRAITS(int, MPI_INT)
        MAKE_TRAITS(unsigned, MPI_UNSIGNED)
        MAKE_TRAITS(long, MPI_LONG)
        MAKE_TRAITS(unsigned long, MPI_UNSIGNED_LONG)
        MAKE_TRAITS(long long, MPI_LONG_LONG)
        MAKE_TRAITS(unsigned long long, MPI_UNSIGNED_LONG_LONG)

        static_assert(std::is_same<std::size_t, unsigned long>::value ||
                          std::is_same<std::size_t, unsigned long long>::value,
                      "size_t is not the same as unsigned long or unsigned long long");

        [[noreturn]] inline void abort_invalid_count(MPI_Comm comm, const char *message) {
            std::cerr << "MARS MPI count validation: " << message << std::endl;
            MPI_Abort(comm, MPI_ERR_COUNT);
            std::abort();
        }

        inline void require_collective_count(bool valid, MPI_Comm comm, const char *message) {
            int all_valid = valid ? 1 : 0;
            MPI_OR_THROW(MPI_Allreduce, MPI_IN_PLACE, &all_valid, 1, MPI_INT, MPI_MIN, comm);
            if (!all_valid) abort_invalid_count(comm, message);
        }

        template <typename Count>
        constexpr bool fits_mpi_count(Count count) {
            static_assert(std::is_integral_v<Count>, "MPI counts must be integral");
            return count >= 0 && static_cast<std::uintmax_t>(count) <= std::numeric_limits<int>::max();
        }

        inline int SizeToInt(unsigned int count) {
            if (!fits_mpi_count(count)) throw std::overflow_error("value exceeds the nonnegative MPI int range");
            return static_cast<int>(count);
        }

        template <typename Count>
        int checked_mpi_count(Count count, MPI_Comm comm) {
            // Per-peer calls must not introduce collectives with different call counts.
            if (!fits_mpi_count(count)) abort_invalid_count(comm, "count exceeds the nonnegative MPI int range");
            return static_cast<int>(count);
        }

        inline int checked_mpi_rank(Integer peer, MPI_Comm comm) {
            // Negative MPI sentinel ranks remain valid for the low-level wrappers.
            if (peer < std::numeric_limits<int>::min() || peer > std::numeric_limits<int>::max())
                abort_invalid_count(comm, "peer rank exceeds the MPI int range");
            return static_cast<int>(peer);
        }

        inline int collective_mpi_count(std::size_t count, std::size_t units, MPI_Comm comm) {
            const auto maximum = static_cast<std::size_t>(std::numeric_limits<int>::max());
            require_collective_count(units > 0 && count <= maximum / units, comm,
                                     "count exceeds the nonnegative MPI int range");
            return static_cast<int>(count * units);
        }

        inline std::vector<int> checked_displacements(const std::vector<int> &counts, MPI_Comm comm) {
            std::vector<int> displacements(counts.size() + 1, 0);
            int total = 0;
            bool valid = true;
            for (std::size_t i = 0; i < counts.size(); ++i) {
                if (counts[i] < 0 || counts[i] > std::numeric_limits<int>::max() - total) {
                    valid = false;
                    break;
                }
                total += counts[i];
                displacements[i + 1] = total;
            }
            require_collective_count(valid, comm, "negative count or MPI displacement overflow");
            return displacements;
        }

        inline bool valid_view_displacements(const Integer *displacements, int ranks,
                                             std::size_t extent, bool repeat_source = false) {
            if (!displacements || displacements[0] < 0) return false;
            for (int i = 0; i < ranks; ++i) {
                if (displacements[i + 1] < displacements[i]) return false;
                const Integer count = displacements[i + 1] - displacements[i];
                if (!fits_mpi_count(count)) return false;
                const auto end = repeat_source ? count : displacements[i + 1];
                if (static_cast<std::uintmax_t>(end) > extent) return false;
            }
            return true;
        }

        inline void require_uniform_count(int count, MPI_Comm comm) {
            int root_count = count;
            MPI_OR_THROW(MPI_Bcast, &root_count, 1, MPI_INT, 0, comm);
            require_collective_count(count == root_count, comm, "collective requires equal counts on every rank");
        }

        template <typename T>
        bool valid_peer_counts(const std::vector<T> &counts, int ranks) {
            return counts.size() == static_cast<std::size_t>(ranks) &&
                   std::all_of(counts.begin(), counts.end(), [](T count) {
                       return count >= 0 && static_cast<long double>(count) <= std::numeric_limits<int>::max() &&
                              count == static_cast<int>(count);
                   });
        }
        // Gather individual values of type T from each rank into a std::vector on
        // the root rank.
        // T must be trivially copyable.
        template <typename T>
        std::vector<T> gather(T value, int root, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(std::is_trivially_copyable_v<T>, "gather requires a trivially copyable type");
            auto buffer_size = (rank(comm) == root) ? size(comm) : 0;
            std::vector<T> buffer(buffer_size);

            MPI_OR_THROW(MPI_Gather,
                         &value,
                         traits::count(),
                         traits::mpi_type(),  // send buffer
                         buffer.data(),
                         traits::count(),
                         traits::mpi_type(),  // receive buffer
                         root,
                         comm);

            return buffer;
        }

        // Gather individual values of type T from each rank into a std::vector on
        // the every rank.
        // T must be trivially copyable
        template <typename T>
        std::vector<T> gather_all(T value, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(std::is_trivially_copyable_v<T>, "gather requires a trivially copyable type");
            std::vector<T> buffer(size(comm));

            MPI_OR_THROW(MPI_Allgather,
                         &value,
                         traits::count(),
                         traits::mpi_type(),  // send buffer
                         buffer.data(),
                         traits::count(),
                         traits::mpi_type(),  // receive buffer
                         comm);

            return buffer;
        }

#ifdef MARS_ENABLE_KOKKOS
        // Scatter sfc values of type T from root into each rank
        template <typename T>
        ViewVectorType<T> scatter(const ViewVectorType<T> global, const ViewVectorType<T> local, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "view scatter requires an MPI native type");
            const int chunk_size = collective_mpi_count(local.extent(0), 1, comm);
            require_uniform_count(chunk_size, comm);
            require_collective_count(rank(comm) != 0 || chunk_size == 0 ||
                                     global.extent(0) / chunk_size >= static_cast<std::size_t>(size(comm)),
                                     comm, "scatter source is smaller than the receive partitions");
            Kokkos::fence("structured MPI scatter producers");
#ifdef MARS_NO_RDMA
            auto host_global = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), global);
            auto host_local = Kokkos::create_mirror_view(local);
            auto send_data = host_global.data();
            auto receive_data = host_local.data();
#else
            auto send_data = global.data();
            auto receive_data = local.data();
#endif
            MPI_OR_THROW(MPI_Scatter,
                         send_data,
                         chunk_size,
                         traits::mpi_type(),  // send buffer
                         receive_data,
                         chunk_size,
                         traits::mpi_type(),  // receive buffer
                         0,
                         comm);
#ifdef MARS_NO_RDMA
            Kokkos::deep_copy(local, host_local);
#endif
            return local;
        }

        template <typename T>
        void scatterv(const ViewVectorType<T> global,
                      const ViewVectorType<T> local,
                      const std::vector<int> &counts,
                      MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "view scatterv requires an MPI native type");
            const int chunk_size = collective_mpi_count(local.extent(0), 1, comm);
            require_collective_count(counts.size() == static_cast<std::size_t>(size(comm)), comm,
                                     "scatterv counts must contain one entry per rank");
            auto displs = checked_displacements(counts, comm);
            require_collective_count(counts[rank(comm)] == chunk_size &&
                                     (rank(comm) != 0 || static_cast<std::size_t>(displs.back()) <= global.extent(0)),
                                     comm, "scatterv counts do not match the source or receive view");

            Kokkos::fence("structured MPI scatter producers");
#ifdef MARS_NO_RDMA
            auto host_global = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), global);
            auto host_local = Kokkos::create_mirror_view(local);
            auto send_data = host_global.data();
            auto receive_data = host_local.data();
#else
            auto send_data = global.data();
            auto receive_data = local.data();
#endif
            MPI_OR_THROW(MPI_Scatterv,
                         send_data,
                         counts.data(),
                         displs.data(),
                         traits::mpi_type(),  // send buffer
                         receive_data,
                         chunk_size,
                         traits::mpi_type(),  // receive buffer
                         0,
                         comm);
#ifdef MARS_NO_RDMA
            Kokkos::deep_copy(local, host_local);
#endif
        }

        template <class View>
        void i_send(const View local,
                    const Integer offset,
                    const Integer count,
                    const Integer remote_proc,
                    int tag,
                    MPI_Comm comm,
                    MPI_Request* request) {
            using traits = mpi_traits<typename View::non_const_value_type>;
            static_assert(traits::is_mpi_native_type(), "view send requires an MPI native type");
            const int mpi_count = checked_mpi_count(count, comm);
            if (offset < 0 || static_cast<std::uintmax_t>(offset) > local.size() ||
                static_cast<std::size_t>(mpi_count) > local.size() - static_cast<std::size_t>(offset))
                abort_invalid_count(comm, "send range exceeds its view");
            MPI_OR_THROW(MPI_Isend, local.data() + offset, mpi_count, traits::mpi_type(),
                         checked_mpi_rank(remote_proc, comm), tag, comm, request);
        }

        template <class View>
        void i_receive(const View local,
                       const Integer offset,
                       const Integer count,
                       const Integer source_proc,
                       int tag,
                       MPI_Comm comm,
                       MPI_Request* request) {
            using traits = mpi_traits<typename View::non_const_value_type>;
            static_assert(traits::is_mpi_native_type(), "view receive requires an MPI native type");
            const int mpi_count = checked_mpi_count(count, comm);
            if (offset < 0 || static_cast<std::uintmax_t>(offset) > local.size() ||
                static_cast<std::size_t>(mpi_count) > local.size() - static_cast<std::size_t>(offset))
                abort_invalid_count(comm, "receive range exceeds its view");
            MPI_OR_THROW(MPI_Irecv, local.data() + offset, mpi_count, traits::mpi_type(),
                         checked_mpi_rank(source_proc, comm), tag, comm, request);
        }

        template <typename T>
        void i_send_v(const std::vector<T> &local,
                      const Integer offset,
                      const Integer count,
                      const Integer remote_proc,
                      int tag,
                      MPI_Comm comm,
                      MPI_Request *request) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "vector send requires an MPI native type");
            const int mpi_count = checked_mpi_count(count, comm);
            if (offset < 0 || static_cast<std::uintmax_t>(offset) > local.size() ||
                static_cast<std::size_t>(mpi_count) > local.size() - static_cast<std::size_t>(offset))
                abort_invalid_count(comm, "send range exceeds its vector");

            MPI_OR_THROW(MPI_Isend,
                         local.data() + offset,
                         mpi_count,
                         traits::mpi_type(),  // send buffer
                         checked_mpi_rank(remote_proc, comm),
                         tag,
                         comm,
                         request);
        }

        template <typename T>
        void i_receive_v(std::vector<T> &local,
                         const Integer offset,
                         const Integer count,
                         const Integer source_proc,
                         int tag,
                         MPI_Comm comm,
                         MPI_Request *request) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "vector receive requires an MPI native type");
            const int mpi_count = checked_mpi_count(count, comm);
            if (offset < 0 || static_cast<std::uintmax_t>(offset) > local.size() ||
                static_cast<std::size_t>(mpi_count) > local.size() - static_cast<std::size_t>(offset))
                abort_invalid_count(comm, "receive range exceeds its vector");

            MPI_OR_THROW(MPI_Irecv,
                         local.data() + offset,
                         mpi_count,
                         traits::mpi_type(),  // send buffer
                         checked_mpi_rank(source_proc, comm),
                         tag,
                         comm,
                         request);
        }

        inline void wait_all_send_recv(const int count, std::vector<MPI_Request> &array_req) {
            MPI_OR_THROW(MPI_Waitall, count, array_req.data(), MPI_STATUSES_IGNORE);
        }

        template <typename T>
        void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count, MPI_Comm comm) {
            if (!valid_peer_counts(send_count, size(comm)) || !valid_peer_counts(receive_count, size(comm)))
                abort_invalid_count(comm, "peer count vectors require one valid MPI count per rank");
            Integer proc_count = std::count_if(send_count.begin(), send_count.end(), [](T count) { return count > 0; });

            Integer proc_count_r =
                std::count_if(receive_count.begin(), receive_count.end(), [](T count) { return count > 0; });

            std::vector<MPI_Request> send_req(proc_count);
            std::vector<MPI_Request> receive_req(proc_count_r);

            int recv_proc = 0;
            for (int i = 0; i < receive_count.size(); ++i) {
                if (receive_count[i] > 0) {
                    i_receive_v(receive_count, i, 1, i, 0, comm, &receive_req[recv_proc++]);
                }
            }

            int send_proc = 0;
            for (int i = 0; i < send_count.size(); ++i) {
                if (send_count[i] > 0) {
                    i_send_v(send_count, i, 1, i, 0, comm, &send_req[send_proc++]);
                }
            }

            if (proc_count > 0) {
                wait_all_send_recv(proc_count, send_req);
            }
            if (proc_count_r > 0) {
                wait_all_send_recv(proc_count_r, receive_req);
            }
        }

        template <typename T>
        void i_send_recv_all_to_all(const std::vector<T> &send_count, std::vector<T> &receive_count, MPI_Comm comm) {
            if (!valid_peer_counts(send_count, size(comm)) || !valid_peer_counts(receive_count, size(comm)))
                abort_invalid_count(comm, "peer count vectors require one valid MPI count per rank");
            Integer proc_count = send_count.size();
            Integer proc_count_r = receive_count.size();

            std::vector<MPI_Request> send_req(proc_count);
            std::vector<MPI_Request> receive_req(proc_count_r);

            int recv_proc = 0;
            for (int i = 0; i < receive_count.size(); ++i) {
                i_receive_v(receive_count, i, 1, i, 0, comm, &receive_req[recv_proc++]);
            }

            int send_proc = 0;
            for (int i = 0; i < send_count.size(); ++i) {
                i_send_v(send_count, i, 1, i, 0, comm, &send_req[send_proc++]);
            }

            if (proc_count > 0) {
                wait_all_send_recv(proc_count, send_req);
            }
            if (proc_count_r > 0) {
                wait_all_send_recv(proc_count_r, receive_req);
            }
        }

        template <typename T>
        void i_send_recv_view_to_all(const ViewVectorType<T> &dest,
                                     const Integer *dest_displ,
                                     const ViewVectorType<T> &src,
                                     const Integer *src_displ,
                                     MPI_Comm comm) {
            auto nranks = size(comm);
            if (!valid_view_displacements(dest_displ, nranks, dest.extent(0)) ||
                !valid_view_displacements(src_displ, nranks, src.extent(0), true))
                abort_invalid_count(comm, "peer counts or displacements exceed the MPI range or view bounds");
            Kokkos::fence("structured MPI producers");
#ifdef MARS_NO_RDMA
            auto host_src = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), src);
            auto host_dest = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dest);
#endif

            int proc_count = 0;
            int proc_count_r = 0;
            for (int i = 0; i < nranks; ++i) {
                Integer count_r = dest_displ[i + 1] - dest_displ[i];
                Integer count = src_displ[i + 1] - src_displ[i];
                if (count > 0) {
                    ++proc_count;
                }
                if (count_r > 0) {
                    ++proc_count_r;
                }
            }

            std::vector<MPI_Request> send_req(proc_count);
            std::vector<MPI_Request> receive_req(proc_count_r);

            int recv_proc = 0;
            for (int i = 0; i < nranks; ++i) {
                Integer count = dest_displ[i + 1] - dest_displ[i];
                if (count > 0) {
#ifdef MARS_NO_RDMA
                    i_receive(host_dest, dest_displ[i], count, i, 0, comm, &receive_req[recv_proc++]);
#else
                    i_receive(dest, dest_displ[i], count, i, 0, comm, &receive_req[recv_proc++]);
#endif
                }
            }

            int send_proc = 0;
            for (int i = 0; i < nranks; ++i) {
                Integer count = src_displ[i + 1] - src_displ[i];
                if (count > 0) {
#ifdef MARS_NO_RDMA
                    i_send(host_src, 0, count, i, 0, comm, &send_req[send_proc++]);
#else
                    i_send(src, 0, count, i, 0, comm, &send_req[send_proc++]);
#endif
                }
            }

            if (proc_count > 0) {
                wait_all_send_recv(proc_count, send_req);
            }
            if (proc_count_r > 0) {
                wait_all_send_recv(proc_count_r, receive_req);
            }
#ifdef MARS_NO_RDMA
            Kokkos::deep_copy(dest, host_dest);
#endif
        }

        template <typename T>
        void i_send_recv_view(const ViewVectorType<T> &dest,
                              const Integer *dest_displ,
                              const ViewVectorType<T> &src,
                              const Integer *src_displ,
                              MPI_Comm comm) {
            auto nranks = size(comm);
            if (!valid_view_displacements(dest_displ, nranks, dest.extent(0)) ||
                !valid_view_displacements(src_displ, nranks, src.extent(0)))
                abort_invalid_count(comm, "peer counts or displacements exceed the MPI range or view bounds");
            Kokkos::fence("structured MPI producers");
#ifdef MARS_NO_RDMA
            auto host_src = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), src);
            auto host_dest = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dest);
#endif

            int proc_count = 0;
            int proc_count_r = 0;
            for (int i = 0; i < nranks; ++i) {
                Integer count_r = dest_displ[i + 1] - dest_displ[i];
                Integer count = src_displ[i + 1] - src_displ[i];
                if (count > 0) {
                    ++proc_count;
                }
                if (count_r > 0) {
                    ++proc_count_r;
                }
            }

            std::vector<MPI_Request> send_req(proc_count);
            std::vector<MPI_Request> receive_req(proc_count_r);

            int recv_proc = 0;
            for (int i = 0; i < nranks; ++i) {
                Integer count = dest_displ[i + 1] - dest_displ[i];

                if (count > 0) {
#ifdef MARS_NO_RDMA
                    i_receive(host_dest, dest_displ[i], count, i, 0, comm, &receive_req[recv_proc++]);
#else
                    i_receive(dest, dest_displ[i], count, i, 0, comm, &receive_req[recv_proc++]);
#endif
                }
            }

            int send_proc = 0;
            for (int i = 0; i < nranks; ++i) {
                Integer count = src_displ[i + 1] - src_displ[i];
                if (count > 0) {
#ifdef MARS_NO_RDMA
                    i_send(host_src, src_displ[i], count, i, 0, comm, &send_req[send_proc++]);
#else
                    i_send(src, src_displ[i], count, i, 0, comm, &send_req[send_proc++]);
#endif
                }
            }

            if (proc_count > 0) {
                wait_all_send_recv(proc_count, send_req);
            }
            if (proc_count_r > 0) {
                wait_all_send_recv(proc_count_r, receive_req);
            }
#ifdef MARS_NO_RDMA
            Kokkos::deep_copy(dest, host_dest);
#endif
        }

        template <typename T>
        void broadcast(const ViewVectorType<T> global, MPI_Comm comm) {
            static_assert(std::is_trivially_copyable<T>::value,
                          "broadcast can only be performed on trivally copyable types");

            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "view broadcast requires an MPI native type");
            const int count = collective_mpi_count(global.extent(0), 1, comm);
            require_uniform_count(count, comm);

            Kokkos::fence("structured MPI broadcast producers");
#ifdef MARS_NO_RDMA
            auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), global);
            MPI_OR_THROW(MPI_Bcast, host.data(), count, traits::mpi_type(), 0, comm);
            Kokkos::deep_copy(global, host);
#else
            MPI_OR_THROW(MPI_Bcast, global.data(), count, traits::mpi_type(), 0, comm);
#endif
        }

        /*
        Gather individual values of type T from each rank
        into a view on  the every rank.
        T must be trivially copyable */
        template <typename T>
        void gather_all_view(T value, const ViewVectorType<T> &buffer, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(std::is_trivially_copyable_v<T>, "gather requires a trivially copyable type");
            require_collective_count(buffer.extent(0) >= static_cast<std::size_t>(size(comm)), comm,
                                     "allgather destination is smaller than the rank count");

            Kokkos::fence("structured MPI gather producers");
#ifdef MARS_NO_RDMA
            auto host_buffer = Kokkos::create_mirror_view(buffer);
            auto receive_data = host_buffer.data();
#else
            auto receive_data = buffer.data();
#endif
            MPI_OR_THROW(MPI_Allgather,
                         &value,
                         traits::count(),
                         traits::mpi_type(),  // send buffer
                         receive_data,
                         traits::count(),
                         traits::mpi_type(),  // receive buffer
                         comm);
#ifdef MARS_NO_RDMA
            Kokkos::deep_copy(buffer, host_buffer);
#endif
        }

        template <typename T>
        void gather_all_view(const ViewVectorType<T> &values, ViewVectorType<T> &buffer, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(std::is_trivially_copyable_v<T>, "gather requires a trivially copyable type");
            auto counts = gather_all(collective_mpi_count(values.extent(0), traits::count(), comm), comm);
            auto displs = checked_displacements(counts, comm);

            buffer = ViewVectorType<T>("Gather View buffer", displs.back() / traits::count());
            Kokkos::fence("structured MPI gather producers");
#ifdef MARS_NO_RDMA
            auto host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
            auto host_buffer = Kokkos::create_mirror_view(buffer);
            auto send_data = host_values.data();
            auto receive_data = host_buffer.data();
#else
            auto send_data = values.data();
            auto receive_data = buffer.data();
#endif
            MPI_OR_THROW(MPI_Allgatherv,
                         // const_cast required for MPI implementations that don't use const* in their interfaces
                         const_cast<T *>(send_data),
                         counts[rank(comm)],
                         traits::mpi_type(),  // send buffer
                         receive_data,
                         counts.data(),
                         displs.data(),
                         traits::mpi_type(),  // receive buffer
                         comm);
#ifdef MARS_NO_RDMA
            Kokkos::deep_copy(buffer, host_buffer);
#endif
        }

#endif

        // Specialize gather for std::string.
        inline std::vector<std::string> gather(std::string str, int root, MPI_Comm comm) {
            using traits = mpi_traits<char>;

            auto counts = gather_all(collective_mpi_count(str.size(), 1, comm), comm);
            auto displs = checked_displacements(counts, comm);

            std::vector<char> buffer(displs.back());

            // const_cast required for MPI implementations that don't use const* in
            // their interfaces.
            std::string::value_type *ptr = const_cast<std::string::value_type *>(str.data());
            MPI_OR_THROW(MPI_Gatherv,
                         ptr,
                         counts[rank(comm)],
                         traits::mpi_type(),  // send
                         buffer.data(),
                         counts.data(),
                         displs.data(),
                         traits::mpi_type(),  // receive
                         root,
                         comm);

            // Unpack the raw string data into a vector of strings.
            std::vector<std::string> result;
            auto nranks = size(comm);
            result.reserve(nranks);
            for (auto i = 0; i < nranks; ++i) {
                result.push_back(std::string(buffer.data() + displs[i], counts[i]));
            }
            return result;
        }

        template <typename T>
        std::vector<T> gather_all(const std::vector<T> &values, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(std::is_trivially_copyable_v<T>, "gather requires a trivially copyable type");
            auto counts = gather_all(collective_mpi_count(values.size(), traits::count(), comm), comm);
            auto displs = checked_displacements(counts, comm);

            std::vector<T> buffer(displs.back() / traits::count());
            MPI_OR_THROW(MPI_Allgatherv,
                         // const_cast required for MPI implementations that don't use const* in their interfaces
                         const_cast<T *>(values.data()),
                         counts[rank(comm)],
                         traits::mpi_type(),  // send buffer
                         buffer.data(),
                         counts.data(),
                         displs.data(),
                         traits::mpi_type(),  // receive buffer
                         comm);

            return buffer;
        }

        /// Gather all of a distributed vector
        /// Retains the meta data (i.e. vector partition)
        template <typename T>
        gathered_vector<T> gather_all_with_partition(const std::vector<T> &values, MPI_Comm comm) {
            using gathered_type = gathered_vector<T>;
            using count_type = typename gathered_vector<T>::count_type;
            using traits = mpi_traits<T>;
            static_assert(std::is_trivially_copyable_v<T>, "gather requires a trivially copyable type");

            // We have to use int for the count and displs vectors instead
            // of count_type because these are used as arguments to MPI_Allgatherv
            // which expects int arguments.
            auto counts = gather_all(collective_mpi_count(values.size(), traits::count(), comm), comm);
            auto displs = checked_displacements(counts, comm);

            std::vector<T> buffer(displs.back() / traits::count());

            MPI_OR_THROW(MPI_Allgatherv,
                         // const_cast required for MPI implementations that don't use const* in their interfaces
                         const_cast<T *>(values.data()),
                         counts[rank(comm)],
                         traits::mpi_type(),  // send buffer
                         buffer.data(),
                         counts.data(),
                         displs.data(),
                         traits::mpi_type(),  // receive buffer
                         comm);

            for (auto &d : displs) {
                d /= traits::count();
            }

            return gathered_type(std::move(buffer), std::vector<count_type>(displs.begin(), displs.end()));
        }

#ifdef MARS_ENABLE_KOKKOS
        template <typename T>
        ViewObject<T> reduce(const ViewObject<T> value, MPI_Op op, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "view reduction requires an MPI native type");

            ViewObject<T> result("result");
            Kokkos::fence("structured MPI reduction producers");
#ifdef MARS_NO_RDMA
            auto host_value = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), value);
            auto host_result = Kokkos::create_mirror_view(result);
            MPI_OR_THROW(MPI_Allreduce, host_value.data(), host_result.data(), 1, traits::mpi_type(), op, comm);
            Kokkos::deep_copy(result, host_result);
#else
            MPI_OR_THROW(MPI_Allreduce, value.data(), result.data(), 1, traits::mpi_type(), op, comm);
#endif
            return result;
        }
#endif

        template <typename T>
        T reduce(T value, MPI_Op op, int root, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "can only perform reductions on MPI native types");

            T result;

            MPI_OR_THROW(MPI_Reduce, &value, &result, 1, traits::mpi_type(), op, root, comm);

            return result;
        }

        template <typename T>
        T reduce(T value, MPI_Op op, MPI_Comm comm) {
            using traits = mpi_traits<T>;
            static_assert(traits::is_mpi_native_type(), "can only perform reductions on MPI native types");

            T result;

            MPI_OR_THROW(MPI_Allreduce, &value, &result, 1, traits::mpi_type(), op, comm);

            return result;
        }

        template <typename T>
        std::pair<T, T> minmax(T value) {
            return {reduce<T>(value, MPI_MIN), reduce<T>(value, MPI_MAX)};
        }

        template <typename T>
        std::pair<T, T> minmax(T value, int root) {
            return {reduce<T>(value, MPI_MIN, root), reduce<T>(value, MPI_MAX, root)};
        }

        template <typename T>
        T broadcast(T value, int root, MPI_Comm comm) {
            static_assert(std::is_trivially_copyable<T>::value,
                          "broadcast can only be performed on trivally copyable types");

            using traits = mpi_traits<T>;

            MPI_OR_THROW(MPI_Bcast, &value, traits::count(), traits::mpi_type(), root, comm);

            return value;
        }

        template <typename T>
        T broadcast(int root, MPI_Comm comm) {
            static_assert(std::is_trivially_copyable<T>::value,
                          "broadcast can only be performed on trivally copyable types");

            using traits = mpi_traits<T>;
            T value;

            MPI_OR_THROW(MPI_Bcast, &value, traits::count(), traits::mpi_type(), root, comm);

            return value;
        }

    }  // namespace mpi
}  // namespace mars
