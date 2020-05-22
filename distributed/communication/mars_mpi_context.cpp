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

// Attempting to acquire an MPI context without MPI enabled will produce
// a link error.

#ifndef MARS_HAVE_MPI
#error "build only if MPI is enabled"
#endif

#include <string>
#include <vector>

#include <mpi.h>

#include "mars_distributed_context.hpp"
#include "mars_mpi.hpp"

namespace mars
{

// Throws mpi::mpi_error if MPI calls fail.
struct mpi_context_impl
{
    int size_;
    int rank_;
    MPI_Comm comm_;

    explicit mpi_context_impl(MPI_Comm comm) : comm_(comm)
    {
        size_ = mpi::size(comm_);
        rank_ = mpi::rank(comm_);
    }

    gathered_vector<Integer>
    gather_gids(const std::vector<Integer> &local_gids) const
    {
        return mpi::gather_all_with_partition(local_gids, comm_);
    }

    ViewVectorType<Integer> scatter_gids(const ViewVectorType<Integer> global,
                                              const ViewVectorType<Integer> local) const
    {
        return mpi::scatter(global, local, comm_);
    }

    void scatterv_gids(const ViewVectorType<Integer> global,
                       const ViewVectorType<Integer> local, const std::vector<int> &counts) const
    {
        mpi::scatterv(global, local, counts, comm_);
    }

    template <typename T>
    void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count,
                         const Integer proc_count) const
    {
        mpi::i_send_recv_vec(send_count, receive_count, proc_count, comm_);
    }

    template <typename T>
    void i_send_recv_view(const ViewVectorType<T> &dest, const Integer* dest_displ,
                          const ViewVectorType<T> &src, const Integer* src_displ, const Integer proc_count) const
    {
        mpi::i_send_recv_view(dest, dest_displ, src, src_displ, proc_count, comm_);
    }

    void broadcast(const ViewVectorType<Integer> global) const
    {
        mpi::broadcast(global, comm_);
    }

    std::string name() const { return "MPI"; }
    int id() const { return rank_; }
    int size() const { return size_; }

    template <typename T>
    T min(T value) const
    {
        return mpi::reduce(value, MPI_MIN, comm_);
    }

    template <typename T>
    T max(T value) const
    {
        return mpi::reduce(value, MPI_MAX, comm_);
    }

    template <typename T>
    T sum(T value) const
    {
        return mpi::reduce(value, MPI_SUM, comm_);
    }

    template <typename T>
    std::vector<T> gather(T value, int root) const
    {
        return mpi::gather(value, root, comm_);
    }

    void barrier() const
    {
        mpi::barrier(comm_);
    }
};

template <>
std::shared_ptr<distributed_context> make_mpi_context(MPI_Comm comm)
{
    return std::make_shared<distributed_context>(mpi_context_impl(comm));
}

} // namespace mars
