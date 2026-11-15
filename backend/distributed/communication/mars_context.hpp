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

#include <memory>

namespace mars {

    // Requested dry-run parameters.
    struct dry_run_info {
        unsigned num_ranks;
        unsigned num_cells_per_rank;
        dry_run_info(unsigned ranks, unsigned cells_per_rank) : num_ranks(ranks), num_cells_per_rank(cells_per_rank) {}
    };

    // Legacy resource inputs retained for source compatibility. The current
    // context ignores these fields; configure threads and devices through Kokkos.

    struct proc_allocation {
        unsigned num_threads;

        // Historical device ordinal; -1 means no device was requested.
        int gpu_id;

        proc_allocation() : proc_allocation(1, -1) {}

        proc_allocation(unsigned threads, int gpu) : num_threads(threads), gpu_id(gpu) {}

        bool has_gpu() const { return gpu_id >= 0; }
    };

    // mars::execution_context holds communication state. Kokkos owns the thread
    // and device execution resources independently of this context.

    // Forward declare execution_context.
    struct execution_context;

    // mars::context is an opaque handle for the execution context for use
    // in the public API, implemented as a unique pointer.
    //
    // As execution_context is an incomplete type, an explicit deleter must be
    // provided.
    struct execution_context_deleter {
        void operator()(execution_context *) const;
    };
    using context = std::unique_ptr<execution_context, execution_context_deleter>;

    // Helpers for creating contexts. These are implemented in the back end.

    // Non-distributed communication context; resources is currently ignored.
    context make_context(const proc_allocation &resources = proc_allocation{});

    // Distributed context borrowing MPI communicator comm, or a dry-run context
    // using dry_run_info. resources is currently ignored.
    template <typename Comm>
    context make_context(const proc_allocation &resources, Comm comm);

    // Queries for properties of execution resources in a context.

    std::string distribution_type(const context &);
    bool has_gpu(const context &);
    unsigned num_threads(const context &);
    bool has_mpi(const context &);
    unsigned num_ranks(const context &);
    unsigned rank(const context &);

}  // namespace mars
