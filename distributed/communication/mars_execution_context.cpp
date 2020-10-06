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

#include <iostream>
#include <memory>

#include "mars_context.hpp"

//#include "mars_gpu_context.hpp"
#include "mars_distributed_context.hpp"
#include "mars_execution_context.hpp"
//#include "mars_threading.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace mars {

void execution_context_deleter::operator()(execution_context* p) const {
  delete p;
}

execution_context::execution_context(const proc_allocation& resources)
    : distributed(
          make_local_context()) /* ,
thread_pool(std::make_shared<threading::task_system>(resources.num_threads)),
gpu(resources.has_gpu()? std::make_shared<gpu_context>(resources.gpu_id)
                      : std::make_shared<gpu_context>()) */
{}

context make_context(const proc_allocation& p) {
  return context(new execution_context(p));
}

#ifdef WITH_MPI
template <>
execution_context::execution_context(const proc_allocation& resources,
                                     MPI_Comm comm)
    : distributed(
          make_mpi_context(comm)) /* ,
thread_pool(std::make_shared<threading::task_system>(resources.num_threads)),
gpu(resources.has_gpu()? std::make_shared<gpu_context>(resources.gpu_id)
                      : std::make_shared<gpu_context>()) */
{}

template <>
context make_context<MPI_Comm>(const proc_allocation& p, MPI_Comm comm) {
  return context(new execution_context(p, comm));
}
#endif
template <>
execution_context::execution_context(const proc_allocation& resources,
                                     dry_run_info d)
    : distributed(
          make_dry_run_context(d.num_ranks, d.num_cells_per_rank)) /* ,
thread_pool(std::make_shared<threading::task_system>(resources.num_threads)),
gpu(resources.has_gpu()? std::make_shared<gpu_context>(resources.gpu_id)
                      : std::make_shared<gpu_context>()) */
{}

template <>
context make_context(const proc_allocation& p, dry_run_info d) {
  return context(new execution_context(p, d));
}

std::string distribution_type(const context& ctx) {
  return ctx->distributed->name();
}

/* bool has_gpu(const context& ctx) {
    return ctx->gpu->has_gpu();
}

unsigned num_threads(const context& ctx) {
    return ctx->thread_pool->get_num_threads();
} */

unsigned num_ranks(const context& ctx) { return ctx->distributed->size(); }

unsigned rank(const context& ctx) { return ctx->distributed->id(); }

bool has_mpi(const context& ctx) { return ctx->distributed->name() == "MPI"; }

}  // namespace mars
