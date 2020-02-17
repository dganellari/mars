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

#include "mars_context.hpp"

#include "mars_distributed_context.hpp"
//#include "mars_threading.hpp"
//#include "mars_gpu_context.hpp"

namespace mars {

// execution_context is a simple container for the state relating to
// execution resources.
// Specifically, it has handles for the distributed context, gpu
// context and thread pool.
//
// Note: the public API uses an opaque handle mars::context for
// execution_context, to hide implementation details of the
// container and its constituent contexts from the public API.

struct execution_context {
    distributed_context_handle distributed;
/*     task_system_handle thread_pool;
    gpu_context_handle gpu; */

    execution_context(const proc_allocation& resources = proc_allocation{});

    // Use a template for constructing with a specific distributed context.
    // Specialised implementations are implemented in execution_context.cpp.
    template <typename Comm>
    execution_context(const proc_allocation& resources, Comm comm);
};

} // namespace mars
