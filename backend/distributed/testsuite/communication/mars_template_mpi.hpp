/* Copyright (c) 2016, Eidgenössische Technische Hochschule Zürich
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

#include "mars.hpp"

namespace mars {

#ifdef MARS_ENABLE_MPI

bool test_mpi(int &argc, char **&argv) {
    using namespace mars;
    
    try {
        // Check if MPI is already initialized
        int initialized;
        MPI_Initialized(&initialized);
        
        if (!initialized) {
            // Use static to make the guard persist
            static marsenv::mpi_guard guard(argc, argv, false);
        }
    } catch (std::exception &e) {
        std::cerr << "error: " << e.what() << "\n";
        return false;
    }
    return true;
}

context test_mpi_context(int &argc, char **&argv) {
    using namespace mars;
    
    try {
        mars::proc_allocation resources;
        
        // Check if MPI is already initialized
        int initialized;
        MPI_Initialized(&initialized);
        
        if (!initialized) {
            // Use static to make the guard persist
            static marsenv::mpi_guard guard(argc, argv, false);
        }
        
        // Create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        
        if(rank(context) == 0) {
            std::cout << "mpi:      " << (has_mpi(context) ? "yes" : "no") << "\n";
            std::cout << "ranks:    " << num_ranks(context) << "\n";
            std::cout << "rank:    " << rank(context) << "\n" << std::endl;
        }
        
        return context;
    } catch (std::exception &e) {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        return {};
    }
}

#endif
}  // namespace mars
