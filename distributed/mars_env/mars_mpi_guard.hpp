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

#include <exception>

#include <mpi.h>

#include "mars_mpi_error.hpp"

namespace marsenv
{

struct mpi_guard
{
    mpi_guard(int &argc, char **&argv, bool fatal_errors = true)
    {
        init(&argc, &argv, fatal_errors);
    }

    explicit mpi_guard(bool fatal_errors = true)
    {
        init(nullptr, nullptr, fatal_errors);
    }

    ~mpi_guard()
    {
        // Test if the stack is being unwound because of an exception.
        // If other ranks have not thrown an exception, there is a very
        // high likelihood that the MPI_Finalize will hang due to the other
        // ranks calling other MPI calls.
        // We don't directly call MPI_Abort in this case because that would
        // force exit the application before the exception that is unwinding
        // the stack has been caught, which would deny the opportunity to print
        // an error message explaining the cause of the exception.
        if (!std::uncaught_exception())
        {
            MPI_Finalize();
        }
    }

private:
    void init(int *argcp, char ***argvp, bool fatal_errors)
    {
        int provided;
        int ev = MPI_Init_thread(argcp, argvp, MPI_THREAD_SERIALIZED, &provided);
        if (ev)
        {
            throw mars::mpi_error(ev, "MPI_Init_thread");
        }
        else if (provided < MPI_THREAD_SERIALIZED)
        {
            throw mars::mpi_error(MPI_ERR_OTHER, "MPI_Init_thread: MPI_THREAD_SERIALIZED unsupported");
        }

        if (!fatal_errors)
        {
            MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
        }
    }
};

} // namespace marsenv
