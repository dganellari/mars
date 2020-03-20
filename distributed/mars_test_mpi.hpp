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

//#include "mars_mesh_kokkos.hpp"

#include "mars_context.hpp"
#include <exception>
#include <iostream>

#ifdef WITH_MPI

#include "mars_mpi_guard.hpp"

#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_generation.hpp"
#endif //WITH_KOKKOS

void test_mpi(int &argc, char **&argv)
{

    using namespace mars;

    try
    {
        // Constructing guard will initialize MPI with a
        // call to MPI_Init_thread()
        marsenv::mpi_guard guard(argc, argv, false);
        printf("test\n");
        // When leaving this scope, the destructor of guard will
        // call MPI_Finalize()
    }
    catch (std::exception &e)
    {
        std::cerr << "error: " << e.what() << "\n";
    }
}

#endif

void test_mpi_context(int &argc, char **&argv)
{

    using namespace mars;

    try
    {
        mars::proc_allocation resources;
        /* 
        // try to detect how many threads can be run on this system
        resources.num_threads = marsenv::thread_concurrency();

        // override thread count if the user set MARS_NUM_THREADS
        if (auto nt = marsenv::get_env_num_threads())
        {
            resources.num_threads = nt;
        } */

#ifdef WITH_MPI
        // initialize MPI
        marsenv::mpi_guard guard(argc, argv, false);

        // assign a unique gpu to this rank if available
        /*  resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD); */

        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int root = mars::rank(context) == 0;
#else
        // resources.gpu_id = marsenv::default_gpu();

        // // create a local context
        // auto context = mars::make_context(resources);
#endif

        // Print a banner with information about hardware configuration
        // std::cout << "gpu:      " << (has_gpu(context) ? "yes" : "no") << "\n";
        //std::cout << "threads:  " << num_threads(context) << "\n";
        std::cout << "mpi:      " << (has_mpi(context) ? "yes" : "no") << "\n";
        std::cout << "ranks:    " << num_ranks(context) << "\n";
        std::cout << "rank:    " << rank(context) << "\n"
                  << std::endl;

        // run some simulations!
    }
    catch (std::exception &e)
    {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
    }
}

void test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(int &argc, char **&argv, const int level)
{

    using namespace mars;

    try
    {
        mars::proc_allocation resources;
        /* 
        // try to detect how many threads can be run on this system
        resources.num_threads = marsenv::thread_concurrency();

        // override thread count if the user set MARS_NUM_THREADS
        if (auto nt = marsenv::get_env_num_threads())
        {
            resources.num_threads = nt;
        } */

#ifdef WITH_MPI
        // initialize MPI
        marsenv::mpi_guard guard(argc, argv, false);

        // assign a unique gpu to this rank if available
        /*  resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD); */

        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        //bool root = mars::rank(context) == 0;
#else
        // resources.gpu_id = marsenv::default_gpu();

        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS
        using namespace Kokkos;

        DistributedQuad4Mesh mesh;
        generate_distributed_cube(context, mesh, level, level, 0);

        /* int proc_num = rank(context);

        ViewMatrixType<Real> poi = mesh.get_view_points();

        parallel_for(
            "print_elem_chunk1", mesh.get_view_points().extent(0), KOKKOS_LAMBDA(const int i) {
                printf(" pt: [(%f, %f) - %i]\n", poi(i, 0), poi(i, 1), proc_num);
            });

        ViewMatrixType<Integer> eeel = mesh.get_view_elements();

        parallel_for(
            "print_elem_chunk", mesh.get_view_elements().extent(0), KOKKOS_LAMBDA(const int i) {
                printf("el 3D: [(  %li, %li, %li, %li)] - %i]\n",
                       eeel(i, 0), eeel(i, 1), eeel(i, 2), eeel(i, 3), proc_num);
            }); */
#endif
    }
    catch (std::exception &e)
    {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
    }
}

void test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(int &argc, char **&argv, const int level)
{

    using namespace mars;

    try
    {
        mars::proc_allocation resources;
        /* 
        // try to detect how many threads can be run on this system
        resources.num_threads = marsenv::thread_concurrency();

        // override thread count if the user set MARS_NUM_THREADS
        if (auto nt = marsenv::get_env_num_threads())
        {
            resources.num_threads = nt;
        } */

#ifdef WITH_MPI
        // initialize MPI
        marsenv::mpi_guard guard(argc, argv, false);

        // assign a unique gpu to this rank if available
        /*  resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD); */

        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        //bool root = mars::rank(context) == 0;
#else
        // resources.gpu_id = marsenv::default_gpu();

        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS
        using namespace Kokkos;

        DistributedHex8Mesh mesh;
        generate_distributed_cube(context, mesh, level, level, level);
/* 
        int proc_num = rank(context);

        ViewMatrixType<Real> poi = mesh.get_view_points();

        parallel_for(
            "print_elem_chunk1", mesh.get_view_points().extent(0), KOKKOS_LAMBDA(const int i) {
                printf(" pt: [(%f, %f, %f) - %i]\n", poi(i, 0), poi(i, 1), poi(i, 2), proc_num);
            });

        ViewMatrixType<Integer> eeel = mesh.get_view_elements();

        parallel_for(
            "print_elem_chunk", mesh.get_view_elements().extent(0), KOKKOS_LAMBDA(const int i) {
                printf("el 3D: [(  %li, %li, %li, %li, %li, %li, %li, %li )] - %i]\n",
                        eeel(i, 0), eeel(i, 1), eeel(i, 2), eeel(i, 3), eeel(i, 4), eeel(i, 5), 
                                    eeel(i, 6), eeel(i, 7), proc_num);
            }); */
#endif
    }
    catch (std::exception &e)
    {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
    }
}