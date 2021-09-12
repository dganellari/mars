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

#ifdef WITH_MPI

    void test_mpi(int &argc, char **&argv) {
        using namespace mars;

        try {
            // Constructing guard will initialize MPI with a
            // call to MPI_Init_thread()
            marsenv::mpi_guard guard(argc, argv, false);
            printf("test\n");
            // When leaving this scope, the destructor of guard will
            // call MPI_Finalize()
        } catch (std::exception &e) {
            std::cerr << "error: " << e.what() << "\n";
        }
    }

#endif

    void test_mpi_context(int &argc, char **&argv) {
        using namespace mars;

        try {
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
            // int root = mars::rank(context) == 0;
#else
            // resources.gpu_id = marsenv::default_gpu();

            // // create a local context
            // auto context = mars::make_context(resources);
#endif

            // Print a banner with information about hardware configuration
            // std::cout << "gpu:      " << (has_gpu(context) ? "yes" : "no") << "\n";
            // std::cout << "threads:  " << num_threads(context) << "\n";
            std::cout << "mpi:      " << (has_mpi(context) ? "yes" : "no") << "\n";
            std::cout << "ranks:    " << num_ranks(context) << "\n";
            std::cout << "rank:    " << rank(context) << "\n" << std::endl;

            // run some simulations!
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        }
    }

    template <typename... T>
    using user_tuple = mars::ViewsTuple<T...>;

    template <typename... T>
    struct functor {
        user_tuple<T...> tuple;

        functor(user_tuple<T...> t) : tuple(t) {}

        MARS_INLINE_FUNCTION
        void operator()(int i) const {
            /* note the use of std get instead */
            std::get<1>(tuple)(i) = 1;
        }
    };

    void test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(const int x, const int y) {
        using namespace mars;
        try {
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
            /* marsenv::mpi_guard guard(argc, argv, false); */

            // assign a unique gpu to this rank if available
            /*  resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD); */

            // create a distributed context
            auto context = mars::make_context(resources, MPI_COMM_WORLD);
            // bool root = mars::rank(context) == 0;
#else
            // resources.gpu_id = marsenv::default_gpu();

            // // create a local context
            // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS

            Kokkos::Timer timer_gen;
            DistributedQuad4Mesh mesh(context);
            generate_distributed_cube(mesh, x, y, 0);

            double time_gen = timer_gen.seconds();
            std::cout << "2D Mesh Generation took: " << time_gen << std::endl;

#endif
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        }
    }

    void test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(const int x, const int y, const int z) {
        using namespace mars;

        try {
            mars::proc_allocation resources;

            /* try to detect how many threads can be run on this system */
            /* resources.num_threads = marsenv::thread_concurrency(); */

            /* override thread count if the user set MARS_NUM_THREADS
            if (auto nt = marsenv::get_env_num_threads())
            {
                resources.num_threads = nt;
            }*/

#ifdef WITH_MPI
            /* initialize MPI */
            /* marsenv::mpi_guard guard(argc, argv, false); */

            /* assign a unique gpu to this rank if available */
            /* resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD);*/

            /* create a distributed context */
            auto context = mars::make_context(resources, MPI_COMM_WORLD);
            /* bool root = mars::rank(context) == 0; */
#else
            /* resources.gpu_id = marsenv::default_gpu();

            // create a local context
            auto context = mars::make_context(resources); */
#endif

#ifdef WITH_KOKKOS
            using namespace Kokkos;

            DistributedHex8Mesh mesh(context);
            Kokkos::Timer timer_gen;
            generate_distributed_cube(mesh, x, y, z);

            double time_gen = timer_gen.seconds();
            std::cout << "3D Mesh Generation took: " << time_gen << std::endl;

#endif
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        }
    }

    template <Integer Type = ElementType::Quad4>
    void test_mars_distributed_vector_valued(const int xDim, const int yDim, const int zDim) {
        using namespace mars;
        mars::proc_allocation resources;

#ifdef WITH_MPI
        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();
        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS

        Kokkos::Timer timer;

        using DMesh = DistributedMesh<Type>;

        // create the quad mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        /* constexpr Integer Dim = DMesh::Dim; */
        /* using Elem = typename DistributedQuad4Mesh::Elem; */
        /* constexpr Integer Type = Elem::ElemType; */

        constexpr Integer Degree = 2;
        constexpr Integer Block = DMesh::Dim;
        /* using DOFHandler = DofHandler<DMesh, Degree>; */
        using DOFHandler = DofHandler<DMesh, Degree, Block>;
        using DMQ = DM<DOFHandler, double, double, double>;
        using FE = FEDofMap<DOFHandler>;
        using SPattern = SparsityPattern<double, Integer, unsigned long, DOFHandler>;

        // use as more readable tuple index to identify the data
        static constexpr int INPUT = 0;
        static constexpr int OUTPUT = 1;
        static constexpr int RHSD = 2;

        DOFHandler dof_handler(&mesh);
        dof_handler.enumerate_dofs();

        /* dof_handler.print_dofs(proc_num); */

        auto fe = build_fe_dof_map(dof_handler);

        // print the global dofs for each element's local dof
        /* print_fe_dof_map(dof_handler, fe); */
        /* print_ghost_dofs(dm); */

        SPattern sp(dof_handler);
        sp.build_pattern(fe);
        /* sp.print_sparsity_pattern(); */

        // create the dm object
        /* DMQ dm(dof_handler); */
        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        //
        // print locally owned dof numbering
        /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */

        // print local dof numbering
        /* print_dof<Type>(dm.get_local_dof_enum(), proc_num); */
        /* print_dofs<Type>(dm, proc_num); */

        ViewVectorType<Integer> x_local("local_data", dof_handler.get_dof_size());
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) { x_local(i) = -1; });

        auto owned_size = dof_handler.get_owned_dof_size();
        ViewVectorType<Integer> owned("owned_data", owned_size);
        dof_handler.owned_dof_iterate(MARS_LAMBDA(const Integer i) { owned(i) = proc_num; });

        set_locally_owned_data(dof_handler, x_local, owned);
        gather_ghost_data(dof_handler, x_local);
        /* scatter_add_ghost_data(dof_handler, x_local); */

        /*print using the index iterate*/
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) {
            auto idata = x_local(i);

            Dof d = dof_handler.local_to_global_dof(i);
            auto base_global = dof_handler.compute_base(d.get_gid());

            printf("lid: %li, u: %li,  global: %li, base: %li, rank: %i\n",
                   i,
                   idata,
                   d.get_gid(),
                   base_global,
                   d.get_proc());
        });

#endif
    }
}  // namespace mars
