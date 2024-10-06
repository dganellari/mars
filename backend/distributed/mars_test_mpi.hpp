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
#include "mars_base.hpp"

namespace mars {

    template <class KeyType = MortonKey<Unsigned>>
    Unsigned test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(const int x, const int y) {
        using namespace mars;
        try {
            mars::proc_allocation resources;

#ifdef MARS_ENABLE_MPI

            // create a distributed context
            auto context = mars::make_context(resources, MPI_COMM_WORLD);
#else
            // resources.gpu_id = marsenv::default_gpu();

            // // create a local context
            // auto context = mars::make_context(resources);
#endif

#ifdef MARS_ENABLE_KOKKOS_KERNELS

            Kokkos::Timer timer_gen;
            DistributedQuad4Mesh<KeyType> mesh(context);
            generate_distributed_cube(mesh, x, y, 0);

            double time_gen = timer_gen.seconds();
            std::cout << "2D Mesh Generation took: " << time_gen << std::endl;

            return get_number_of_elements(mesh);

#endif
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
            return 0;
        }
    }

    template <class KeyType = MortonKey<Unsigned>>
    Unsigned test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(const int x, const int y, const int z) {
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

#ifdef MARS_ENABLE_MPI
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

#ifdef MARS_ENABLE_KOKKOS_KERNELS
            using namespace Kokkos;

            DistributedHex8Mesh<KeyType> mesh(context);
            Kokkos::Timer timer_gen;
            generate_distributed_cube(mesh, x, y, z);

            double time_gen = timer_gen.seconds();
            std::cout << "3D Mesh Generation took: " << time_gen << std::endl;
            return get_number_of_elements(mesh);

#endif
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
            return 0;
        }
    }

    template <Integer Type = ElementType::Quad4,
              Integer Degree = 1,
              bool Overlap = true,
              class KeyType = MortonKey<Unsigned>>
    void test_mars_distributed_vector_valued(const int xDim, const int yDim, const int zDim, const int block) {
        using namespace mars;
        mars::proc_allocation resources;

#ifdef MARS_ENABLE_MPI
        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();
        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef MARS_ENABLE_KOKKOS_KERNELS

        Kokkos::Timer timer;

        using DMesh = DistributedMesh<KeyType, Type>;

        // create the quad mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        constexpr Integer Block = 0;
        using DOFHandler = DofHandler<DMesh, Degree, Block>;
        using DMQ = DM<DOFHandler, double, double, double>;
        using FE = FEDofMap<DOFHandler>;
        using SPattern = SparsityPattern<double, Integer, unsigned long, DOFHandler>;
        using SMatrix = SparsityMatrix<SPattern>;

        // use as more readable tuple index to identify the data
        /* static constexpr int INPUT = 0;
        static constexpr int OUTPUT = 1;
        static constexpr int RHSD = 2; */

        double time_gen = timer.seconds();
        std::cout << "Mesh Generation took: " << time_gen << std::endl;

        Kokkos::Timer timer_dof;

        DOFHandler dof_handler(mesh);
        dof_handler.enumerate_dofs();
        dof_handler.set_block(block);

        double time_dh = timer_dof.seconds();
        std::cout << "DOFHandler enum took: " << time_dh << std::endl;

        Kokkos::Timer timer_map;

        // dof_handler.print_dofs(proc_num);

        auto fe = build_fe_dof_map<DOFHandler, Overlap>(dof_handler);

        double time_map = timer_map.seconds();
        std::cout << "DOFMAP took: " << time_map << std::endl;

        Kokkos::Timer timer_sp;
        // print the global dofs for each element's local dof
        /* fe.print(); */

        SPattern sp(dof_handler);
        sp.build_pattern(fe);

        SMatrix sm(sp);
        sm.build_crs_matrix();

        double time_sp = timer_sp.seconds();
        std::cout << "SP build took: " << time_sp << std::endl;

        double time_all = timer.seconds();
        std::cout << "Discretization: " << time_all << std::endl;

        Kokkos::Timer timer_as;
        ViewVectorType<Integer> x_local("local_data", dof_handler.get_dof_size());
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) {
            x_local(i) = -1;
            if (dof_handler.is_owned(i)) sm.set_value(i, i, 9);
        });

        double time_aall = timer_as.seconds();
        std::cout << "Assembly: " << time_aall << std::endl;

        /* sp.print_sparsity_pattern(); */
        /* sm.write("overlap_sp.txt"); */

        auto owned_size = dof_handler.get_owned_dof_size();
        ViewVectorType<Integer> owned("owned_data", owned_size);
        dof_handler.owned_dof_iterate(MARS_LAMBDA(const Integer i) { owned(i) = proc_num; });

        set_locally_owned_data(dof_handler, x_local, owned);

        Kokkos::Timer timer_c;
        gather_ghost_data(dof_handler, x_local);
        /* scatter_add_ghost_data(dof_handler, x_local); */

        double time_call = timer_c.seconds();
        std::cout << "Collect: " << time_call << std::endl;

        dof_handler.boundary_dof_iterate(
            MARS_LAMBDA(const Integer local_dof) { x_local(local_dof) = 7; }, 0);

        dof_handler.boundary_dof_iterate(
            MARS_LAMBDA(const Integer local_dof) { x_local(local_dof) = 8; }, 1);

        /*print using the index iterate*/
        /* dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) {
            auto idata = x_local(i);

            Dof d = dof_handler.local_to_global_dof(i);
            auto base_global = dof_handler.compute_base(d.get_gid());

            printf("lid: %li, u: %li,  global: %li, base: %li, rank: %i\n",
                   i,
                   idata,
                   d.get_gid(),
                   base_global,
                   d.get_proc());
        }); */

#endif
    }
}  // namespace mars
