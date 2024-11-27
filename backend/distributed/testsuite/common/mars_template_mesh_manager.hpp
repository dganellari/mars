
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

#pragma once

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

#ifdef MARS_ENABLE_KOKKOS

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

#ifdef MARS_ENABLE_KOKKOS
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

#ifdef MARS_ENABLE_CUDA
    template <Integer Type, class KeyType = MortonKey<Unsigned>>
    bool test_mars_distributed_nonsimplex_mesh_ghost_layer(DistributedMesh<KeyType, Type> &mesh_thrust_boundary,
                                                           DistributedMesh<KeyType, Type> &mesh_predicate_boundary) {
        using namespace mars;

        partition_mesh(mesh_thrust_boundary);
        mesh_thrust_boundary.template build_boundary_element_sets<Type, 0>();

        partition_mesh(mesh_predicate_boundary);
        mesh_predicate_boundary.template build_boundary_element_sets<Type, 1>();

        bool equal_scan_boundary = thrust::equal(
            thrust::device,
            mesh_thrust_boundary.get_view_scan_boundary().data(),
            mesh_thrust_boundary.get_view_scan_boundary().data() + mesh_thrust_boundary.get_view_scan_boundary().size(),
            mesh_predicate_boundary.get_view_scan_boundary().data());

        bool equal_boundary = thrust::equal(
            thrust::device,
            mesh_thrust_boundary.get_view_boundary().data(),
            mesh_thrust_boundary.get_view_boundary().data() + mesh_thrust_boundary.get_view_boundary().size(),
            mesh_predicate_boundary.get_view_boundary().data());

        bool equal_boundary_lsfc = thrust::equal(thrust::device,
                                                 mesh_thrust_boundary.get_view_boundary_sfc_index().data(),
                                                 mesh_thrust_boundary.get_view_boundary_sfc_index().data() +
                                                     mesh_thrust_boundary.get_view_boundary_sfc_index().size(),
                                                 mesh_predicate_boundary.get_view_boundary_sfc_index().data());

        return equal_scan_boundary && equal_boundary && equal_boundary_lsfc;
    }

    template <class KeyType = MortonKey<Unsigned>>
    bool test_mars_distributed_nonsimplex_mesh_ghost_layer_2D(const int x, const int y) {
        using namespace mars;
        try {
            mars::proc_allocation resources;

#ifdef MARS_ENABLE_MPI
            // create a distributed context
            auto context = mars::make_context(resources, MPI_COMM_WORLD);
#endif

#ifdef MARS_ENABLE_KOKKOS
            DistributedQuad4Mesh<KeyType> mesh1(context);
            generate_distributed_cube(mesh1, x, y, 0);

            DistributedQuad4Mesh<KeyType> mesh2(context);
            generate_distributed_cube(mesh2, x, y, 0);

            return test_mars_distributed_nonsimplex_mesh_ghost_layer(mesh1, mesh2);
#endif
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
            return 0;
        }
    }

    template <class KeyType = MortonKey<Unsigned>>
    bool test_mars_distributed_nonsimplex_mesh_ghost_layer_3D(const int x, const int y, const int z) {
        using namespace mars;
        try {
            mars::proc_allocation resources;

#ifdef MARS_ENABLE_MPI
            // create a distributed context
            auto context = mars::make_context(resources, MPI_COMM_WORLD);
#endif

#ifdef MARS_ENABLE_KOKKOS
            DistributedHex8Mesh<KeyType> mesh1(context);
            generate_distributed_cube(mesh1, x, y, z);

            DistributedHex8Mesh<KeyType> mesh2(context);
            generate_distributed_cube(mesh2, x, y, z);

            return test_mars_distributed_nonsimplex_mesh_ghost_layer(mesh1, mesh2);
#endif
        } catch (std::exception &e) {
            std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
            return 0;
        }
    }
#endif

    template <class DofHandler>
    MARS_INLINE_FUNCTION bool check_unique_dofs(DofHandler &handler) {
        auto size = handler.get_dof_size();
        ViewVectorType<bool> unique_ids("unique_ids", 1);
        Kokkos::deep_copy(unique_ids, true);

        // Check for duplicates
        Kokkos::parallel_for(
            "check_duplicates", size - 1, KOKKOS_LAMBDA(const size_t i) {
                auto sfc_i = handler.template local_to_sfc<DofHandler::Block>(i);
                auto sfc_j = handler.template local_to_sfc<DofHandler::Block>(i + 1);

                // Get the components of the block local dof in order to compute the local dof
                auto component_i = handler.get_component(i);
                auto component_j = handler.get_component(i + 1);

                // in this way we test component and sfc to local map to get to the same local dof value
                auto local_i = handler.sfc_to_local(sfc_i, component_i);
                auto local_j = handler.sfc_to_local(sfc_j, component_j);
                assert(local_i == i && local_j == i + 1);

                // it is unique only if one of the two is not equal
                if (sfc_i == sfc_j && local_i == local_j) {
                    bool expected = true;
                    Kokkos::atomic_compare_exchange(&unique_ids(0), expected, false);
                }
            });

        // Copy the result back to the host
        auto host_unique_ids = Kokkos::create_mirror_view(unique_ids);
        Kokkos::deep_copy(host_unique_ids, unique_ids);
        bool result = host_unique_ids(0);

        return result;
    }

    template <Integer Type, Integer Degree = 1, bool Overlap = true, class KeyType = MortonKey<Unsigned>>
    bool test_mars_distributed_dof_handler(const int xDim, const int yDim, const int zDim, const int block) {
        using namespace mars;
        mars::proc_allocation resources;

        bool result = false;

#ifdef MARS_ENABLE_MPI
        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();
        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef MARS_ENABLE_KOKKOS

        Kokkos::Timer timer;

        using DMesh = DistributedMesh<KeyType, Type>;

        // create the quad mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        constexpr Integer Block = 0;
        using DOFHandler = DofHandler<DMesh, Degree, Block>;
        double time_gen = timer.seconds();
        std::cout << "Mesh Generation took: " << time_gen << std::endl;

        Kokkos::Timer timer_dof;

        DOFHandler dof_handler(mesh);
        dof_handler.enumerate_dofs();
        dof_handler.set_block(block);

        double time_dh = timer_dof.seconds();
        std::cout << "DOFHandler enum took: " << time_dh << std::endl;

        result = check_unique_dofs(dof_handler);
#endif
        return result;
    }

    // Helper function to compute the expected value for a given element and DOF index
    KOKKOS_INLINE_FUNCTION
    int expected_value(int elem_index, int dof_index, int block) {
        // Compute the expected value based on the element index, DOF index, and block size
        return elem_index * block + dof_index;
    }

    template <Integer Type = ElementType::Quad4,
              Integer Degree = 1,
              bool Overlap = true,
              class KeyType = MortonKey<Unsigned>>
    bool test_mars_distributed_dof_map(const int xDim, const int yDim, const int zDim, const int block) {
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

#ifdef MARS_ENABLE_KOKKOS

        using DMesh = DistributedMesh<KeyType, Type>;

        // create the quad mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        constexpr Integer Block = 0;
        using DOFHandler = DofHandler<DMesh, Degree, Block>;
        using DMQ = DM<DOFHandler, double, double, double>;
        using FE = FEDofMap<DOFHandler>;

        DOFHandler dof_handler(mesh);
        dof_handler.enumerate_dofs();
        dof_handler.set_block(block);

        auto fe_dof_map = build_fe_dof_map<DOFHandler, Overlap>(dof_handler);
        auto elem_dof_enum = fe_dof_map.get_elem_dof_map();

        bool success = true;

        Kokkos::parallel_reduce(
            "VerifyTopologicalOrder2DLarge",
            elem_dof_enum.extent(0),
            KOKKOS_LAMBDA(const int elem_index, bool &local_success) {
                for (int dof_index = 0; dof_index < elem_dof_enum.extent(1); ++dof_index) {
                    // Check that the DOFs are in the expected topological order
                    if (elem_dof_enum(elem_index, dof_index) != expected_value(elem_index, dof_index, block)) {
                        local_success = false;
                    }
                }
            },
            Kokkos::LAnd<bool>(success));

        return success;
#endif
        return false;
    }

/* #ifdef MARS_ENABLE_KOKKOS_KERNELS

    template <Integer Type, Integer Degree = 1, bool Overlap = true, class KeyType = MortonKey<Unsigned>>
    bool test_mars_distributed_sparsity_pattern(const int xDim, const int yDim, const int zDim, const int block) {
        // Create a DofHandler
        mars::proc_allocation resources;

#ifdef MARS_ENABLE_MPI
        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
#else
        // resources.gpu_id = marsenv::default_gpu();
        // // create a local context
        // auto context = mars::make_context(resources);
#endif

        using DMesh = DistributedMesh<KeyType, Type>;

        // create the mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        constexpr Integer Block = 0;
        using DOFHandler = DofHandler<DMesh, Degree, Block>;
        DOFHandler dof_handler(mesh);
        dof_handler.enumerate_dofs();
        dof_handler.set_block(block);

        auto fe_dof_map = build_fe_dof_map<DOFHandler, Overlap>(dof_handler);
        // Initialize SparsityPattern
        using SPattern = SparsityPattern<double, Integer, unsigned long, DOFHandler>;
        SPattern sparsity_pattern(dof_handler);

        // Build the sparsity pattern
        sparsity_pattern.build_pattern(fe_dof_map);

        // Verify the sparsity pattern
        auto row_map = sparsity_pattern.get_row_map();
        auto entries = sparsity_pattern.get_col();

        bool success = true;

        // Verify specific values for a known scenario
        Kokkos::parallel_reduce(
            "VerifySpecificValues",
            row_map.extent(0) - 1,
            KOKKOS_LAMBDA(const int row, bool &local_success) {
                // Example check: Ensure that the diagonal element is present
                bool has_diagonal = false;
                for (int i = row_map(row); i < row_map(row + 1); ++i) {
                    if (entries(i) == row) {
                        has_diagonal = true;
                        break;
                    }
                }
                if (!has_diagonal) {
                    local_success = false;
                }
            },
            Kokkos::LAnd<bool>(success));

        SparsityMatrix<SPattern> sparsity_matrix(sparsity_pattern);
        sparsity_matrix.build_crs_matrix();

        auto crs_matrix = sparsity_matrix.get_crs_matrix();

        // Verify the CRS matrix
        Kokkos::parallel_reduce(
            "VerifyCRSMatrix",
            crs_matrix.numRows(),
            KOKKOS_LAMBDA(const int row, bool &local_success) {
                if (crs_matrix.graph.row_map(row + 1) < crs_matrix.graph.row_map(row)) {
                    local_success = false;
                }
            },
            Kokkos::LAnd<bool>(success));

        return success;
    }

#endif */

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

#ifdef MARS_ENABLE_KOKKOS

        Kokkos::Timer timer;

        using DMesh = DistributedMesh<KeyType, Type>;

        // create the quad mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        constexpr Integer Block = 0;
        using DOFHandler = DofHandler<DMesh, Degree, Block>;
        using DMQ = DM<DOFHandler, double, double, double>;
        using FE = FEDofMap<DOFHandler>;

#ifdef MARS_ENABLE_KOKKOS_KERNELS
        using SPattern = SparsityPattern<double, Integer, unsigned long, DOFHandler>;
        using SMatrix = SparsityMatrix<SPattern>;
#endif

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
#ifdef MARS_ENABLE_KOKKOS_KERNELS

        SPattern sp(dof_handler);
        sp.build_pattern(fe);

        SMatrix sm(sp);
        sm.build_crs_matrix();

        double time_sp = timer_sp.seconds();
        std::cout << "SP build took: " << time_sp << std::endl;
#endif
        double time_all = timer.seconds();
        std::cout << "Discretization: " << time_all << std::endl;

        Kokkos::Timer timer_as;
        ViewVectorType<Integer> x_local("local_data", dof_handler.get_dof_size());
#ifdef MARS_ENABLE_KOKKOS_KERNELS
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) {
            x_local(i) = -1;
            if (dof_handler.is_owned(i)) sm.set_value(i, i, 9);
        });

        double time_aall = timer_as.seconds();
        std::cout << "Assembly: " << time_aall << std::endl;
#endif
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

        dof_handler.boundary_dof_iterate(MARS_LAMBDA(const Integer local_dof) { x_local(local_dof) = 7; }, 0);

        dof_handler.boundary_dof_iterate(MARS_LAMBDA(const Integer local_dof) { x_local(local_dof) = 8; }, 1);

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
