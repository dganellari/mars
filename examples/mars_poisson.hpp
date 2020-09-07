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

#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Bitset.hpp"
#include "Kokkos_Parallel_Reduce.hpp"
#include "mars_context.hpp"
#include "mars_globals.hpp"
#include <bits/c++config.h>
#include <exception>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#ifdef WITH_MPI

#include "mars_mpi_guard.hpp"

#ifdef WITH_KOKKOS
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_mesh_generation.hpp"
#endif //WITH_KOKKOS
#endif

namespace mars
{

using DMQ2 = DM<DistributedQuad4Mesh, 2>;
/* enum DataDesc : int
{
    u = 0,
    du_0 = 1,
    du_1 = 2,
    dudt = 3
};
template <Integer idx>
using DataType = typename Data::UserDataType<idx>;

template <typename... T>
using user_tuple = mars::ViewsTuple<T...>;
 */

template<Integer Type>
MARS_INLINE_FUNCTION void print_dof(const SFC<Type> &dof, const int rank)
{
    Kokkos::parallel_for("for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
        Integer sfc_elem = dof.get_view_elements()(i);

        double point[3];
        get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

        printf("dof: %li - (%lf, %lf) - rank: %i\n", i, point[0], point[1], rank);
    });
}

//print thlocal and the global number of the dof within each element.
MARS_INLINE_FUNCTION void print_global_dof_enumeration(const DMQ2 dm, const int rank)
{
    dm.get_data().elem_iterate(MARS_LAMBDA(const Integer elem_index) {
        for (int i = 0; i < DMQ2::elem_nodes; i++)
        {
            const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
            Dof d = dm.local_to_global_dof(local_dof);

            //do something. In this case we are printing.
            printf("lgm: local: %li, global: %li, proc: %i, rank:%i\n",  local_dof, d.get_gid(), d.get_proc(), rank);
        }
    });
}

void poisson(int &argc, char **&argv, const int level)
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
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();

        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS

        DistributedQuad4Mesh mesh;
        generate_distributed_cube(context, mesh, level, level, 0);
        /* mesh.set_periodic(); //set the domain to be periodic */
/*
        const Integer xDim = mesh.get_XDim();
        const Integer yDim = mesh.get_YDim();
        const Integer zDim = mesh.get_ZDim();
 */
        constexpr Integer Dim = DistributedQuad4Mesh::Dim;

        using Elem = typename DistributedQuad4Mesh::Elem;
        constexpr Integer Type = Elem::ElemType;

        std::cout << "Type: " << Type << std::endl;

        DMQ2 dm(&mesh, context);
        dm.enumerate_dofs(context);
        /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */
        print_dof<Type>(dm.get_local_dof_enum(), proc_num);

        print_global_dof_enumeration(dm, proc_num);

        /* ViewVectorType<Integer> sfc = mesh.get_view_sfc(); */
        Kokkos::Timer timer;

        double time = timer.seconds();
        std::cout << "face iterate took: " << time << " seconds." << std::endl;

        /* print_derivatives<Type, DataDesc::du_0, DataDesc::du_1>(data); */
#endif
    }
    catch (std::exception &e)
    {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
    }
}
} // namespace mars
