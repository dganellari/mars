/* ---------------------
 * Matrix free implementation of the Poisson equation
 *
//global to local
vk = P^T v

w_j := quadrature weight
x_j := quadrature point

// Q is the number of quad points
j=1,...,Q
// L is the number of basis functions
l,i=1,...,L

g_glob_lj :=  (J^-T grad \phi_l (x_j) ) \in R^2

// linear number of basis functions
g_j = sum_l g_glob_lj * v_l


dX_j = det(J) * ref_vol * w_j


// linear number of quadrature points
ui = sum_j dot(g_j, (J^-T grad \phi_i (x_j) )) * dX_j

// number of basis functions times number of quad points
L * Q

uk =[ u1, ... uL ], L = 4

---------------------
//local to global
u = P uk */

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

using DMQ2 = DM<DistributedQuad4Mesh, 2, double, double>;
/*
enum class DMDataDesc
{
    v = 0,
    u = 1
};
 */

template <Integer idx>
using DMDataType = typename DMQ2::UserDataType<idx>;

template <typename... T>
using tuple = mars::ViewsTuple<T...>;

using dm_tuple = typename DMQ2::user_tuple;

template<Integer Type>
void print_dof(const SFC<Type> &dof, const int rank)
{
    Kokkos::parallel_for("for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
        Integer sfc_elem = dof.get_view_elements()(i);

        double point[3];
        get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

        printf("dof: %li - (%lf, %lf) - rank: %i\n", i, point[0], point[1], rank);
    });
}

//print thlocal and the global number of the dof within each element.
void print_global_dof_enumeration(const DMQ2 dm, const int rank)
{
    dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
        for (int i = 0; i < DMQ2::elem_nodes; i++)
        {
            const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
            Dof d = dm.local_to_global_dof(local_dof);

            //do something. In this case we are printing.
            printf("lgm: local: %li, global: %li, proc: %i, rank:%i\n",  local_dof, d.get_gid(), d.get_proc(), rank);
        }
    });
}

void form_operator(const DMQ2 dm, const ViewVectorType<double> &u, const
        ViewVectorType<double> &v, const int rank)
{
    dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
        for (int i = 0; i < DMQ2::elem_nodes; i++)
        {
            const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
            /* Dof d = dm.local_to_global_dof(local_dof); */

            double ui = v(local_dof) + 1;
            //TODO: apply the operator

            Kokkos::atomic_add(&u(local_dof), ui);
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

        /* ViewVectorType<double> u("u", dm.get_local_dof_enum().get_elem_size());
        ViewVectorType<double> v("v", dm.get_local_dof_enum().get_elem_size());

 */
        //use as more readable tuple index to identify the data
        constexpr int u = 0;
        constexpr int v = 1;

        //specify the tuple indices of the tuplelements that are needed to gather.
        //if no index specified it gathers all views of the tuple. All data.
        dm.gather_ghost_data<v>(context);

        /* form_operator(dm, u, v, proc_num); */

        dm_tuple boundary_data = dm.scatter_ghost_data<u>(context);
        dm.scatter_add<u>(boundary_data);
        /* dm.scatter_max(boundary_data); */

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
