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
#include "Kokkos_Parallel.hpp"
#include "Kokkos_Parallel_Reduce.hpp"
#include "mars_context.hpp"
#include "mars_globals.hpp"
// #include <bits/c++config.h>
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
#endif // WITH_KOKKOS
#endif

namespace mars {

using DMQ2 = DM<DistributedQuad4Mesh, 2, double, double>;
/*
enum class DMDataDesc
{
    v = 0,
    u = 1
};
 */

template <Integer idx> using DMDataType = typename DMQ2::UserDataType<idx>;

template <typename... T> using tuple = mars::ViewsTuple<T...>;

using dm_tuple = typename DMQ2::user_tuple;

template <Integer Type> void print_dof(const SFC<Type> &dof, const int rank) {
  Kokkos::parallel_for(
      "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
        Integer sfc_elem = dof.get_view_elements()(i);

        double point[3];
        get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(),
                                              dof.get_YDim(), dof.get_ZDim());

        printf("dof: %li - (%lf, %lf) - rank: %i\n", i, point[0], point[1],
               rank);
      });
}

template <Integer Type> void print_dofs(const DMQ2 &dm, const int rank) {
  SFC<Type> dof = dm.get_local_dof_enum();
  Kokkos::parallel_for(
      "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
        const Integer sfc_elem = dm.local_to_sfc(i);
        const Integer global_dof = dm.local_to_global(i);

        double point[3];
        get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(),
                                              dof.get_YDim(), dof.get_ZDim());

        printf("dof: %li - gdof: %li --- (%lf, %lf) - rank: %i\n", i,
               global_dof, point[0], point[1], rank);
      });
}

// print thlocal and the global number of the dof within each element.
// the dof enumeration within eachlement is topological
void print_elem_global_dof(const DMQ2 dm, const int rank) {
  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
    // go through all the dofs of the elem_index element
    for (int i = 0; i < DMQ2::elem_nodes; i++) {
      // get the local dof of the i-th index within thelement
      const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
      // convert the local dof number to global dof number
      Dof d = dm.local_to_global_dof(local_dof);

      // do something. In this case we are printing.
      printf("lgm: local: %li, global: %li, proc: %i, rank:%i\n", local_dof,
             d.get_gid(), d.get_proc(), rank);
    }
  });
}

// print the local and global numbering of the ghost dofs per process
void print_ghost_dofs(const DMQ2 dm, const int rank) {
  Kokkos::parallel_for(
      dm.get_ghost_lg_map().capacity(), KOKKOS_LAMBDA(Integer i) {
        if (dm.get_ghost_lg_map().valid_at(i)) {
          auto sfc = dm.get_ghost_lg_map().key_at(i);
          auto global_dof = dm.get_ghost_lg_map().value_at(i);
          printf("local: %li, global: %li - proc: %li - rank: %li\n",
                 dm.sfc_to_local(sfc), global_dof.get_gid(),
                 global_dof.get_proc(), rank);
        }
      });
}

// form the matrix free operator
template <Integer u, Integer v>
void form_operator(const DMQ2 dm, const int rank) {
  // go through each element and for each element iterate through its dofs
  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
    for (int i = 0; i < DMQ2::elem_nodes; i++) {
      // forach dof get the local number
      const Integer local_dof = dm.get_elem_local_dof(elem_index, i);

      // use the local number to read the corresponding user data
      double ui = dm.get_dof_data<v>(local_dof) + 1;
      // TODO: apply the operator

      // atomically updated the contributions to the same dof
      Kokkos::atomic_add(&dm.get_dof_data<u>(local_dof), ui);
    }
  });
}

template <Integer... dataidx>
void scatter_add_ghost_data(DMQ2 &dm, const context &context) {
  // scatter the data to the procs and keep them in a boundary data tuple
  // again if no template argument is specified all the data is scattered.
  // if not all of them then be careful since the tuple is not initialized on
  // the others example: dm_tuple boundary_data =
  // dm.scatter_ghost_data<1>(context);
  dm_tuple boundary_data = dm.scatter_ghost_data<dataidx...>(context);

  // use the scattered data "boundary_data" to do ops like max, add or min in
  // the dof contributions. Otherwise you can use predifined features like
  // scatter_add as following. careful to use the same template argument for
  // specifing the data as in the scatter_ghost_data since otherwise you might
  // try to access uninitialized tuplelement and get seg faults. example::
  // dm.scatter_add<1>(boundary_data); If: dm.scatter_add<0>(boundary_data) then
  // seg faults.
  dm.scatter_add<dataidx...>(boundary_data);
  /* dm.scatter_max<u>(boundary_data); */
  /* dm.scatter_min<u>(boundary_data); */
}

void poisson(int &argc, char **&argv, const int level) {

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
    int proc_num = mars::rank(context);
#else
    // resources.gpu_id = marsenv::default_gpu();

    // // create a local context
    // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS

    Kokkos::Timer timer;
    // create the quad mesh distributed through the mpi procs.
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
    // the type of the mesh elements. In this case quad4 (Type=4)
    constexpr Integer Type = Elem::ElemType;

    // create the dm object
    DMQ2 dm(&mesh, context);
    // enumerate the dofs locally and globally. The ghost dofs structures
    // are now created and ready to use for the gather and scatter ops.
    dm.enumerate_dofs(context);
    // print locally owned dof numbering
    /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */

    // print local dof numbering
    /* print_dof<Type>(dm.get_local_dof_enum(), proc_num); */
    /* print_dofs<Type>(dm, proc_num); */

    // print the global dofs for each element's local dof
    /* print_elem_global_dof(dm, proc_num); */

    // use as more readable tuple index to identify the data
    constexpr int u = 0;
    constexpr int v = 1;
    // local dof enum size
    const Integer dof_size = dm.get_dof_size();

    // initialize the values by iterating through local dofs
    Kokkos::parallel_for(
        "initdatavalues", dof_size, MARS_LAMBDA(const Integer i) {
          dm.get_dof_data<u>(i) = 0;
          dm.get_dof_data<v>(i) = i;
        });

    // specify the tuple indices of the tuplelements that are needed to gather.
    // if no index specified it gathers all views of the tuple. All data.
    dm.gather_ghost_data<v>(context);

    form_operator<u, v>(dm, proc_num);

    // iterate through the local dofs and print the local number and the data
    /* dm.dof_iterate(
        MARS_LAMBDA(const Integer i) {
            printf("lid: %li, u: %lf, v: %lf, rank: %i\n", i,
                   dm.get_dof_data<u>(i), dm.get_dof_data<v>(i), proc_num);
        }); */

    scatter_add_ghost_data<u>(dm, context);

    dm.dof_iterate(MARS_LAMBDA(const Integer i) {
      printf("ggid: %li, u: %lf, v: %lf, rank: %i\n", i, dm.get_dof_data<u>(i),
             dm.get_dof_data<v>(i), proc_num);
    });

    double time = timer.seconds();
    std::cout << "Matrix free took: " << time << " seconds." << std::endl;

#endif
  } catch (std::exception &e) {
    std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
  }
}
} // namespace mars
