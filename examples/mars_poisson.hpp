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
#include "mars_quad4.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_boundary_conditions.hpp"
#include  "mars_poisson_operator.hpp"
#endif // WITH_KOKKOS
#endif


namespace mars {

/* using DMQ2 = DM<DistributedQuad4Mesh, 2, double, double>; */
using DMQ2 = DM<DistributedQuad4Mesh, 1, double, double>;
/*
enum class DMDataDesc
{
    v = 0,
    u = 1
};
 */

// use as more readable tuple index to identify the data
static constexpr int INPUT = 0;
static constexpr int OUTPUT = 1;

template <Integer idx>
using DMDataType = typename DMQ2::UserDataType<idx>;

template <typename... T>
using tuple = mars::ViewsTuple<T...>;

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
void print_elem_global_dof(const DMQ2 dm) {
  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
    // go through all the dofs of the elem_index element
    for (int i = 0; i < DMQ2::elem_nodes; i++) {
      // get the local dof of the i-th index within thelement
      const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
      // convert the local dof number to global dof number
      Dof d = dm.local_to_global_dof(local_dof);

      // do something. In this case we are printing.
      printf("lgm: i: %li, local: %li, global: %li, proc: %li\n",i, local_dof,
             d.get_gid(), d.get_proc());
    }
  });
}

// print the local and global numbering of the ghost dofs per process
void print_ghost_dofs(const DMQ2 dm) {
  Kokkos::parallel_for(
      dm.get_ghost_lg_map().capacity(), KOKKOS_LAMBDA(Integer i) {
        if (dm.get_ghost_lg_map().valid_at(i)) {
          auto sfc = dm.get_ghost_lg_map().key_at(i);
          auto global_dof = dm.get_ghost_lg_map().value_at(i);
          printf("local: %li, global: %li - proc: %li \n",
                 dm.sfc_to_local(sfc), global_dof.get_gid(),
                 global_dof.get_proc());
        }
      });
}

template <Integer Type>
MARS_INLINE_FUNCTION void
get_elem_coordinates(const DMQ2 &dm, const Integer elem_index, double *points) {

  SFC<Type> dof = dm.get_local_dof_enum();

  for (int i = 0; i < DMQ2::elem_nodes; i++) {
    const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
    const Integer sfc_elem = dm.local_to_sfc(local_dof);
    get_vertex_coordinates_from_sfc<Type>(sfc_elem, points + 3 * i,
                                          dof.get_XDim(), dof.get_YDim(),
                                          dof.get_ZDim());
  }
}

// template <typename T, typename Scalar>
// MARS_INLINE_FUNCTION bool invert3(const T *mat, T *mat_inv, const Scalar
// &det) {
//   assert(det != 0.);

//   if (det == 0.) {
//     return false;
//   }

//   mat_inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / det;
//   mat_inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / det;
//   mat_inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / det;
//   mat_inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / det;
//   mat_inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / det;
//   mat_inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / det;
//   mat_inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / det;
//   mat_inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / det;
//   mat_inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / det;
//   return true;
// }

// q.points = {{0.5, 0.5},
// {0.98304589153964795245728880523899, 0.5},
// {0.72780186391809642112479237299488, 0.074042673347699754349082179816666},
// {0.72780186391809642112479237299488, 0.92595732665230024565091782018333},
// {0.13418502421343273531598225407969, 0.18454360551162298687829339850317},
// {0.13418502421343273531598225407969, 0.81545639448837701312170660149683}};

// q.weights = {0.28571428571428571428571428571428,
// 0.10989010989010989010989010989011,
// 0.14151805175188302631601261486295,
// 0.14151805175188302631601261486295,
// 0.16067975044591917148618518733485,
// 0.16067975044591917148618518733485};

/*
template <typename T, typename Scalar>
MARS_INLINE_FUNCTION bool invert2(const T *mat, T *mat_inv, const Scalar &det) {
  mat_inv[0] = mat[3] / det;
  mat_inv[1] = -mat[1] / det;
  mat_inv[2] = -mat[2] / det;
  mat_inv[3] = mat[0] / det;
  return true;
} */

template<Integer Type>
void compute_invJ_and_detJ(const DMQ2 &dm, ViewVectorType<double> detJ,
                           ViewMatrixType<double> invJ) {

  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {

    double J[4];
    double point_ref[2];
    double next_point[2];

    Integer local_dof = dm.get_elem_local_dof(elem_index, 0);
    dm.get_dof_coordinates_from_local<Type>(local_dof, point_ref);

    local_dof = dm.get_elem_local_dof(elem_index, 1);
    dm.get_dof_coordinates_from_local<Type>(local_dof, next_point);

    // col 0, p1
    J[0] = next_point[0] - point_ref[0];
    J[2] = next_point[1] - point_ref[1];

    // we skip p2

    local_dof = dm.get_elem_local_dof(elem_index, 3);
    dm.get_dof_coordinates_from_local<Type>(local_dof, next_point);

    // col 1, p3
    J[1] = next_point[0] - point_ref[0];
    J[3] = next_point[1] - point_ref[1];

    // determinant
    const double det_J = J[0] * J[3] - J[2] * J[1];

    //fill out the views
    invert2(J, &invJ(elem_index, 0), det_J);
    detJ(elem_index) = det_J;
  });
}


template<Integer INPUT>
MARS_INLINE_FUNCTION
void gather_elem_data(const DMQ2& dm, const Integer elem_index,
        DMDataType<INPUT>* sol)
{
  for (int i = 0; i < DMQ2::elem_nodes; i++)
  {
    // forach dof get the local number
    const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
    // use the local number to read the corresponding user data
    sol[i] = dm.get_dof_data<INPUT>(local_dof);
  }
}

//if non-linear than the quad rule should be computed per quad point and not
//anymore for each element so the coalescing of the Jinv will not matter.
//In that case maybe a better way to go is parallel through the quad points.
template <Integer INPUT>
void integrate(const DMQ2 &dm, const FEQuad4<double>::Quadrature& quad,
        ViewVectorType<double> det_J, ViewMatrixType<double> J_inv,
        ViewMatrixType<double> res) {

  constexpr int n_qp = FEQuad4<double>::Quadrature::n_points();
  constexpr int dim = FEQuad4<double>::Quadrature::dim();

  ViewVectorTextureC<double, n_qp> q_weights = quad.q_w;
  ViewMatrixTextureC<double, n_qp, dim> q_points = quad.q_p;

  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
    double gi[dim], g[dim];

    double sol[DMQ2::elem_nodes];
    gather_elem_data<INPUT>(dm, elem_index, sol);

    for (int k = 0; k < n_qp; ++k) {
      // we need 2nd order quadrature rule for quads!
      /* double q[2] = {q_points[k][0], q_points[k][1]}; */
      /* double w = q_weights[k]; */

      FEQuad4<double>::Grad::ref(&q_points(k,0), sol, gi); // compute it once

      g[0] = J_inv(elem_index, 0) * gi[0] + J_inv(elem_index, 2) * gi[1];
      g[1] = J_inv(elem_index, 1) * gi[0] + J_inv(elem_index, 3) * gi[1];

      const double detj = det_J(elem_index);

      for (int i = 0; i < DMQ2::elem_nodes; i++) {
        // forach dof get the local number
        const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
        FEQuad4<double>::Grad::affine_f(i, &J_inv(elem_index, 0), &q_points(k, 0), gi);

        // dot(grad sol, grad phi_i)
        const double dot_grads = (gi[0] * g[0] + gi[1] * g[1]);
        const double integr = q_weights(k) * detj * dot_grads;

        res(elem_index, i) += integr;
        assert(res(elem_index, i) == res(elem_index, i));
      }
    }
  });
}

template <Integer OUTPUT>
void add_dof_contributions(const DMQ2& dm, const ViewMatrixType<double> &res) {
  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
    // update output
    for (int i = 0; i < DMQ2::elem_nodes; i++) {
      const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
      /* atomically updated the contributions to the same dof */
      Kokkos::atomic_add(&dm.get_dof_data<OUTPUT>(local_dof), res(elem_index, i));
    }
  });
}

// form the matrix free operator
void form_operator(const DMQ2 dm, const FEQuad4<double>::Quadrature& quad) {

  using Elem = DistributedQuad4Mesh::Elem;

  ViewVectorType<double> detJ("detJ", dm.get_elem_size());
  ViewMatrixType<double> invJ("J_inv", dm.get_elem_size(), 4);
  compute_invJ_and_detJ<Elem::ElemType>(dm, detJ, invJ);

  ViewMatrixType<double> res("res", dm.get_elem_size(), DMQ2::elem_nodes);
  integrate<INPUT>(dm, quad, detJ, invJ, res);

  add_dof_contributions<OUTPUT>(dm, res);
}

/*
// form the matrix free operator
template <Integer u, Integer v>
void form_operator(const DMQ2 dm, const int rank) {
  // go through each element and for each element iterate through its dofs
  using Elem = DistributedQuad4Mesh::Elem;

  dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
    double points[DMQ2::elem_nodes * 3];
    double J[2 * 2], J_inv[2 * 2];
    double sol[DMQ2::elem_nodes];
    double res[DMQ2::elem_nodes];
    double gi[2], g[2];

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Uniform data ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    static const int n_qp = 6;
    double q_points[n_qp][2] = {{0.5, 0.5},
                                {0.98304589153964795245728880523899, 0.5},
                                {0.72780186391809642112479237299488,
                                 0.074042673347699754349082179816666},
                                {0.72780186391809642112479237299488,
                                 0.92595732665230024565091782018333},
                                {0.13418502421343273531598225407969,
                                 0.18454360551162298687829339850317},
                                {0.13418502421343273531598225407969,
                                 0.81545639448837701312170660149683}};

    double q_weights[n_qp] = {
        0.28571428571428571428571428571428, 0.10989010989010989010989010989011,
        0.14151805175188302631601261486295, 0.14151805175188302631601261486295,
        0.16067975044591917148618518733485, 0.16067975044591917148618518733485};

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Geometric data ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    get_elem_coordinates<Elem::ElemType>(dm, elem_index, points);

    // p0
    const double x0 = points[0];
    const double y0 = points[1];

    // col 0, p1
    J[0] = points[3] - x0;
    J[2] = points[4] - y0;

    // we skip p2

    // col 1, p3
    J[1] = points[9] - x0;
    J[3] = points[10] - y0;

    // determinant
    const double det_J = J[0] * J[3] - J[2] * J[1];

    invert2(J, J_inv, det_J);

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Read solution (global 2 local) ////////////////
    ////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < DMQ2::elem_nodes; i++) {
      // forach dof get the local number
      const Integer local_dof = dm.get_elem_local_dof(elem_index, i);

      // use the local number to read the corresponding user data
      sol[i] = dm.get_dof_data<INPUT>(local_dof);
      res[i] = 0.0;
    }

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Integrate /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    for (int k = 0; k < n_qp; ++k) {
      // we need 2nd order quadrature rule for quads!
      double q[2] = {q_points[k][0], q_points[k][1]};
      double w = q_weights[k];

      FEQuad4<double>::Grad::ref(q, sol, gi); //compute it once

      g[0] = J_inv[0] * gi[0] + J_inv[2] * gi[1];
      g[1] = J_inv[1] * gi[0] + J_inv[3] * gi[1];

      for (int i = 0; i < DMQ2::elem_nodes; i++) {
        // forach dof get the local number
        const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
        FEQuad4<double>::Grad::affine_f(i, J_inv, q, gi);

        // dot(grad sol, grad phi_i)
        const double dot_grads = (gi[0] * g[0] + gi[1] * g[1]);
        const double integr = w * det_J * dot_grads;

        res[i] += integr;

        assert(res[i] == res[i]);
      }
    }

    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Local 2 global ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < DMQ2::elem_nodes; i++) {
      const Integer local_dof = dm.get_elem_local_dof(elem_index, i);

      // atomically updated the contributions to the same dof
      Kokkos::atomic_add(&dm.get_dof_data<OUTPUT>(local_dof), res[i]);
    }
  });
} */

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

// template <Integer Dim> struct BC {
//   template <Integer Type, Integer Dir>
//   MARS_INLINE_FUNCTION void operator()(Face<Type, Dir> &face) const {
//     std::cout << face.get_sides()[0].elem_id << " "
//               << face.get_sides()[0].is_boundary() << std::endl;
//   }
// };

template<class BC, class RHS, class AnalyticalFun>
void poisson_2D(int &argc, char **&argv, const int level) {

  using namespace mars;
  try {
    mars::proc_allocation resources;
#ifdef WITH_MPI
    // initialize MPI
    marsenv::mpi_guard guard(argc, argv, false);

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
    /* print_elem_global_dof(dm); */

    /* using SMesh = mars::Mesh<Dim, PMesh::ManifoldDim>; */
    using VectorReal = mars::ViewVectorType<Real>;
    using VectorInt = mars::ViewVectorType<Integer>;
    using VectorBool = mars::ViewVectorType<bool>;

    BC bc_fun;
    RHS rhs_fun;
    AnalyticalFun an_fun;

    const Integer dof_size = dm.get_dof_size();

    VectorReal x("X", dof_size);
    VectorReal rhs("rhs", dof_size);

    if (proc_num == 0) {
        std::cout << "form_operator..." << std::flush;
    }

    PoissonOperator<DMQ2, INPUT, OUTPUT> pop(context, dm);
    pop.init();

    if (proc_num == 0) {
        std::cout << "DONE" << std::endl;
    }

    pop.assemble_rhs(rhs, rhs_fun);

    /* print_ghost_dofs(dm); */
    double time = timer.seconds();
    std::cout << "Setup took: " << time << " seconds." << std::endl;

    //if no facenr specified all the boundary is processed. If more than one and less than all
    //is needed than choose all (do not provide face number) and check manually within the lambda
    dm.boundary_dof_iterate<INPUT>(
        MARS_LAMBDA(const Integer local_dof, DMDataType<INPUT> &value) {
            /* do something with the local dof number if needed.
            For example: If no face nr is specified at the boundary dof iterate: */
            /*if (dm.is_boundary<Type, left>(local_dof) || dm.is_boundary<Type, up>(local_dof))*/
            double point[2];
            dm.get_dof_coordinates_from_local<Type>(local_dof, point);
            /* bc_fun(point, value); */
            bc_fun(point, x(local_dof));
            bc_fun(point, rhs(local_dof));
        });

    /* Integer num_iter = 0;
    bcg_stab(pop, *prec_ptr, rhs, 10 * rhs.extent(0), x, num_iter); */

    std::cout << "[" << proc_num
              << "] ndofs : " << dm.get_local_dof_enum().get_elem_size()
              << std::endl;

    /* dm.dof_iterate(MARS_LAMBDA(const Integer i) {
      printf("ggid: %li, INPUT: %lf, OUTPUT: %lf, rank: %i\n", i,
             dm.get_dof_data<INPUT>(i), dm.get_dof_data<OUTPUT>(i),
             proc_num);
    });
 */
    time = timer.seconds();
    std::cout << "Matrix free took: " << time << " seconds." << std::endl;

#endif
  } catch (std::exception &e) {
    std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
  }
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
    /* print_elem_global_dof(dm); */

    /* print_ghost_dofs(dm); */
    double time = timer.seconds();
    std::cout << "Setup took: " << time << " seconds." << std::endl;

    // use as more readable tuple index to identify the data
    constexpr int u = 0;
    constexpr int v = 1;

    /* const Integer locall_owned_dof_size = dm.get_locally_owned_dof_size();
    ViewVectorType<double> x("x", locall_owned_dof_size);

    Kokkos::parallel_for(
        "initglobaldatavalues", locall_owned_dof_size,
        MARS_LAMBDA(const Integer i) { x(i) = 2 * i; });

    dm.set_locally_owned_data<INPUT>(x); */

    constexpr Integer left = 0; // 0 for the left face, 1 for the right, 2 and 3 for down and up.
    constexpr Integer up = 3;
    //if no facenr specified all the boundary is processed. If more than one and less than all
    //is needed than choose all (do not provide face number) and check manually within the lambda
    dm.boundary_dof_iterate<INPUT>(
        MARS_LAMBDA(const Integer local_dof, DMDataType<INPUT> &value) {
            //do something with the local dof number if needed.
            //For example: If no face nr is specified at the boundary dof iterate:
            /* if (dm.is_boundary<Type, left>(local_dof) || dm.is_boundary<Type, up>(local_dof)) */
            Example2Dirichlet func;

            double point[2];
            dm.get_dof_coordinates_from_local<Type>(local_dof, point);
            func(point, value);
        });

    /* dm.dof_iterate(MARS_LAMBDA(const Integer i) {
        const Integer gd= dm.local_to_global(i);
        printf("llid: %li, ggid: %li, INPUT: %lf, OUTPUT: %lf, rank: %i\n", i, gd,
               dm.get_dof_data<INPUT>(i), dm.get_dof_data<OUTPUT>(i), proc_num);
    }); */

    // local dof enum size
    const Integer dof_size = dm.get_dof_size();

    // initialize the values by iterating through local dofs
    Kokkos::parallel_for(
        "initdatavalues", dof_size, MARS_LAMBDA(const Integer i) {
          dm.get_dof_data<INPUT>(i) = 1.0;
          // dm.get_dof_data<INPUT>(i) = i;
          dm.get_dof_data<OUTPUT>(i) = 0.0;
        });


    /* dm.get_locally_owned_data<INPUT>(x); */

    /* Kokkos::parallel_for(
        "printglobaldatavalues", locall_owned_dof_size,
        MARS_LAMBDA(const Integer i) {
          printf("i: %li, gdata: %lf - rank: %i\n", i, x(i), proc_num);
        }); */

    // specify the tuple indices of the tuplelements that are needed to gather.
    // if no index specified it gathers all views of the tuple. All data.
    dm.gather_ghost_data<INPUT>(context);

    if (proc_num == 0) {
      std::cout << "form_operator..." << std::flush;
    }

    const FEQuad4<double>::Quadrature quad = FEQuad4<double>::Quadrature::make();
    form_operator(dm, quad);

    if (proc_num == 0) {
      std::cout << "DONE" << std::endl;
    }

    // BC<2> bc;
    // dm.face_iterate(bc);

    // iterate through the local dofs and print the local number and the data
    /* dm.dof_iterate(
        MARS_LAMBDA(const Integer i) {
            printf("lid: %li, u: %lf, v: %lf, rank: %i\n", i,
                   dm.get_dof_data<u>(i), dm.get_dof_data<v>(i), proc_num);
        }); */

    scatter_add_ghost_data<OUTPUT>(dm, context);

    std::cout << "[" << proc_num
              << "] ndofs : " << dm.get_local_dof_enum().get_elem_size()
              << std::endl;

    /* dm.dof_iterate(MARS_LAMBDA(const Integer i) {
      printf("ggid: %li, INPUT: %lf, OUTPUT: %lf, rank: %i\n", i,
             dm.get_dof_data<INPUT>(i), dm.get_dof_data<OUTPUT>(i),
             proc_num);
    });
 */
    time = timer.seconds();
    std::cout << "Matrix free took: " << time << " seconds." << std::endl;

#endif
  } catch (std::exception &e) {
    std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
  }
}
} // namespace mars
