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

#include "mars_base.hpp"

#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include <KokkosBlas1_sum.hpp>

#include "mars.hpp"
#include "mars_boundary_conditions.hpp"
#include "mars_dm_interpolate.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_poisson_operator.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_quad4.hpp"

// #include "mars_pvtu_writer.hpp"  // VTK

namespace mars {

    template <Integer Type>
    void print_dof(const SFC<Type> &dof, const int rank) {
        Kokkos::parallel_for(
            "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
                Integer sfc_elem = dof.get_view_elements()(i);

                double point[3];
                get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

                printf("dof: %li - (%lf, %lf) - rank: %i\n", i, point[0], point[1], rank);
            });
    }

    /* template <typename DH>
    void print_dofs(const DH &dofhandler, const int rank) {
        SFC<DH::ElemType> dof = dofhandler.get_local_dof_enum();
        Kokkos::parallel_for(
            "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
                const Integer sfc_elem = dofhandler.local_to_sfc(i);
                const Integer global_dof = dofhandler.local_to_global(i);

                double point[3];
                dofhandler.get_dof_coordinates_from_sfc<DH::ElemType>(sfc_elem, point);
                printf("dof: %li - gdof: %li --- (%lf, %lf) - rank: %i\n", i, global_dof, point[0], point[1], rank);
            });
    } */

    // print thlocal and the global number of the dof within each element.
    // the dof enumeration within eachlement is topological
    template <typename DH, typename FE>
    void print_fe_dof_map(const DH &dof_handler, const FE &fe) {
        dof_handler.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
            // go through all the dofs of the elem_index element
            for (int i = 0; i < fe.get_elem_nodes(); i++) {
                // get the local dof of the i-th index within thelement
                const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                // convert the local dof number to global dof number
                Dof d = dof_handler.local_to_global_dof(local_dof);
                auto octant = dof_handler.get_octant_from_local(local_dof);

                auto base_global = dof_handler.compute_base(d.get_gid());
                // do something. In this case we are printing.
                printf("lgm: i: %li, local: %li, octant: [%li, %li, %li], global: %li, base_global: %li, proc: %li\n",
                       i,
                       local_dof,
                       octant.x,
                       octant.y,
                       octant.z,
                       d.get_gid(),
                       base_global,
                       d.get_proc());
            }
        });
    }

    // print the local and global numbering of the ghost dofs per process
    template <typename DH>
    void print_ghost_dofs(const DH &dm) {
        Kokkos::parallel_for(
            dm.get_ghost_lg_map().capacity(), KOKKOS_LAMBDA(Integer i) {
                if (dm.get_ghost_lg_map().valid_at(i)) {
                    auto sfc = dm.get_ghost_lg_map().key_at(i);
                    auto global_dof = dm.get_ghost_lg_map().value_at(i);
                    printf("local: %li, global: %li - proc: %li \n",
                           dm.sfc_to_local(sfc),
                           global_dof.get_gid(),
                           global_dof.get_proc());
                }
            });
    }

    template <class BC, class RHS, class AnalyticalFun, Integer Type = ElementType::Quad4>
    void matrix_free_poisson(const int xDim, const int yDim, const int zDim) {
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

        using DMesh = DistributedMesh<Type>;

        // create the quad mesh distributed through the mpi procs.
        DMesh mesh(context);
        generate_distributed_cube(mesh, xDim, yDim, zDim);

        /* constexpr Integer Dim = DMesh::Dim; */
        /* using Elem = typename DistributedQuad4Mesh::Elem; */
        /* constexpr Integer Type = Elem::ElemType; */

        constexpr Integer Degree = 1;
        using DOFHandler = DofHandler<DMesh, Degree>;
        using DMQ = DM<DOFHandler, double, double, double>;
        using FE = FEDofMap<DOFHandler, false>;
        /* using SPattern = SparsityPattern<double, Integer, unsigned long, DOFHandler>; */

        static_assert(DOFHandler::Block == 1, "Poisson example does not support yet vector valued block structures.");

        // use as more readable tuple index to identify the data
        static constexpr int INPUT = 0;
        static constexpr int OUTPUT = 1;
        static constexpr int RHSD = 2;

        DOFHandler dof_handler(mesh);
        dof_handler.enumerate_dofs();

        /* dof_handler.print_dofs(proc_num); */

        auto fe = build_fe_dof_map<DOFHandler, false>(dof_handler);

        // print the global dofs for each element's local dof
        /* print_fe_dof_map(dof_handler, fe); */
        /* print_ghost_dofs(dm); */

        /* SPattern sp(dof_handler);
        sp.build_pattern(fe);
        sp.print_sparsity_pattern(); */

        // create the dm object
        DMQ dm(dof_handler);
        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        //
        // print locally owned dof numbering
        /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */

        // print local dof numbering
        /* print_dof<Type>(dm.get_local_dof_enum(), proc_num); */
        /* print_dofs<Type>(dm, proc_num); */

        using VectorReal = mars::ViewVectorType<Real>;
        using VectorInt = mars::ViewVectorType<Integer>;
        using VectorBool = mars::ViewVectorType<bool>;

        BC bc_fun;
        RHS rhs_fun;
        // AnalyticalFun an_fun;

        const Integer locally_owned_dof_size = dof_handler.get_owned_dof_size();
        printf("Locally owned dof size: %li\n", locally_owned_dof_size);

        VectorReal x("X", locally_owned_dof_size);
        VectorReal rhs("rhs", locally_owned_dof_size);

        if (proc_num == 0) {
            std::cout << "Init PoissonOperator..." << std::endl;
        }

        PoissonOperator<INPUT, OUTPUT, RHSD, DMQ, FE> pop(context, dm, fe);
        pop.init();

        if (proc_num == 0) {
            std::cout << "DONE" << std::endl;
        }

        pop.assemble_rhs(rhs_fun, rhs);

        // if no facenr specified all the boundary is processed. If more than one and less than all
        // is needed than choose all (do not provide face number) and check manually within the lambda
        dof_handler.boundary_owned_dof_iterate(MARS_LAMBDA(const Integer owned_dof_index, const Integer sfc) {
            /* do something with the local dof number if needed.
            For example: If no face nr is specified at the boundary dof iterate: */
            /*if (dm.is_boundary<Type, left>(local_dof) || dm.is_boundary<Type, up>(local_dof))*/
            double point[2];
            dm.get_dof_handler().template get_dof_coordinates_from_sfc<Type>(sfc, point);
            /* bc_fun(point, value); */
            bc_fun(point, x(owned_dof_index));
            bc_fun(point, rhs(owned_dof_index));
        });
        /* dm.boundary_dof_iterate<INPUT> does the same except that provides the local_dof num and the reference
         * on the INPUT data to be updated. Suitable when working with the internal data tuple directly */

        /* setting boundary BoundaryConditions */
        if (proc_num == 0) {
            std::cout << "Boundary conditions set." << std::flush;
        }

        double time = timer.seconds();
        std::cout << "Setup took: " << time << " seconds." << std::endl;

        CopyOperator preconditioner;
        Integer num_iter = 0;
        Integer max_iter = pop.comm().sum(rhs.extent(0));
        std::cout << "Starting bcg on proc:" << proc_num << std::endl;

        Kokkos::Timer timer_bcg;
        bcg_stab(pop, preconditioner, rhs, max_iter, x, num_iter);

        /* print_ghost_dofs(dm); */
        double time_bcg = timer_bcg.seconds();
        std::cout << "Bicgstab took: " << time_bcg << " seconds on proc: " << proc_num << std::endl;

        std::cout << "[" << proc_num << "] ndofs : " << dm.get_dof_handler().get_local_dof_enum().get_elem_size()
                  << std::endl;

        /* dm.dof_iterate(MARS_LAMBDA(const Integer i) {
          printf("ggid: %li, INPUT: %lf, OUTPUT: %lf, rank: %i\n", i,
                 dm.get_dof_data<INPUT>(i), dm.get_dof_data<OUTPUT>(i),
                 proc_num);
        });  */

        time = timer.seconds();
        std::cout << "Matrix free took: " << time << " seconds on proc: " << proc_num << std::endl;

        /* VectorReal x_exact("X_E", locally_owned_dof_size); */

        /* Real sum_rhs = KokkosBlas::sum(rhs);
        Real tot = context->distributed->sum(sum_rhs);
        printf("sum_rh****: %lf\n", tot);
 */
        /* DMInterpolate<DMQ> ip(dm);
        ip.apply(an_fun, x_exact); */

        // PVTUMeshWriter<DMQ, Type> w;                   // VTK
        // w.write_vtu("poisson_exact.vtu", dm, x_exact);  // VTK
        // w.write_vtu("poisson_rhs.vtu", dm, rhs);        // VTK
        /*
                VectorReal diff("diff", locally_owned_dof_size);
                pop.apply(x_exact, diff);

                Real sum_diff = KokkosBlas::sum(diff);
                Real diffT = context->distributed->sum(sum_diff);
                printf("diff****: %lf\n", diffT);
         */

#endif
    }

}  // namespace mars
#endif
