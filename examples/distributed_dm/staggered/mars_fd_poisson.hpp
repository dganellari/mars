#ifndef MARS_FD_POISSON_
#define MARS_FD_POISSON_

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
#include <KokkosBlas1_sum.hpp>
#include "mars_distributed_staggered_data_management.hpp"
#endif  // WITH_KOKKOS
#endif

// #include "mars_pvtu_writer.hpp"  // VTK

namespace mars {

    /* using SDM = FDDM<DistributedQuad4Mesh, 2, double, double>; */
    /* using SDM = StagDM<DistributedQuad4Mesh, true, double, double>; */
    using SDM = StagDM<DistributedQuad4Mesh, false, double, double>;
    /*
    enum class DMDataDesc
    {
        v = 0,
        u = 1
    };
     */

    // use as more readable tuple index to identify the data
    static constexpr int IN = 0;
    static constexpr int OU = 1;

    template <Integer idx>
    using SDMDataType = typename SDM::UserDataType<idx>;

    /* template <typename... T>
    using tuple = mars::ViewsTuple<T...>;

    using dm_tuple = typename SDM::user_tuple; */

    template <Integer Type>
    void print_dofs(const SDM &dm, const int rank) {
        SFC<Type> dof = dm.get_local_dof_enum();
        Kokkos::parallel_for(
            "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
                const Integer sfc_elem = dm.local_to_sfc(i);
                const Integer global_dof = dm.local_to_global(i);

                double point[3];
                dm.get_dof_coordinates_from_sfc<Type>(sfc_elem, point);
                printf("dof: %li - gdof: %li --- (%lf, %lf) - rank: %i\n", i, global_dof, point[0], point[1], rank);
            });
    }

    void print_face_locally_owned(const SDM dm) {
        const Integer size = dm.get_locally_owned_face_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_face_dof(index);
                // convert the local dof number to global dof number
                Dof d = dm.local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf("lofd: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    void print_volume_locally_owned(const SDM dm) {
        const Integer size = dm.get_locally_owned_volume_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_volume_dof(index);
                // convert the local dof number to global dof number
                Dof d = dm.local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "lovd: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    void print_corner_locally_owned(const SDM dm) {
        const Integer size = dm.get_locally_owned_corner_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_corner_dof(index);
                // convert the local dof number to global dof number
                Dof d = dm.local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "locd: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    // print thlocal and the global number of the dof within each element.
    // the dof enumeration within eachlement is topological
    void print_elem_global_dof(const SDM dm) {
        dm.elem_iterate(MARS_LAMBDA(const Integer elem_index) {
            // go through all the dofs of the elem_index element
            for (int i = 0; i < SDM::elem_nodes; i++) {
                // get the local dof of the i-th index within thelement
                const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
                // convert the local dof number to global dof number
                Dof d = dm.local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf("lgm: i: %li, local: %li, global: %li, proc: %li\n", i, local_dof, d.get_gid(), d.get_proc());
            }
        });
    }

    // print the local and global numbering of the ghost dofs per process
    void print_ghost_dofs(const SDM dm) {
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

    // print thlocal and the global number of the dof within each element.
    // the dof enumeration within eachlement is topological

    template <typename Stencil>
    void print_stencil(const SDM dm, const Stencil stencil) {
        stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            // go through all the dofs of the elem_index element
            for (int i = 0; i < stencil.get_length(); i++) {
                // get the local dof of the i-th index within thelement
                const Integer local_dof = stencil.get_value(stencil_index, i);

                if (local_dof > -1) {
                    // convert the local dof number to global dof number
                    Dof d = dm.local_to_global_dof(local_dof);

                    // do something. In this case we are printing.
                    printf(
                        "Stencil: i: %li, local: %li, global: %li, proc: %li\n", i, local_dof, d.get_gid(), d.get_proc());
                } else {
                    printf("Stencil: i: %li, local: %li\n", i, local_dof);
                }
            }
            printf("\n\n");
        });
    }

    template <class BC, class RHS, class AnalyticalFun>
    void staggered_poisson_2D(const int level) {
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
        // create the quad mesh distributed through the mpi procs.
        DistributedQuad4Mesh mesh;
        generate_distributed_cube(context, mesh, level, level, 0);

        constexpr Integer Dim = DistributedQuad4Mesh::Dim;

        using Elem = typename DistributedQuad4Mesh::Elem;
        // the type of the mesh elements. In this case quad4 (Type=4)
        constexpr Integer Type = Elem::ElemType;

        // create the dm object
        SDM dm(&mesh, context);
        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        dm.enumerate_dofs(context);
        // print locally owned dof numbering
        /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */

        // print local dof numbering
        /* print_dof<Type>(dm.get_local_dof_enum(), proc_num); */
        /* print_dofs<Type>(dm, proc_num); */

        // print the global dofs for each element's local dof
        /* print_elem_global_dof(dm);

        print_face_locally_owned(dm);
        print_volume_locally_owned(dm);
        print_corner_locally_owned(dm); */
        /* print_ghost_dofs(dm); */

        /* classic width 1 stencil on volume nodes. */
        /* auto volume_stencil = mars::build_volume_stencil(dm);
        print_stencil(dm, volume_stencil); */

        /* classic width 2 stencil on face nodes. */
        /* auto face_stencil = mars::build_face_stencil<2>(dm);
        print_stencil(dm, face_stencil); */

        /* classic width 2 stencil on face nodes. */
        /* auto corner_stencil = mars::build_corner_stencil<1>(dm);
        print_stencil(dm, corner_stencil); */


        /* print_stencil(dm, dm.get_volume_stencil()); */
        print_stencil(dm, dm.get_face_stencil());


        const Integer dof_size = dm.get_dof_size();

           // initialize the values by iterating through local dofs
        Kokkos::parallel_for(
            "initdatavalues", dof_size, MARS_LAMBDA(const Integer i) {
                dm.get_dof_data<IN>(i) = 1.0;
                // dm.get_dof_data<INPUT>(i) = i;
                dm.get_dof_data<OU>(i) = i;
            });

        dm.gather_ghost_data<OU>(context);
        // iterate through the local dofs and print the local number and the data
        dm.dof_iterate(MARS_LAMBDA(const Integer i) {
            Dof d = dm.local_to_global_dof(i);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n",
                    i, dm.get_dof_data<IN>(i), dm.get_dof_data<OU>(i), d.get_gid(), d.get_proc());
        });

        /* using VectorReal = mars::ViewVectorType<Real>;
        using VectorInt = mars::ViewVectorType<Integer>;
        using VectorBool = mars::ViewVectorType<bool>; */

        // if no facenr specified all the boundary is processed. If more than one and less than all
        // is needed than choose all (do not provide face number) and check manually within the lambda
        /* dm.boundary_owned_dof_iterate(MARS_LAMBDA(const Integer owned_dof_index, const Integer sfc) { */
        /* do something with the local dof number if needed.
        For example: If no face nr is specified at the boundary dof iterate: */
        /*if (dm.is_boundary<Type, left>(local_dof) || dm.is_boundary<Type, up>(local_dof))*/
        /* double point[2];
        dm.get_dof_coordinates_from_sfc<Type>(sfc, point);
        [>bc_fun(point, value);<]
        bc_fun(point, x(owned_dof_index));
        bc_fun(point, rhs(owned_dof_index));
    }); */
        /* dm.boundary_dof_iterate<INPUT> does the same except that provides the local_dof num and the reference
         * on the INPUT data to be updated. Suitable when working with the internal data tuple directly */

        /* setting boundary BoundaryConditions */
        /* if (proc_num == 0) {
            std::cout << "Boundary conditions set." << std::flush;
        } */

        double time = timer.seconds();
        std::cout << "Stag DM Setup took: " << time << " seconds." << std::endl;

#endif
    }

}  // namespace mars

#endif
