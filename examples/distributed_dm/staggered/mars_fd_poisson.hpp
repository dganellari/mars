#ifndef MARS_FD_POISSON_
#define MARS_FD_POISSON_

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
#include "Kokkos_ArithTraits.hpp"
#include <KokkosBlas1_sum.hpp>
#include "mars_distributed_staggered_data_management.hpp"
#endif  // WITH_KOKKOS
#endif

// #include "mars_pvtu_writer.hpp"  // VTK

namespace mars {

    static constexpr  Integer DIM = DistributedQuad4Mesh::Dim;
    static constexpr  bool Orient = true;

    using SStencil = StokesStencil<DIM>;
    using FSStencil = FullStokesStencil<DIM>;
    //general width 2 stencil used as constant viscosity stokes.
    using FStencil = Stencil<DIM, 2>;
    //general width 1 stencil used as pressure stencil.
    using VCStencil = Stencil<DIM, 1>;

    //use the FDDM (finite difference DM object). In this case build your stencil manually!
    /* using SDM = FDDM<DistributedQuad4Mesh, 2, double, double>; */

    //Or use the staggered DM object which builds needed staggered stencils for you!
    /* using SDM = StagDM<DistributedQuad4Mesh, FSStencil, double, double>; */

    /*
    enum class DMDataDesc
    {
        v = 0,
        u = 1
    };
     */

    using VolumeDM = VDM<DistributedQuad4Mesh, 2, double, double>;
    using DofH = DofHandler<DistributedQuad4Mesh, 2>;

    // use as more readable tuple index to identify the data
    static constexpr int IN = 0;
    static constexpr int OUT = 1;

    template <Integer idx>
    using VDMDataType = typename VolumeDM::UserDataType<idx>;

    /* template <typename... T>
    using tuple = mars::ViewsTuple<T...>;

    using dm_tuple = typename SDM::user_tuple; */

    template <class SDM, Integer Type>
    void print_dofs(const SDM &data, const int rank) {
        auto dm = data.get_dof_handler();
        SFC<Type> dof = dm.get_local_dof_enum();
        Kokkos::parallel_for(
            "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
                const Integer sfc_elem = dm.local_to_sfc(i);
                const Integer label = dof.get_label(i);
                Dof d = dm.local_to_global_dof(i);

                double point[3];
                dm.get_dof_coordinates_from_sfc<Type>(sfc_elem, point);
                printf("dof: %li - gdof: %li - label: %li --- (%lf, %lf) - rank: %i\n", i, d.get_gid(), label, point[0], point[1], d.get_proc());
            });
    }

    template<typename SDM>
    void print_partition_boundary_dofs(const SDM dm) {
        const Integer size = dm.get_boundary_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof= dm.get_boundary_dof(index);
                /* const Integer local_dof = dm.sfc_to_local(sfc); */
                // convert the local dof number to global dof number
                Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "Bounary Volume dof: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    template <typename SDM>
    void print_ghost_dofs(const SDM dm) {
        const Integer size = dm.get_ghost_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_ghost_dof(index);
                /* const Integer local_dof = dm.sfc_to_local(sfc); */
                // convert the local dof number to global dof number
                Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "Ghost Volume dof: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    template <typename SDM>
    void print_local_dofs(const SDM dm) {
        const Integer size = dm.get_local_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_local_dof(index);
                // convert the local dof number to global dof number
                Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "Dof: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    template <typename SDM>
    void print_owned_dofs(const SDM dm) {
        const Integer size = dm.get_owned_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_owned_dof(index);
                // convert the local dof number to global dof number
                Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "Owned Dof: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(), d.get_proc());
            });
    }

    /* void print_face_locally_owned(const SDM dm) {
        const Integer size = dm.get_locally_owned_face_dofs().extent(0);
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_face_dof(index);
                // convert the local dof number to global dof number
                Dof d = dm.local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf("lofd: i: %li, local: %li, global: %li, proc: %li\n", index, local_dof, d.get_gid(),
    d.get_proc());
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
    } */

    // print thlocal and the global number of the dof within each element.
    // the dof enumeration within eachlement is topological

    template <typename SDM>
    void print_elem_global_dof(const SDM data) {
        auto dm = data.get_dof_handler();
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
    template <typename SDM>
    void print_ghost_map(const SDM data) {
        auto dm = data.get_dof_handler();
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

    template <typename SDM, typename Stencil>
    void print_stencil(const SDM data, const Stencil stencil) {
        auto dm = data.get_dof_handler();
        stencil.dof_iterate(MARS_LAMBDA(const Integer stencil_index, const Integer local_dof) {
            if (local_dof > -1) {
                // convert the local dof number to global dof number
                Dof d = dm.local_to_global_dof(local_dof);

                printf(
                    "Stencil: i: %li, local: %li, global: %li, proc: %li\n", stencil_index, local_dof, d.get_gid(), d.get_proc());
            } else {
                printf("Stencil: i: %li, local: %li\n", stencil_index, local_dof);
            }
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
        DofH dof_handler(&mesh, context);
        dof_handler.enumerate_dofs(context);

        VolumeDM vdm(dof_handler);
        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        // print locally owned dof numbering
        /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */

        // print local dof numbering
        /* print_dof<Type>(dm.get_local_dof_enum(), proc_num); */
        /* print_dofs<Type>(dm, proc_num); */
        print_local_dofs(vdm);
        print_owned_dofs(vdm);

        print_partition_boundary_dofs(vdm);
        print_ghost_dofs(vdm);

        /* classic width 1 stencil on volume nodes. */
        //If finite difference DM fd_dm object is used then:
        /* auto volume_stencil = mars::build_volume_stencil<VCStencil>(vdm); */

        //Recommended: if separated volume dm is used then use as follow:
        auto volume_stencil = vdm.build_stencil<VCStencil>();
        print_stencil(vdm, volume_stencil);

        /* classic width 2 stencil on face nodes. */
        /* auto face_stencil = mars::build_face_stencil<FStencil, Orient>(dm);
        print_stencil(dm, face_stencil); */

        /* classic width 2 stencil on face nodes. */
        /* auto corner_stencil = mars::build_corner_stencil<VCStencil>(dm);
        print_stencil(dm, corner_stencil); */


        /* dm.build_pressure_stencil();
        print_stencil(dm, dm.get_pressure_stencil()); */
        /* dm.build_stokes_stencil<Orient>(); */
        /* print_stencil(dm, dm.get_stokes_stencil()); */

        //it gives the size of the local dofs of the dm. If volume then only volume dofs.
        const Integer dof_size = vdm.get_dof_size();

           // initialize the values by iterating through local dofs
        /* Kokkos::parallel_for(
            "initdatavalues", dof_size, MARS_LAMBDA(const Integer i) {
                vdm.get_volume_data<IN>(i) = 1.0;
                vdm.get_volume_data<OUT>(i) = proc_num;
            });
 */
        vdm.iterate(MARS_LAMBDA(const Integer i) {
            vdm.get_data<IN>(i) = 1.0;
            vdm.get_data<OUT>(i) = proc_num;
        });

        vdm.gather_ghost_data<OUT>(context);
        scatter_add_ghost_data<VolumeDM, OUT>(vdm, context);


        //print using the dof iterate
        /* vdm.dof_iterate(MARS_LAMBDA(const Integer local_dof) {
            const auto idata = vdm.get_dof_data<IN>(local_dof);
            const auto odata = vdm.get_dof_data<OUT>(local_dof);

            Dof d = vdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n",
                    i, idata, odata, d.get_gid(), d.get_proc());
        }); */

        //print using the data iterate
        vdm.iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = vdm.get_local_dof(i);

            const auto idata = vdm.get_data<IN>(i);
            const auto odata = vdm.get_data<OUT>(i);

            Dof d = vdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n",
                    i, idata, odata, d.get_gid(), d.get_proc());
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
