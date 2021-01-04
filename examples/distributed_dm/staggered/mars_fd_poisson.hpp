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

    /*
    enum class DMDataDesc
    {
        v = 0,
        u = 1
    };
     */

    using FEDM = DM<DistributedQuad4Mesh, 2, double, double>;

    using VolumeDM = VDM<DistributedQuad4Mesh, 2, double, double>;
    /* using FaceDM = FDM<DistributedQuad4Mesh, 2, double, double>; */
    /* using CornerDM = CDM<DistributedQuad4Mesh, 1, double>; */

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
                const Integer dir = dm.get_dof_handler().get_orientation(local_dof);
                // convert the local dof number to global dof number
                Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf(
                    "Dof: i: %li, local: %li, Dir: %li, global: %li, proc: %li\n", index, local_dof, dir, d.get_gid(), d.get_proc());
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

    // print thlocal and the global number of the dof within each element.
    // the dof enumeration within eachlement is topological

    template <typename H, typename FE>
    void print_elem_global_dof(const H handler, const FE fe) {
        fe.iterate(MARS_LAMBDA(const Integer elem_index) {
            // go through all the dofs of the elem_index element
            for (int i = 0; i < fe.get_fe_size(); i++) {
                // get the local dof of the i-th index within thelement
                const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                // convert the local dof number to global dof number
                Dof d = handler.local_to_global_dof(local_dof);

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

        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        DofH dof_handler(&mesh, context);
        dof_handler.enumerate_dofs();

        //create the DM object from the dof handler
        VolumeDM vdm(dof_handler);
        /* FaceDM fdm(dof_handler); */

        FEDM fedm(dof_handler);

        auto fe = build_fe_dof_map(dof_handler);
        print_elem_global_dof(dof_handler, fe);


        // it gives the size of the local dofs of the dm. If volume then only volume dofs.
        const Integer dof_size = dof_handler.get_dof_size();

        //if manually managed the data view should have the size of local dof size.
        ViewVectorType<double> data("IN", dof_size);
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) { data(i) = proc_num; });
        gather_ghost_data(dof_handler, data);

        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) {
            const auto idata = data(i);

            Dof d = dof_handler.local_to_global_dof(i);

            printf("vlid: %li, u: %lf, global: %li, rank: %i\n",
                    i, idata, d.get_gid(), d.get_proc());
        });

        print_local_dofs(vdm);
        print_owned_dofs(vdm);

        print_partition_boundary_dofs(vdm);
        print_ghost_dofs(vdm);

        /* If you Orient then the same order is applied through all stencil
         * based on the orientation. Otherwise no order but just a normal stencil. */
        /* auto volume_stencil = vdm.build_stencil<SStencil, Orient>(); */
        auto volume_stencil = build_stencil<SStencil>(vdm);
        print_stencil(vdm, volume_stencil);

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

        vdm.gather_ghost_data<OUT>();
        scatter_add_ghost_data<VolumeDM, OUT>(vdm);


        //print using the dof iterate
        /* vdm.dof_iterate(MARS_LAMBDA(const Integer local_dof) {
            const auto idata = vdm.get_dof_data<IN>(local_dof);
            const auto odata = vdm.get_dof_data<OUT>(local_dof);

            Dof d = vdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n",
                    i, idata, odata, d.get_gid(), d.get_proc());
        }); */

        //print using the index iterate
        vdm.iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = vdm.get_local_dof(i);

            const auto idata = vdm.get_data<IN>(i);
            const auto odata = vdm.get_data<OUT>(i);

            Dof d = vdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n",
                    i, idata, odata, d.get_gid(), d.get_proc());
        });

        double time = timer.seconds();
        std::cout << "Stag DM Setup took: " << time << " seconds." << std::endl;

#endif
    }

}  // namespace mars

#endif
