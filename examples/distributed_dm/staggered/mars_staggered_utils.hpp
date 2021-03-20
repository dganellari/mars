#ifndef MARS_SU_STOKES_
#define MARS_SU_STOKES_

#include "mars.hpp"

namespace mars {

    namespace stag {

        static constexpr Integer Degree = 2;

        using DHandler = DofHandler<DistributedQuad4Mesh, Degree>;

        // Finite Element DM
        using FEDM = DM<DHandler, double, double>;

        // Staggered DofHandlers
        using VDH = VolumeDofHandler<DHandler>;
        using CDH = CornerDofHandler<DHandler>;
        using FDH = FaceDofHandler<DHandler>;

        using FVDH = FaceVolumeDofHandler<DHandler>;
        using CVDH = CornerVolumeDofHandler<DHandler>;

        // Staggered DMs
        using VolumeDM = SDM<VDH, double, double>;
        using CornerDM = SDM<CDH, double>;
        using FaceDM = SDM<FDH, double, double>;

        using FaceVolumeDM = SDM<FVDH, double, double>;
        using CornerVolumeDM = SDM<CVDH, double>;

        /*
            using VolumeDM = VDM<DHandler, double, double>;
            using CornerDM = CDM<DHandler, double>;
            using FaceDM = FDM<DHandler, double, double, double>;
         */

        static constexpr bool Orient = true;

        using SStencil = StokesStencil<DofLabel::lFace>;
        using FSStencil = FullStokesStencil<DofLabel::lFace>;
        // general width 2 stencil used as constant viscosity stokes.
        using VStencil = Stencil<DofLabel::lVolume>;
        // general width 1 stencil used as pressure stencil.
        using CStencil = Stencil<DofLabel::lCorner>;

        using SPattern = SparsityPattern<double, default_lno_t, unsigned long, FVDH, VStencil, SStencil>;
        /* using SPattern = SparsityPattern<double, VStencil, FSStencil>; */
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
                    printf("dof: %li - gdof: %li - label: %li --- (%lf, %lf) - rank: %i\n",
                           i,
                           d.get_gid(),
                           label,
                           point[0],
                           point[1],
                           d.get_proc());
                });
        }

        template <typename SDM>
        void print_partition_boundary_dofs(const SDM dm) {
            const Integer size = dm.get_dof_handler().get_boundary_dof_size();
            Kokkos::parallel_for(
                "for", size, MARS_LAMBDA(const int index) {
                    // go through all the dofs of the elem_index element
                    const Integer local_dof = dm.get_dof_handler().get_boundary_dof(index);
                    /* const Integer local_dof = dm.sfc_to_local(sfc); */
                    // convert the local dof number to global dof number
                    Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                    // do something. In this case we are printing.
                    printf("Bounary Volume dof: i: %li, local: %li, global: %li, proc: %li\n",
                           index,
                           local_dof,
                           d.get_gid(),
                           d.get_proc());
                });
        }

        template <typename SDM>
        void print_ghost_dofs(const SDM dm) {
            const Integer size = dm.get_dof_handler().get_ghost_dof_size();
            Kokkos::parallel_for(
                "for", size, MARS_LAMBDA(const int index) {
                    // go through all the dofs of the elem_index element
                    const Integer local_dof = dm.get_dof_handler().get_ghost_dof(index);
                    /* const Integer local_dof = dm.sfc_to_local(sfc); */
                    // convert the local dof number to global dof number
                    Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                    // do something. In this case we are printing.
                    printf("Ghost Volume dof: i: %li, local: %li, global: %li, proc: %li\n",
                           index,
                           local_dof,
                           d.get_gid(),
                           d.get_proc());
                });
        }

        template <typename S>
        void print_sparsity_pattern(S &sp) {
            const Integer size = sp.get_num_rows();
            Kokkos::parallel_for(
                "for", size, MARS_LAMBDA(const int row) {
                    const Integer start = sp.get_row_map(row);
                    const Integer end = sp.get_row_map(row + 1);

                    // print only if end - start > 0. Otherwise segfaults.
                    // The row index is not a global index of the current process!
                    for (int i = start; i < end; ++i) {
                        auto value = sp.get_value(i);
                        auto col = sp.get_col(i);
                        // do something. In this case we are printing.

                        const Integer local_dof = sp.get_dof_handler().get_owned_dof(row);
                        const Integer global_row = sp.get_dof_handler().get_dof_handler().local_to_global(local_dof);

                        const Integer local_col = sp.get_dof_handler().global_to_local(col);
                        const Integer global_col = sp.get_dof_handler().get_dof_handler().local_to_global(local_col);

                        printf("row_dof: %li - %li, col_dof: %li - %li, value: %lf\n",
                               row,
                               global_row,
                               col,
                               global_col,
                               value);
                    }
                });
        }

        template <typename SDM>
        void print_local_dofs(const SDM dm) {
            const Integer size = dm.get_dof_size();
            Kokkos::parallel_for(
                "for", size, MARS_LAMBDA(const int index) {
                    // go through all the dofs of the elem_index element
                    const Integer local_dof = dm.get_local_dof(index);
                    const Integer sfc = dm.local_to_sfc(local_dof);
                    /* const Integer dir = dm.get_orientation(local_dof); */
                    // convert the local dof number to global dof number
                    Dof d = dm.local_to_global_dof(local_dof);
                    auto o = dm.get_octant_from_local(local_dof);

                    // do something. In this case we are printing.
                    printf("Dof: i: %li, local: %li, Dir: %li,  o: [%li, %li], global: %li, proc: %li\n",
                           index,
                           local_dof,
                           sfc,
                           o.x,
                           o.y,
                           d.get_gid(),
                           d.get_proc());
                });
        }

        template <typename SDM>
        void print_owned_dofs(const SDM dm) {
            const Integer size = dm.get_dof_handler().get_owned_dof_size();
            Kokkos::parallel_for(
                "for", size, MARS_LAMBDA(const int index) {
                    // go through all the dofs of the elem_index element
                    const Integer local_dof = dm.get_dof_handler().get_owned_dof(index);
                    // convert the local dof number to global dof number
                    Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                    // do something. In this case we are printing.
                    printf("Owned Dof: i: %li, local: %li, global: %li, proc: %li\n",
                           index,
                           local_dof,
                           d.get_gid(),
                           d.get_proc());
                });
        }

        // print thlocal and the global dof nbh of each owned dofs.
        template <typename H, typename FE>
        void print_dof_to_dof_map(const H handler, const FE fe) {
            fe.dof_to_dof_iterate(MARS_LAMBDA(const Integer owned_index) {
                // go through all the dofs of the elem_index element
                for (int i = 0; i < fe.get_nbh_dof_size(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = fe.get_nbh_dof(owned_index, i);
                    // convert the local dof number to global dof number
                    Dof d = handler.local_to_global_dof(local_dof);

                    // do something. In this case we are printing.
                    printf(
                        "DTD: i: %li, local: %li, global: %li, proc: %li\n", i, local_dof, d.get_gid(), d.get_proc());
                }
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
                    printf("lgm: %li - i: %li, local: %li, global: %li, proc: %li\n",
                           H::dofLabel,
                           i,
                           local_dof,
                           d.get_gid(),
                           d.get_proc());
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
        void print_stencil(const SDM dm, const Stencil stencil) {
            printf("Length: %li\n", stencil.get_length());
            stencil.dof_iterate(MARS_LAMBDA(const Integer stencil_index, const Integer local_dof) {
                if (local_dof > -1) {
                    // convert the local dof number to global dof number
                    Dof d = dm.local_to_global_dof(local_dof);

                    printf("Stencil: i: %li, local: %li, global: %li, proc: %li\n",
                           stencil_index,
                           local_dof,
                           d.get_gid(),
                           d.get_proc());
                } else {
                    printf("Stencil: i: %li, local: %li\n", stencil_index, local_dof);
                }
            });
        }

    //shows how to use most of the mars utilities.
    template <class BC, class RHS, class AnalyticalFun>
    void staggered_dm_multi_test(const int level) {
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

        double mtime = timer.seconds();
        std::cout << "Mesh generation took: " << mtime << " seconds." << std::endl;

        /* constexpr Integer Dim = DistributedQuad4Mesh::Dim; */

        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        DHandler dof_handler(&mesh, context);
        dof_handler.enumerate_dofs();

        double dtime = timer.seconds();
        std::cout << "Dof Handler enumeration took: " << dtime << " seconds." << std::endl;

        FEDM fedm(dof_handler);

        /* print_ghost_dofs(fedm); */

        auto fe = build_fe_dof_map(dof_handler);
        /* print_elem_global_dof(dof_handler, fe); */

        /* auto fed = build_dof_to_dof_map(dof_handler); */
        /* print_dof_to_dof_map(dof_handler, fed); */

        // it gives the size of the local dofs of the dm. If volume then only volume dofs.
        const Integer dof_size = dof_handler.get_dof_size();

        // if manually managed the data view should have the size of local dof size.
        ViewVectorType<double> data("IN", dof_size);
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) { data(i) = proc_num; });

        gather_ghost_data(dof_handler, data);
        scatter_add_ghost_data(dof_handler, data);

        /* dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) {
            const auto idata = data(i);

            Dof d = dof_handler.local_to_global_dof(i);

            printf("vlid: %li, u: %lf, global: %li, rank: %i\n", i, idata, d.get_gid(), d.get_proc());
        }); */

        // create the DM object from the dof handler
        /* VolumeDofHandler vdh(dof_handler); */
        /* FaceDofHandler fdh(dof_handler); */

        // use a part of the dof handler (separate it and use only volume for example)
        /* VolumeDM vdm(vdh); */
        /* VolumeDM vdm(dof_handler); */
        FaceDM fdm(dof_handler);
        /* CornerDM cdm(dof_handler); */

        /* auto fe = build_fe_dof_map(vdm.get_dof_handler()); */
        /* print_elem_global_dof(vdm.get_dof_handler(), fe); */

        auto ffe = build_fe_dof_map(fdm.get_dof_handler());
        print_elem_global_dof(fdm.get_dof_handler(), ffe);

        /* auto cfe = build_fe_dof_map(cdm.get_dof_handler());
        print_elem_global_dof(cdm.get_dof_handler(), cfe);<] */

        /* print_local_dofs(fdm); */
        /* print_owned_dofs(fdm); */

        /* print_partition_boundary_dofs(vdm); */
        /* print_ghost_dofs(vdm); */

        /* If you Orient then the same order is applied through all stencil
         * based on the orientation. Otherwise no order but just a normal stencil. */
        /* auto volume_stencil = vdm.build_stencil<VStencil, Orient>(); */
        /* auto volume_stencil = build_stencil<VStencil>(vdm.get_dof_handler()); */
        /* print_stencil(vdm.get_dof_handler(), volume_stencil); */

        auto face_stencil = build_stencil<SStencil>(fdm.get_dof_handler());
        /* print_stencil(fdm.get_dof_handler(), face_stencil); */

        double time = timer.seconds();
        std::cout << "Stag DM test took: " << time << " seconds." << std::endl;

#endif
    }

    }  // namespace stag
}  // namespace mars

#endif