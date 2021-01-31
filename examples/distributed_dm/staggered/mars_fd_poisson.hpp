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
#include <KokkosBlas1_sum.hpp>
#include "Kokkos_ArithTraits.hpp"
#include "mars_distributed_sparsity_pattern.hpp"
#include "mars_distributed_staggered_data_management.hpp"
#include "mars_distributed_staggered_dof_management.hpp"
#endif  // WITH_KOKKOS
#endif

// #include "mars_pvtu_writer.hpp"  // VTK

namespace mars {

    using DHandler = DofHandler<DistributedQuad4Mesh, 2>;

    // Finite Element DM
    using FEDM = DM<DHandler, double, double>;

    // Staggered DofHandlers
    using VDH = VolumeDofHandler<DHandler>;
    using CDH = CornerDofHandler<DHandler>;
    using FDH = FaceDofHandler<DHandler>;

    // Staggered DMs
    using VolumeDM = SDM<VDH, double, double>;
    using CornerDM = SDM<CDH, double, double, double>;
    using FaceDM = SDM<FDH, double, double>;

    /*
        using VolumeDM = VDM<DHandler, double, double>;
        using CornerDM = CDM<DHandler, double>;
        using FaceDM = FDM<DHandler, double, double, double>;
     */

    static constexpr bool Orient = true;

    using SStencil = StokesStencil<FDH>;
    using FSStencil = FullStokesStencil<FDH>;
    // general width 2 stencil used as constant viscosity stokes.
    using VStencil = Stencil<VDH, 1>;
    // general width 1 stencil used as pressure stencil.
    using CStencil = Stencil<CDH, 1>;

    using SPattern = SparsityPattern<double, DHandler, VStencil, SStencil>;
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

                //print only if end - start > 0. Otherwise segfaults.
                //The row index is not a global index of the current process!
                for (int i = start; i < end; ++i) {
                    auto value = sp.get_value(i);
                    auto col = sp.get_col(i);
                    // do something. In this case we are printing.

                    const Integer local_dof = sp.get_owned_map().global_to_local(row);
                    const Integer global_row = sp.get_owned_map().get_dof_handler().local_to_global(local_dof);

                    const Integer local_col = sp.get_owned_map().global_to_local(col);
                    const Integer global_col = sp.get_owned_map().get_dof_handler().local_to_global(local_col);

                    printf("row_dof: %li - %li, col_dof: %li - %li, value: %lf\n", row, global_row, col, global_col, value);
                }
            });
    }

    template <typename SDM>
    void print_local_dofs(const SDM dm) {
        const Integer size = dm.get_dof_handler().get_dof_size();
        Kokkos::parallel_for(
            "for", size, MARS_LAMBDA(const int index) {
                // go through all the dofs of the elem_index element
                const Integer local_dof = dm.get_dof_handler().get_local_dof(index);
                const Integer dir = dm.get_dof_handler().get_orientation(local_dof);
                // convert the local dof number to global dof number
                Dof d = dm.get_dof_handler().local_to_global_dof(local_dof);

                // do something. In this case we are printing.
                printf("Dof: i: %li, local: %li, Dir: %li, global: %li, proc: %li\n",
                       index,
                       local_dof,
                       dir,
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
                printf("DTD: i: %li, local: %li, global: %li, proc: %li\n", i, local_dof, d.get_gid(), d.get_proc());
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
    void print_stencil(const SDM data, const Stencil stencil) {
        auto dm = data.get_dof_handler();
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

        /* constexpr Integer Dim = DistributedQuad4Mesh::Dim; */

        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        DHandler dof_handler(&mesh, context);
        dof_handler.enumerate_dofs();

        /* FEDM fedm(dof_handler); */

        /* print_ghost_dofs(fedm); */

        /* auto fe = build_fe_dof_map(dof_handler); */
        /* print_elem_global_dof(dof_handler, fe); */

        /* auto fed = build_dof_to_dof_map(dof_handler); */
        /* print_dof_to_dof_map(dof_handler, fed); */

        // it gives the size of the local dofs of the dm. If volume then only volume dofs.
        const Integer dof_size = dof_handler.get_dof_size();

        // if manually managed the data view should have the size of local dof size.
        /* ViewVectorType<double> data("IN", dof_size);
        dof_handler.dof_iterate(MARS_LAMBDA(const Integer i) { data(i) = proc_num; });

        gather_ghost_data(dof_handler, data);
        scatter_add_ghost_data(dof_handler, data); */

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
        VolumeDM vdm(dof_handler);
        FaceDM fdm(dof_handler);
        CornerDM cdm(dof_handler);
        /*
                auto fe = build_fe_dof_map(vdm.get_dof_handler());
                print_elem_global_dof(vdm.get_dof_handler(), fe);

                auto ffe = build_fe_dof_map(fdm.get_dof_handler());
                print_elem_global_dof(fdm.get_dof_handler(), ffe);

                auto cfe = build_fe_dof_map(cdm.get_dof_handler());
                print_elem_global_dof(cdm.get_dof_handler(), cfe);
         */
        /* print_local_dofs(fdm); */
        /* print_owned_dofs(fdm); */

        /* print_partition_boundary_dofs(vdm); */
        /* print_ghost_dofs(vdm); */

        /* If you Orient then the same order is applied through all stencil
         * based on the orientation. Otherwise no order but just a normal stencil. */
        /* auto volume_stencil = vdm.build_stencil<VStencil, Orient>(); */

        auto volume_stencil = build_stencil<VStencil>(vdm.get_dof_handler());
        /* print_stencil(vdm, volume_stencil); */

        auto face_stencil = build_stencil<SStencil>(fdm.get_dof_handler());

        /* print_stencil(fdm, face_stencil); */

        /* using Pattern = SparsityPattern<DofLabel::lCorner, VStencil>;
        Pattern sp(volume_stencil); */
        VFDofMap<DHandler> map(dof_handler);
        SPattern sp(map, volume_stencil, face_stencil);

        // get the first owned dof for the first process.
        const Integer first_volume_dof = volume_stencil.get_dof_handler().get_owned_dof(0);
        printf(
            "First volume Dof: %li - %li\n", first_volume_dof, vdm.get_dof_handler().local_to_global(first_volume_dof));

        // TODO:: optimization idea. Iterate through the colidx instead of the stencil for better coalesing.
        // for each col idx (global dof) find the row pointer from the scan
        // in this way you would know directly the value index to be add and you need to figure out the labels.
        volume_stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            const Integer diag_dof = volume_stencil.get_value(stencil_index, SLabel::Diagonal);

            if (proc_num == 0 && diag_dof == first_volume_dof) {
                sp.set_value(diag_dof, diag_dof, 1);
            } else {
                const Integer right_dof = volume_stencil.get_value(stencil_index, SLabel::Right);
                sp.set_value(diag_dof, right_dof, 1);

                const Integer left_dof = volume_stencil.get_value(stencil_index, SLabel::Left);
                sp.set_value(diag_dof, left_dof, -1);

                const Integer up_dof = volume_stencil.get_value(stencil_index, SLabel::Up);
                sp.set_value(diag_dof, up_dof, 1);

                const Integer down_dof = volume_stencil.get_value(stencil_index, SLabel::Down);
                sp.set_value(diag_dof, down_dof, -1);
            }
        });

        face_stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            const Integer diag_dof = face_stencil.get_value(stencil_index, SLabel::Diagonal);

            /* printf("Diag_dof: %li global: %li\n", diag_dof, dof_handler.local_to_global(diag_dof)); */
            if (!face_stencil.get_dof_handler().is_boundary_dof(diag_dof)) {
                if (face_stencil.get_dof_handler().get_orientation(diag_dof) == DofOrient::xDir) {
                    sp.set_value(diag_dof, diag_dof, -4);

                    const Integer pr = face_stencil.get_value(stencil_index, SSXLabel::PRight);
                    sp.set_value(diag_dof, pr, -1);

                    const Integer pl = face_stencil.get_value(stencil_index, SSXLabel::PLeft);
                    sp.set_value(diag_dof, pl, 1);

                    const Integer vxr = face_stencil.get_value(stencil_index, SSXLabel::VXRight);
                    sp.set_value(diag_dof, vxr, 2);

                    const Integer vxl = face_stencil.get_value(stencil_index, SSXLabel::VXLeft);
                    sp.set_value(diag_dof, vxl, 2);

                    const Integer vxu = face_stencil.get_value(stencil_index, SSXLabel::VXUp);
                    if (vxu == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vxu, 1);
                    }
                    const Integer vxd = face_stencil.get_value(stencil_index, SSXLabel::VXDown);
                    if (vxd == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vxd, 1);
                    }

                    const Integer vyur = face_stencil.get_face_value(stencil_index, SSXFLabel::VYUpRight);
                    sp.set_value(diag_dof, vyur, 1);

                    const Integer vyul = face_stencil.get_face_value(stencil_index, SSXFLabel::VYUpLeft);
                    sp.set_value(diag_dof, vyul, -1);

                    const Integer vydr = face_stencil.get_face_value(stencil_index, SSXFLabel::VYDownRight);
                    sp.set_value(diag_dof, vydr, -1);

                    const Integer vydl = face_stencil.get_face_value(stencil_index, SSXFLabel::VYDownLeft);
                    sp.set_value(diag_dof, vydl, 1);

                } else if (face_stencil.get_dof_handler().get_orientation(diag_dof) == DofOrient::yDir) {
                    sp.set_value(diag_dof, diag_dof, -4);

                    const Integer pu = face_stencil.get_value(stencil_index, SSYLabel::PUp);
                    sp.set_value(diag_dof, pu, -1);

                    const Integer pd = face_stencil.get_value(stencil_index, SSYLabel::PDown);
                    sp.set_value(diag_dof, pd, 1);

                    const Integer vyr = face_stencil.get_value(stencil_index, SSYLabel::VYRight);
                    if (vyr == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vyr, 1);
                    }
                    const Integer vyl = face_stencil.get_value(stencil_index, SSYLabel::VYLeft);
                    if (vyl == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vyl, 1);
                    }

                    const Integer vyu = face_stencil.get_value(stencil_index, SSYLabel::VYUp);
                    sp.set_value(diag_dof, vyu, 2);

                    const Integer vyd = face_stencil.get_value(stencil_index, SSYLabel::VYDown);
                    sp.set_value(diag_dof, vyd, 2);

                    const Integer vxur = face_stencil.get_face_value(stencil_index, SSYFLabel::VXUpRight);
                    sp.set_value(diag_dof, vxur, 1);

                    const Integer vxul = face_stencil.get_face_value(stencil_index, SSYFLabel::VXUpLeft);
                    sp.set_value(diag_dof, vxul, -1);

                    const Integer vxdr = face_stencil.get_face_value(stencil_index, SSYFLabel::VXDownRight);
                    sp.set_value(diag_dof, vxdr, -1);

                    const Integer vxdl = face_stencil.get_face_value(stencil_index, SSYFLabel::VXDownLeft);
                    sp.set_value(diag_dof, vxdl, 1);
                }
            }
        });

        dof_handler.boundary_dof_iterate(
            MARS_LAMBDA(const Integer local_dof) { sp.set_value(local_dof, local_dof, 1); });

        print_sparsity_pattern(sp);

        vdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            vdm.get_data<IN>(i) = 0.0;
            vdm.get_data<OUT>(i) = 0.0;
        });

        vdm.gather_ghost_data<OUT>();
        /* scatter_add_ghost_data<VolumeDM, OUT>(vdm); */

        /* ViewVectorType<double> rhs("rhs", sp.get_num_rows()); */

        /* auto volume_handler = vdm.get_dof_handler();
        volume_handler.owned_dof_iterate(MARS_LAMBDA(const Integer local_dof) {
                double point[3];
                volume_handler.get_vertex_coordinates_from_local(local_dof, point);

                const Integer index = volume_handler.local_to_global(local_dof);
                rhs(index) = x;


                });
 */
/*
        auto face_handler = fdm.get_dof_handler();
        face_handler.owned_dof_iterate(MARS_LAMBDA(const Integer local_dof) {
            double point[2];
            face_handler.get_local_dof_coordinates(local_dof, point);

            double rval = 0;
            const Integer index = map.local_to_global(local_dof);
            rhs(index) = rval;
        }); */

        // print using the dof iterate
        /* vdm.get_dof_handler().dof_iterate(MARS_LAMBDA(const Integer local_dof) {
            const auto idata = vdm.get_dof_data<IN>(local_dof);
            const auto odata = vdm.get_dof_data<OUT>(local_dof);

            Dof d = vdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n",
                    i, idata, odata, d.get_gid(), d.get_proc());
        }); */

        // print using the index iterate
        /* vdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = vdm.get_dof_handler().get_local_dof(i);

            const auto idata = vdm.get_data<IN>(i);
            const auto odata = vdm.get_data<OUT>(i);

            Dof d = vdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n", i, idata, odata, d.get_gid(), d.get_proc());
        }); */

        double time = timer.seconds();
        std::cout << "Stag DM Setup took: " << time << " seconds." << std::endl;

#endif
    }

}  // namespace mars

#endif
