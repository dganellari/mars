#ifndef MARS_CV_STOKES_
#define MARS_CV_STOKES_

#include "mars.hpp"
#include "mars_staggered_utils.hpp"
#include "mars_stokes_common.hpp"

namespace mars {

    using namespace stag;

    template <typename S, typename SP>
    void assemble_oriented_face(S face_stencil, SP sp) {
        auto fv_dof_handler = sp.get_dof_handler();

        face_stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            const Integer diag_dof = face_stencil.get_value(stencil_index, SLabel::Diagonal);

            if (!fv_dof_handler.is_boundary_dof(diag_dof)) {
                sp.set_value(diag_dof, diag_dof, -4);

                const Integer pu = face_stencil.get_value(stencil_index, SSOLabel::VolumeUp);
                sp.set_value(diag_dof, pu, -1);

                const Integer pd = face_stencil.get_value(stencil_index, SSOLabel::VolumeDown);
                sp.set_value(diag_dof, pd, 1);

                const Integer vyr = face_stencil.get_value(stencil_index, SSOLabel::FaceRight);
                if (vyr == -1) {
                    sp.set_value(diag_dof, diag_dof, -3);
                } else {
                    sp.set_value(diag_dof, vyr, 1);
                }
                const Integer vyl = face_stencil.get_value(stencil_index, SSOLabel::FaceLeft);
                if (vyl == -1) {
                    sp.set_value(diag_dof, diag_dof, -3);
                } else {
                    sp.set_value(diag_dof, vyl, 1);
                }

                const Integer vyu = face_stencil.get_value(stencil_index, SSOLabel::FaceUp);
                sp.set_value(diag_dof, vyu, 2);

                const Integer vyd = face_stencil.get_value(stencil_index, SSOLabel::FaceDown);
                sp.set_value(diag_dof, vyd, 2);

                const Integer vxur = face_stencil.get_value(stencil_index, SSOLabel::FaceUpRight);
                sp.set_value(diag_dof, vxur, 1);

                const Integer vxul = face_stencil.get_value(stencil_index, SSOLabel::FaceUpLeft);
                sp.set_value(diag_dof, vxul, -1);

                const Integer vxdr = face_stencil.get_value(stencil_index, SSOLabel::FaceDownRight);
                sp.set_value(diag_dof, vxdr, -1);

                const Integer vxdl = face_stencil.get_value(stencil_index, SSOLabel::FaceDownLeft);
                sp.set_value(diag_dof, vxdl, 1);
            }
        });
    }

    template <typename S, typename SP>
    void assemble_face(S face_stencil, SP sp) {
        auto fv_dof_handler = sp.get_dof_handler();

        face_stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            const Integer diag_dof = face_stencil.get_value(stencil_index, SLabel::Diagonal);

            if (!fv_dof_handler.is_boundary_dof(diag_dof)) {
                if (fv_dof_handler.get_orientation(diag_dof) == DofOrient::xDir) {
                    sp.set_value(diag_dof, diag_dof, -4);

                    const Integer pr = face_stencil.get_value(stencil_index, SSXLabel::VolumeXRight);
                    sp.set_value(diag_dof, pr, -1);

                    const Integer pl = face_stencil.get_value(stencil_index, SSXLabel::VolumeXLeft);
                    sp.set_value(diag_dof, pl, 1);

                    const Integer vxr = face_stencil.get_value(stencil_index, SSXLabel::FaceXRight);
                    sp.set_value(diag_dof, vxr, 2);

                    const Integer vxl = face_stencil.get_value(stencil_index, SSXLabel::FaceXLeft);
                    sp.set_value(diag_dof, vxl, 2);

                    const Integer vxu = face_stencil.get_value(stencil_index, SSXLabel::FaceXUp);
                    if (vxu == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vxu, 1);
                    }
                    const Integer vxd = face_stencil.get_value(stencil_index, SSXLabel::FaceXDown);
                    if (vxd == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vxd, 1);
                    }

                    const Integer vyur = face_stencil.get_value(stencil_index, SSXLabel::FaceYUpRight);
                    sp.set_value(diag_dof, vyur, 1);

                    const Integer vyul = face_stencil.get_value(stencil_index, SSXLabel::FaceYUpLeft);
                    sp.set_value(diag_dof, vyul, -1);

                    const Integer vydr = face_stencil.get_value(stencil_index, SSXLabel::FaceYDownRight);
                    sp.set_value(diag_dof, vydr, -1);

                    const Integer vydl = face_stencil.get_value(stencil_index, SSXLabel::FaceYDownLeft);
                    sp.set_value(diag_dof, vydl, 1);

                } else if (fv_dof_handler.get_orientation(diag_dof) == DofOrient::yDir) {
                    sp.set_value(diag_dof, diag_dof, -4);

                    const Integer pu = face_stencil.get_value(stencil_index, SSYLabel::VolumeYUp);
                    sp.set_value(diag_dof, pu, -1);

                    const Integer pd = face_stencil.get_value(stencil_index, SSYLabel::VolumeYDown);
                    sp.set_value(diag_dof, pd, 1);

                    const Integer vyr = face_stencil.get_value(stencil_index, SSYLabel::FaceYRight);
                    if (vyr == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vyr, 1);
                    }
                    const Integer vyl = face_stencil.get_value(stencil_index, SSYLabel::FaceYLeft);
                    if (vyl == -1) {
                        sp.set_value(diag_dof, diag_dof, -3);
                    } else {
                        sp.set_value(diag_dof, vyl, 1);
                    }

                    const Integer vyu = face_stencil.get_value(stencil_index, SSYLabel::FaceYUp);
                    sp.set_value(diag_dof, vyu, 2);

                    const Integer vyd = face_stencil.get_value(stencil_index, SSYLabel::FaceYDown);
                    sp.set_value(diag_dof, vyd, 2);

                    const Integer vxur = face_stencil.get_value(stencil_index, SSYLabel::FaceXUpRight);
                    sp.set_value(diag_dof, vxur, 1);

                    const Integer vxul = face_stencil.get_value(stencil_index, SSYLabel::FaceXUpLeft);
                    sp.set_value(diag_dof, vxul, -1);

                    const Integer vxdr = face_stencil.get_value(stencil_index, SSYLabel::FaceXDownRight);
                    sp.set_value(diag_dof, vxdr, -1);

                    const Integer vxdl = face_stencil.get_value(stencil_index, SSYLabel::FaceXDownLeft);
                    sp.set_value(diag_dof, vxdl, 1);
                }
            }
        });
    }

    template <class BC, class RHS, class AnalyticalFun>
    void staggered_constant_viscosty_stokes_2D(const int level) {
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

        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        DHandler dof_handler(&mesh, context);
        dof_handler.enumerate_dofs();

        FVDH fv_dof_handler(dof_handler);

        auto volume_stencil = build_stencil<VStencil>(fv_dof_handler);
        /* print_stencil(fv_dof_handler, volume_stencil); */

        /* auto face_stencil = build_stencil<SStencil, Orient>(fv_dof_handler); */
        auto face_stencil = build_stencil<SStencil>(fv_dof_handler);
        /* print_stencil(fv_dof_handler, face_stencil); */

        SPattern sp(fv_dof_handler);
        sp.build_pattern(volume_stencil, face_stencil);

        assemble_volume(volume_stencil, sp, proc_num);
        assemble_face(face_stencil, sp);
        /* assemble_oriented_face(face_stencil, sp); */

        fv_dof_handler.boundary_dof_iterate(
            MARS_LAMBDA(const Integer local_dof) { sp.set_value(local_dof, local_dof, 1); });

        /* print_sparsity_pattern(sp); */
        sp.write("Spattern");

        CornerDM cdm(dof_handler);
        set_data_in_circle(cdm, 0, 1);
        cdm.gather_ghost_data<IN>();

        cdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = cdm.get_dof_handler().get_local_dof(i);

            const auto idata = cdm.get_data<IN>(i);
            Dof d = cdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, global: %li, rank: %i\n", i, idata, d.get_gid(), d.get_proc());
        });

        auto rhs = assemble_rhs(fv_dof_handler, cdm);

        /* ********************************gather scatter ghost data**************************************** */

        /* FaceVolumeDM fvdm(fv_dof_handler);
        fvdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            fvdm.get_data<IN>(i) = 3.0;
            fvdm.get_data<OUT>(i) = proc_num;
        });

        fvdm.gather_ghost_data<OUT>();
        scatter_add_ghost_data<FaceVolumeDM, OUT>(fvdm); */

        // print using the index iterate
        /* fvdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = fvdm.get_dof_handler().get_local_dof(i);

            const auto idata = fvdm.get_data<IN>(i);
            const auto odata = fvdm.get_data<OUT>(i);

            Dof d = fvdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, v: %lf, global: %li, rank: %i\n", i, idata, odata, d.get_gid(), d.get_proc());
        }); */

        /* ********************************end ghost data gather scatter************************************ */

        double time = timer.seconds();
        std::cout << "Stag DM Setup took: " << time << " seconds." << std::endl;

#endif
    }

}  // namespace mars

#endif
