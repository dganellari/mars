#ifndef MARS_CV_STOKES_
#define MARS_CV_STOKES_

#include "mars.hpp"
#include "mars_staggered_utils.hpp"

namespace mars {

    using namespace stag;

    template <typename S, typename SP>
    void assemble_volume(S volume_stencil, SP sp, const Integer proc_num) {
        // TODO:: optimization idea. Iterate through the colidx instead of the stencil for better coalesing.
        // for each col idx (global dof) find the row pointer from the scan
        // in this way you would know directly the value index to be add and you need to figure out the labels.
        volume_stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            const Integer diag_dof = volume_stencil.get_value(stencil_index, SLabel::Diagonal);

            /* first volume dof of the first process */
            if (proc_num == 0 && stencil_index == 0) {
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
    }

    template <typename S, typename SP>
    void assemble_face(S face_stencil, SP sp) {
        auto fv_dof_handler = sp.get_dof_handler();

        face_stencil.iterate(MARS_LAMBDA(const Integer stencil_index) {
            const Integer diag_dof = face_stencil.get_value(stencil_index, SLabel::Diagonal);

            if (!fv_dof_handler.is_boundary_dof(diag_dof)) {
                if (fv_dof_handler.get_orientation(diag_dof) == DofOrient::xDir) {
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

                    const Integer vyur = face_stencil.get_value(stencil_index, SSXLabel::VYUpRight);
                    sp.set_value(diag_dof, vyur, 1);

                    const Integer vyul = face_stencil.get_value(stencil_index, SSXLabel::VYUpLeft);
                    sp.set_value(diag_dof, vyul, -1);

                    const Integer vydr = face_stencil.get_value(stencil_index, SSXLabel::VYDownRight);
                    sp.set_value(diag_dof, vydr, -1);

                    const Integer vydl = face_stencil.get_value(stencil_index, SSXLabel::VYDownLeft);
                    sp.set_value(diag_dof, vydl, 1);

                } else if (fv_dof_handler.get_orientation(diag_dof) == DofOrient::yDir) {
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

                    const Integer vxur = face_stencil.get_value(stencil_index, SSYLabel::VXUpRight);
                    sp.set_value(diag_dof, vxur, 1);

                    const Integer vxul = face_stencil.get_value(stencil_index, SSYLabel::VXUpLeft);
                    sp.set_value(diag_dof, vxul, -1);

                    const Integer vxdr = face_stencil.get_value(stencil_index, SSYLabel::VXDownRight);
                    sp.set_value(diag_dof, vxdr, -1);

                    const Integer vxdl = face_stencil.get_value(stencil_index, SSYLabel::VXDownLeft);
                    sp.set_value(diag_dof, vxdl, 1);
                }
            }
        });
    }

    template <typename H>
    ViewVectorType<double> assemble_rhs(const H &fv_dof_handler) {
        ViewVectorType<double> rhs("rhs", fv_dof_handler.get_owned_dof_size());

        fv_dof_handler.template owned_dof_iterate<DofLabel::lFace>(MARS_LAMBDA(const Integer local_dof) {
            double point[2];
            fv_dof_handler.get_local_dof_coordinates(local_dof, point);

            /*auto global = fv_dof_handler.get_dof_handler().local_to_global(local_dof);
            auto o = fv_dof_handler.get_octant_from_local(local_dof);
            printf("local_dof: %li, x, y, z %li, %li, %li, global: %li\n", local_dof, o.x, o.y, o.z, global);*/

            if (!fv_dof_handler.is_boundary_dof(local_dof)) {
                double rval = 0;
                if (fv_dof_handler.get_orientation(local_dof) == DofOrient::yDir) {
                    rval = 1;
                }
                const Integer index = fv_dof_handler.local_to_owned_index(local_dof);
                rhs(index) = rval;
            }
        });

        return rhs;
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

        auto face_stencil = build_stencil<SStencil>(fv_dof_handler);
        /* print_stencil(fv_dof_handler, face_stencil); */

        SPattern sp(fv_dof_handler);
        sp.build_pattern(volume_stencil, face_stencil);

        assemble_volume(volume_stencil, sp, proc_num);
        assemble_face(face_stencil, sp);

        fv_dof_handler.boundary_dof_iterate(
            MARS_LAMBDA(const Integer local_dof) { sp.set_value(local_dof, local_dof, 1); });

        print_sparsity_pattern(sp);
        sp.write("Spattern");

        auto rhs = assemble_rhs(fv_dof_handler);

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
