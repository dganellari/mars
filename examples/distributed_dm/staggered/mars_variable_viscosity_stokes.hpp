#ifndef MARS_VV_STOKES_
#define MARS_VV_STOKES_

#include "mars.hpp"
#include "mars_staggered_utils.hpp"
#include "mars_stokes_common.hpp"

namespace mars {

    using namespace stag;

    template <class BC, class RHS, class AnalyticalFun>
    void staggered_variable_viscosty_stokes_2D(const int level) {
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


        CornerDM cdm(dof_handler);
        set_data_in_circle(cdm, 0, 1);
        cdm.gather_ghost_data<IN>();


        CornerVolumeDM cvdm(dof_handler);
        set_data_in_circle(cvdm, 1, 10);
        cvdm.gather_ghost_data<IN>();

        /* cdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = cdm.get_dof_handler().get_local_dof(i);

            const auto idata = cdm.get_data<IN>(i);
            Dof d = cdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, global: %li, rank: %i\n", i, idata, d.get_gid(), d.get_proc());
        });
 */

        cvdm.get_dof_handler().iterate(MARS_LAMBDA(const Integer i) {
            const Integer local_dof = cvdm.get_dof_handler().get_local_dof(i);

            const auto idata = cvdm.get_data<IN>(i);
            Dof d = cvdm.get_dof_handler().local_to_global_dof(local_dof);

            printf("lid: %li, u: %lf, global: %li, rank: %i\n", i, idata, d.get_gid(), d.get_proc());
        });

        /* assemble_volume(volume_stencil, sp, proc_num); */
        /* assemble_face(face_stencil, sp); */

        fv_dof_handler.boundary_dof_iterate(
            MARS_LAMBDA(const Integer local_dof) { sp.set_value(local_dof, local_dof, 1); });

        /* print_sparsity_pattern(sp); */
        sp.write("Spattern");

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
