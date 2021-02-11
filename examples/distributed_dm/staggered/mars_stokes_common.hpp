#ifndef MARS_SC_STOKES_
#define MARS_SC_STOKES_

#include "mars.hpp"
#include "mars_staggered_utils.hpp"

namespace mars {

    using namespace stag;

    template <typename H, typename DM>
    ViewVectorType<double> assemble_rhs(const H &fv_dof_handler, const DM &dm) {
        ViewVectorType<double> rhs("rhs", fv_dof_handler.get_owned_dof_size());

        fv_dof_handler.template owned_dof_iterate<DofLabel::lFace>(MARS_LAMBDA(const Integer local_dof) {
            /* double point[2];
            fv_dof_handler.get_local_dof_coordinates(local_dof, point);
            auto global = fv_dof_handler.get_dof_handler().local_to_global(local_dof);
            auto o = fv_dof_handler.get_octant_from_local(local_dof);
            printf("local_dof: %li, x, y, z %li, %li, %li, global: %li\n", local_dof, o.x, o.y, o.z, global); */

            auto o = fv_dof_handler.get_octant_from_local(local_dof);

            if (!fv_dof_handler.is_boundary_dof(local_dof)) {
                double rval = 0;
                if (fv_dof_handler.get_orientation(local_dof) == DofOrient::yDir) {
                    auto o1 = Octant(o.x - 1, o.y, o.z);
                    auto local_dof1 = fv_dof_handler.get_local_from_octant(o1);

                    auto o2 = Octant(o.x + 1, o.y, o.z);
                    auto local_dof2 = fv_dof_handler.get_local_from_octant(o2);

                    auto avg =
                        (dm.template get_dof_data<IN>(local_dof1) + dm.template get_dof_data<IN>(local_dof2)) / 2;
                    rval = -10 * avg;
                }

                const Integer index = fv_dof_handler.local_to_owned_index(local_dof);
                rhs(index) = rval;

                auto global = fv_dof_handler.get_dof_handler().local_to_global(local_dof);
                printf("local_dof: %li, %li, global: %li - rval: %lf\n", local_dof, global, rval);
            }
        });

        return rhs;
    }

    template <typename DM>
    void set_data_in_circle(const DM &dm, const double a, double b) {
        dm.template owned_data_iterate<IN>(MARS_LAMBDA(const Integer local_dof, double &data) {
            auto octant = dm.get_dof_handler().get_octant_from_local(local_dof);
            double xMax = dm.get_dof_handler().get_XMax();
            double yMax = dm.get_dof_handler().get_YMax();

            double x = octant.x / xMax;
            double y = octant.y / yMax;

            auto r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);

            if (r2 > 0.25 * 0.25)
                data = a;
            else
                data = b;
        });
    }

}  // namespace mars

#endif
