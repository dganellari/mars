#ifndef MARS_DMINTERPOLATE_HPP
#define MARS_DMINTERPOLATE_HPP

#include "mars_base.hpp"

#include "mars_globals.hpp"

namespace mars {

    template <class DM>
    class DMInterpolate {
    public:
        DMInterpolate(DM &dm) : dm_(dm) {}

        using Elem = typename DM::simplex_type;

        template <typename F, typename T>
        void apply(F fun, ViewVectorType<T> &v) {
            auto dm = dm_;

            dm.owned_dof_iterate(MARS_LAMBDA(const Integer index) {
                double p[2];

                const Integer sfc = dm.get_global_dof_enum().get_view_elements()(index);
                dm.template get_dof_coordinates_from_sfc<Elem::ElemType>(sfc, p);

                v(index) = fun(p);
            });
        }

    private:
        DM &dm_;
    };

}  // namespace mars

#endif  // MARS_INTERPOLATE_HPP
