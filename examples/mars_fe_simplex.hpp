#ifndef MARS_FE_SIMPLEX_HPP
#define MARS_FE_SIMPLEX_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"

namespace mars {
    template <int Dim>
    class FESimplex {
    public:
        MARS_INLINE_FUNCTION static Real fun(const int i, const Real *p) {
            Real ret = 0.0;

            if (i == 0) {
                ret = 1.0;
                for (int d = 0; d < Dim; ++d) {
                    ret -= p[d];
                }
            } else {
                ret = p[i - 1];
            }

            return ret;
        }
    };
}  // namespace mars

#endif  // MARS_FE_SIMPLEX_HPP
