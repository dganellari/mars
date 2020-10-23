#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP


#include "mars_base.hpp"

namespace mars {

    /* template <Integer Dim, Integer degree, Integer Type = Star> */
    template <Integer Dim, Integer degree>
    class Stencil
    {
        public:
            ViewMatrixTypeRC<Integer> x_stencil;
            ViewMatrixTypeRC<Integer> y_stencil;
            ViewMatrixTypeRC<Integer> z_stencil;
    };

    template<>
    class Stencil<2, 2>
    {
        public:
        /* constexpr Integer length = 4^Dim + 1; */
        constexpr Integer Length = 17;
        ViewMatrixTypeRC<Integer, Length> x_stencil;
        ViewMatrixTypeRC<Integer, Length> y_stencil;

        void reserve_stencil(const Integer size) {
            x_stencil = ViewMatrixType<Integer,Length>(size);
            y_stencil = ViewMatrixType<Integer,Length>(size);
        }

        void build_stencil() {

        }
}

#endif //mars_stencil
