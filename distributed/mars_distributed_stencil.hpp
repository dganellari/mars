#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP


#include "mars_base.hpp"

namespace mars {

    /* template <Integer Dim, Integer degree, Integer Type = Star> */
    template <Integer Dim, Integer Degree, Integer Width, bool Fill = false>
    class Stencil {
    public:
        ViewMatrixTypeRC<Integer> stencil;

        constexpr Integer Length = 2 * Dim * Width + 1;

        void reserve_stencil(const Integer size) {
            stencil = ViewMatrixTypeRC<Integer, Length>(size);
        }
    };

    /* template <>
    class Stencil<2, 2, 1, true> {
    public:
    };
 */
    //stokes fd stencil used ex: variable viscosity
    template<>
    class Stencil<2, 2, 2, true>
    {
        public:
        /* constexpr Integer length = 4^Dim + 1; */
        constexpr Integer Length = 2 * Dim * Width + 1;
        constexpr Integer Face_Length = 2 * Width;
        constexpr Integer Corner_Length = 2 * Width;

        /* constexpr Integer Length = core_Length + face_Length + corner_Length; */

        ViewMatrixTypeRC<Integer, Length> stencil;
        ViewMatrixTypeRC<Integer, Face_Length> face_extension;
        ViewMatrixTypeRC<Integer, Corner_Length> corner_extension;

        void reserve_stencil(const Integer size) {
            stencil = ViewMatrixTypeRC<Integer,Length>(size);
            face_extension = ViewMatrixTypeRC<Integer,Face_Length>(size);
            corner_extension = ViewMatrixTypeRC<Integer,Corner_Length>(size);
        }
    };


    //TODO: specialize without the corners just facextension
}

#endif //mars_stencil
