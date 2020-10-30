#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP

#include "mars_base.hpp"

namespace mars {

    /* template <Integer Dim, Integer degree, Integer Type = Star> */
    template <Integer Dim, Integer Degree, Integer Width, bool Fill = false>
    class Stencil {
    public:
        static constexpr Integer Length = 2 * Dim * Width + 1;
        Integer stencil_size;

        MARS_INLINE_FUNCTION
        Stencil() = default;

        MARS_INLINE_FUNCTION
        ViewMatrixTypeRC<Integer, Length> get_stencil() const { return stencil; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_width() const { return Width; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_length() const { return Length; }

        MARS_INLINE_FUNCTION
        Integer get_value(const Integer row, const Integer col) const { return stencil(row, col); }

        MARS_INLINE_FUNCTION
        Integer set_value(const Integer row, const Integer col, const Integer value) const { return stencil(row, col) = value; }

        template <typename F>
        void iterate(F f) const {
            Kokkos::parallel_for("stencil_dof_iter", stencil_size, f);
        }

        void reserve_stencil(const Integer size) {
            stencil_size = size;
            stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size);
        }

    private:
        ViewMatrixTypeRC<Integer, Length> stencil;
    };

    /* template <>
    class Stencil<2, 2, 1, true> {
    public:
    };
 */
    // stokes fd stencil used ex: variable viscosity
    template <>
    class Stencil<2, 2, 2, true> {
    public:
        /* constexpr Integer length = 4^Dim + 1; */
        static constexpr Integer Dim = 2;
        static constexpr Integer Width = 2;
        static constexpr Integer Length = 2 * Dim * Width + 1;
        static constexpr Integer Face_Length = 2 * Width;
        static constexpr Integer Corner_Length = 2 * Width;

        /* constexpr Integer Length = core_Length + face_Length + corner_Length; */

        ViewMatrixTypeRC<Integer, Length> stencil;
        ViewMatrixTypeRC<Integer, Face_Length> face_extension;
        ViewMatrixTypeRC<Integer, Corner_Length> corner_extension;

        void reserve_stencil(const Integer size) { stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size); }

        void reserve_stencil_face_extension(const Integer size) {
            face_extension = ViewMatrixTypeRC<Integer, Face_Length>("stencil_face_ext", size);
        }

        void reserve_stencil_corner_extension(const Integer size) {
            corner_extension = ViewMatrixTypeRC<Integer, Corner_Length>("stencil_corner_ext", size);
        }
    };

    // TODO: specialize without the corners just facextension
}  // namespace mars

#endif  // mars_stencil
