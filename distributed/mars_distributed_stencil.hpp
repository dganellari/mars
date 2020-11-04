#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP

#include "mars_base.hpp"

namespace mars {

    enum StencilType : int {
        // 1D
        Classic = 0,
        // 2D
        FaceComplete = 1,
        CornerComplete = 2,
        InvalidType = -1
    };

    /* template <Integer Dim, Integer Degree, Integer Width, Integer Type = StencilType::Classic> */
    template <Integer Dim, Integer Degree, Integer Width>
    class Stencil {
    public:
        static constexpr Integer Length = 2 * Dim * Width + 1;

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
        Integer set_value(const Integer row, const Integer col, const Integer value) const {
            return stencil(row, col) = value;
        }

        template <typename F>
        void iterate(F f) const {
            Kokkos::parallel_for("stencil_dof_iter", get_stencil_size(), f);
        }

        virtual void reserve_stencil(const Integer size) { stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size); }

        Integer get_stencil_size() const { return stencil.extent(0); }

        template <typename DM>
        MARS_INLINE_FUNCTION void build_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index,
                                                const Integer Orientation = 1) const {
            set_value(stencil_index, 0, localid);
            const Integer sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

            Integer face_nr;
            for (int dir = 0; dir < 2; ++dir) {
                for (int side = 0; side < 2; ++side) {
                    if (side == 0)
                        face_nr = 2 * dir + 1;
                    else
                        face_nr = 2 * dir;

                    // this gives the index for different face orientation. Corner and volume have no
                    // extra orientation and the default  is 1. Orientation=0 -> x orientation).
                    const Integer dir_dim = !(Orientation ^ dir);
                    Integer index = 2 * dir_dim + side + 1;

                    Octant o = oc;
                    for (int w = 0; w < get_width(); ++w) {
                        /* const Integer nbh_sfc = dm.get_sfc_face_nbh(oc, face_nr); */

                        o = o.sfc_face_nbh<DM::simplex_type::ElemType>(face_nr);
                        const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                        Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                        Integer offset = w * 2 * Dim;
                        set_value(stencil_index, index + offset, nbh_id);
                    }
                }
            }
        }

    private:
        ViewMatrixTypeRC<Integer, Length> stencil;
    };

    /* template <Integer Dim, Integer degree, Integer Type = Star> */
    template <Integer Dim, Integer Degree, Integer Width>
    class StokesStencil: public Stencil<Dim, Degree, 2> {
    public:
        static constexpr Integer Face_Length = 2 * Width;

        using SuperStencil = Stencil<Dim, Degree, 2>;

        virtual void reserve_stencil_face_extension(const Integer size) override {
            SuperStencil::reserve_stencil(size);
            face_extension = ViewMatrixTypeRC<Integer, Face_Length>("stencil_face_ext", size);
        }

        template <typename DM>
        MARS_INLINE_FUNCTION void build_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index,
                                                const Integer Orientation = 1) const {
            const Integer sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

            for (int corner = 0; corner < power_of_2(DM::Mesh::ManifoldDim); ++corner) {
                // this gives the index for different face orientation. Corner and volume have no
                // extra orientation and the default  is 1. Orientation=0 -> x orientation).
                const Integer dir_dim = !(Orientation ^ dir);
                Integer index = 2 * dir_dim + side + 1;

                Octant o = oc.sfc_corner_nbh<DM::simplex_type::ElemType>(corner);
                const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                set_value(stencil_index, index, nbh_id);
            }
        }

    private:
        ViewMatrixTypeRC<Integer, Face_Length> face_extension;
    };

    /* template <>
    class Stencil<2, 2, 1, true> {
    public:
    };
 */
    // stokes fd stencil used ex: variable viscosity
    /* template <>
    class Stencil<2, 2, 2, StencilType::FaceComplete> {
    public:
        [>constexpr Integer length = 4^Dim + 1;<]
        static constexpr Integer Dim = 2;
        static constexpr Integer Width = 2;
        static constexpr Integer Length = 2 * Dim * Width + 1;
        static constexpr Integer Face_Length = 2 * Width;
        static constexpr Integer Corner_Length = 2 * Width;

        [>constexpr Integer Length = core_Length + face_Length + corner_Length;<]

        void reserve_stencil(const Integer size) { stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size); }

        void reserve_stencil_face_extension(const Integer size) {
            face_extension = ViewMatrixTypeRC<Integer, Face_Length>("stencil_face_ext", size);
        }

        void reserve_stencil_corner_extension(const Integer size) {
            corner_extension = ViewMatrixTypeRC<Integer, Corner_Length>("stencil_corner_ext", size);
        }


    private:
        ViewMatrixTypeRC<Integer, Length> stencil;
        ViewMatrixTypeRC<Integer, Face_Length> face_extension;
        ViewMatrixTypeRC<Integer, Corner_Length> corner_extension;

    }; */

    // TODO: specialize without the corners just facextension
}  // namespace mars

#endif  // mars_stencil
