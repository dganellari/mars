#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP

#include "mars_base.hpp"

namespace mars {

    /* enum StencilType : int {
        // 1D
        Classic = 0,
        // 2D
        FaceComplete = 1,
        CornerComplete = 2,
        InvalidType = -1
    }; */

    template <Integer Dim, Integer Width = 1>
    class Stencil {
    public:
        static constexpr Integer Length = 2 * Dim * Width + 1;

        MARS_INLINE_FUNCTION
        Stencil() = default;

        MARS_INLINE_FUNCTION
        Stencil(const Integer size) { stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size); }

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

        Integer get_stencil_size() const { return stencil.extent(0); }

        template <typename F>
        void dof_iterate(F f) const {
            auto st = get_stencil();
            Kokkos::parallel_for(
                "stencil_dof_iter", get_stencil_size(), MARS_LAMBDA(const Integer stencil_index) {
                    for (int i = 0; i < get_length(); i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = st(stencil_index, i);
                        f(stencil_index, local_dof);
                    }
                });
        }

        template <typename F>
        void iterate(F f) const {
            Kokkos::parallel_for("stencil_dof_iter", get_stencil_size(), f);
        }

        virtual void reserve_stencil(const Integer size) {
            stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size);
        }
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

    template <Integer Dim>
    class StokesStencil : public Stencil<Dim, 2> {
    public:
        static constexpr Integer Face_Length = 2 * Dim;

        using SuperStencil = Stencil<Dim, 2>;

        MARS_INLINE_FUNCTION
        StokesStencil() = default;

        MARS_INLINE_FUNCTION
        StokesStencil(const Integer size) : SuperStencil(size) {
            face_extension = ViewMatrixTypeRC<Integer, Face_Length>("face_ext", size);
        }

        virtual void reserve_stencil(const Integer size) override {
            SuperStencil::reserve_stencil(size);
            face_extension = ViewMatrixTypeRC<Integer, Face_Length>("stencil_face_ext", size);
        }

        MARS_INLINE_FUNCTION
        Integer get_face_value(const Integer row, const Integer col) const { return face_extension(row, col); }

        MARS_INLINE_FUNCTION
        Integer set_face_value(const Integer row, const Integer col, const Integer value) const {
            return face_extension(row, col) = value;
        }

        template <typename F>
        void dof_iterate(F f) const {
            auto st = SuperStencil::get_stencil();
            auto fe = get_face_stencil();
            auto length = SuperStencil::get_length();
            Kokkos::parallel_for(
                "stencil_dof_iter", SuperStencil::get_stencil_size(), MARS_LAMBDA(const Integer stencil_index) {
                    for (int i = 0; i < length; i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = st(stencil_index, i);
                        f(stencil_index, local_dof);
                    }

                    for (int i = 0; i < get_face_length(); i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = fe(stencil_index, i);
                        f(stencil_index, local_dof);
                    }
                });
        }

        template <typename DM>
        MARS_INLINE_FUNCTION void build_diagonal_stencil(const DM &dm,
                                                         const Integer localid,
                                                         const Integer stencil_index,
                                                         const Integer Orientation = 1) const {
            const Integer sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

            for (int corner = power_of_2(DM::ManifoldDim) - 1; corner != -1; --corner) {
                Octant o = oc.sfc_corner_nbh<DM::simplex_type::ElemType>(corner);
                const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                set_face_value(stencil_index, corner, nbh_id);
            }

            /*FIXME for 3D.*/
            if (Orientation == 0) {
                const Integer tmp = face_extension(stencil_index, 1);
                face_extension(stencil_index, 1) = face_extension(stencil_index, 2);
                face_extension(stencil_index, 2) = tmp;
            }
        }

        template <typename DM>
        MARS_INLINE_FUNCTION void build_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index,
                                                const Integer Orientation = 1) const {
            SuperStencil::build_stencil(dm, localid, stencil_index, Orientation);

            build_diagonal_stencil(dm, localid, stencil_index, Orientation);
        }

        MARS_INLINE_FUNCTION
        ViewMatrixTypeRC<Integer, Face_Length> get_face_stencil() const { return face_extension; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_face_length() const { return Face_Length; }


    private:
        ViewMatrixTypeRC<Integer, Face_Length> face_extension;
    };

}  // namespace mars

#endif  // mars_stencil
