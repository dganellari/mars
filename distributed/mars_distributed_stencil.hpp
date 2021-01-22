#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP

#include "mars_base.hpp"

namespace mars {

    /* Face numbering on the stencil => ordering in the stencil stencil[1,0,3,2]
           ----3----
           |       |
           0   x   1
           |       |
           ----2---- */
    // building the stencil is the responsibility of the specialized DofHandler.
    template <typename ST, bool Orient = false, typename DofHandler>
    ST build_stencil(const DofHandler &handler) {
        ST stencil(handler);

        handler.owned_iterate(MARS_LAMBDA(const Integer i) {
            const Integer localid = handler.get_owned_dof(i);
            stencil.template build_stencil<Orient>(localid, i);
        });
        return stencil;
    }


    template <typename DofHandler, Integer Width = 1>
    class Stencil {
    public:
        static constexpr Integer Dim = DofHandler::ManifoldDim;

        static constexpr Integer Length = 2 * Dim * Width + 1;

        MARS_INLINE_FUNCTION
        Stencil() = default;

        MARS_INLINE_FUNCTION
        Stencil(DofHandler handler) : dof_handler(handler) {
            const Integer size = get_dof_handler().get_owned_dof_size();
            stencil = ViewMatrixTypeRC<Integer, Length>("stencil", size);
        }

        MARS_INLINE_FUNCTION
        ViewMatrixTypeRC<Integer, Length> get_stencil() const { return stencil; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_width() const { return Width; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_length() const { return Length; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_face_length() const { return 0; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_corner_length() const { return 0; }

        MARS_INLINE_FUNCTION
        Integer get_value(const Integer row, const Integer col) const { return stencil(row, col); }

        MARS_INLINE_FUNCTION
        Integer set_value(const Integer row, const Integer col, const Integer value) const {
            return stencil(row, col) = value;
        }

        const Integer get_stencil_size() const { return stencil.extent(0); }

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

        template <bool Orient = false>
        MARS_INLINE_FUNCTION void build_stencil(const Integer localid, const Integer stencil_index) const {
            const Integer orientation = get_dof_handler().get_orientation(localid);
            generate_stencil<Orient>(get_dof_handler(), localid, stencil_index, orientation);
        }
/*
        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index) const {
            auto handler = dm.get_dof_handler();
            build_stencil(handler, localid, stencil_index);
        } */

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void generate_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index,
                                                Integer orientation) const {
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

                    // if the user chooses no orientation then orientation 1
                    // means discard the orientation values and do
                    // not use them. !(x ^ y) leaves y the same when x == 1.
                    if (!Orient) orientation = DofOrient::yDir;
                    // this gives the index for different face orientation. Corner and volume have no
                    // extra orientation and the default  is 1. Orientation=0 -> x orientation).
                    const Integer dir_dim = !(orientation ^ dir);
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

        MARS_INLINE_FUNCTION
        const DofHandler &get_dof_handler() const { return dof_handler; }

    private:
        DofHandler dof_handler;
        ViewMatrixTypeRC<Integer, Length> stencil;
    };

    template <typename DofHandler>
    class StokesStencil : public Stencil<DofHandler, 2> {
    public:

        static constexpr Integer Dim = DofHandler::ManifoldDim;

        static constexpr Integer Face_Length = 2 * Dim;

        using SuperStencil = Stencil<DofHandler, 2>;

        MARS_INLINE_FUNCTION
        StokesStencil() = default;

        MARS_INLINE_FUNCTION
        StokesStencil(DofHandler handler) : SuperStencil(handler) {
            const Integer size = SuperStencil::get_dof_handler().get_owned_dof_size();
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

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_diagonal_stencil(const DM &dm,
                                                         const Integer localid,
                                                         const Integer stencil_index,
                                                         const Integer orientation) const {
            const Integer sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

            for (int corner = power_of_2(DM::ManifoldDim) - 1; corner != -1; --corner) {
                Octant o = oc.sfc_corner_nbh<DM::simplex_type::ElemType>(corner);
                const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                set_face_value(stencil_index, corner, nbh_id);
            }

            /*FIXME for 3D.*/
            if (Orient && orientation == DofOrient::xDir) {
                const Integer tmp = face_extension(stencil_index, 1);
                face_extension(stencil_index, 1) = face_extension(stencil_index, 2);
                face_extension(stencil_index, 2) = tmp;
            }
        }

        template <bool Orient = false>
        MARS_INLINE_FUNCTION void build_stencil(const Integer localid, const Integer stencil_index) const {
            const Integer orientation = SuperStencil::get_dof_handler().get_orientation(localid);
            generate_stencil<Orient>(SuperStencil::get_dof_handler(), localid, stencil_index, orientation);
        }
/*
        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index) const {
            auto handler = dm.get_dof_handler();
            build_stencil(handler, localid, stencil_index);
        } */

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void generate_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index,
                                                const Integer orientation) const {
            SuperStencil::template generate_stencil<Orient>(dm, localid, stencil_index, orientation);
            build_diagonal_stencil<Orient>(dm, localid, stencil_index, orientation);
        }

        MARS_INLINE_FUNCTION
        ViewMatrixTypeRC<Integer, Face_Length> get_face_stencil() const { return face_extension; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_face_length() const { return Face_Length; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_corner_length() const { return 0; }

    private:
        ViewMatrixTypeRC<Integer, Face_Length> face_extension;
    };

    template <typename DofHandler>
    class FullStokesStencil : StokesStencil<DofHandler> {
    public:
        static constexpr Integer Dim = DofHandler::ManifoldDim;
        static constexpr Integer Corner_Length = 2 * Dim;
        static constexpr Integer Width = 2;

        using SuperStencil = Stencil<DofHandler, Width>;
        using SuperStokesStencil = StokesStencil<DofHandler>;

        MARS_INLINE_FUNCTION
        FullStokesStencil() = default;

        MARS_INLINE_FUNCTION
        FullStokesStencil(DofHandler handler) : SuperStokesStencil(handler) {
            const Integer size = SuperStencil::get_dof_handler().get_owned_dof_size();
            corner_extension = ViewMatrixTypeRC<Integer, Corner_Length>("corner_ext", size);
        }

        virtual void reserve_stencil(const Integer size) override {
            SuperStokesStencil::reserve_stencil(size);
            corner_extension = ViewMatrixTypeRC<Integer, Corner_Length>("stencil_corner_ext", size);
        }

        MARS_INLINE_FUNCTION
        Integer get_corner_value(const Integer row, const Integer col) const { return corner_extension(row, col); }

        MARS_INLINE_FUNCTION
        Integer set_corner_value(const Integer row, const Integer col, const Integer value) const {
            return corner_extension(row, col) = value;
        }

        template <typename F>
        void dof_iterate(F f) const {
            auto st = SuperStencil::get_stencil();
            auto fe = SuperStokesStencil::get_face_stencil();
            auto ce = get_corner_stencil();
            auto length = SuperStencil::get_length();
            auto flength = SuperStokesStencil::get_face_length();
            Kokkos::parallel_for(
                "stencil_dof_iter", SuperStencil::get_stencil_size(), MARS_LAMBDA(const Integer stencil_index) {
                    for (int i = 0; i < length; i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = st(stencil_index, i);
                        f(stencil_index, local_dof);
                    }

                    for (int i = 0; i < flength; i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = fe(stencil_index, i);
                        f(stencil_index, local_dof);
                    }

                    for (int i = 0; i < get_corner_length(); i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = ce(stencil_index, i);
                        f(stencil_index, local_dof);
                    }
                });
        }

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_corner_face_stencil(const DM &dm,
                                                            const Integer localid,
                                                            const Integer stencil_index,
                                                            const Integer orientation) const {
            const Integer sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

            for (int corner = power_of_2(DM::ManifoldDim) - 1; corner != -1; --corner) {
                Octant co = oc.sfc_corner_nbh<DM::simplex_type::ElemType>(corner);

                Integer face_nr;
                // FIXME 3D
                /* The orientation needed here anyway. However the reorder only if Orient true.
                This needs the orientation unlike the previous stencil impl.which can be computed
                without the orientation as only dependend on sfc. ex: face or corner nbh. */
                if (orientation == DofOrient::yDir) {
                    if (corner < 2)
                        face_nr = 2;  // sfc face nbh downwards direction
                    else
                        face_nr = 3;  // sfc face nbh upwards direction
                } else {
                    face_nr = corner % 2;
                }

                // after locating the corner octant co, use it to find the stokes corner.
                Octant o = co.sfc_face_nbh<DM::simplex_type::ElemType>(face_nr);
                const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                set_corner_value(stencil_index, corner, nbh_id);
            }

            /*FIXME for 3D.*/
            // Reorder only if Orient set.
            if (Orient && orientation == DofOrient::xDir) {
                const Integer tmp = corner_extension(stencil_index, 1);
                corner_extension(stencil_index, 1) = corner_extension(stencil_index, 2);
                corner_extension(stencil_index, 2) = tmp;
            }
        }

        template <bool Orient = false>
        MARS_INLINE_FUNCTION void build_stencil(const Integer localid, const Integer stencil_index) const {
            const Integer orientation = SuperStencil::get_dof_handler().get_orientation(localid);
            SuperStokesStencil::template generate_stencil<Orient>(SuperStencil::get_dof_handler(), localid, stencil_index, orientation);
            build_corner_face_stencil<Orient>(SuperStencil::get_dof_handler(), localid, stencil_index, orientation);
        }

        /* template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_stencil(const DM &dm,
                                                const Integer localid,
                                                const Integer stencil_index) const {
            auto handler = dm.get_dof_handler();
            build_stencil(handler, localid, stencil_index);
        } */

        MARS_INLINE_FUNCTION
        ViewMatrixTypeRC<Integer, Corner_Length> get_corner_stencil() const { return corner_extension; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_corner_length() const { return Corner_Length; }

    private:
        ViewMatrixTypeRC<Integer, Corner_Length> corner_extension;
    };
}  // namespace mars

#endif  // mars_stencil
