#ifndef MARS_STENCIL_HPP
#define MARS_STENCIL_HPP

#include "mars_base.hpp"

#ifdef MARS_ENABLE_KOKKOS_KERNELS
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
        static_assert(DofHandler::Block == 1, "Stencil does not support yet vector valued block structures.");

        ViewVectorType<Integer> locally_owned_dofs;
        compact_owned_dofs<ST::dofLabel>(handler, locally_owned_dofs);

        const Integer size = locally_owned_dofs.extent(0);
        ST stencil(size);

        Kokkos::parallel_for(
            "owned_dofs_label_iterate", size, MARS_LAMBDA(const Integer i) {
                const Integer localid = locally_owned_dofs(i);
                stencil.template build_stencil<Orient>(handler, localid, i);
            });
        return stencil;
    }

    template <Integer Label, Integer Dim = 2, Integer Width = 1, Integer L = 2 * Dim *Width + 1>
    class Stencil {
    public:
        /* static constexpr Integer Dim = DofHandler::ManifoldDim; */
        static constexpr Integer dofLabel = Label;
        static constexpr Integer Length = L;

        Stencil() = default;

        Stencil(const Integer size) { reserve_stencil(size); }

        void reserve_stencil(const Integer size) { stencil = ViewMatrixTypeRC<Integer, L>("stencil", size); }

        ViewMatrixTypeRC<Integer, L> get_stencil() const { return stencil; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_label() const { return dofLabel; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_width() const { return Width; }

        MARS_INLINE_FUNCTION
        constexpr Integer get_length() const { return Length; }

        MARS_INLINE_FUNCTION
        Integer get_value(const Integer row, const Integer col) const { return stencil(row, col); }

        MARS_INLINE_FUNCTION
        void set_value(const Integer row, const Integer col, const Integer value) const {
            stencil(row, col) = value;
        }

        Integer get_stencil_size() const { return stencil.extent(0); }

        template <typename F>
        void dof_iterate(F f) const {
            auto stencil = *this;
            auto st = stencil.get_stencil();
            auto size = stencil.get_stencil_size();
            auto length = stencil.get_length();
            Kokkos::parallel_for(
                "stencil_dof_iter", size, MARS_LAMBDA(const Integer stencil_index) {
                    for (int i = 0; i < length; i++) {
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

        template <bool Orient = false, typename DofHandler>
        MARS_INLINE_FUNCTION void build_stencil(const DofHandler &dof_handler,
                                                const Integer localid,
                                                const Integer stencil_index) const {
            const Integer orientation = dof_handler.get_orientation(localid);
            generate_stencil<Orient>(dof_handler, localid, stencil_index, orientation);
        }

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void generate_stencil(const DM &dm,
                                                   const Integer localid,
                                                   const Integer stencil_index,
                                                   Integer orientation) const {
            set_value(stencil_index, 0, localid);
            const auto sfc = dm.local_to_sfc(localid);

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
                        const auto nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                        Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                        Integer offset = w * 2 * Dim;
                        set_value(stencil_index, index + offset, nbh_id);
                    }
                }
            }
        }

    private:
        ViewMatrixTypeRC<Integer, L> stencil;
    };

    template <Integer Label, Integer Dim = 2, Integer L = 6 * Dim + 1>
    class StokesStencil : public Stencil<Label, Dim, 2, L> {
    public:
        static constexpr Integer dofLabel = Label;
        /* static constexpr Integer Face_Length = 2 * Dim; */
        using SuperStencil = Stencil<Label, Dim, 2, L>;

        StokesStencil() = default;

        StokesStencil(const Integer size) : SuperStencil(size) {}

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_diagonal_stencil(const DM &dm,
                                                         const Integer localid,
                                                         const Integer stencil_index,
                                                         const Integer orientation) const {
            const auto sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

            for (int corner = power_of_2(DM::ManifoldDim) - 1; corner != -1; --corner) {
                Octant o = oc.sfc_corner_nbh<DM::simplex_type::ElemType>(corner);
                const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                set_face_value(stencil_index, corner, nbh_id);
            }

            /*FIXME for 3D.*/
            if (Orient && orientation == DofOrient::xDir) {
                const Integer tmp = get_face_value(stencil_index, 1);
                set_face_value(stencil_index, 1, get_face_value(stencil_index, 2));
                set_face_value(stencil_index, 2, tmp);
            }
        }

        template <bool Orient = false, typename DofHandler>
        MARS_INLINE_FUNCTION void build_stencil(const DofHandler &dof_handler,
                                                const Integer localid,
                                                const Integer stencil_index) const {
            const Integer orientation = dof_handler.get_orientation(localid);
            generate_stencil<Orient>(dof_handler, localid, stencil_index, orientation);
        }

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void generate_stencil(const DM &dm,
                                                   const Integer localid,
                                                   const Integer stencil_index,
                                                   const Integer orientation) const {
            SuperStencil::template generate_stencil<Orient>(dm, localid, stencil_index, orientation);
            build_diagonal_stencil<Orient>(dm, localid, stencil_index, orientation);
        }

    private:
        static constexpr Integer Face_Start = 4 * Dim + 1;

        MARS_INLINE_FUNCTION
        Integer get_face_value(const Integer row, const Integer col) const {
            return SuperStencil::get_value(row, Face_Start + col);
        }

        MARS_INLINE_FUNCTION
        void set_face_value(const Integer row, const Integer col, const Integer value) const {
            SuperStencil::set_value(row, Face_Start + col, value);
        }
    };

    template <Integer Label, Integer Dim = 2, Integer L = 8 * Dim + 1>
    class FullStokesStencil : public StokesStencil<Label, Dim, L> {
    public:
        static constexpr Integer dofLabel = Label;
        /* static constexpr Integer Corner_Length = 2 * Dim; */
        static constexpr Integer Width = 2;

        using SuperStencil = Stencil<Label, Dim, Width, L>;
        using SuperStokesStencil = StokesStencil<Label, Dim, L>;

        FullStokesStencil() = default;

        FullStokesStencil(const Integer size) : SuperStokesStencil(size) {}

        template <bool Orient = false, typename DM>
        MARS_INLINE_FUNCTION void build_corner_face_stencil(const DM &dm,
                                                            const Integer localid,
                                                            const Integer stencil_index,
                                                            const Integer orientation) const {
            const auto sfc = dm.local_to_sfc(localid);

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
                const auto nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                set_corner_value(stencil_index, corner, nbh_id);
            }

            /*FIXME for 3D.*/
            // Reorder only if Orient set.
            if (Orient && orientation == DofOrient::xDir) {
                const Integer tmp = get_corner_value(stencil_index, 1);
                set_corner_value(stencil_index, 1, get_corner_value(stencil_index, 2));
                set_corner_value(stencil_index, 2, tmp);
            }
        }

        template <bool Orient = false, typename DofHandler>
        MARS_INLINE_FUNCTION void build_stencil(const DofHandler &dof_handler,
                                                const Integer localid,
                                                const Integer stencil_index) const {
            const Integer orientation = dof_handler.get_orientation(localid);
            SuperStokesStencil::template generate_stencil<Orient>(dof_handler, localid, stencil_index, orientation);
            build_corner_face_stencil<Orient>(dof_handler, localid, stencil_index, orientation);
        }

    private:
        static constexpr Integer Corner_Start = 6 * Dim + 1;

        MARS_INLINE_FUNCTION
        Integer get_corner_value(const Integer row, const Integer col) const {
            return SuperStencil::get_value(row, Corner_Start + col);
        }

        MARS_INLINE_FUNCTION
        void set_corner_value(const Integer row, const Integer col, const Integer value) const {
            SuperStencil::set_value(row, Corner_Start + col, value);
        }
    };
}  // namespace mars

#endif
#endif  // mars_stencil
