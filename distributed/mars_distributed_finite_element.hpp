#ifndef GENERATION_MARS_DISTRIBUTED_FE_HPP_
#define GENERATION_MARS_DISTRIBUTED_FE_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    template <class DofHandler>
    class FEDofMap {
    public:
        static constexpr Integer degree = DofHandler::Degree;
        static constexpr Integer Dim = DofHandler::ManifoldDim;

        static constexpr Integer ElemType = DofHandler::ElemType;

        static constexpr Integer volume_nodes = DofHandler::volume_dofs;
        static constexpr Integer face_nodes = DofHandler::face_dofs;
        static constexpr Integer corner_nodes = DofHandler::corner_dofs;
        static constexpr Integer elem_nodes = DofHandler::elem_dofs;
        /* static constexpr Integer elem_nodes = (degree + 1) ^ Dim; Should be in general. */

        using DHandler = DofHandler;
        //! regular grids only. This can be a specialized version of the  non uniform impl.
        // Each node has max 2^Dim neighboring elements and  3^DIM neighboring nodes.
        // To account for the dofs we multiply with the degree.
        /* static constexpr Integer max_dof_to_dof_size = power((degree * 2 + 1), Dim); */

        MARS_INLINE_FUNCTION
        FEDofMap() = default;

        MARS_INLINE_FUNCTION
        FEDofMap(DofHandler handler) : dof_handler(handler) {}

        struct EnumLocalDofs {
            MARS_INLINE_FUNCTION void enumerate_volume(const Octant &o, const Integer sfc_index, Integer &index) const {
                Integer sfc = get_sfc_from_octant<ElemType>(o);
                Integer localid = dofHandler.sfc_to_local(sfc);
                elem_dof_enum(sfc_index, index++) = localid;
            }

            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void> volume_iterate(
                const Integer sfc_index,
                Integer &index) const {
                Octant o;
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);
                // go through all the inside dofs for the current element
                // Currently no topological order within the volume dofs if more than 1.
                for (int j = 0; j < degree - 1; j++) {
                    for (int i = 0; i < degree - 1; i++) {
                        o.x = degree * oc.x + i + 1;
                        o.y = degree * oc.y + j + 1;
                        enumerate_volume(o, sfc_index, index);
                    }
                }
            }

            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void> volume_iterate(const Integer sfc_index,
                                                                                                Integer &index) const {
                Octant o;
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);
                // go through all the inside dofs for the current element
                // Currently no topological order within the volume dofs if more than 1.
                for (int j = 0; j < degree - 1; j++) {
                    for (int i = 0; i < degree - 1; i++) {
                        for (int k = 0; k < degree - 1; k++) {
                            o.x = degree * oc.x + i + 1;
                            o.y = degree * oc.y + j + 1;
                            o.z = degree * oc.z + k + 1;
                            enumerate_volume(o, sfc_index, index);
                        }
                    }
                }
            }

            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void> corner_iterate(
                const Integer sfc_index,
                Integer &index) const {
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);

                // go through all the corner dofs for the current element counterclockwise
                /* for (const auto &x : {{0,0}, {1, 0}, {1, 1}, {0, 1}}) */
                /* for (i,j from 0 to 2) first = (i + j) % 2; second = i;*/

                Integer localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 1);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 1);
                elem_dof_enum(sfc_index, index++) = localid;
            }

            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void> corner_iterate(const Integer sfc_index,
                                                                                                Integer &index) const {
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);

                // go through all the corner dofs for the current element counterclockwise

                Integer localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 0, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 0, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 1, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 1, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 0, 1);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 0, 1);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 1, 1);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 1, 1);
                elem_dof_enum(sfc_index, index++) = localid;
            }

            template <Integer part>
            MARS_INLINE_FUNCTION void face_iterate(const Integer sfc_index, Integer &index) const {
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);
                // side  0 means origin side and 1 destination side.
                for (int dir = 0; dir < 2; ++dir) {
                    Octant face_cornerA = DofHandler::template enum_face_corner<part, ElemType>(oc, dir);

                    for (int j = 0; j < face_nodes; j++) {
                        Integer localid =
                            DofHandler::template enum_face_node<part, ElemType>(sfc_to_local, face_cornerA, j, dir);
                        elem_dof_enum(sfc_index, index++) = localid;
                    }
                }
            }

            // counter clockwise enumeration of the element dofs.
            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void> ordered_dof_enumeration(
                const Integer i,
                Integer &index) const {
                if (DofHandler::dofLabel & DofLabel::lCorner) {
                    corner_iterate(i, index);
                }
                if (DofHandler::dofLabel & DofLabel::lFace) {
                    face_iterate<0>(i, index);
                    face_iterate<1>(i, index);
                }
                if (DofHandler::dofLabel & DofLabel::lVolume) {
                    volume_iterate(i, index);
                }
            }

            /* Handling cases for DofLabel::lAll (classical dof handler)  or separated handler like: lCorner + lFace */
            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void> ordered_dof_enumeration(
                const Integer i,
                Integer &index) const {
                if (DofHandler::dofLabel & DofLabel::lCorner) {
                    corner_iterate(i, index);
                }
                if (DofHandler::dofLabel & DofLabel::lEdge) {
                }
                if (DofHandler::dofLabel & DofLabel::lVolume) {
                    volume_iterate(i, index);
                }
                if (DofHandler::dofLabel & DofLabel::lFace) {
                    face_iterate<2>(i, index);  // z direction first for the counterclockwise
                    face_iterate<0>(i, index);  // x dir
                    face_iterate<1>(i, index);  // y dir
                }
            }

            MARS_INLINE_FUNCTION void operator()(const Integer i) const {
                Integer index = 0;
                // topological order within the element
                ordered_dof_enumeration(i, index);
            }

            EnumLocalDofs(DofHandler d, ViewMatrixType<Integer> ede, ViewVectorType<Integer> stl)
                : dofHandler(d), elem_dof_enum(ede), sfc_to_local(stl) {}

            DofHandler dofHandler;
            ViewMatrixType<Integer> elem_dof_enum;
            ViewVectorType<Integer> sfc_to_local;
        };

        void enumerate_local_dofs() {
            auto handler = get_dof_handler();
            const Integer size = handler.get_mesh_manager().get_host_mesh()->get_chunk_size();

            elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", size, elem_nodes);
            /* enumerates the dofs within each element topologically */
            Kokkos::parallel_for(
                "enum_local_dofs",
                size,
                EnumLocalDofs(handler, elem_dof_enum, handler.get_local_dof_enum().get_view_sfc_to_local()));
        }

        /*
                template <typename DM, typename H>
                MARS_INLINE_FUNCTION void add_directional_neighbors(const DM &dm,
                                                                    H map,
                                                                    const Integer index,
                                                                    const Integer dir) const {
                    const Integer end = (2 * degree + 1) ^ dir;

                    Integer col_index = end;
                    for (int i = 0; i < end; ++i) {
                        const Integer sfc = dm.local_to_sfc(map(index, i));
                        Octant oc = get_octant_from_sfc<ElemType>(sfc);

                        Integer face_nr;
                        for (int side = 0; side < 2; ++side) {
                            if (side == 0)
                                face_nr = 2 * dir + 1;
                            else
                                face_nr = 2 * dir;

                            Octant o = oc;
                            for (int w = 0; w < degree; ++w) {
                                o = o.sfc_face_nbh<ElemType>(face_nr);
                                const Integer nbh_sfc = get_sfc_from_octant<ElemType>(o);
                                Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                                map(index, col_index) = nbh_id;
                                ++col_index;
                            }
                        }
                    }
                }

                //TODO: Check for 3D. Should work as well for 3D.
                template <typename DM, typename H>
                MARS_INLINE_FUNCTION void generate_dof_to_dof_map(const DM &dm,
                                                                  H map,
                                                                  const Integer localid,
                                                                  const Integer index) const {
                    map(index, 0) = localid;

                    for (int dir = 0; dir < Dim; ++dir) {
                        add_directional_neighbors(dm, map, index, dir);
                    }
                }

                template <typename DM>
                void build_dof_to_dof_map(const DM &dm) {
                    dof_to_dof_map =
                        ViewMatrixTypeRC<Integer, max_dof_to_dof_size>("node_to_node_dof", dm.get_owned_dof_size());

                    //due to the mars lambda. Not needed if c++17 is used.
                    auto map = dof_to_dof_map;

                    dm.owned_iterate(MARS_LAMBDA(const Integer i) {
                        const Integer localid = dm.get_owned_dof(i);
                        generate_dof_to_dof_map(dm, map, localid, i);
                    });
                } */

        MARS_INLINE_FUNCTION
        const ViewMatrixType<Integer> get_elem_dof_map() const { return elem_dof_enum; }

        /* MARS_INLINE_FUNCTION
        const ViewMatrixTypeRC<Integer, max_dof_to_dof_size> get_dof_to_dof_map() const { return dof_to_dof_map; } */

        MARS_INLINE_FUNCTION
        const Integer get_elem_local_dof(const Integer elem_index, const Integer i) const {
            return elem_dof_enum(elem_index, i);
        }

        /* MARS_INLINE_FUNCTION
        const Integer get_nbh_dof(const Integer elem_index, const Integer i) const {
            return dof_to_dof_map(elem_index, i);
        } */

        /* template <typename F>
        void dof_to_dof_iterate(F f) const {
            Kokkos::parallel_for("fedomap iterate", get_dof_to_dof_map_size(), f);
        } */

        template <typename F>
        void iterate(F f) const {
            Kokkos::parallel_for("fedomap iterate", get_fe_dof_map_size(), f);
        }

        const Integer get_fe_dof_map_size() const { return elem_dof_enum.extent(0); }
        const Integer get_fe_size() const { return elem_dof_enum.extent(1); }

        /* const Integer get_dof_to_dof_map_size() const { return dof_to_dof_map.extent(0); } */
        /* const Integer get_nbh_dof_size() const { return dof_to_dof_map.extent(1); } */
        MARS_INLINE_FUNCTION
        const DofHandler &get_dof_handler() const { return dof_handler; }

    private:
        DofHandler dof_handler;
        // local enumeration of the dofs topologically foreach element
        ViewMatrixType<Integer> elem_dof_enum;
        // owned dof to local dofs map (including ghosts). owned dof rows, (>), local dof columns. Sparse.
        /* ViewMatrixTypeRC<Integer, max_dof_to_dof_size> dof_to_dof_map; */
    };

    /* template <typename Mesh, Integer Degree>
    auto build_dof_to_dof_map(const DofHandler<Mesh, Degree> &handler) {
        FEDofMap<Degree> fe;
        using H = DofHandler<Mesh, Degree>;
        fe.template build_dof_to_dof_map<H>(handler);
        return fe;
    }

 */
    template <class DofHandler>
    auto build_fe_dof_map(const DofHandler &handler) {
        FEDofMap<DofHandler> fe(handler);
        fe.enumerate_local_dofs();
        return fe;
    }

}  // namespace mars

#endif
#endif

#endif
