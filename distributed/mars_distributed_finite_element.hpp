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

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);
        /* static constexpr Integer elem_nodes = (degree + 1) ^ Dim; Should be in general. */

        using DHandler = DofHandler;
        //! regular grids only. This can be a specialized version of the  non uniform impl.
        // Each node has max 2^Dim neighboring elements and  3^DIM neighboring nodes.
        // To account for the dofs we multiply with the degree.
        static constexpr Integer max_dof_to_dof_size = power((degree * 2 + 1), Dim);

        MARS_INLINE_FUNCTION
        FEDofMap() = default;

        MARS_INLINE_FUNCTION
        FEDofMap(DofHandler handler) : dof_handler(handler) {}

        struct EnumLocalDofs {
            using SimplexType = typename DofHandler::simplex_type;
            static constexpr Integer Type = SimplexType::ElemType;

            MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc_index, Integer &index) const {
                Octant o;
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);
                // go through all the inside dofs for the current element
                // Currently no topological order within the volume dofs if more than 1.
                for (int j = 0; j < degree - 1; j++) {
                    for (int i = 0; i < degree - 1; i++) {
                        o.x = degree * oc.x + i + 1;
                        o.y = degree * oc.y + j + 1;
                        Integer sfc = get_sfc_from_octant<Type>(o);
                        Integer localid = dofHandler.sfc_to_local(sfc);
                        elem_dof_enum(sfc_index, index++) = localid;
                    }
                }
            }

            MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc_index, Integer &index) const {
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);

                // go through all the corner dofs for the current element counterclockwise
                /* for (const auto &x : {{0,0}, {1, 0}, {1, 1}, {0, 1}}) */
                /* for (i,j from 0 to 2) first = (i + j) % 2; second = i;*/

                Integer localid = DofHandler::template enum_corner<Type>(sfc_to_local, oc, 0, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<Type>(sfc_to_local, oc, 1, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<Type>(sfc_to_local, oc, 1, 1);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = DofHandler::template enum_corner<Type>(sfc_to_local, oc, 0, 1);
                elem_dof_enum(sfc_index, index++) = localid;
            }

            template <Integer part>
            MARS_INLINE_FUNCTION void face_iterate(const Integer sfc_index, Integer &index) const {
                Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index);
                // side  0 means origin side and 1 destination side.
                for (int dir = 0; dir < 2; ++dir) {
                    Octant face_cornerA = DofHandler::template enum_face_corner<part>(oc, dir);

                    for (int j = 0; j < face_nodes; j++) {
                        Integer localid =
                            DofHandler::template enum_face_node<part, Type>(sfc_to_local, face_cornerA, j, dir);
                        elem_dof_enum(sfc_index, index++) = localid;
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                Integer index = 0;
                // topological order within the element
                switch (DofHandler::dofLabel) {
                    case DofLabel::lAll: {
                        corner_iterate(i, index);
                        face_iterate<0>(i, index);
                        face_iterate<1>(i, index);
                        volume_iterate(i, index);
                        break;
                        // TODO: 3D part
                    }
                    case DofLabel::lCorner: {
                        corner_iterate(i, index);
                        break;
                    }
                    case DofLabel::lFace: {
                        face_iterate<0>(i, index);
                        face_iterate<1>(i, index);
                        break;
                    }

                    case DofLabel::lVolume: {
                        volume_iterate(i, index);
                        break;
                    }
                }
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

            Integer elem_size = elem_nodes;

            switch (DofHandler::dofLabel) {
                case DofLabel::lCorner: {
                    elem_size = corner_nodes * power_of_2(Dim);
                    break;
                }
                case DofLabel::lFace: {
                    elem_size = face_nodes * 2 * Dim;
                    break;
                }

                case DofLabel::lVolume: {
                    elem_size = volume_nodes;
                    break;
                }
            }

            elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", size, elem_size);

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
                        Octant oc = get_octant_from_sfc<DM::simplex_type::ElemType>(sfc);

                        Integer face_nr;
                        for (int side = 0; side < 2; ++side) {
                            if (side == 0)
                                face_nr = 2 * dir + 1;
                            else
                                face_nr = 2 * dir;

                            Octant o = oc;
                            for (int w = 0; w < degree; ++w) {
                                o = o.sfc_face_nbh<DM::simplex_type::ElemType>(face_nr);
                                const Integer nbh_sfc = get_sfc_from_octant<DM::simplex_type::ElemType>(o);
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
