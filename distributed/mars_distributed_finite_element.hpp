#ifndef GENERATION_MARS_DISTRIBUTED_FE_HPP_
#define GENERATION_MARS_DISTRIBUTED_FE_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    template <class DofHandler, Integer Label = DofHandler::dofLabel>
    class FEDofMap {
    public:
        static constexpr Integer degree = DofHandler::Degree;
        static constexpr Integer Dim = DofHandler::ManifoldDim;

        static constexpr Integer ElemType = DofHandler::ElemType;

        static constexpr Integer volume_nodes = DofHandler::volume_dofs;
        static constexpr Integer face_nodes = DofHandler::face_dofs;
        static constexpr Integer corner_nodes = DofHandler::corner_dofs;
        /* static constexpr Integer elem_nodes = DofHandler::elem_dofs; */

        using NDofs = NumDofs<degree, Label, ElemType>;
        static constexpr Integer elem_nodes = NDofs::elem_dofs();
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

        template <typename F, bool G>
        struct EnumLocalDofs {
            MARS_INLINE_FUNCTION void enumerate_volume(const Octant &o, const Integer sfc_index, Integer &index) const {
                Integer sfc = get_sfc_from_octant<ElemType>(o);
                Integer localid = dofHandler.sfc_to_local(sfc);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */
            }

            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void> volume_iterate(
                const Integer sfc_index,
                Integer &index) const {
                Octant o;
                /* Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index); */
                Octant oc =
                    DofHandler::template get_octant_ghost_or_local<G>(dofHandler.get_mesh_manager().get_mesh(), sfc_index);
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
                /* Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index); */
                Octant oc =
                    DofHandler::template get_octant_ghost_or_local<G>(dofHandler.get_mesh_manager().get_mesh(), sfc_index);
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
                /* Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index); */
                Octant oc = DofHandler::template get_octant_ghost_or_local<G>(dofHandler.get_mesh_manager().get_mesh(),
                                                                             sfc_index);

                // go through all the corner dofs for the current element counterclockwise
                /* for (const auto &x : {{0,0}, {1, 0}, {1, 1}, {0, 1}}) */
                /* for (i,j from 0 to 2) first = (i + j) % 2; second = i;*/

                Integer localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 0);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 0);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 1);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 1);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */
            }

            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void> corner_iterate(const Integer sfc_index,
                                                                                                Integer &index) const {
                /* Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index); */
                Octant oc = DofHandler::template get_octant_ghost_or_local<G>(dofHandler.get_mesh_manager().get_mesh(),
                                                                             sfc_index);

                // go through all the corner dofs for the current element counterclockwise

                Integer localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 0, 0);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 0, 0);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 1, 0);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 1, 0);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 0, 1);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 0, 1);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 1, 1, 1);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */

                localid = DofHandler::template enum_corner<ElemType>(sfc_to_local, oc, 0, 1, 1);
                f(sfc_index, index, localid);
                /* elem_dof_enum(sfc_index, index++) = localid; */
            }

            template <Integer dir>
            MARS_INLINE_FUNCTION void face_iterate(const Integer sfc_index, Integer &index, const int side) const {
                /* Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index); */
                Octant oc = DofHandler::template get_octant_ghost_or_local<G>(dofHandler.get_mesh_manager().get_mesh(),
                                                                             sfc_index);
                Octant face_cornerA = DofHandler::template enum_face_corner<dir, ElemType>(oc, side);

                for (int j = 0; j < face_nodes; j++) {
                    Integer localid =
                        DofHandler::template enum_face_node<dir, ElemType>(sfc_to_local, face_cornerA, j, side);
                    f(sfc_index, index, localid);
                    /* elem_dof_enum(sfc_index, index++) = localid; */
                }
            }

            template <Integer T = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<T == ElementType::Hex8, void> edge_iterate(const Integer sfc_index,
                                                                                             Integer &index,
                                                                                             const int edge) const {
                /* const Octant oc = dofHandler.get_mesh_manager().get_mesh()->get_octant(sfc_index); */
                const Octant oc = DofHandler::template get_octant_ghost_or_local<G>(
                    dofHandler.get_mesh_manager().get_mesh(), sfc_index);

                const Octant start = dofHandler.get_mesh_manager().get_mesh()->get_octant_edge_start(oc, edge);
                const Integer direction = oc.get_edge_direction(edge);

                for (int j = 0; j < DofHandler::edge_dofs; j++) {
                    Integer dof_sfc = DofHandler::template process_edge_node<ElemType>(
                        mars::get_sfc_from_octant<ElemType>, start, direction, j);
                    auto localid = sfc_to_local(dof_sfc);
                    f(sfc_index, index, localid);
                    /* elem_dof_enum(sfc_index, index++) = sfc_to_local(dof_sfc); */
                }
            }

            // counter clockwise enumeration of the element dofs.
            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void> ordered_dof_enumeration(
                const Integer i,
                Integer &index) const {
                if (Label & DofLabel::lCorner) {
                    corner_iterate(i, index);
                }
                if (Label & DofLabel::lFace) {
                    face_iterate<0>(i, index, 0);
                    face_iterate<0>(i, index, 1);
                    face_iterate<1>(i, index, 0);
                    face_iterate<1>(i, index, 1);
                }
                if (Label & DofLabel::lVolume) {
                    volume_iterate(i, index);
                }
            }

            /* Handling cases for DofLabel::lAll (classical dof handler)  & separated handler like: lCorner + lFace */
            template <Integer ET = ElemType>
            MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void> ordered_dof_enumeration(
                const Integer i,
                Integer &index) const {
                if (Label & DofLabel::lCorner) {
                    corner_iterate(i, index);
                }
                if (Label & DofLabel::lEdge) {
                    edge_iterate(i, index, 0);
                    edge_iterate(i, index, 5);
                    edge_iterate(i, index, 1);
                    edge_iterate(i, index, 4);

                    edge_iterate(i, index, 8);
                    edge_iterate(i, index, 9);
                    edge_iterate(i, index, 11);
                    edge_iterate(i, index, 10);

                    edge_iterate(i, index, 2);
                    edge_iterate(i, index, 7);
                    edge_iterate(i, index, 3);
                    edge_iterate(i, index, 6);
                }
                if (Label & DofLabel::lVolume) {
                    volume_iterate(i, index);
                }
                if (Label & DofLabel::lFace) {
                    face_iterate<2>(i, index, 1);  // z direction first for the counterclockwise
                    face_iterate<2>(i, index, 0);  // z direction first for the counterclockwise
                    face_iterate<0>(i, index, 1);  // x dir
                    face_iterate<0>(i, index, 0);  // x dir
                    face_iterate<1>(i, index, 1);  // y dir
                    face_iterate<1>(i, index, 0);  // y dir
                }
            }

            MARS_INLINE_FUNCTION void operator()(const Integer i) const {
                Integer index = 0;
                // topological order within the element
                ordered_dof_enumeration(i, index);
            }

            /* EnumLocalDofs(DofHandler d, ViewMatrixType<Integer> ede, ViewVectorType<Integer> stl) */
                /* : dofHandler(d), elem_dof_enum(ede), sfc_to_local(stl) {} */

            MARS_INLINE_FUNCTION
            EnumLocalDofs(DofHandler d, F fun, ViewVectorType<Integer> stl)
                : dofHandler(d), f(fun), sfc_to_local(stl) {}

            DofHandler dofHandler;
            F f;
            /* ViewMatrixType<Integer> elem_dof_enum; */
            ViewVectorType<Integer> sfc_to_local;
        };

        struct DofMap {
            ViewMatrixType<Integer> dof_enum;
            Integer size;

            MARS_INLINE_FUNCTION
            DofMap(ViewMatrixType<Integer> ede, Integer s) : dof_enum(ede), size(s) {}

            MARS_INLINE_FUNCTION void operator()(const Integer sfc_index, Integer &index, const Integer localid) const {
                dof_enum(sfc_index + size, index++) = localid;
            }
        };

        void enumerate_local_dofs() {
            auto handler = get_dof_handler();
            const Integer size = handler.get_mesh_manager().get_host_mesh()->get_chunk_size();
            const Integer ghost_size = handler.get_mesh_manager().get_host_mesh()->get_ghost_size();

            // This is the way to go to avoid atomic on critical code parts when building the unique sorted sparsity
            // pattern. Further on, it can be optimized by storing it to shared memory.
            elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", size + ghost_size, elem_nodes);
            /* enumerates the dofs within each element topologically */
            Kokkos::parallel_for(
                "enum_local_dofs",
                size,
                EnumLocalDofs<DofMap, false>(
                    handler, DofMap(elem_dof_enum, 0), handler.get_local_dof_enum().get_view_sfc_to_local()));

            //go through the ghost layer
            Kokkos::parallel_for(
                "enum_local_dofs",
                ghost_size,
                EnumLocalDofs<DofMap, true>(
                    handler, DofMap(elem_dof_enum, size), handler.get_local_dof_enum().get_view_sfc_to_local()));
        }

        /* void enumerate_ghost_local_dofs() {
            auto handler = get_dof_handler();
            const Integer ghost_size = handler.get_mesh_manager().get_host_mesh()->get_ghost_size();

            ghost_elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", ghost_size, elem_nodes);
            [>enumerates the dofs within each element topologically<]
            Kokkos::parallel_for(
                "enum_local_dofs",
                ghost_size,
                EnumLocalDofs<DofMap, true>(
                    handler, DofMap(ghost_elem_dof_enum), handler.get_local_dof_enum().get_view_sfc_to_local()));
        } */

        struct NodeElementDofMap {
            DofHandler handler;
            ViewMatrixType<Integer> dof_enum;
            ViewVectorType<Integer> owned_index;
            Integer size;

            MARS_INLINE_FUNCTION
            NodeElementDofMap(DofHandler h, ViewMatrixType<Integer> ede, ViewVectorType<Integer> o, Integer s)
                : handler(h), dof_enum(ede), owned_index(o), size(s) {}

            MARS_INLINE_FUNCTION void operator()(const Integer sfc_index, Integer &index, const Integer localid) const {
                auto id = handler.local_to_owned_index(localid);

                if (id > INVALID_INDEX) {
                    auto aindex = Kokkos::atomic_fetch_add(&owned_index(id), 1);
                    dof_enum(id, aindex) = sfc_index + size;
                }
            }
        };

        ViewMatrixType<Integer> build_node_element_dof_map() {
            auto handler = get_dof_handler();

            const Integer size = handler.get_mesh_manager().get_host_mesh()->get_chunk_size();
            const Integer ghost_size = handler.get_mesh_manager().get_host_mesh()->get_ghost_size();

            auto owned_size = handler.get_owned_dof_size();

            ViewMatrixType<Integer> dof_enum("build_node_element_dof_map", owned_size, ElemType);
            Kokkos::parallel_for(
                "init", owned_size, MARS_LAMBDA(const Integer i) {
                    for (int j = 0; j < ElemType; j++) {
                        dof_enum(i, j) = -1;
                    }
                });

            ViewVectorType<Integer> owned_index("owned_index", owned_size);
            /* enumerates the dofs within each element topologically */
            Kokkos::parallel_for(
                "enum_local_dofs",
                size,
                EnumLocalDofs<NodeElementDofMap, false>(handler,
                                                        NodeElementDofMap(handler, dof_enum, owned_index, 0),
                                                        handler.get_local_dof_enum().get_view_sfc_to_local()));
            Kokkos::parallel_for(
                "enum_local_dofs",
                ghost_size,
                EnumLocalDofs<NodeElementDofMap, true>(handler,
                                                       NodeElementDofMap(handler, dof_enum, owned_index, size),
                                                       handler.get_local_dof_enum().get_view_sfc_to_local()));
            return dof_enum;
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


        MARS_INLINE_FUNCTION
        bool is_valid(const Integer sfc_index) const { return sfc_index > INVALID_INDEX; }

        MARS_INLINE_FUNCTION
        bool is_ghost(const Integer sfc_index) const {
            const Integer size = get_dof_handler().get_mesh_manager().get_host_mesh()->get_chunk_size();
            return (is_valid(sfc_index) && sfc_index >= size);
        }

        MARS_INLINE_FUNCTION
        bool is_owned(const Integer sfc_index) const {
            const Integer size = get_dof_handler().get_mesh_manager().get_host_mesh()->get_chunk_size();
            return (is_valid(sfc_index) && sfc_index < size);
        }


        MARS_INLINE_FUNCTION
        Integer get_elem_sfc(const Integer sfc_index) const {
            const auto m = get_dof_handler().get_mesh_manager().get_mesh();
            const auto size = m->get_chunk_size();

            if (is_owned(sfc_index)) {
                return DofHandler::template get_sfc_ghost_or_local<false>(m, sfc_index);
            } else if (is_ghost(sfc_index)) {
                return DofHandler::template get_sfc_ghost_or_local<true>(m, sfc_index - size);
            } else {
                return INVALID_INDEX;
            }
        }

        template <typename F>
        void owned_iterate(F f) const {
            const Integer size = get_dof_handler().get_mesh_manager().get_host_mesh()->get_chunk_size();
            Kokkos::parallel_for("fedomap iterate", size, f);
        }

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
        /* ViewMatrixType<Integer> ghost_elem_dof_enum; */
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
    template <class DofHandler, Integer Label = DofHandler::dofLabel>
    auto build_fe_dof_map(const DofHandler &handler) {
        FEDofMap<DofHandler, Label> fe(handler);
        fe.enumerate_local_dofs();
        /* fe.enumerate_ghost_local_dofs(); */
        return fe;
    }

}  // namespace mars

#endif
#endif

#endif
