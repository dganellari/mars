#ifndef GENERATION_MARS_DISTRIBUTED_FE_HPP_
#define GENERATION_MARS_DISTRIBUTED_FE_HPP_

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    class IFEDofMap {
    public:
        virtual MARS_INLINE_FUNCTION ~IFEDofMap() {}
    };

    template <class DofHandler, bool Overlap = true, Integer Label = DofHandler::dofLabel>
    class FEDofMap : public IFEDofMap {
    public:
        using DHandler = DofHandler;
        static constexpr Integer degree = DofHandler::Degree;
        static constexpr Integer Dim = DofHandler::ManifoldDim;

        static constexpr Integer Block = DofHandler::Block;

        static constexpr Integer ElemType = DofHandler::ElemType;

        static constexpr Integer volume_nodes = DofHandler::volume_dofs;
        static constexpr Integer face_nodes = DofHandler::face_dofs;
        static constexpr Integer corner_nodes = DofHandler::corner_dofs;

        using NDofs = NumDofs<degree, Label, ElemType>;
        /* elem_nodes = (degree + 1) ^ Dim; Should be in general. */
        static constexpr Integer elem_nodes = Block * NDofs::elem_dofs();

        MARS_INLINE_FUNCTION
        Integer get_elem_nodes() const { return get_dof_handler().get_block() * NDofs::elem_dofs(); }

        //! regular grids only. This can be a specialized version of the  non uniform impl.
        // Each node has max 2^Dim neighboring elements and  3^DIM neighboring nodes.
        // To account for the dofs we multiply with the degree.
        /* static constexpr Integer max_dof_to_dof_size = power((degree * 2 + 1), Dim); */

        MARS_INLINE_FUNCTION
        FEDofMap() = default;

        MARS_INLINE_FUNCTION
        FEDofMap(DofHandler handler) : dof_handler(handler) {}

        static MARS_INLINE_FUNCTION Integer enumerate_volume(const DofHandler &handler, const Octant &o) {
            Integer sfc = get_sfc_from_octant<ElemType>(o);
            Integer localid = handler.sfc_to_local(sfc);
            return localid;
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        volume_iterate(const DofHandler &handler, F f, const Octant &oc, Integer &index) {
            Octant o;
            // go through all the inside dofs for the current element
            // Currently no topological order within the volume dofs if more than 1.
            for (int j = 0; j < degree - 1; j++) {
                for (int i = 0; i < degree - 1; i++) {
                    o.x = degree * oc.x + i + 1;
                    o.y = degree * oc.y + j + 1;
                    f(index, enumerate_volume(handler, o));
                }
            }
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        volume_iterate(const DofHandler &handler, F f, const Octant &oc, Integer &index) {
            Octant o;
            // go through all the inside dofs for the current element
            // Currently no topological order within the volume dofs if more than 1.
            for (int j = 0; j < degree - 1; j++) {
                for (int i = 0; i < degree - 1; i++) {
                    for (int k = 0; k < degree - 1; k++) {
                        o.x = degree * oc.x + i + 1;
                        o.y = degree * oc.y + j + 1;
                        o.z = degree * oc.z + k + 1;
                        f(index, enumerate_volume(handler, o));
                    }
                }
            }
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        corner_iterate(const DofHandler &dofHandler, F f, const Octant &oc, Integer &index) {
            // go through all the corner dofs for the current element counterclockwise
            /* for (const auto &x : {{0,0}, {1, 0}, {1, 1}, {0, 1}}) */
            /* for (i,j from 0 to 2) first = (i + j) % 2; second = i;*/
            Integer localid = dofHandler.template enum_corner<ElemType>(oc, 0, 0);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 1, 0);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 1, 1);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 0, 1);
            f(index, localid);
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        corner_iterate(const DofHandler &dofHandler, F f, const Octant &oc, Integer &index) {
            // go through all the corner dofs for the current element counterclockwise
            Integer localid = dofHandler.template enum_corner<ElemType>(oc, 0, 0, 0);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 1, 0, 0);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 1, 1, 0);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 0, 1, 0);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 0, 0, 1);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 1, 0, 1);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 1, 1, 1);
            f(index, localid);

            localid = dofHandler.template enum_corner<ElemType>(oc, 0, 1, 1);
            f(index, localid);
        }

        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const DofHandler &dofHandler,
                                                      F f,
                                                      const Octant &oc,
                                                      Integer &index,
                                                      const int side) {
            Octant face_cornerA = DofHandler::template enum_face_corner<dir, ElemType>(oc, side);

            for (int j = 0; j < face_nodes; j++) {
                Integer localid = dofHandler.template enum_face_node<dir, ElemType>(face_cornerA, j, side);
                f(index, localid);
            }
        }

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T == ElementType::Hex8, void>
        edge_iterate(const DofHandler &dofHandler, F f, const Octant &oc, Integer &index, const int edge) {
            const Octant start = dofHandler.get_mesh().get_octant_edge_start(oc, edge);
            const Integer direction = oc.get_edge_direction(edge);

            for (int j = 0; j < DofHandler::edge_dofs; j++) {
                Integer dof_sfc = DofHandler::template process_edge_node<ElemType>(
                    mars::get_sfc_from_octant<ElemType>, start, direction, j);
                auto localid = dofHandler.sfc_to_local(dof_sfc);
                f(index, localid);
            }
        }

        // counter clockwise enumeration of the element dofs.
        template <typename F, bool G, Integer L = Label, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        ordered_dof_enumeration(const DofHandler &handler, F f, const Integer i, Integer &index) {
            Octant oc = DofHandler::template get_octant_ghost_or_local<G>(handler.get_mesh(), i);

            if (L & DofLabel::lCorner) {
                corner_iterate(handler, f, oc, index);
            }
            if (L & DofLabel::lFace) {
                face_iterate<0>(handler, f, oc, index, 0);
                face_iterate<0>(handler, f, oc, index, 1);
                face_iterate<1>(handler, f, oc, index, 0);
                face_iterate<1>(handler, f, oc, index, 1);
            }
            if (L & DofLabel::lVolume) {
                volume_iterate(handler, f, oc, index);
            }
        }

        /* Handling cases for DofLabel::lAll (classical dof handler)  & separated handler like: lCorner + lFace */
        template <typename F, bool G, Integer L = Label, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        ordered_dof_enumeration(const DofHandler &handler, F f, const Integer i, Integer &index) {
            Octant oc = DofHandler::template get_octant_ghost_or_local<G>(handler.get_mesh(), i);

            if (L & DofLabel::lCorner) {
                corner_iterate(handler, f, oc, index);
            }
            if (L & DofLabel::lEdge) {
                edge_iterate(handler, f, oc, index, 0);
                edge_iterate(handler, f, oc, index, 5);
                edge_iterate(handler, f, oc, index, 1);
                edge_iterate(handler, f, oc, index, 4);

                edge_iterate(handler, f, oc, index, 8);
                edge_iterate(handler, f, oc, index, 9);
                edge_iterate(handler, f, oc, index, 11);
                edge_iterate(handler, f, oc, index, 10);

                edge_iterate(handler, f, oc, index, 2);
                edge_iterate(handler, f, oc, index, 7);
                edge_iterate(handler, f, oc, index, 3);
                edge_iterate(handler, f, oc, index, 6);
            }
            if (L & DofLabel::lVolume) {
                volume_iterate(handler, f, oc, index);
            }
            if (L & DofLabel::lFace) {
                face_iterate<2>(handler, f, oc, index, 1);  // z direction first for the counterclockwise
                face_iterate<2>(handler, f, oc, index, 0);  // z direction first for the counterclockwise
                face_iterate<0>(handler, f, oc, index, 1);  // x dir
                face_iterate<0>(handler, f, oc, index, 0);  // x dir
                face_iterate<1>(handler, f, oc, index, 1);  // y dir
                face_iterate<1>(handler, f, oc, index, 0);  // y dir
            }
        }

        template <typename F, bool G, Integer L = Label>
        struct EnumLocalDofs {
            MARS_INLINE_FUNCTION void operator()(const Integer i) const {
                Integer index = 0;
                // topological order within the element
                ordered_dof_enumeration<F, G, L>(dofHandler, f, i, index);
            }

            MARS_INLINE_FUNCTION
            EnumLocalDofs(DofHandler d, F fun) : dofHandler(d), f(fun) {}

            DofHandler dofHandler;
            F f;
        };

        /* template <typename F, bool G, Integer L = Label, typename O = std::nullopt_t>
        struct EnumLocalDofs {
            MARS_INLINE_FUNCTION void operator()(const Integer i) const {
                Integer index = 0;
                topological order within the element
                ordered_dof_enumeration<F, G, L>(dofHandler, f, i, index);

                optinal functor to be used in special cases.
                if (o) {
                    auto callable = o.value();
                    callable(i, index);
                }
            }

            MARS_INLINE_FUNCTION
            EnumLocalDofs(DofHandler d, F fun) : dofHandler(d), f(fun) {}

            MARS_INLINE_FUNCTION
            EnumLocalDofs(DofHandler d, F fun, std::optional<O> opt) : dofHandler(d), f(fun), o(opt) {}

            DofHandler dofHandler;
            F f;
            std::optional<O> o;
        }; */

        template <typename P, typename S>
        void invert_predicate_and_reset_scan(const P &t, const S &v) {
            Kokkos::parallel_for(
                "invert_predicate_and_reset_scan", t.extent(0), MARS_LAMBDA(const Integer i) {
                    t(i) ^= 1;  // invert from 0 to 1 and viceversa.
                    v(i + 1) = 0;
                });
        }

        template <typename P, typename S>
        void build_map_from_scan(const P &predicate, const S &scan, const Integer offset = 0) {
            /* auto map = sfc_to_elem_index; */
            auto view = elem_index;
            Kokkos::parallel_for(
                "build_map_from_scan", predicate.extent(0), MARS_LAMBDA(const Integer i) {
                    if (predicate(i) == 1) {
                        auto index = scan(i) + offset;
                        /* map.insert(i, index); */
                        view(index) = i;
                    }
                });
        }

        // Build predicate for elements that contain only owned dofs and for all the others.
        struct CountOwned {
            DofHandler handler;

            MARS_INLINE_FUNCTION
            CountOwned(DofHandler h) : handler(h) {}

            MARS_INLINE_FUNCTION void operator()(Integer &index, const Integer localid) const {
                // Set to 1 only elements that contain at least one non-owned dof.
                if (localid > INVALID_INDEX && !handler.template is_owned<1>(localid)) {
                    ++index;
                }
            }
        };

        // Handling block structures (vector valued) FE
        struct DofMap {
            ViewMatrixType<Integer> dof_enum;
            Integer sfc_index;
            Integer block;

            MARS_INLINE_FUNCTION
            DofMap(ViewMatrixType<Integer> ede, Integer s, Integer b) : dof_enum(ede), sfc_index(s), block(b) {}

            MARS_INLINE_FUNCTION void operator()(Integer &index, const Integer localid) const {
                for (Integer i = 0; i < block; ++i) {
                    dof_enum(sfc_index, index++) = block * localid + i;
                }
            }
        };

        template <typename V>
        Integer get_size_from_scan(const V &scan) {
            auto ps = Kokkos::subview(scan, scan.extent(0) - 1);
            auto h_ps = Kokkos::create_mirror_view(ps);
            // Deep copy device view to host view.
            Kokkos::deep_copy(h_ps, ps);

            return h_ps();
        }

        // Store the elements with all owned dofs first (first matrix rows) and then elements with non-owned dofs.
        // Finally add ghost elements. This to overlap assembly with communication of ghost layer data.
        template <bool O = Overlap>
        std::enable_if_t<O == true, void> build_elem_index() {
            auto handler = get_dof_handler();
            const Integer size = handler.get_mesh().get_chunk_size();

            elem_index = ViewVectorType<Integer>("elem_index", size);

            ViewVectorType<bool> predicate("predicate_elements", size);
            ViewVectorType<Integer> scan("predicate_elements", size + 1);

            // possibly when we switch to C++17
            /* Kokkos::parallel_for("predicate_local_dofs",
                                 size,
                                 EnumLocalDofs<CountOwned, false, Label, PredicateOwned>(
                                     handler, CountOwned(handler), std::make_optional<PredicateOwned>(predicate))); */

            Kokkos::parallel_for(
                "predicate_local_dofs", size, MARS_LAMBDA(const Integer i) {
                    Integer index = 0;
                    ordered_dof_enumeration<CountOwned, false>(handler, CountOwned(handler), i, index);
                    // If index remained 0 (no unowned dof were found in the element) then set the predicate.
                    if (index == 0) {
                        predicate(i) = 1;
                    }
                });

            incl_excl_scan(0, size, predicate, scan);
            owned_size = get_size_from_scan(scan);
            build_map_from_scan(predicate, scan);

            // First predicate is used to filter the elements that contain the non-owned.
            // After the invertion the predicate filters the elements with only owned dofs.
            invert_predicate_and_reset_scan(predicate, scan);

            incl_excl_scan(0, size, predicate, scan);
            non_owned_size = get_size_from_scan(scan);
            build_map_from_scan(predicate, scan, owned_size);

            assert(size == non_owned_size + owned_size);
        }

        template <bool O = Overlap>
        std::enable_if_t<O == false, void> build_elem_index() {}

        // Store the elements with all owned dofs first (first matrix rows) and then elements with non-owned dofs.
        // Finally add ghost elements. This to overlap assembly with communication of ghost layer data.
        void enumerate_local_dofs() {
            auto handler = get_dof_handler();
            const Integer size = handler.get_mesh().get_chunk_size();
            const Integer ghost_size = handler.get_mesh().get_ghost_size();
            const Integer block = handler.get_block();

            build_elem_index();

            // This is the way to go to avoid atomic on critical code parts when building the unique sorted sparsity
            // pattern. Further on, it can be optimized by storing it to shared memory.
            elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", size + ghost_size, get_elem_nodes());
            /* enumerates the dofs within each element topologically */

            // This view copy is not needed if C++17 is used!
            auto elem_dof_enum_view = elem_dof_enum;
            auto elem_index_view = elem_index;

            Kokkos::parallel_for(
                "enum_local_dofs", size, MARS_LAMBDA(const Integer i) {
                    Integer index = 0;
                    // write coalesced by going through the new order of ids and mapping it to the old order.
                    auto sfc_index = get_element_index(elem_index_view, i);
                    // topological order within the element
                    ordered_dof_enumeration<DofMap, false>(
                        handler, DofMap(elem_dof_enum_view, i, block), sfc_index, index);
                });

            // go through the ghost layer
            Kokkos::parallel_for(
                "enum_ghost_dofs", ghost_size, MARS_LAMBDA(const Integer i) {
                    Integer index = 0;
                    ordered_dof_enumeration<DofMap, true>(
                        handler, DofMap(elem_dof_enum_view, i + size, block), i, index);
                });
        }

        // Unique number of elemnts that share a dof
        template <Integer LL>
        constexpr Integer label_based_element_count() const {
            switch (LL) {
                case DofLabel::lVolume: {
                    return 1;
                }
                case DofLabel::lCorner: {
                    return DofHandler::num_corners;
                }
                case DofLabel::lFace: {
                    return 2;
                }
                case DofLabel::lEdge: {
                    return 4;
                }
                default: {
                    printf("Invalid Label!\n");
                    return 0;
                }
            }
        }

        struct NodeElementDofMap {
            DofHandler handler;
            ViewMatrixType<Integer> dof_enum;
            ViewVectorType<Integer> owned_index;
            ViewVectorType<Integer> owned_map;
            Integer sfc_index;

            MARS_INLINE_FUNCTION
            NodeElementDofMap(DofHandler h,
                              ViewMatrixType<Integer> ede,
                              ViewVectorType<Integer> o,
                              ViewVectorType<Integer> om,
                              Integer s)
                : handler(h), dof_enum(ede), owned_index(o), owned_map(om), sfc_index(s) {}

            MARS_INLINE_FUNCTION void operator()(Integer &index, const Integer localid) const {
                // base local ID is received in all cases, i.e Block=1. No need for vector valued.
                if (localid > INVALID_INDEX) {
                    auto lid = handler.template local_to_owned_index<1>(localid);

                    if (lid > INVALID_INDEX) {
                        auto id = owned_map(lid);
                        for (Integer bi = 0; bi < handler.get_block(); ++bi) {
                            auto bid = handler.get_block() * id + bi;
                            auto aindex = Kokkos::atomic_fetch_add(&owned_index(bid), 1);
                            dof_enum(bid, aindex) = sfc_index;
                        }
                    }
                }
            }
        };

        template <bool O = Overlap>
        static MARS_INLINE_FUNCTION std::enable_if_t<O == true, Integer> get_element_index(
            const ViewVectorType<Integer> &view,
            const Integer i) {
            return view(i);
        }

        template <bool O = Overlap>
        static MARS_INLINE_FUNCTION std::enable_if_t<O == false, Integer> get_element_index(
            const ViewVectorType<Integer> &view,
            const Integer i) {
            return i;
        }

        template <Integer L = Label>
        ViewMatrixType<Integer> build_node_element_dof_map(ViewVectorType<Integer> &locally_owned_dofs) const {
            auto handler = get_dof_handler();

            const Integer size = handler.get_mesh().get_chunk_size();
            const Integer ghost_size = handler.get_mesh().get_ghost_size();

            /* ViewVectorType<Integer> locally_owned_dofs; */
            auto owned_dof_map = compact_owned_dofs<L>(get_dof_handler(), locally_owned_dofs);
            const Integer owned_size = handler.get_block() * locally_owned_dofs.extent(0);

            auto node_max_size = label_based_element_count<L>();
            ViewMatrixType<Integer> dof_enum("build_node_element_dof_map", owned_size, node_max_size);
            Kokkos::parallel_for(
                "init", owned_size, MARS_LAMBDA(const Integer i) {
                    for (int j = 0; j < node_max_size; j++) {
                        dof_enum(i, j) = -1;
                    }
                });

            ViewVectorType<Integer> owned_index("owned_index", owned_size);
            /* enumerates the dofs within each element topologically */
            auto elem_index_view = elem_index;
            Kokkos::parallel_for(
                "enum_local_dofs", size, MARS_LAMBDA(const Integer i) {
                    Integer index = 0;
                    auto sfc_index = get_element_index(elem_index_view, i);
                    ordered_dof_enumeration<NodeElementDofMap, false, L>(
                        handler, NodeElementDofMap(handler, dof_enum, owned_index, owned_dof_map, i), sfc_index, index);
                });

            Kokkos::parallel_for(
                "ghost_enum_local_dofs", ghost_size, MARS_LAMBDA(const Integer i) {
                    Integer index = 0;
                    ordered_dof_enumeration<NodeElementDofMap, true, L>(
                        handler, NodeElementDofMap(handler, dof_enum, owned_index, owned_dof_map, i + size), i, index);
                });

            return dof_enum;
        }

        /*      template <typename DM, typename H>
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
        const ViewMatrixTypeRC<Integer, max_dof_to_dof_size> get_dof_to_dof_map() const { return dof_to_dof_map; }
      */

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
            const Integer size = get_dof_handler().get_mesh().get_chunk_size();
            return (is_valid(sfc_index) && sfc_index >= size);
        }

        MARS_INLINE_FUNCTION
        bool is_owned(const Integer sfc_index) const {
            const Integer size = get_dof_handler().get_mesh().get_chunk_size();
            return (is_valid(sfc_index) && sfc_index < size);
        }

        /* MARS_INLINE_FUNCTION
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
        } */

        template <typename F>
        void owned_dof_element_iterate(F f) const {
            Kokkos::parallel_for("fedomap owned iterate", owned_size, f);
            /* Kokkos::parallel_for("fedomap owned iterate", owned_size, f(elem_index(i)); */
        }

        template <typename F>
        void ghost_dof_element_iterate(F f) const {
            const Integer size = get_fe_dof_map_size();
            Kokkos::parallel_for("fedomap ghost iterate", Kokkos::RangePolicy<>(owned_size, size), f);
            /* Kokkos::parallel_for("fedomap ghost iterate", Kokkos::RangePolicy<>(owned_size, size), f(elem_index(i));
             */
        }

        // iterate on the elements that contain only owned local dofs.
        template <typename F>
        void non_owned_dof_element_iterate(F f) const {
            const Integer local_size = get_dof_handler().get_mesh().get_chunk_size();
            Kokkos::parallel_for("fedomap non owned dof iterate", Kokkos::RangePolicy<>(owned_size, local_size), f);
        }

        // iterate on the owned elements
        template <typename F>
        void owned_element_iterate(F f) const {
            const Integer local_size = get_dof_handler().get_mesh().get_chunk_size();
            Kokkos::parallel_for("fedomap ghost iterate", local_size, f);
        }

        // iterate on the elements that containt ghost dofs.
        template <typename F>
        void ghost_element_iterate(F f) const {
            const Integer local_size = get_dof_handler().get_mesh().get_chunk_size();
            const Integer size = get_fe_dof_map_size();
            Kokkos::parallel_for("fedomap ghost iterate", Kokkos::RangePolicy<>(local_size, size), f);
        }

        // iterate on all elements, owned + ghost
        template <typename F>
        void iterate(F f) const {
            Kokkos::parallel_for("fedomap iterate", get_fe_dof_map_size(), f);
        }

        // print thlocal and the global number of the dof within each element.
        // the dof enumeration within eachlement is topological
        void print() {
            auto fe = *this;
            iterate(MARS_LAMBDA(const Integer elem_index) {
                // go through all the dofs of the elem_index element
                for (int i = 0; i < fe.get_elem_nodes(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                    // convert the local dof number to global dof number
                    Dof d = fe.get_dof_handler().local_to_global_dof(local_dof);
                    auto octant = fe.get_dof_handler().get_octant_from_local(local_dof);

                    auto base_global = fe.get_dof_handler().compute_base(d.get_gid());

                    printf(
                        "fe_dof: i: %li, local: %li, octant: [%li, %li, %li], global: %li, base_global: %li, proc: "
                        "%li\n",
                        i,
                        local_dof,
                        octant.x,
                        octant.y,
                        octant.z,
                        d.get_gid(),
                        base_global,
                        d.get_proc());
                }
            });
        }

        const Integer get_fe_dof_map_size() const { return elem_dof_enum.extent(0); }
        const Integer get_fe_size() const { return elem_dof_enum.extent(1); }
        const Integer get_owned_dof_elements_size() const { return owned_size; }
        const Integer get_ghost_dof_elements_size() const {
            const Integer size = get_fe_dof_map_size();
            return size - owned_size;
        }
        const Integer get_non_owned_dof_elements_size() const { return non_owned_size; }

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
        Integer owned_size;
        Integer non_owned_size;

        ViewVectorType<Integer> elem_index;
        // A map serving as an sfc to local for the dof map element indices. Maybe needed in the future.
        /* UnorderedMap<Integer, Integer> sfc_to_elem_index; */
    };

    /* template <typename Mesh, Integer Degree>
    auto build_dof_to_dof_map(const DofHandler<Mesh, Degree> &handler) {
        FEDofMap<Degree> fe;
        using H = DofHandler<Mesh, Degree>;
        fe.template build_dof_to_dof_map<H>(handler);
        return fe;
    }

 */
    template <class DofHandler, bool Overlap = true, Integer Label = DofHandler::dofLabel>
    auto build_fe_dof_map(const DofHandler &handler) {
        FEDofMap<DofHandler, Overlap, Label> fe(handler);
        fe.enumerate_local_dofs();
        return fe;
    }

    template <class DofHandler, class Mesh, bool Overlap = true, Integer Label = DofHandler::dofLabel>
    auto build_fe_dof_map(Mesh &mesh, const Integer block = 1) {
        DofHandler handler(mesh);
        handler.enumerate_dofs();
        handler.set_block(block);
        return build_fe_dof_map<DofHandler, Overlap, Label>(handler);
    }

}  // namespace mars

#endif
#endif

#endif
