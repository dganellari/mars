#ifndef GENERATION_MARS_DISTRIBUTED_FE_HPP_
#define GENERATION_MARS_DISTRIBUTED_FE_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"

namespace mars {

    template <Integer degree>
    class FEDofMap {
    public:

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        FEDofMap() = default;

        template <class DofHandler>
        struct EnumLocalDofs {
            using SimplexType = typename DofHandler::simplex_type;
            static constexpr Integer Type = SimplexType::ElemType;

            MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc_index, Integer &index) const {
                Octant o;
                Octant oc = dofHandler.get_data().get_mesh()->get_octant(sfc_index);
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
                Octant oc = dofHandler.get_data().get_mesh()->get_octant(sfc_index);

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
                Octant oc = dofHandler.get_data().get_mesh()->get_octant(sfc_index);
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
                corner_iterate(i, index);
                face_iterate<0>(i, index);
                face_iterate<1>(i, index);
                volume_iterate(i, index);
                // TODO: 3D part
            }

            EnumLocalDofs(DofHandler d, ViewMatrixType<Integer> ede, ViewVectorType<Integer> stl)
                : dofHandler(d), elem_dof_enum(ede), sfc_to_local(stl) {}

            DofHandler dofHandler;
            ViewMatrixType<Integer> elem_dof_enum;
            ViewVectorType<Integer> sfc_to_local;
        };

        template <class DofHandler>
        void enumerate_local_dofs(DofHandler &handler) {
            const Integer size = handler.get_data().get_host_mesh()->get_chunk_size();
            elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", size, elem_nodes);

            /* enumerates the dofs within each element topologically */
            Kokkos::parallel_for("enum_local_dofs",
                                 size,
                                 EnumLocalDofs<DofHandler>(
                                     handler, elem_dof_enum,
                                     handler.get_local_dof_enum().get_view_sfc_to_local()));
        }

        MARS_INLINE_FUNCTION
        const ViewMatrixType<Integer> get_elem_dof_map() const { return elem_dof_enum; }

        MARS_INLINE_FUNCTION
        const Integer get_elem_local_dof(const Integer elem_index, const Integer i) const {
            return elem_dof_enum(elem_index, i);
        }

    private:
        // local enumeration of the dofs topologically foreach element
        ViewMatrixType<Integer> elem_dof_enum;
    };

    template <typename DofHandler>
    auto build_fe_dof_map(const DofHandler &handler) {
        FEDofMap<DofHandler::Degree> fe;
        fe.enumerate_local_dofs(handler);
        return fe;
    }

}  // namespace mars

#endif
#endif

#endif
