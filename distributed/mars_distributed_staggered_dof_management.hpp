#ifndef GENERATION_MARS_DISTRIBUTED_SDH_HPP_
#define GENERATION_MARS_DISTRIBUTED_SDH_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <Integer Label, class Mesh_, Integer degree>
    class SDofHandler {
    public:
        using Mesh = Mesh_;

        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        static constexpr Integer dofLabel = Label;
        static constexpr Integer Degree = degree;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        MARS_INLINE_FUNCTION
        SDofHandler(DofHandler<Mesh, degree> d) : dof_handler(d) { prepare_separated_dofs(); }

        ViewVectorType<bool> build_label_dof_predicate(const ViewVectorType<Integer> element_labels) {
            const Integer local_size = element_labels.extent(0);
            ViewVectorType<bool> dof_predicate("label_dof_predicate", local_size);
            Kokkos::parallel_for(
                "separatedoflabelss", local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (element_labels(i) == Label) {
                        dof_predicate(i) = 1;
                    }
                });

            return dof_predicate;
        }

        template <typename V>
        void compact_owned_dofs(V &locally_owned_dofs) {
            using namespace Kokkos;

            const Integer local_size = get_dof_handler().get_global_dof_enum().get_elem_size();
            auto dof_predicate =
                build_label_dof_predicate(get_dof_handler().get_global_dof_enum().get_view_element_labels());

            assert(local_size == get_dof_handler().get_global_dof_enum().get_view_element_labels().extent(0));

            /* perform a scan on the dof predicate*/
            ViewVectorType<Integer> owned_dof_scan("owned_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, dof_predicate, owned_dof_scan);

            auto vol_subview = subview(owned_dof_scan, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            locally_owned_dofs = ViewVectorType<Integer>("locally_owned_dofs", h_vs());
            const ViewVectorType<Integer> global_to_sfc = get_dof_handler().get_global_dof_enum().get_view_elements();

            auto dofhandler = get_dof_handler();

            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (dof_predicate(i) == 1) {
                        Integer vindex = owned_dof_scan(i);
                        const Integer local = dofhandler.sfc_to_local(global_to_sfc(i));
                        locally_owned_dofs(vindex) = local;
                    }
                });
        }

        template <typename V>
        void compact_local_dofs(V &local_dof_map, V &local_dofs) {
            using namespace Kokkos;

            const Integer local_size = get_dof_handler().get_local_dof_enum().get_elem_size();
            auto dof_predicate =
                build_label_dof_predicate(get_dof_handler().get_local_dof_enum().get_view_element_labels());

            assert(local_size == get_dof_handler().get_local_dof_enum().get_view_element_labels().extent(0));

            /* perform a scan on the dof predicate*/
            local_dof_map = ViewVectorType<Integer>("local_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, dof_predicate, local_dof_map);

            auto vol_subview = subview(local_dof_map, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            local_dofs = ViewVectorType<Integer>("local_dofs", h_vs());

            /* Compact the predicate into the volume and face dofs views */
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (dof_predicate(i) == 1) {
                        Integer vindex = local_dof_map(i);
                        local_dofs(vindex) = i;
                    }
                });
        }
        struct IsSeparatedDof {
            ViewVectorType<Integer> local_separated_dof_map;

            MARS_INLINE_FUNCTION
            IsSeparatedDof(ViewVectorType<Integer> map) : local_separated_dof_map(map) {}

            MARS_INLINE_FUNCTION
            bool operator()(const Integer local_dof) const {
                if ((local_dof + 1) >= local_separated_dof_map.extent(0)) return false;
                /*use the map which is the scan of the predicate.
                 * To get the predicate value the difference with the successive index is needed.*/
                const Integer pred_value = local_separated_dof_map(local_dof + 1) - local_separated_dof_map(local_dof);
                return (pred_value > 0);
            }
        };

        void prepare_separated_dofs() {
            compact_local_dofs(local_dof_map, local_dofs);
            compact_owned_dofs(locally_owned_dofs);

            auto is_separated = IsSeparatedDof(local_dof_map);
            // building the counts for boundary and ghost separations to use for gather and scatter separated data only!
            auto boundary_predicate = compact_sfc_to_local(
                get_dof_handler(), is_separated, get_dof_handler().get_boundary_dofs(), boundary_dofs);
            auto boundary_scan = count_sfc_to_local(get_dof_handler().get_view_scan_send(), boundary_predicate);
            scan_send_mirror = create_mirror_view(boundary_scan);
            Kokkos::deep_copy(scan_send_mirror, boundary_scan);

            auto ghost_predicate =
                compact_sfc_to_local(get_dof_handler(), is_separated, get_dof_handler().get_ghost_dofs(), ghost_dofs);
            auto ghost_scan = count_sfc_to_local(get_dof_handler().get_view_scan_recv(), ghost_predicate);
            scan_recv_mirror = create_mirror_view(ghost_scan);
            Kokkos::deep_copy(scan_recv_mirror, ghost_scan);
        }

        MARS_INLINE_FUNCTION Integer get_dof_index(const Integer local_dof) const {
            return get_local_dof_map(local_dof);
        }

        template <typename F>
        void owned_iterate(F f) const {
            Kokkos::parallel_for("owned_separated_dof_iter", get_owned_dof_size(), f);
        }

        template <typename F>
        void iterate(F f) const {
            Kokkos::parallel_for("separated_dof_iter", get_dof_size(), f);
        }

        template <typename F>
        void ghost_iterate(F f) const {
            Kokkos::parallel_for("ghost_dof_iter", get_ghost_dofs().extent(0), f);
        }

        template <typename F>
        void owned_dof_iterate(F f) const {
            auto dofs = get_owned_dofs();
            Kokkos::parallel_for(
                "separated_dof_iter", dofs.extent(0), MARS_LAMBDA(const Integer i) { f(dofs(i)); });
        }

        template <typename F>
        void dof_iterate(F f) const {
            auto dofs = get_local_dofs();
            Kokkos::parallel_for(
                "separated_dof_iter", dofs.extent(0), MARS_LAMBDA(const Integer i) { f(dofs(i)); });
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_local_dofs() const { return local_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_local_dof(const Integer i) const { return local_dofs(i); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_owned_dofs() const { return locally_owned_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_owned_dof(const Integer i) const { return locally_owned_dofs(i); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_local_dof_map() const { return local_dof_map; }

        MARS_INLINE_FUNCTION
        const Integer get_local_dof_map(const Integer local_dof) const { return local_dof_map(local_dof); }

        MARS_INLINE_FUNCTION
        Integer get_boundary_dof(const Integer index) const { return boundary_dofs(index); }

        MARS_INLINE_FUNCTION
        Integer get_ghost_dof(const Integer index) const { return ghost_dofs(index); }

        MARS_INLINE_FUNCTION
        Integer get_boundary_dof_size() const { return boundary_dofs.extent(0); }

        MARS_INLINE_FUNCTION
        Integer get_ghost_dof_size() const { return ghost_dofs.extent(0); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_recv_mirror() const { return scan_recv_mirror; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_send_mirror() const { return scan_send_mirror; }

        MARS_INLINE_FUNCTION
        Integer get_dof_size() const { return get_local_dofs().extent(0); }

        MARS_INLINE_FUNCTION
        Integer get_owned_dof_size() const { return get_owned_dofs().extent(0); }

        MARS_INLINE_FUNCTION
        const DofHandler<Mesh, degree> &get_dof_handler() const { return dof_handler; }

        const context &get_context() const { return get_dof_handler().get_context(); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_boundary_dofs() const { return boundary_dofs; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_ghost_dofs() const { return ghost_dofs; }


        /* *******dof handler related functionalities for completing the handler.******* */
        /* chose this way to hide the full interface of the general handler. Inheritance is the other way*/

        MARS_INLINE_FUNCTION
        const Integer get_orientation(const Integer local_dof) const {
            return get_dof_handler().get_orientation(local_dof);
        }

        MARS_INLINE_FUNCTION
        Dof local_to_global_dof(const Integer local) const {
            return get_dof_handler().local_to_global_dof(local);
        }

        MARS_INLINE_FUNCTION
        UD get_data() const { return get_dof_handler().get_data(); }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_local_dof_enum() const { return get_dof_handler().get_local_dof_enum(); }

        MARS_INLINE_FUNCTION
        Integer local_to_sfc(const Integer local) const { return get_local_dof_enum().get_view_elements()(local); }


        MARS_INLINE_FUNCTION
        Integer sfc_to_local(const Integer sfc) const { return get_local_dof_enum().get_view_sfc_to_local()(sfc); }

        template <Integer Type>
        static MARS_INLINE_FUNCTION Integer
        enum_corner(const ViewVectorType<Integer> &sfc_to_local, const Octant &oc, const int i, const int j) {
            return DofHandler<Mesh, degree>::template enum_corner<Type>(sfc_to_local, oc, i, j);
        }

        MARS_INLINE_FUNCTION
        bool is_local(const Integer sfc) const { return get_dof_handler().is_local(sfc); }

        template <Integer part>
        static MARS_INLINE_FUNCTION Octant enum_face_corner(Octant &oc, const int dir) {
            return DofHandler<Mesh, degree>::template enum_face_corner<part>(oc, dir);
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION Integer enum_face_node(const ViewVectorType<Integer> &sfc_to_local,
                                                           const Octant &face_cornerA,
                                                           const int j,
                                                           const int dir) {
            return DofHandler<Mesh, degree>::template enum_face_node<part, Type>(sfc_to_local, face_cornerA, j, dir);
        }

        /* *************************************************************************** */

    private:
        DofHandler<Mesh, degree> dof_handler;
        // local dofs vector of locally owned dofs. Needed to build the stencils.
        ViewVectorType<Integer> locally_owned_dofs;

        // needed to assign the data to each local dof (including ghosts)
        ViewVectorType<Integer> local_dofs;
        ViewVectorType<Integer> local_dof_map;

        // boundary and ghost local dofs predicated for the separated dm.
        ViewVectorType<Integer> boundary_dofs;
        ViewVectorType<Integer>::HostMirror scan_send_mirror;

        ViewVectorType<Integer> ghost_dofs;
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;
        };

    template <class DofHandler>
    using VolumeDofHandler = SDofHandler<DofLabel::lVolume, typename DofHandler::Mesh, DofHandler::Degree>;

    template <class DofHandler>
    using FaceDofHandler = SDofHandler<DofLabel::lFace, typename DofHandler::Mesh, DofHandler::Degree>;

    template <class DofHandler>
    using CornerDofHandler = SDofHandler<DofLabel::lCorner, typename DofHandler::Mesh, DofHandler::Degree>;

}  // namespace mars

#endif
#endif

#endif
