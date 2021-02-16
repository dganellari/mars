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

        static constexpr Integer ElemType = simplex_type::ElemType;

        static constexpr Integer dofLabel = Label;

        static constexpr Integer Degree = degree;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        MARS_INLINE_FUNCTION
        SDofHandler(DofHandler<Mesh, degree> d) : dof_handler(d) { prepare_separated_dofs(); }

        struct IsSeparatedDof {
            ViewVectorType<Integer> local_separated_dof_map;
            DofHandler<Mesh, degree> handler;

            MARS_INLINE_FUNCTION
            IsSeparatedDof(ViewVectorType<Integer> map, DofHandler<Mesh, degree> d)
                : local_separated_dof_map(map), handler(d) {}

            MARS_INLINE_FUNCTION
            bool operator()(const Integer sfc) const {
                const Integer local_dof = handler.sfc_to_local(sfc);
                if ((local_dof + 1) >= local_separated_dof_map.extent(0)) return false;
                /*use the map which is the scan of the predicate.
                 * To get the predicate value the difference with the successive index is needed.*/
                const Integer pred_value = local_separated_dof_map(local_dof + 1) - local_separated_dof_map(local_dof);
                return (pred_value > 0);
            }
        };

        /* struct IsSeparatedOwnedDof {
            ViewVectorType<Integer> local_separated_dof_map;
            DofHandler handler;

            MARS_INLINE_FUNCTION
            IsSeparatedOwnedDof(ViewVectorType<Integer> map, DofHandler d) : local_separated_dof_map(map), handler(d) {}

            MARS_INLINE_FUNCTION
            bool operator()(const Integer sfc) const {
                const Integer dof = handler.sfc_to_owned(sfc);
                if ((dof + 1) >= local_separated_dof_map.extent(0)) return false;
                const Integer pred_value = local_separated_dof_map(dof + 1) - local_separated_dof_map(dof);
                return (pred_value > 0);
            }
        }; */

        void compute_global_offset() {
            const int rank_size = num_ranks(get_dof_handler().get_context());

            ViewVectorType<Integer> global_dof_size_per_rank("global_dof_size_per_rank", rank_size);
            global_dof_offset = ViewVectorType<Integer>("global_dof_offset", rank_size + 1);

            const Integer global_size = get_owned_dof_size();
            get_dof_handler().get_context()->distributed->gather_all_view(global_size, global_dof_size_per_rank);

            incl_excl_scan(0, rank_size, global_dof_size_per_rank, global_dof_offset);
        }

        void prepare_separated_dofs() {
            local_dof_map = compact_local_dofs<Label>(get_dof_handler(), local_dofs);
            owned_dof_map = compact_owned_dofs<Label>(get_dof_handler(), locally_owned_dofs);

            compute_global_offset();

            auto is_separated = IsSeparatedDof(local_dof_map, get_dof_handler());
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

            ViewVectorType<Integer> boundary_dofs_index =
                ViewVectorType<Integer>("bounary_dofs_index", get_boundary_dof_size());
            build_boundary_dof_sets(boundary_dofs_index);

            ViewVectorType<Integer> ghost_dofs_index("ghost_index_dof", get_ghost_dof_size());
            get_dof_handler().get_context()->distributed->i_send_recv_view(ghost_dofs_index,
                                                                           get_view_scan_recv_mirror().data(),
                                                                           boundary_dofs_index,
                                                                           get_view_scan_send_mirror().data());

            std::cout << "Separated DofHandler:MPI send receive for the ghost dofs layer done." << std::endl;

            build_ghost_local_global_map(ghost_dofs_index);
        }

        MARS_INLINE_FUNCTION Integer get_dof_index(const Integer local_dof) const {
            return get_local_dof_map(local_dof);
        }

        MARS_INLINE_FUNCTION Integer get_owned_dof_index(const Integer owned_dof) const {
            return get_owned_dof_map(owned_dof);
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

        template <Integer FLabel, typename F>
        void owned_dof_iterate(F f) const {
            ViewVectorType<Integer> lowned_dofs;
            compact_owned_dofs<FLabel>(*this, lowned_dofs);

            Kokkos::parallel_for(
                "separated_dof_iter", lowned_dofs.extent(0), MARS_LAMBDA(const Integer i) { f(lowned_dofs(i)); });
        }

        template <typename F>
        void dof_iterate(F f) const {
            auto dofs = get_local_dofs();
            Kokkos::parallel_for(
                "separated_dof_iter", dofs.extent(0), MARS_LAMBDA(const Integer i) { f(dofs(i)); });
        }

        template <Integer FLabel, typename F>
        void dof_iterate(F f) const {
            ViewVectorType<Integer> dofs;
            compact_local_dofs<FLabel>(*this, dofs);

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
        const ViewVectorType<Integer> get_owned_dof_map() const { return owned_dof_map; }

        MARS_INLINE_FUNCTION
        const Integer get_owned_dof_map(const Integer local_dof) const { return owned_dof_map(local_dof); }

        MARS_INLINE_FUNCTION
        Integer get_boundary_dof(const Integer index) const { return boundary_dofs(index); }

        MARS_INLINE_FUNCTION
        Integer get_ghost_dof(const Integer index) const { return ghost_dofs(index); }

        MARS_INLINE_FUNCTION
        Integer get_boundary_dof_size() const { return boundary_dofs.extent(0); }

        MARS_INLINE_FUNCTION
        Integer get_ghost_dof_size() const { return ghost_dofs.extent(0); }

        MARS_INLINE_FUNCTION
        const Integer get_scan_recv_mirror_size() const { return scan_recv_mirror.extent(0); }

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

        /* *******************************local_to_global********************************** */

        MARS_INLINE_FUNCTION
        Integer locally_owned_dof(const Integer local_dof) const {
            const Integer owned = get_dof_handler().local_to_owned_dof(local_dof);
            const Integer pred_value = owned_dof_map(owned + 1) - owned_dof_map(owned);
            return (pred_value > 0);
        }

        MARS_INLINE_FUNCTION
        Integer local_to_owned_index(const Integer local_dof) const {
            const Integer owned = get_dof_handler().local_to_owned_dof(local_dof);
            const Integer pred_value = owned_dof_map(owned + 1) - owned_dof_map(owned);
            if (pred_value > 0)
                return owned_dof_map(owned);
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Integer local_to_owned_dof(const Integer local_dof) const {
            const Integer owned = get_dof_handler().local_to_owned_dof(local_dof);
            const Integer pred_value = owned_dof_map(owned + 1) - owned_dof_map(owned);
            if (pred_value > 0)
                return owned;
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        const Integer global_to_local(const Integer global_index) const {
            if (global_index == INVALID_INDEX) return INVALID_INDEX;

            const Integer proc = get_dof_handler().get_proc();
            const auto it = ghost_global_to_local_map.find(global_index);

            if (ghost_global_to_local_map.valid_at(it)) {
                return ghost_global_to_local_map.value_at(it);
            } else {
                const Integer owned_index = global_index - global_dof_offset(proc);
                return get_owned_dof(owned_index);
            }
        }

        MARS_INLINE_FUNCTION
        const Dof local_to_separated_global_dof(const Integer local_dof) const {
            Dof dof;
            const Integer proc = get_dof_handler().get_proc();
            const Integer owned = local_to_owned_dof(local_dof);

            if (owned != INVALID_INDEX) {
                const Integer index = owned_dof_map(owned) + global_dof_offset(proc);
                dof.set_gid(index);
                dof.set_proc(proc);
            } else {
                const Integer sfc = get_dof_handler().local_to_sfc(local_dof);
                const auto it = ghost_local_to_global_map.find(sfc);
                if (!ghost_local_to_global_map.valid_at(it)) {
                    dof.set_invalid();
                } else {
                    dof = ghost_local_to_global_map.value_at(it);
                }
            }
            return dof;
        }

        MARS_INLINE_FUNCTION
        Integer local_to_global_proc(const Integer local) const {
            Dof dof = local_to_separated_global_dof(local);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        const Integer local_to_global(const Integer local_dof) const {
            Dof dof = local_to_separated_global_dof(local_dof);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        const Integer get_global_dof_size() const {
            const Integer rank_size = num_ranks(get_dof_handler().get_context());

            auto ss = subview(global_dof_offset, rank_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            return h_ss();
        }

        void build_ghost_local_global_map(ViewVectorType<Integer> ghost_dofs_index) {
            auto handler = *this;

            int rank_size = num_ranks(handler.get_context());
            Integer size = get_view_scan_recv_mirror()(rank_size);

            ViewVectorType<Integer> scan_recv_proc("scan_recv_device", get_scan_recv_mirror_size());
            Kokkos::deep_copy(scan_recv_proc, get_view_scan_recv_mirror());

            ghost_local_to_global_map = UnorderedMap<Integer, Dof>(size);
            ghost_global_to_local_map = UnorderedMap<Integer, Integer>(size);
            auto glgm = ghost_local_to_global_map;
            auto gglm = ghost_global_to_local_map;
            /* iterate through the unique ghost dofs and build the map */
            Kokkos::parallel_for(
                "BuildLocalGlobalmap", size, MARS_LAMBDA(const Integer i) {
                    const Integer ghost_dof = handler.get_ghost_dof(i);
                    const Integer ghost_sfc = handler.local_to_sfc(ghost_dof);
                    const Integer gid = ghost_dofs_index(i);

                    /* find the process by binary search in the scan_recv_proc view of size rank_size
                     * and calculate the global id by adding the ghost local id to the global offset for
                     * ghost process*/
                    const Integer owner_proc = find_owner_processor(scan_recv_proc, i, 1, 0);

                    // build the ghost dof object and insert into the map
                    const auto result1 = glgm.insert(ghost_sfc, Dof(gid, owner_proc));
                    assert(result1.success());
                    const auto result2 = gglm.insert(gid, ghost_dof);
                    assert(result2.success());
                });

            /* In the end the size of the map should be as the size of the ghost_dofs.
             * Careful map size  is not capacity */
            assert(size == ghost_local_to_global_map.size());
        }

        // build boundary dof sets as global indexes so that the ghost indices to be global and poissible
        // therefor to go from a global index to a local one also for the ghosts.
        void build_boundary_dof_sets(ViewVectorType<Integer> boundary_dofs_index) {
            auto handler = *this;
            Kokkos::parallel_for(
                "boundary_iterate", get_boundary_dof_size(), MARS_LAMBDA(const Integer i) {
                    const Integer local_dof = handler.get_boundary_dof(i);
                    const Integer index = handler.local_to_owned_index(local_dof);
                    assert(index != INVALID_INDEX);
                    auto proc = handler.get_dof_handler().get_proc();
                    boundary_dofs_index(i) = index + handler.get_global_dof_offset(proc);
                });
        }

        /* *******dof handler related functionalities for completing the handler.******* */
        /* chose this way to hide the full interface of the general handler. Inheritance is the other way*/

        MARS_INLINE_FUNCTION Integer get_local_from_octant(const Octant &o) const {
            return get_dof_handler().get_local_from_octant(o);
        }

        MARS_INLINE_FUNCTION Octant get_octant_from_local(const Integer local) const {
            return get_dof_handler().get_octant_from_local(local);
        }

        MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const Integer sfc) const {
            return get_dof_handler().get_octant_from_sfc(sfc);
        }

        MARS_INLINE_FUNCTION void get_local_dof_coordinates(const Integer local, double *point) const {
            get_dof_handler().template get_dof_coordinates_from_local<ElemType>(local, point);
        }

        template <Integer Type, Integer FaceNr = -1>
        MARS_INLINE_FUNCTION bool is_boundary(const Integer local) const {
            return get_dof_handler().template is_boundary<Type, FaceNr>(local);
        }

        template <Integer FaceNr = -1>
        MARS_INLINE_FUNCTION bool is_boundary_dof(const Integer local) const {
            return is_boundary<ElemType, FaceNr>(local);
        }

        template <Integer face_nr = -1, typename F>
        void boundary_dof_iterate(F f) {
            auto handler = *this;
            Kokkos::parallel_for(
                "boundary_iterate", get_owned_dof_size(), MARS_LAMBDA(const Integer i) {
                    const Integer local = handler.get_owned_dof(i);
                    if (handler.template is_boundary_dof<face_nr>(local)) {
                        f(local);
                    }
                });
        }

        MARS_INLINE_FUNCTION
        constexpr Integer get_elem_type() { return simplex_type::ElemType; }

        MARS_INLINE_FUNCTION
        const Integer get_global_dof_offset(const Integer proc) const { return global_dof_offset(proc); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_global_dof_offset() const { return global_dof_offset; }

        /* MARS_INLINE_FUNCTION
        const Integer get_global_dof_size() const { return get_dof_handler().get_global_dof_size(); }
 */

        MARS_INLINE_FUNCTION
        const Integer get_orientation(const Integer local_dof) const {
            return get_dof_handler().get_orientation(local_dof);
        }

        MARS_INLINE_FUNCTION
        const Integer get_label(const Integer local_dof) const { return get_dof_handler().get_label(local_dof); }

        MARS_INLINE_FUNCTION
        const Integer get_owned_label(const Integer owned_dof) const {
            return get_dof_handler().get_owned_label(owned_dof);
        }

        /* MARS_INLINE_FUNCTION
        const Integer local_to_global(const Integer local) const { return get_dof_handler().local_to_global(local); } */

        MARS_INLINE_FUNCTION
        const Dof local_to_global_dof(const Integer local) const {
            return get_dof_handler().local_to_global_dof(local);
        }

        MARS_INLINE_FUNCTION
        UD get_data() const { return get_dof_handler().get_data(); }

        MARS_INLINE_FUNCTION
        const Integer get_proc() const { return get_dof_handler().get_proc(); }

        /* MARS_INLINE_FUNCTION
        const Integer local_to_owned_dof(const Integer local) const { return
        get_dof_handler().local_to_owned_dof(local); }
 */
        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_global_dof_enum() const {
            return get_dof_handler().get_global_dof_enum();
        }

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

        MARS_INLINE_FUNCTION
        const Integer get_XMax() const { return get_dof_handler().get_XMax(); }

        MARS_INLINE_FUNCTION
        const Integer get_YMax() const { return get_dof_handler().get_YMax(); }

        MARS_INLINE_FUNCTION
        const Integer get_ZMax() const { return get_dof_handler().get_ZMax(); }

        /* *************************************************************************** */

    private:
        DofHandler<Mesh, degree> dof_handler;

        // local dofs vector of locally owned dofs. Needed to build the stencils.
        ViewVectorType<Integer> locally_owned_dofs;
        ViewVectorType<Integer> owned_dof_map;

        UnorderedMap<Integer, Dof> ghost_local_to_global_map;
        UnorderedMap<Integer, Integer> ghost_global_to_local_map;
        ViewVectorType<Integer> global_dof_offset;

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
    using CornerVolumeDofHandler =
        SDofHandler<DofLabel::lVolume + DofLabel::lCorner, typename DofHandler::Mesh, DofHandler::Degree>;

    template <class DofHandler>
    using FaceVolumeDofHandler =
        SDofHandler<DofLabel::lVolume + DofLabel::lFace, typename DofHandler::Mesh, DofHandler::Degree>;

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
