#ifndef GENERATION_MARS_DISTRIBUTED_SDMAP_HPP_
#define GENERATION_MARS_DISTRIBUTED_SDMAP_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <Integer Label1, Integer Label2, typename DofHandler>
    class SDofMap {
    public:
        MARS_INLINE_FUNCTION
        SDofMap(DofHandler d) : dof_handler(d) { prepare_separated_dofs(); }

        ViewVectorType<bool> build_label_dof_predicate(const ViewVectorType<Integer> element_labels) {
            const Integer local_size = element_labels.extent(0);
            ViewVectorType<bool> dof_predicate("label_dof_predicate", local_size);
            Kokkos::parallel_for(
                "separatedoflabelss", local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (element_labels(i) == Label1 || element_labels(i) == Label2) {
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
            owned_dof_map = ViewVectorType<Integer>("owned_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, dof_predicate, owned_dof_map);

            auto vol_subview = subview(owned_dof_map, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            locally_owned_dofs = ViewVectorType<Integer>("locally_owned_dofs", h_vs());
            const ViewVectorType<Integer> global_to_sfc = get_dof_handler().get_global_dof_enum().get_view_elements();

            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (dof_predicate(i) == 1) {
                        Integer vindex = owned_dof_map(i);
                        const Integer local = get_dof_handler().sfc_to_local(global_to_sfc(i));
                        locally_owned_dofs(vindex) = local;
                    }
                });

            const int rank_size = num_ranks(get_dof_handler().get_context());

            ViewVectorType<Integer> global_dof_size_per_rank("global_dof_size_per_rank", rank_size);
            global_dof_offset = ViewVectorType<Integer>("global_dof_offset", rank_size + 1);

            const Integer global_size = get_owned_dof_size();
            get_dof_handler().get_context()->distributed->gather_all_view(global_size, global_dof_size_per_rank);

            incl_excl_scan(0, rank_size, global_dof_size_per_rank, global_dof_offset);
        }

        MARS_INLINE_FUNCTION Integer get_owned_dof_index(const Integer owned_dof) const {
            return get_owned_dof_map(owned_dof);
        }

        template <typename F>
        void owned_iterate(F f) const {
            Kokkos::parallel_for("owned_separated_dof_iter", get_owned_dof_size(), f);
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

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_owned_dofs() const { return locally_owned_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_owned_dof(const Integer i) const { return locally_owned_dofs(i); }

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
        const ViewVectorType<Integer> get_boundary_dofs() const { return boundary_dofs; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_ghost_dofs() const { return ghost_dofs; }


        MARS_INLINE_FUNCTION
        Integer get_owned_dof_size() const { return get_owned_dofs().extent(0); }

        MARS_INLINE_FUNCTION
        const Integer global_to_local(const Integer global_index) const {
            const Integer proc = get_dof_handler().get_proc();
            const auto it = ghost_global_to_local_map.find(global_index);

            if (!ghost_global_to_local_map.valid_at(it)) {
                const Integer owned_index = global_index - global_dof_offset(proc);
                return get_owned_dof(owned_index);
            } else {
                return ghost_global_to_local_map.value_at(it);
            }
        }

        MARS_INLINE_FUNCTION
        const Dof local_to_separated_global_dof(const Integer local_dof) const {
            Dof dof;
            const Integer proc = get_dof_handler().get_proc();
            const Integer owned = get_dof_handler().local_to_owned_dof(local_dof);

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

        void build_ghost_local_global_map(const ViewVectorType<Integer> &ghost_dofs_index) {
            auto handler = get_dof_handler();

            int rank_size = num_ranks(handler.get_context());
            Integer size = get_view_scan_recv_mirror()(rank_size);

            // TODO: when integrated into the staggered dof handler remove the handler because it is the handler itself.
            ViewVectorType<Integer> scan_recv_proc("scan_recv_device", get_scan_recv_mirror_size());
            Kokkos::deep_copy(scan_recv_proc, get_view_scan_recv_mirror());

            auto gdoffset = global_dof_offset;

            ghost_local_to_global_map = UnorderedMap<Integer, Dof>(size);
            ghost_global_to_local_map = UnorderedMap<Integer, Integer>(size);
            auto glgm = ghost_local_to_global_map;
            auto gglm = ghost_global_to_local_map;
            /* iterate through the unique ghost dofs and build the map */
            Kokkos::parallel_for(
                "BuildLocalGlobalmap", size, MARS_LAMBDA(const Integer i) {
                    const Integer ghost_dof = get_ghost_dof(i);
                    const Integer ghost_sfc = handler.local_to_sfc(ghost_dof);
                    const Integer gid = ghost_dofs_index(i);

                    /* find the process by binary search in the scan_recv_proc view of size rank_size
                     * and calculate the global id by adding the ghost local id to the global offset for
                     * ghost process*/
                    const Integer owner_proc = find_owner_processor(scan_recv_proc, i, 1, 0);
                    /* const Integer gid = ghost_gid + gdoffset(owner_proc); */

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


        //build boundary dof sets as global indexes so that the ghost indices to be global and poissible
        //therefor to go from a global index to a local one also for the ghosts.
        void build_boundary_dof_sets(ViewVectorType<Integer> &boundary_dofs_index) {
            auto handler = get_dof_handler();

            auto size = get_boundary_dof_size();
            boundary_dofs_index = ViewVectorType<Integer>("bounary_dofs_index", size);
            auto map = owned_dof_map;

            auto gdoffset = global_dof_offset;

            Kokkos::parallel_for(
                "boundary_iterate", size, MARS_LAMBDA(const Integer i) {
                    const Integer local_dof = get_boundary_dof(i);
                    const Integer owned_dof = handler.local_to_owned_dof(local_dof);
                    boundary_dofs_index(i) = map(owned_dof) + gdoffset(handler.get_proc());
                });
        }

        struct IsSeparatedOwnedDof {
            ViewVectorType<Integer> local_separated_dof_map;
            DofHandler handler;

            MARS_INLINE_FUNCTION
            IsSeparatedOwnedDof(ViewVectorType<Integer> map, DofHandler d) : local_separated_dof_map(map), handler(d) {}

            MARS_INLINE_FUNCTION
            bool operator()(const Integer sfc) const {
                const Integer dof = handler.sfc_to_owned(sfc);
                if ((dof + 1) >= local_separated_dof_map.extent(0)) return false;
                /*use the map which is the scan of the predicate.
                 * To get the predicate value the difference with the successive index is needed.*/
                const Integer pred_value = local_separated_dof_map(dof + 1) - local_separated_dof_map(dof);
                return (pred_value > 0);
            }
        };


        template <typename V>
        void compact_local_dofs(V &local_dof_map) {
            using namespace Kokkos;

            const Integer local_size = get_dof_handler().get_local_dof_enum().get_elem_size();
            auto dof_predicate =
                build_label_dof_predicate(get_dof_handler().get_local_dof_enum().get_view_element_labels());

            assert(local_size == get_dof_handler().get_local_dof_enum().get_view_element_labels().extent(0));

            /* perform a scan on the dof predicate*/
            local_dof_map = ViewVectorType<Integer>("local_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, dof_predicate, local_dof_map);
/*
            auto vol_subview = subview(local_dof_map, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            local_dofs = ViewVectorType<Integer>("local_dofs", h_vs());

            [>Compact the predicate into the volume and face dofs views<]
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (dof_predicate(i) == 1) {
                        Integer vindex = local_dof_map(i);
                        local_dofs(vindex) = i;
                    }
                }); */
        }

        struct IsSeparatedDof {
            ViewVectorType<Integer> local_separated_dof_map;
            DofHandler handler;

            MARS_INLINE_FUNCTION
            IsSeparatedDof(ViewVectorType<Integer> map, DofHandler d)
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

        void prepare_separated_dofs() {
            compact_owned_dofs(locally_owned_dofs);
            compact_local_dofs(local_dof_map);

            /* auto is_separated = IsSeparatedOwnedDof(owned_dof_map, get_dof_handler()); */
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

            ViewVectorType<Integer> boundary_dofs_index;
            build_boundary_dof_sets(boundary_dofs_index);
            auto owned_dofs = locally_owned_dofs;
            printf("boundary size: %li - %li\n", boundary_dofs_index.extent(0), get_boundary_dof_size());

            auto handler = get_dof_handler();
            auto gdoffset = global_dof_offset;

            Kokkos::parallel_for(
                "print", boundary_dofs_index.extent(0), MARS_LAMBDA(const Integer i) {
                    const Integer index = boundary_dofs_index(i) - gdoffset(handler.get_proc());
                    const Integer local = owned_dofs(index);
                    const Integer global = handler.local_to_global(local);
                    printf("i: %li, boundary_dofs_index(i) %li, global: %li\n", i, index, global);
                });
            auto ss = get_dof_handler().get_view_scan_send();
            auto sr = get_dof_handler().get_view_scan_recv();

            const Integer proc_size = num_ranks(get_dof_handler().get_context());
            Kokkos::parallel_for(
                "print", proc_size+1, MARS_LAMBDA(const Integer i) {
                    printf("device rank: %li, recv %li, send: %li\n", i, ghost_scan(i), boundary_scan(i));
                    printf("device orig.rank: %li, recv %li, send: %li\n", i, sr(i), ss(i));
                });

            for (int i = 0; i < proc_size; ++i) {
                std::cout << "proc: " << i << "recv(i): " << scan_recv_mirror(i) << "send(i): " << scan_send_mirror(i)
                          << std::endl;
            }
                ViewVectorType<Integer> ghost_dofs_index("ghost_index_dof", get_ghost_dof_size());
                get_dof_handler().get_context()->distributed->i_send_recv_view(ghost_dofs_index,
                                                                               get_view_scan_recv_mirror().data(),
                                                                               boundary_dofs_index,
                                                                               get_view_scan_send_mirror().data());

                std::cout << "Separated DofHandler:MPI send receive for the ghost dofs layer done." << std::endl;

                build_ghost_local_global_map(ghost_dofs_index);
        }

        MARS_INLINE_FUNCTION
        const DofHandler &get_dof_handler() const { return dof_handler; }

    private:
        DofHandler dof_handler;

        ViewVectorType<Integer> local_dof_map;

        ViewVectorType<Integer> locally_owned_dofs;
        ViewVectorType<Integer> owned_dof_map;

        // boundary and ghost local dofs predicated for the separated dm.
        ViewVectorType<Integer> boundary_dofs;
        ViewVectorType<Integer>::HostMirror scan_send_mirror;

        ViewVectorType<Integer> ghost_dofs;
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;

        UnorderedMap<Integer, Dof> ghost_local_to_global_map;
        UnorderedMap<Integer, Integer> ghost_global_to_local_map;
        ViewVectorType<Integer> global_dof_offset;
    };

    template <class DofHandler>
    using VFDofMap = SDofMap<DofLabel::lVolume, DofLabel::lFace, DofHandler>;

}  // namespace mars

#endif
#endif

#endif
