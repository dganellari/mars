#ifndef GENERATION_MARS_DISTRIBUTED_RDofMap_HPP_
#define GENERATION_MARS_DISTRIBUTED_RDofMap_HPP_

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_dof.hpp"

namespace mars {

    template <typename DofHandler>
    class RemoteDofMap {
    public:
        MARS_INLINE_FUNCTION
        RemoteDofMap(DofHandler handler) : dof_handler(handler) {}

        static_assert(DofHandler::Block == 1, "RemoteDofMap does not support yet vector valued block structures.");

        template <typename V>
        void build_global_dof_per_rank(const V &view,
                                       V &global_dof,
                                       V &global_sfc,
                                       const V &scan_recv_view,
                                       const Integer rank_size) {
            using namespace Kokkos;

            const Integer size = view.extent(0);
            ViewMatrixTypeLeft<bool> rank_view("count_per_proc", size, rank_size);

            auto handler = dof_handler;
            // build predicate
            Kokkos::parallel_for(
                "build global dofs", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = view(i);
                    if (handler.is_local(sfc)) {
                        const Integer proc = find_owner_processor(scan_recv_view, i, 1, 0);
                        rank_view(i, proc) = 1;
                    }
                });

            ViewMatrixTypeLeft<Integer> rank_scan("rank_scan", size + 1, rank_size);
            for (int i = 0; i < rank_size; ++i) {
                auto subpredicate = subview(rank_view, ALL, i);
                auto subscan = subview(rank_scan, ALL, i);
                incl_excl_scan(0, size, subpredicate, subscan);
            }

            scan_view_rank = V("scan_missing", rank_size + 1);
            auto scan_view = scan_view_rank;
            // perform a scan on the last column to get the total sum.
            row_scan(rank_size, size, rank_scan, scan_view);

            auto index_subview = subview(scan_view, rank_size);
            auto h_ic = create_mirror_view(index_subview);
            // Deep copy device view to host view.
            deep_copy(h_ic, index_subview);

            global_dof = V("global_dofs", h_ic());
            global_sfc = V("global_sfcs", h_ic());

            parallel_for(
                MDRangePolicy<Rank<2>>({0, 0}, {size, rank_size}), KOKKOS_LAMBDA(const Integer i, const Integer j) {
                    if (rank_view(i, j) == 1) {
                        Integer index = scan_view(j) + rank_scan(i, j);

                        const Integer sfc = view(i);
                        const Integer local = handler.sfc_to_local(sfc);
                        global_dof(index) = handler.local_to_global(local);
                        global_sfc(index) = sfc;
                    }
                });
        }

        void exchange_ghost_counts(const mars::context &context,
                                   std::vector<mars::Integer> &send_count,
                                   std::vector<mars::Integer> &receive_count)

        {
            using namespace Kokkos;
            using namespace mars;

            Kokkos::Timer timer;

            int proc_num = rank(context);
            int rank_size = num_ranks(context);

            auto scan_snd = create_mirror_view(scan_view_rank);
            Kokkos::deep_copy(scan_snd, scan_view_rank);

            /* Integer proc_count = 0; */
            for (int i = 0; i < rank_size; ++i) {
                Integer count = scan_snd(i + 1) - scan_snd(i);
                if (count > 0) {
                    send_count[i] = count;
                }
            }

            context->distributed->i_send_recv_all_to_all(send_count, receive_count);

            /* for (int i = 0; i < rank_size; ++i) {
                if (receive_count[i] > 0) {
                    std::cout << "Second count:-----FromProc: " << i << " count:" << receive_count[i]
                              << " Proc: " << proc_num << std::endl;
                }
            } */
        }

        template <typename V>
        void create_scan_mirrors(const mars::context &context, const V &send_count, const V &receive_count) {
            using namespace Kokkos;

            int proc_num = rank(context);
            int rank_size = num_ranks(context);

            /* create the scan recv mirror view from the receive count */
            scan_recv = ViewVectorType<Integer>("scan_recv_", rank_size + 1);
            scan_recv_mirror = create_mirror_view(scan_recv);
            make_scan_index_mirror(scan_recv_mirror, receive_count);
            Kokkos::deep_copy(scan_recv, scan_recv_mirror);

            /*create the scan send mirror view from the send count*/
            scan_send = ViewVectorType<Integer>("scan_send_", rank_size + 1);
            scan_send_mirror = create_mirror_view(scan_send);
            make_scan_index_mirror(scan_send_mirror, send_count);
            Kokkos::deep_copy(scan_send, scan_send_mirror);
            /*
                        std::cout << "Missing Send_scan: " << std::endl;
                        for (int i = 0; i < rank_size + 1; ++i) {
                            std::cout << scan_send_mirror(i) << " ";
                        }
                        std::cout << std::endl;

                        std::cout << "Missing recv_scan:" << std::endl;
                        for (int i = 0; i < rank_size + 1; ++i) {
                            std::cout << scan_recv_mirror(i) << " ";
                        }
                        std::cout << std::endl; */
        }

        template <typename VW>
        void exchange_ghost_globals(const mars::context &context, const VW &owned_dofs, const VW &owned_sfcs) {
            using namespace Kokkos;

            int proc_num = rank(context);
            int rank_size = num_ranks(context);

            Integer ghost_size = scan_recv_mirror(rank_size);
            global_ghost_dofs = ViewVectorType<Integer>("ghost_dofs", ghost_size);
            global_ghost_sfcs = ViewVectorType<Integer>("ghost_dofs", ghost_size);

            context->distributed->i_send_recv_view(
                global_ghost_dofs, scan_recv_mirror.data(), owned_dofs, scan_send_mirror.data());

            context->distributed->i_send_recv_view(
                global_ghost_sfcs, scan_recv_mirror.data(), owned_sfcs, scan_send_mirror.data());
        }

        template <typename V>
        void build_ghost_global_map(const mars::Integer rank_size, const V &ghost_dofs, const V &ghost_sfcs) {
            using namespace mars;
            using namespace Kokkos;

            Integer size = scan_recv_mirror(rank_size);

            ghost_sfc_to_global_map = UnorderedMap<Integer, Dof>(size);
            auto scan_recv_proc = scan_recv;
            auto gsgm = ghost_sfc_to_global_map;

            /* iterate through the unique ghost dofs and build the map */
            Kokkos::parallel_for(
                "BuildLocalGlobalmap", size, MARS_LAMBDA(const Integer i) {
                    /* read the sfc and the local dof id from the ghost */
                    const Integer gid = ghost_dofs(i);
                    const Integer ghost_sfc = ghost_sfcs(i);
                    /* find the process by binary search in the scan_recv_proc view of size rank_size
                    and calculate the global id by adding the ghost local id to the global offset for ghost process */
                    const Integer owner_proc = find_owner_processor(scan_recv_proc, i, 1, 0);

                    // build the ghost dof object and insert into the map
                    const auto result = gsgm.insert(ghost_sfc, Dof(gid, owner_proc));
                    assert(!result.failed());
                });
        }

        void build_remote_dof_map(ViewVectorType<Integer> missing_dofs_sfc) {
            using namespace Kokkos;

            const auto &context = dof_handler.get_context();

            int proc_num = rank(context);
            int rank_size = num_ranks(context);

            auto missing_dofs_sfc_size = missing_dofs_sfc.extent(0);

            std::vector<Integer> send_count(rank_size, 0);
            std::vector<Integer> receive_count(rank_size, 0);

            for (int i = 0; i < rank_size; ++i) {
                if (i != proc_num) send_count[i] = missing_dofs_sfc_size;
            }

            context->distributed->i_send_recv_all_to_all(send_count, receive_count);

            /* create the scan recv mirror view from the receive count */
            auto scan_rcv = ViewVectorType<Integer>("scan_rcv_", rank_size + 1);
            auto scan_rcv_mirror = create_mirror_view(scan_rcv);
            make_scan_index_mirror(scan_rcv_mirror, receive_count);
            Kokkos::deep_copy(scan_rcv, scan_rcv_mirror);

            /*create the scan send mirror view from the send count*/
            auto scn_send = ViewVectorType<Integer>("scn_send_", rank_size + 1);
            auto scn_send_mirror = create_mirror_view(scn_send);
            make_scan_index_mirror(scn_send_mirror, send_count);
            Kokkos::deep_copy(scn_send, scn_send_mirror);

            Integer ghost_size = scan_rcv_mirror(rank_size);
            auto ghost_dofs_sfc = ViewVectorType<Integer>("ghost_dofs", ghost_size);

            context->distributed->i_send_recv_view_to_all(
                ghost_dofs_sfc, scan_rcv_mirror.data(), missing_dofs_sfc, scn_send_mirror.data());

            ViewVectorType<Integer> global_dofs, global_sfcs;
            build_global_dof_per_rank(ghost_dofs_sfc, global_dofs, global_sfcs, scan_rcv, rank_size);

            Kokkos::fence();

            std::vector<Integer> snd_count(rank_size, 0);
            std::vector<Integer> rcv_count(rank_size, 0);
            exchange_ghost_counts(context, snd_count, rcv_count);

            create_scan_mirrors(context, snd_count, rcv_count);
            exchange_ghost_globals(context, global_dofs, global_sfcs);

            build_ghost_global_map(rank_size, global_ghost_dofs, global_ghost_sfcs);
            std::cout << "Remote DOf Map: Ending mpi send receive for the ghost sfc dofs " << std::endl;
        }

        MARS_INLINE_FUNCTION
        const Dof sfc_to_global_dof(const Integer sfc) const {
            Dof dof;
            const auto it = ghost_sfc_to_global_map.find(sfc);
            if (ghost_sfc_to_global_map.valid_at(it)) {
                dof = ghost_sfc_to_global_map.value_at(it);
            } else {
                dof.set_invalid();
            }
            return dof;
        }

        MARS_INLINE_FUNCTION
        const Integer sfc_to_global(const Integer sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        const Integer sfc_to_proc(const Integer sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

    private:
        DofHandler dof_handler;

        // remote sfc to global_map.
        UnorderedMap<Integer, Dof> ghost_sfc_to_global_map;

        ViewVectorType<Integer> scan_view_rank;
        // mirror view on the dof scan boundary view used as offset for the mpi send dofs.
        ViewVectorType<Integer> scan_send;
        ViewVectorType<Integer>::HostMirror scan_send_mirror;

        // ghost dofs received from other procs's boundary sfcs.
        ViewVectorType<Integer> global_ghost_dofs;
        ViewVectorType<Integer> global_ghost_sfcs;
        // mirror view on the scan_ghost view used as an offset to receive ghost dofs.
        ViewVectorType<Integer> scan_recv;
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;
    };

}  // namespace mars

#endif
#endif

#endif
