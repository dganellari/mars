#ifndef GENERATION_MARS_DISTRIBUTED_VDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_VDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_base_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class VDM : public BDM<Mesh, degree, T...> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = BDM<Mesh, degree, T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        MARS_INLINE_FUNCTION
        VDM(DofHandler<Mesh, degree> d) : SuperDM(d) { prepare_volume_dofs(); }

        ViewVectorType<bool> build_volume_predicate(const Integer local_size,
                                                    const ViewVectorType<Integer> element_labels) {
            ViewVectorType<bool> volume_dof_predicate("volume_predicate", local_size);
            Kokkos::parallel_for(
                "separatedoflabelss", local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (element_labels(i) == DofLabel::lVolume) {
                        volume_dof_predicate(i) = 1;
                    }
                });

            return volume_dof_predicate;
        }

        void compact_owned_volume_dofs() {
            using namespace Kokkos;

            const Integer local_size = SuperDM::get_dof_handler().get_global_dof_enum().get_elem_size();
            auto volume_dof_predicate = build_volume_predicate(
                local_size, SuperDM::get_dof_handler().get_global_dof_enum().get_view_element_labels());

            /* perform a scan on the volume dof predicate*/
            ViewVectorType<Integer> owned_volume_dof_scan("owned_volume_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, volume_dof_predicate, owned_volume_dof_scan);

            auto vol_subview = subview(owned_volume_dof_scan, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            locally_owned_volume_dofs = ViewVectorType<Integer>("local_volume_dofs", h_vs());
            ViewVectorType<Integer> lovd = locally_owned_volume_dofs;

            const ViewVectorType<Integer> global_to_sfc =
                SuperDM::get_dof_handler().get_global_dof_enum().get_view_elements();
            auto dofhandler = *this;
            /* Compact the predicate into the volume and face dofs views */
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (volume_dof_predicate(i) == 1) {
                        Integer vindex = owned_volume_dof_scan(i);
                        const Integer local = dofhandler.sfc_to_local(global_to_sfc(i));
                        lovd(vindex) = local;
                    }
                });
        }

        void compact_volume_dofs() {
            using namespace Kokkos;

            const Integer local_size = SuperDM::get_dof_handler().get_local_dof_enum().get_elem_size();
            auto volume_dof_predicate = build_volume_predicate(
                local_size, SuperDM::get_dof_handler().get_local_dof_enum().get_view_element_labels());

            /* perform a scan on the volume dof predicate*/
            local_volume_dof_map = ViewVectorType<Integer>("volume_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, volume_dof_predicate, local_volume_dof_map);

            auto vol_subview = subview(local_volume_dof_map, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            local_volume_dofs = ViewVectorType<Integer>("local_volume_dofs", h_vs());
            ViewVectorType<Integer> lovd = local_volume_dofs;
            ViewVectorType<Integer> map = local_volume_dof_map;

            /* Compact the predicate into the volume and face dofs views */
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (volume_dof_predicate(i) == 1) {
                        Integer vindex = map(i);
                        lovd(vindex) = i;
                    }
                });
        }

        struct IsVolumeDof {
            ViewVectorType<Integer> local_volume_dof_map;

            MARS_INLINE_FUNCTION
            IsVolumeDof(ViewVectorType<Integer> map) : local_volume_dof_map(map) {}

            MARS_INLINE_FUNCTION
            bool operator()(const Integer local_dof) const {
                if ((local_dof + 1) >= local_volume_dof_map.extent(0)) return false;
                /*use the map which is the scan of the predicate.
                 * To get the predicate value the difference with the successive index is needed.*/
                const Integer pred_value = local_volume_dof_map(local_dof + 1) - local_volume_dof_map(local_dof);
                return (pred_value > 0);
            }
        };

        void prepare_volume_dofs() compact_volume_dofs();
        compact_owned_volume_dofs();

        reserve_user_data(vdata, "volume_user_data tuple", get_volume_dofs().extent(0));

        auto is_volume = IsVolumeDof(local_volume_dof_map);
        // building the counts for boundary and ghost separations to use for gather and scatter volume data only!
        auto boundary_predicate = compact_sfc_to_local(*this,
                                                       is_volume,
                                                       SuperDM::get_dof_handler().get_boundary_dofs(),
                                                       boundary_volume_dofs_sfc);
        auto boundary_scan = count_sfc_to_local(SuperDM::get_dof_handler().get_view_scan_send(), boundary_predicate);
        volume_scan_send_mirror = create_mirror_view(boundary_scan);
        Kokkos::deep_copy(volume_scan_send_mirror, boundary_scan);

        auto ghost_predicate =
            compact_sfc_to_local(*this, is_volume, SuperDM::get_dof_handler().get_ghost_dofs(), ghost_volume_dofs_sfc);
        auto ghost_scan = count_sfc_to_local(SuperDM::get_dof_handler().get_view_scan_recv(), ghost_predicate);
        volume_scan_recv_mirror = create_mirror_view(ghost_scan);
        Kokkos::deep_copy(volume_scan_recv_mirror, ghost_scan);
    }

    template <typename F>
    void owned_volume_dof_iterate(F f) const {
        Kokkos::parallel_for("volume_dof_iter", locally_owned_volume_dofs.extent(0), f);
    }

    template <typename F>
    void volume_dof_iterate(F f) const {
        Kokkos::parallel_for("volume_dof_iter", local_volume_dofs.extent(0), f);
    }

    /* building the stencil is the responsibility of the specialized DM. */
    template <typename ST>
    ST build_stencil() {
        return build_volume_stencil<ST>(*this);
    }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> get_volume_dofs() const { return local_volume_dofs; }

    MARS_INLINE_FUNCTION
    const Integer get_volume_dof(const Integer i) const { return local_volume_dofs(i); }

    MARS_INLINE_FUNCTION const Integer get_volume_dof_size() const { return local_volume_dofs.extent(0); }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> get_owned_volume_dofs() const { return locally_owned_volume_dofs; }

    MARS_INLINE_FUNCTION
    const Integer get_owned_volume_dof(const Integer i) const { return locally_owned_volume_dofs(i); }

    MARS_INLINE_FUNCTION const Integer get_owned_volume_dof_size() const { return locally_owned_volume_dofs.extent(0); }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> get_boundary_volume_dofs() const { return boundary_volume_dofs_sfc; }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> get_ghost_volume_dofs() const { return ghost_volume_dofs_sfc; }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer>::HostMirror &get_volume_scan_recv_mirror() const { return volume_scan_recv_mirror; }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer>::HostMirror &get_volume_scan_send_mirror() const { return volume_scan_send_mirror; }

    MARS_INLINE_FUNCTION
    const user_tuple &get_volume_data() const { return vdata; }

    template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
    MARS_INLINE_FUNCTION const H get_dof_data() const {
        return std::get<idx>(vdata);
    }

    template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
    MARS_INLINE_FUNCTION H &get_dof_data(const Integer i) const {
        return std::get<idx>(vdata)(i);
    }

    // gather operation: fill the data from the received ghost data
    template <Integer... dataidx>
    void gather_ghost_data(const context &context) {
        using namespace Kokkos;

        Kokkos::Timer timer;

        // exchange the ghost dofs first since it will be used to find the address
        // of the userdata based on the sfc code.

        int proc_num = rank(context);
        int size = num_ranks(context);

        Integer ghost_size = get_volume_scan_recv_mirror()(size);
        user_tuple ghost_user_data;
        reserve_user_data<dataidx...>(ghost_user_data, "ghost_user_data", ghost_size);

        // prepare the buffer to send the boundary data
        const Integer buffer_size = get_boundary_volume_dofs().extent(0);
        user_tuple buffer_data;
        reserve_user_data<dataidx...>(buffer_data, "buffer_data", buffer_size);

        /* TODO: create a view with local dofs corrensponding to the boundary_volume dofs sfc then do
        auto boundary_volume_local_dofs = translate_to_local(get_boundary_volume_dofs());
        fill_buffer_data<0, dataidx...>(buffer_data, boundary_volume_local_dofs, local_volume_dof_map); */

        //OR MAYBE USE DOFS instead of SFC for boundary and ghosts!!!!!
        /* fill_buffer_data<0, dataidx...>(buffer_data, get_boundary_volume_dofs(), local_volume_dof_map); */

        expand_tuple<ExchangeGhostDofsData, dataidx...>(
            ExchangeGhostDofsData(context, get_volume_scan_recv_mirror().data(), get_volume_scan_send_mirror().data()),
            ghost_user_data,
            buffer_data);

        // use the received ghost data and the sfc to put them to the unified local data
        fill_user_data<0, dataidx...>(ghost_user_data, get_ghost_volume_dofs());

        /* print_nth_tuple<1>(proc_num); */
    }

    template <Integer... dataidx>
    user_tuple scatter_ghost_data(const context &context) {
        using namespace Kokkos;

        Kokkos::Timer timer;

        // exchange the ghost dofs first since it will be used to find the address
        // of the userdata based on the sfc code.

        int proc_num = rank(context);
        int size = num_ranks(context);

        Integer ghost_size = get_volume_scan_recv_mirror()(size);
        user_tuple ghost_buffer_data;
        reserve_user_data<dataidx...>(ghost_buffer_data, "ghost_user_data", ghost_size);

        fill_user_data<1, dataidx...>(ghost_buffer_data, get_ghost_volume_dofs());

        const Integer boundary_size = get_boundary_volume_dofs().extent(0);
        user_tuple boundary_user_data;
        reserve_user_data<dataidx...>(boundary_user_data, "boundary_user_data", boundary_size);

        // prepare the buffer to send the boundary data
        expand_tuple<ExchangeGhostDofsData, dataidx...>(
            ExchangeGhostDofsData(context, get_volume_scan_send_mirror().data(), get_volume_scan_recv_mirror().data()),
            boundary_user_data,
            ghost_buffer_data);
        /* print_nth_tuple<1>(proc_num); */

        return boundary_user_data;
    }

private:
    // needed to build the stencils (only on the owned local dofs).
    ViewVectorType<Integer> locally_owned_volume_dofs;

    // needed to assign the data to each local dof (including ghosts)
    ViewVectorType<Integer> local_volume_dofs;
    ViewVectorType<Integer> local_volume_dof_map;
    // data assigned to each volume local dof
    user_tuple vdata;

    // boundary and ghost sfc predicated for the volume dm.
    ViewVectorType<Integer> boundary_volume_dofs_sfc;
    ViewVectorType<Integer>::HostMirror volume_scan_send_mirror;

    ViewVectorType<Integer> ghost_volume_dofs_sfc;
    ViewVectorType<Integer>::HostMirror volume_scan_recv_mirror;

};  // namespace mars

}  // namespace mars

#endif
#endif

#endif
