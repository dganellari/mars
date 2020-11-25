#ifndef GENERATION_MARS_DISTRIBUTED_VDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_VDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class VDM : public DOFM<Mesh, degree> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        using SuperDM = DOFM<Mesh, degree>;
        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;


        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        static constexpr Integer Degree = degree;

        MARS_INLINE_FUNCTION
        VDM(Mesh *mesh, const context &c) : SuperDM(mesh, c) {}

        ViewVectorType<bool> build_volume_predicate(const Integer local_size, const ViewVectorType<Integer> element_labels) {
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

            const Integer local_size = SuperDM::get_global_dof_enum().get_elem_size();
            auto volume_dof_predicate =
                build_volume_predicate(local_size, SuperDM::get_global_dof_enum().get_view_element_labels());

            /* perform a scan on the volume dof predicate*/
            ViewVectorType<Integer> owned_volume_dof_scan("owned_volume_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, volume_dof_predicate, owned_volume_dof_scan);

            auto vol_subview = subview(owned_volume_dof_scan, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            locally_owned_volume_dofs = ViewVectorType<Integer>("local_volume_dofs", h_vs());
            ViewVectorType<Integer> lovd = locally_owned_volume_dofs;

            const ViewVectorType<Integer> global_to_sfc = SuperDM::get_global_dof_enum().get_view_elements();
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

            const Integer local_size = SuperDM::get_local_dof_enum().get_elem_size();
            auto volume_dof_predicate =
                build_volume_predicate(local_size, SuperDM::get_local_dof_enum().get_view_element_labels());

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

        virtual void enumerate_dofs(const context &context) override {
            SuperDM::enumerate_dofs(context);
            compact_volume_dofs();
            compact_owned_volume_dofs();
            //reserve TODO
            auto is_volume = IsVolumeDof(local_volume_dof_map);
            compact_sfc_to_local(*this, is_volume, SuperDM::get_boundary_dofs(), boundary_volume_dofs_sfc);
            compact_sfc_to_local(*this, is_volume, SuperDM::get_ghost_dofs(), ghost_volume_dofs_sfc);
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

    private:
        //needed to build the stencils (only on the owned local dofs).
        ViewVectorType<Integer> locally_owned_volume_dofs;

        //needed to assign the data to each local dof (including ghosts)
        ViewVectorType<Integer> local_volume_dofs;
        ViewVectorType<Integer> local_volume_dof_map;

        //boundary and ghost sfc predicated for the volume dm.
        ViewVectorType<Integer> boundary_volume_dofs_sfc;
        ViewVectorType<Integer> ghost_volume_dofs_sfc;

        user_tuple vdata;
   };

}  // namespace mars

#endif
#endif

#endif
