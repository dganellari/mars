#ifndef GENERATION_MARS_DISTRIBUTED_FDDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_FDDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    // building the stencil is the responsibility of the specialized DM.
    template <typename ST, typename DM>
    ST build_volume_stencil(const DM &dm) {
        ST vstencil(dm.get_volume_dof_size());

        dm.volume_dof_iterate(MARS_LAMBDA(const Integer i) {
            const Integer localid = dm.get_volume_dof(i);
            vstencil.build_stencil(dm, localid, i);
        });
        return vstencil;
    }

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

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        VDM(Mesh *mesh, const context &c) : SuperDM(mesh, c) {}

        template <bool Ghost>
        struct VolumeOwnedDof {
            ViewVectorType<bool> predicate;
            ViewVectorType<bool> l_predicate;
            ViewVectorType<Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            VolumeOwnedDof(ViewVectorType<bool> rp, ViewVectorType<bool> lp, ViewVectorType<Integer> l, Integer p)
                : predicate(rp), g_predicate(lp), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            void volume_owned_dof(const Mesh *mesh, const Integer i, const Integer sfc, std::true_type) const {
                Integer index = sfc_to_local(sfc);
                l_predicate(index) = 1;
            }
            MARS_INLINE_FUNCTION
            void volume_owned_dof(const Mesh *mesh, const Integer i, const Integer sfc, std::false_type) const {
                Integer index = sfc_to_local(sfc);
                l_predicate(index) = 1;
                predicate(index) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh, const Integer i, const Integer dof_sfc) const {
                volume_owned_dof(mesh, i, dof_sfc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct SeperateVolumeDofs {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = SuperDM::get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh->get_proc();

                if (volume_nodes > 0) {
                    SuperDM::template volume_iterate(
                        sfc,
                        mesh,
                        i,
                        VolumeOwnedDof<Ghost>(volume_predicate, local_volume_predicate, sfc_to_local, proc));
                }

                // TODO: 3D part
            }

            SeperateVolumeDofs(Mesh *m,
                         ViewVectorType<bool> vp,
                         ViewVectorType<bool> lvp,
                         ViewVectorType<Integer> sl)
                : mesh(m),
                  volume_predicate(vp),
                  local_volume_predicate(lvp),
                  sfc_to_local(sl) {}

            Mesh *mesh;
            ViewVectorType<bool> volume_predicate;
            ViewVectorType<bool> local_volume_predicate;
            ViewVectorType<Integer> sfc_to_local;
        };

        void build_volume_dofs() {
            using namespace Kokkos;

            const Integer size = SuperDM::get_data().get_host_mesh()->get_chunk_size();

            Integer xDim = SuperDM::get_data().get_host_mesh()->get_XDim();
            Integer yDim = SuperDM::get_data().get_host_mesh()->get_YDim();
            Integer zDim = SuperDM::get_data().get_host_mesh()->get_ZDim();

            const Integer local_size = SuperDM::get_local_dof_enum().get_elem_size();

            ViewVectorType<bool> volume_dof_predicate("volume_predicate", local_size);
            ViewVectorType<bool> local_volume_dof_predicate("volume_predicate", local_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("separatedofs",
                                 size,
                                 SeperateVolumeDofs<false>(SuperDM::get_data().get_mesh(),
                                              volume_dof_predicate, local_volume_predicate,
                                              SuperDM::get_local_dof_enum().get_view_sfc_to_local()));

            Kokkos::parallel_for("separatedofs",
                                 size,
                                 SeperateVolumeDofs<true>(SuperDM::get_data().get_mesh(),
                                                    volume_dof_predicate, local_volume_predicate,
                                                    SuperDM::get_local_dof_enum().get_view_sfc_to_local()));

            /* perform a scan on the volume dof predicate*/
            ViewVectorType<Integer> volume_dof_scan("volume_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, volume_dof_predicate, volume_dof_scan);

            auto vol_subview = subview(volume_dof_scan, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            locally_owned_volume_dofs = ViewVectorType<Integer>("locally_owned_volume_dofs", h_vs());

            ViewVectorType<Integer> lovd = locally_owned_volume_dofs;

            /* Compact the predicate into the volume and face dofs views */
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (volume_dof_predicate(i) == 1) {
                        Integer vindex = volume_dof_scan(i);
                        lovd(vindex) = i;
                    }
                });
        }

        virtual void enumerate_dofs(const context &context) override {
            SuperDM::enumerate_dofs(context);
            build_volume_dofs();
            //reserve TODO
        }

        template <typename F>
        void volume_dof_iterate(F f) const {
            Kokkos::parallel_for("volume_dof_iter", locally_owned_volume_dofs.extent(0), f);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_locally_owned_volume_dofs() const { return locally_owned_volume_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_volume_dof(const Integer i) const { return locally_owned_volume_dofs(i); }

        MARS_INLINE_FUNCTION const Integer get_volume_dof_size() const { return locally_owned_volume_dofs.extent(0); }

    private:
        ViewVectorType<Integer> locally_owned_volume_dofs;
        ViewVectorType<Integer> local_volume_dofs;
        user_tuple vdata;
   };

}  // namespace mars

#endif
#endif

#endif
