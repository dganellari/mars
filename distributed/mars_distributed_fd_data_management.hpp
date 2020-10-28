#ifndef GENERATION_MARS_DISTRIBUTED_FDDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_FDDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_data_management.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class FDDM : public DM<Mesh, degree, T...> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        using SuperDM = DM<Mesh, degree, T...>;

        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        FDDM(Mesh *mesh, const context &c) : DM<Mesh, degree, T...>(mesh, c) {}

        /* template <bool Ghost>
        struct VolumeRankBoundary {
            ViewMatrixType<bool> rank_boundary;
            ViewVectorType<Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            VolumeRankBoundary(ViewMatrixType<bool> rb,
                               ViewVectorType<Integer> sl,
                               ViewVectorType<Integer> m,
                               Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            MARS_INLINE_FUNCTION
            void volume_rank_boundary(const Mesh *mesh,
                                const Integer i,
                                const Integer sfc,
                                std::false_type) const {
                const Integer ghost_proc =
                    find_owner_processor(mesh->get_view_scan_boundary(), i, 1, proc);
                    Integer index = sfc_to_locally_owned(sfc);
                    rank_boundary(index, map(ghost_proc)) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh* mesh, const Integer i, const Integer dof_sfc) const {
                volume_rank_boundary(mesh, i, dof_sfc, std::integral_constant<bool, Ghost>{});
            }
        };
 */

        template <bool Ghost>
        struct VolumeOwnedDof {
            ViewVectorType<bool> predicate;
            ViewVectorType<Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            VolumeOwnedDof(ViewVectorType<bool> rp, ViewVectorType<Integer> l, Integer p)
                : predicate(rp), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            void volume_owned_dof(const Mesh *mesh, const Integer i, const Integer sfc, std::false_type) const {
                Integer index = sfc_to_local(sfc);
                predicate(index) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh, const Integer i, const Integer dof_sfc) const {
                volume_owned_dof(mesh, i, dof_sfc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct FaceOwnedDof {
            ViewVectorType<bool> predicate;
            ViewVectorType<Integer> dir;
            ViewVectorType<Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            FaceOwnedDof(ViewVectorType<bool> rp, ViewVectorType<Integer> d, ViewVectorType<Integer> l, Integer p)
                : predicate(rp), dir(d), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            void face_owned_dof(const Mesh *mesh,
                                const Integer i,
                                const Integer sfc,
                                const Integer owner_proc,
                                const Integer d,
                                std::false_type) const {
                if (proc >= owner_proc) {
                    Integer index = sfc_to_local(sfc);
                    predicate(index) = 1;
                    dir(index) = d;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                face_owned_dof(mesh, i, dof_sfc, owner_proc, dir, std::integral_constant<bool, Ghost>{});
            }
        };
        /*
                template <bool Ghost, Integer part, typename F>
                static MARS_INLINE_FUNCTION void face_dof_iterate(const Octant &oc) const {
                    [>Octant oc = mesh->get_octant(sfc_index);<]
                    // side  0 means origin side and 1 destination side.
                    for (int dir = 0; dir < 2; ++dir) {
                        Octant face_cornerA = enum_face_corner<part>(oc, dir);

                        for (int j = 0; j < face_nodes; j++) {
                            Integer dof_sfc = sfc_face_node<part, simplex_type::ElemType>(face_cornerA, j, dir);
                            f(dof_sfc, dir);
                        }
                    }
                } */

        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_dof_iterate(const Integer sfc,
                                                          const Mesh *mesh,
                                                          const Integer index,
                                                          F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh->octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc =
                    SuperDM::template process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_nodes; j++) {
                    Integer dof_sfc = SuperDM::template process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(mesh, index, dof_sfc, owner_proc, dir);
                }
            }
        }

        struct IdentifyFaceDofs {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = mesh->get_sfc(i);
                const Integer proc = mesh->get_proc();

                constexpr bool BoundaryIter = false;
                if (face_nodes > 0) {
                    FaceOwnedDof<BoundaryIter> fo =
                        FaceOwnedDof<BoundaryIter>(face_predicate, face_dir, sfc_to_local, proc);
                    face_dof_iterate<0>(sfc, mesh, i, fo);
                    face_dof_iterate<1>(sfc, mesh, i, fo);
                }

                if (volume_nodes > 0) {
                    SuperDM::template volume_iterate(
                        sfc, mesh, i, VolumeOwnedDof<BoundaryIter>(volume_predicate, sfc_to_local, proc));
                }
                // TODO: 3D part
            }

            IdentifyFaceDofs(Mesh *m, ViewVectorType<bool> sp, ViewVectorType<bool> vp, ViewVectorType<Integer> sd, ViewVectorType<Integer> sl)
                : mesh(m), face_predicate(sp), volume_predicate(vp), face_dir(sd), sfc_to_local(sl) {}

            Mesh *mesh;
            ViewVectorType<bool> face_predicate;
            ViewVectorType<bool> volume_predicate;
            ViewVectorType<Integer> face_dir;
            ViewVectorType<Integer> sfc_to_local;
        };

        void build_locally_owned_face_dofs() {
            using namespace Kokkos;

            const Integer size = SuperDM::get_data().get_host_mesh()->get_chunk_size();

            Integer xDim = SuperDM::get_data().get_host_mesh()->get_XDim();
            Integer yDim = SuperDM::get_data().get_host_mesh()->get_YDim();
            Integer zDim = SuperDM::get_data().get_host_mesh()->get_ZDim();

            const Integer local_size = SuperDM::get_local_dof_enum().get_elem_size();

            ViewVectorType<bool> volume_dof_predicate("volume_predicate", local_size);

            ViewVectorType<bool> face_dof_predicate("face_predicate", local_size);
            ViewVectorType<Integer> face_dof_dir("face_dir_predicate", local_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("identify_face_dofs",
                                 size,
                                 IdentifyFaceDofs(SuperDM::get_data().get_mesh(),
                                                  face_dof_predicate,
                                                  volume_dof_predicate,
                                                  face_dof_dir,
                                                  SuperDM::get_local_dof_enum().get_view_sfc_to_local()));

            /* perform a scan on the face_dof_predicate*/
            ViewVectorType<Integer> face_dof_scan("face_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, face_dof_predicate, face_dof_scan);

            auto face_subview = subview(face_dof_scan, local_size);
            auto h_fs = create_mirror_view(face_subview);
            // Deep copy device view to host view.
            deep_copy(h_fs, face_subview);

            locally_owned_face_dofs = ViewMatrixTypeRC<Integer, 2>("locally_owned_face_dofs", h_fs());

            /* perform a scan on the volume dof predicate*/
            ViewVectorType<Integer> volume_dof_scan("volume_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, volume_dof_predicate, volume_dof_scan);

            auto vol_subview = subview(volume_dof_scan, local_size);
            auto h_vs = create_mirror_view(vol_subview);
            // Deep copy device view to host view.
            deep_copy(h_vs, vol_subview);

            locally_owned_volume_dofs = ViewVectorType<Integer>("locally_owned_volume_dofs", h_vs());

            ViewMatrixTypeRC<Integer, 2> lofd = locally_owned_face_dofs;
            ViewVectorType<Integer> lovd = locally_owned_volume_dofs;

            /* Compact the predicate into the volume and face dofs views */
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (face_dof_predicate(i) == 1) {
                        Integer index = face_dof_scan(i);
                        lofd(index, 0) = i;
                        lofd(index, 1) = face_dof_dir(i);
                    }

                    if(volume_dof_predicate(i) == 1 ) {
                        Integer vindex = volume_dof_scan(i);
                        lovd(vindex) = i;
                    }
                });
        }

        virtual void enumerate_dofs(const context &context) override {
            SuperDM::enumerate_dofs(context);
            build_locally_owned_face_dofs();
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_locally_owned_volume_dofs() const { return locally_owned_volume_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_locally_owned_volume_dof(const Integer i) const { return locally_owned_volume_dofs(i); }

        MARS_INLINE_FUNCTION
        const ViewMatrixTypeRC<Integer, 2> get_locally_owned_face_dofs() const { return locally_owned_face_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_locally_owned_face_dof(const Integer i) const { return locally_owned_face_dofs(i, 0); }

        MARS_INLINE_FUNCTION
        const Integer get_locally_owned_face_dof_dir(const Integer i) const { return locally_owned_face_dofs(i, 1); }

        template <typename F>
        void volume_dof_iterate(F f) {
            Kokkos::parallel_for("volume_dof_iter", locally_owned_volume_dofs.extent(0), f);
        }

    private:
        /* Stencil<Dim, degree> stencil; */
        ViewVectorType<Integer> locally_owned_volume_dofs;
        ViewMatrixTypeRC<Integer, 2> locally_owned_face_dofs;
        //build corner only when degree = 1
        /* Stencil<Dim, degree, 1, false> corner_stencil; */
        //build if face node and volume node > 0
        /* Stencil<Dim, degree, 1, false> volume_stencil;
        Stencil<Dim, degree, 1, false> face_stencil; */

    };

}  // namespace mars

#endif
#endif

#endif
