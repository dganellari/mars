#ifndef GENERATION_MARS_DISTRIBUTED_FDDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_FDDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    /*
         * Face numbering on the stencil => ordering in the stencil stencil[1,0,3,2]
                ----3----
                |       |
                0   x   1
                |       |
                ----2---- */
    // building the stencil is the responsibility of the specialized DM.
    template <typename ST, typename DM>
    ST build_volume_stencil(const DM &dm) {
        ST vstencil(dm.get_owned_volume_dof_size());

        dm.owned_volume_dof_iterate(MARS_LAMBDA(const Integer i) {
            const Integer localid = dm.get_owned_volume_dof(i);
            vstencil.build_stencil(dm.get_dof_handler(), localid, i);
        });
        return vstencil;
    }

    template <typename ST, bool Orient = false, typename DM>
    ST build_face_stencil(const DM &dm) {
        ST fstencil(dm.get_owned_face_dof_size());

        dm.owned_face_dof_iterate(MARS_LAMBDA(const Integer i) {
            const Integer localid = dm.get_owned_face_dof(i);
            const Integer dir = dm.get_owned_face_dof_dir(i);
            fstencil.template build_stencil<Orient>(dm.get_dof_handler(), localid, i, dir);
        });
        return fstencil;
    }

    /* template <bool Orient, typename ST, typename DM>
    typename std::enable_if<Orient == true, ST>::type build_face_stencils(const DM &dm) {
        ST fstencil(dm.get_face_dof_size());

        dm.face_dof_iterate(MARS_LAMBDA(const Integer i) {
            const Integer localid = dm.get_face_dof(i);
            const Integer dir = dm.get_face_dof_dir(i);
            [>DM::fill_stencil(dm, fstencil, localid, i, dir);<]
            fstencil.build_stencil(dm, localid, i, dir);
        });
        return fstencil;
    }*/

    template <typename ST, typename DM>
    ST build_corner_stencil(const DM &dm) {
        ST cstencil(dm.get_owned_corner_dof_size());

        dm.owned_corner_dof_iterate(MARS_LAMBDA(const Integer i) {
            const Integer localid = dm.get_owned_corner_dof(i);
            cstencil.build_stencil(dm.get_dof_handler(), localid, i);
        });
        return cstencil;
    }

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

        static constexpr Integer Degree = degree;

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        FDDM(Mesh *mesh, const context &c) : SuperDM(mesh, c) {}

        template <bool Ghost>
        struct CornerOwnedDof {
            ViewVectorType<bool> predicate;
            ViewVectorType<Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            CornerOwnedDof(ViewVectorType<bool> rp, ViewVectorType<Integer> l, Integer p)
                : predicate(rp), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            void corner_owned_dof(const Integer sfc, const Integer max_proc, std::false_type) const {
                if (proc >= max_proc) {
                    Integer index = sfc_to_local(sfc);
                    predicate(index) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer max_proc,
                            Integer oro[simplex_type::ElemType]) const {
                corner_owned_dof(dof_sfc, max_proc, std::integral_constant<bool, Ghost>{});
            }
        };

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
            void face_owned_dof(const Integer sfc, const Integer owner_proc, const Integer d, std::false_type) const {
                if (proc >= owner_proc) {
                    Integer index = sfc_to_local(sfc);
                    predicate(index) = 1;
                    dir(index) = d;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Integer dof_sfc, const Integer owner_proc, const Integer dir) const {
                face_owned_dof(dof_sfc, owner_proc, dir, std::integral_constant<bool, Ghost>{});
            }
        };
/*
        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_dof_iterate(const Integer sfc, const Mesh *mesh, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh->octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc =
                    SuperDM::template process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_nodes; j++) {
                    Integer dof_sfc = SuperDM::template process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(dof_sfc, owner_proc, dir);
                }
            }
        } */

        struct SeperateDofs {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = mesh->get_sfc(i);
                const Integer proc = mesh->get_proc();

                constexpr bool BoundaryIter = false;
                if (face_nodes > 0) {
                    FaceOwnedDof<BoundaryIter> fo =
                        FaceOwnedDof<BoundaryIter>(face_predicate, face_dir, sfc_to_local, proc);
                    face_iterate<0>(sfc, mesh, fo);
                    face_iterate<1>(sfc, mesh, fo);
                }

                if (volume_nodes > 0) {
                    SuperDM::template volume_iterate(
                        sfc, mesh, i, VolumeOwnedDof<BoundaryIter>(volume_predicate, sfc_to_local, proc));
                }

                SuperDM::template corner_iterate(
                    sfc, mesh, i, CornerOwnedDof<BoundaryIter>(corner_predicate, sfc_to_local, proc));
                // TODO: 3D part
            }

            SeperateDofs(Mesh *m,
                         ViewVectorType<bool> sp,
                         ViewVectorType<bool> vp,
                         ViewVectorType<bool> cp,
                         ViewVectorType<Integer> sd,
                         ViewVectorType<Integer> sl)
                : mesh(m),
                  face_predicate(sp),
                  volume_predicate(vp),
                  corner_predicate(cp),
                  face_dir(sd),
                  sfc_to_local(sl) {}

            Mesh *mesh;
            ViewVectorType<bool> face_predicate;
            ViewVectorType<bool> volume_predicate;
            ViewVectorType<bool> corner_predicate;
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
            ViewVectorType<bool> corner_dof_predicate("corner_predicate", local_size);
            ViewVectorType<bool> face_dof_predicate("face_predicate", local_size);
            ViewVectorType<Integer> face_dof_dir("face_dir_predicate", local_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("separatedofs",
                                 size,
                                 SeperateDofs(SuperDM::get_data().get_mesh(),
                                              face_dof_predicate,
                                              volume_dof_predicate,
                                              corner_dof_predicate,
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

            /* perform a scan on the corner dof predicate*/
            ViewVectorType<Integer> corner_dof_scan("corner_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, corner_dof_predicate, corner_dof_scan);

            auto cor_subview = subview(corner_dof_scan, local_size);
            auto h_cs = create_mirror_view(cor_subview);
            // Deep copy device view to host view.
            deep_copy(h_cs, cor_subview);

            locally_owned_corner_dofs = ViewVectorType<Integer>("locally_owned_corner_dofs", h_cs());

            ViewMatrixTypeRC<Integer, 2> lofd = locally_owned_face_dofs;
            ViewVectorType<Integer> lovd = locally_owned_volume_dofs;
            ViewVectorType<Integer> locd = locally_owned_corner_dofs;

            /* Compact the predicate into the volume and face dofs views */
            parallel_for(
                local_size, KOKKOS_LAMBDA(const Integer i) {
                    if (face_dof_predicate(i) == 1) {
                        Integer index = face_dof_scan(i);
                        lofd(index, 0) = i;
                        lofd(index, 1) = face_dof_dir(i);
                    }

                    if (volume_dof_predicate(i) == 1) {
                        Integer vindex = volume_dof_scan(i);
                        lovd(vindex) = i;
                    }

                    if (corner_dof_predicate(i) == 1) {
                        Integer cindex = corner_dof_scan(i);
                        locd(cindex) = i;
                    }
                });
        }

        virtual void enumerate_dofs(const context &context) override {
            SuperDM::enumerate_dofs(context);
            build_locally_owned_face_dofs();
            //reserve TODO
        }

        template <typename Stencil>
        static MARS_INLINE_FUNCTION void fill_stencil(const SuperDM &dm,
                                                      Stencil stencil,
                                                      const Integer localid,
                                                      const Integer stencil_index,
                                                      const Integer Orientation = 1) {
            /* stencil(stencil_index, 0) = localid; */
            stencil.set_value(stencil_index, 0, localid);
            const Integer sfc = dm.local_to_sfc(localid);

            Octant oc = get_octant_from_sfc<simplex_type::ElemType>(sfc);

            Integer face_nr;
            for (int dir = 0; dir < 2; ++dir) {
                for (int side = 0; side < 2; ++side) {
                    if (side == 0)
                        face_nr = 2 * dir + 1;
                    else
                        face_nr = 2 * dir;

                    // this gives the index for different face orientation. Corner and volume have no
                    // extra orientation and the default  is 1. Orientation=0 -> x orientation).
                    const Integer dir_dim = !(Orientation ^ dir);
                    Integer index = 2 * dir_dim + side + 1;

                    Octant o = oc;
                    for (int w = 0; w < stencil.get_width(); ++w) {
                        /* const Integer nbh_sfc = dm.get_sfc_face_nbh(oc, face_nr); */

                        o = o.sfc_face_nbh<simplex_type::ElemType>(face_nr);
                        const Integer nbh_sfc = get_sfc_from_octant<simplex_type::ElemType>(o);
                        Integer nbh_id = dm.is_local(nbh_sfc) ? dm.sfc_to_local(nbh_sfc) : -1;
                        Integer offset = w * 2 * Dim;
                        /* stencil(stencil_index, index + offset) = nbh_id; */
                        stencil.set_value(stencil_index, index + offset, nbh_id);
                    }
                }
            }
        }

        template <typename F>
        void owned_volume_dof_iterate(F f) const {
            Kokkos::parallel_for("volume_dof_iter", locally_owned_volume_dofs.extent(0), f);
        }

        template <typename F>
        void owned_face_dof_iterate(F f) const {
            Kokkos::parallel_for("face_dof_iter", locally_owned_face_dofs.extent(0), f);
        }

        template <typename F>
        void owned_corner_dof_iterate(F f) const {
            Kokkos::parallel_for("corner_dof_iter", locally_owned_corner_dofs.extent(0), f);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_locally_owned_corner_dofs() const { return locally_owned_corner_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_owned_corner_dof(const Integer i) const { return locally_owned_corner_dofs(i); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_locally_owned_volume_dofs() const { return locally_owned_volume_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_owned_volume_dof(const Integer i) const { return locally_owned_volume_dofs(i); }

        MARS_INLINE_FUNCTION
        const ViewMatrixTypeRC<Integer, 2> get_locally_owned_face_dofs() const { return locally_owned_face_dofs; }

        MARS_INLINE_FUNCTION
        const Integer get_owned_face_dof(const Integer i) const { return locally_owned_face_dofs(i, 0); }

        MARS_INLINE_FUNCTION
        const Integer get_owned_face_dof_dir(const Integer i) const { return locally_owned_face_dofs(i, 1); }

        MARS_INLINE_FUNCTION const Integer get_owned_volume_dof_size() const { return locally_owned_volume_dofs.extent(0); }
        MARS_INLINE_FUNCTION const Integer get_owned_corner_dof_size() const { return locally_owned_corner_dofs.extent(0); }

        MARS_INLINE_FUNCTION const Integer get_owned_face_dof_size() const { return locally_owned_face_dofs.extent(0); }

    private:
        /* Stencil<Dim, degree> stencil; */
        ViewVectorType<Integer> locally_owned_corner_dofs;
        ViewVectorType<Integer> locally_owned_volume_dofs;
        ViewMatrixTypeRC<Integer, 2> locally_owned_face_dofs;
        /* user_tuple vdata; */
    };

}  // namespace mars

#endif
#endif

#endif
