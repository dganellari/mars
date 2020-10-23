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

        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        static constexpr Integer volume_nodes = (degree - 1) * (degree - 1);
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;
        static constexpr Integer elem_nodes = (degree + 1) * (degree + 1);

        MARS_INLINE_FUNCTION
        FDDM(Mesh *mesh, const context &c) : DM(mesh, c) {}


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
        struct FaceOwnedDof{
            ViewVectorType<bool> predicate;
            ViewVectorType<Integer> dir;
            ViewVectorType<Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            FaceOwnedDof(ViewVectorType<bool> rp,
                               ViewVectorType<Integer> d,
                               ViewVectorType<Integer> l,
                               Integer p)
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

        template <bool Ghost, Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_dof_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh->octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>
                    (face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_nodes; j++) {
                    Integer dof_sfc = process_face_node<simplex_type::ElemType, dir>
                        (face_cornerA, j);
                    f(mesh, index, dof_sfc, owner_proc, dir);
                }
            }
        }

       struct IdentifyFaceDofs{
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = mesh->get_sfc(i);
                const Integer proc = mesh->get_proc();

                constexpr bool BoundaryIter = false;
                if (face_nodes > 0) {
                    FaceOwnedDof<BoundaryIter> fp =
                        FaceOwnedDof<BoundaryIter>(face_predicate, face_dir, sfc_to_local, proc);
                    face_iterate<BoundaryIter, 0>(sfc, mesh, i, fp);
                    face_iterate<BoundaryIter, 1>(sfc, mesh, i, fp);
                }

                /* if (volume_nodes > 0) {
                    volume_iterate<BoundaryIter>(
                        sfc, mesh, i, VolumeRankBoundary<BoundaryIter>(rank_boundary, sfc_to_locally_owned, map, proc));
                } */
                // TODO: 3D part
            }

            IdentifyFaceDofs(Mesh *m, ViewVectorType<bool> sp, ViewVectorType<Integer> sd, ViewVectorType<Integer> sl)
                : mesh(m), face_predicate(sp), face_dir(sd), sfc_to_local(sl) {}

            Mesh *mesh;
            ViewVectorType<bool> face_predicate;
            ViewVectorType<Integer> face_dir;
            ViewVectorType<Integer> sfc_to_local;
        };


        void build_locally_owned_face_dofs(ViewVectorType<Integer> &scan_boundary,
                                     ViewVectorType<Integer> &boundary_lsfc_index,
                                     ViewVectorType<Integer> &sender_ranks_scan,
                                     const Integer nbh_rank_size) {
            using namespace Kokkos;

            const Integer size = data.get_host_mesh()->get_chunk_size();

            Integer xDim = data.get_host_mesh()->get_XDim();
            Integer yDim = data.get_host_mesh()->get_YDim();
            Integer zDim = data.get_host_mesh()->get_ZDim();

            const Integer local_size = local_dof_enum.get_elem_size();

            ViewMatrixType<bool> face_dof_predicate("stencil_predicate", local_size);
            ViewMatrixType<Integer> face_dof_dir("stencil_predicate", local_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for(
                "identify_face_dofs",
                size,
                IdentifyFaceDofs(
                    data.get_mesh(), face_dof_predicate, face_dof_dir, local_dof_enum.get_view_sfc_to_local()));

            /* perform a scan for each row with the sum at the end for each rank */
            ViewVectorType<Integer> face_dof_scan("face_dof_scan", local_size + 1);
            incl_excl_scan(0, local_size, face_dof_predicate, face_dof_scan);

            auto index_subview = subview(face_dof_scan, local_size);
            auto h_ic = create_mirror_view(index_subview);
            // Deep copy device view to host view.
            deep_copy(h_ic, index_subview);

            locally_owned_face_dofs = ViewMatrixTypeRC<Integer, 2>("locally_owned_face_dofs", h_ic());
            /*   parallel_for(
                "print scan", rank_size + 1, KOKKOS_LAMBDA(const int i) {
                    printf(" scan boundary: %i-%li\n", i, scan_boundary(i));
                }); */

            ViewMatrixTypeRC<Integer, 2> lofd = locally_owned_face_dofs;;
            /* We use this strategy so that the compacted elements from the local_sfc
            would still be sorted and unique. */
            parallel_for(local_size,
                KOKKOS_LAMBDA(const Integer i) {
                    if (face_dof_predicate(i) == 1) {
                        Integer index = face_dof_scan(i);
                        lofd(index, 0) = i;
                        lofd(index, 1) = face_dof_dir(i);
                    }
                });

            /* parallel_for(
                "print set", h_ic(), KOKKOS_LAMBDA(const Integer i) {
                    Integer proc = data.get_host_mesh()->get_proc();
                    const Integer rank = find_owner_processor(scan_boundary, i, 1, proc);

                    printf("i:%li - boundary_ : %i - %li (%li) - proc: %li - rank: %li\n", i, boundary_dofs_sfc(i),
               boundary_lsfc_index(i), get_octant_from_sfc<simplex_type::ElemType>(boundary_dofs_sfc(i)).template
               get_global_index<simplex_type::ElemType>(xDim, yDim), rank, proc);
                }); */
        }

        void enumerate_dofs(const context &context) {
            this->enumerate_dofs(context);
            build_locally_owned_face_dofs();
        }


        MARS_INLINE_FUNCTION
        const ViewMatrixTypeRC<Integer, 2> get_locally_owned_face_dofs() const { return locally_owned_face_dofs; }

    private:

        // local enumeration of the dofs topologically foreach element
        ViewMatrixType<Integer> elem_dof_enum;
        Stencil<Dim, degree> stencil;
        ViewMatrixTypeRC<Integer, 2> locally_owned_face_dofs;
    };

}  // namespace mars

#endif
#endif

#endif
