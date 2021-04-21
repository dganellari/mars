#ifndef GENERATION_MARS_DISTRIBUTED_DofHandler_HPP_
#define GENERATION_MARS_DISTRIBUTED_DofHandler_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof.hpp"
#include "mars_distributed_mesh_management.hpp"

namespace mars {

    template <Integer degree, Integer Type>
    struct NumDofs {};

    template <Integer degree>
    struct NumDofs<degree, ElementType::Quad4> {
        static constexpr Integer corner_dofs = 1;
        static constexpr Integer face_dofs = (degree - 1);
        static constexpr Integer volume_dofs = (degree - 1) * (degree - 1);
        static constexpr Integer edge_dofs = 0;
        static constexpr Integer elem_dofs = (degree + 1) * (degree + 1);
    };

    template <Integer degree>
    struct NumDofs<degree, ElementType::Hex8> {
        static constexpr Integer corner_dofs = 1;
        static constexpr Integer face_dofs = (degree - 1) * (degree - 1);
        static constexpr Integer volume_dofs = (degree - 1) * (degree - 1) * (degree - 1);
        static constexpr Integer edge_dofs = (degree - 1);
        static constexpr Integer elem_dofs = (degree + 1) * (degree + 1) * (degree + 1);
    };

    template <class Mesh_, Integer degree>
    class DofHandler {
    public:
        using Mesh = Mesh_;

        using MM = MeshManager<Mesh>;
        using simplex_type = typename Mesh::Elem;

        static constexpr Integer ElemType = simplex_type::ElemType;

        static constexpr Integer dofLabel = DofLabel::lAll;

        static constexpr Integer Degree = degree;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        using NDofs = NumDofs<degree, ElemType>;

        static constexpr Integer corner_dofs = NDofs::corner_dofs;
        static constexpr Integer face_dofs = NDofs::face_dofs;
        static constexpr Integer volume_dofs = NDofs::volume_dofs;
        static constexpr Integer edge_dofs = NDofs::edge_dofs;
        static constexpr Integer elem_dofs = NDofs::elem_dofs;

        MARS_INLINE_FUNCTION
        constexpr Integer get_elem_type() { return simplex_type::ElemType; }

        MARS_INLINE_FUNCTION
        DofHandler(Mesh *mesh, const context &c) : mesh_manager(MM(mesh)), ctx(c) {}

        template <typename H>
        MARS_INLINE_FUNCTION void owned_iterate(H f) const {
            owned_dof_iterate(f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void owned_dof_iterate(H f) const {
            const Integer size = global_dof_enum.get_elem_size();
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void dof_iterate(H f) const {
            const Integer size = local_dof_enum.get_elem_size();
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void elem_iterate(H f) const {
            get_mesh_manager().elem_iterate(f);
        }

        MARS_INLINE_FUNCTION Integer get_elem_size() const {
            return get_mesh_manager().get_host_mesh()->get_chunk_size();
        }

        template <typename H>
        MARS_INLINE_FUNCTION void face_iterate(H f) const {
            get_mesh_manager().face_iterate(f);
        }

        MARS_INLINE_FUNCTION
        Integer get_sfc_face_nbh(const Octant &oc, const Integer face_nr) const {
            Octant o = oc.sfc_face_nbh<simplex_type::ElemType>(face_nr);
            return get_sfc_from_octant(o);
        }

        template <Integer Type>
        static MARS_INLINE_FUNCTION Integer
        enum_corner(const ViewVectorType<Integer> &sfc_to_local, const Octant &oc, const int i, const int j) {
            Octant o;
            // get the next corner using the elem sfc
            o.x = oc.x + i;
            o.y = oc.y + j;
            // convert the octant value into the new nodal sfc system
            o.x *= degree;
            o.y *= degree;
            Integer sfc = mars::get_sfc_from_octant<Type>(o);

            return sfc_to_local(sfc);
        }

        template <Integer part>
        static MARS_INLINE_FUNCTION Octant enum_face_corner(Octant &oc, const int dir) {
            Octant face_cornerA;
            /*This is how the face nr is computed to get a topological order of the faces
            Integer face_nr;
            if (side == 0)
                face_nr = 2 + dir;
            else
                face_nr = 1 - dir; */

            /* find the starting corner "face_cornerA" of the face and use 0 as lower right part and
            1 as upper right part of the quad to order topologically:dir 1 y and dir 0 x direction. */
            if (part == 0) {
                face_cornerA.x = oc.x + dir;
                face_cornerA.y = oc.y;
            }
            if (part == 1) {
                face_cornerA.x = oc.x + (dir ^ 1);
                face_cornerA.y = oc.y + 1;
            }

            return face_cornerA;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION Integer sfc_face_node(const Octant &face_cornerA, const int j, const int dir) {
            int sign = (part == 1) ? -1 : 1;
            // dir is 1 for y direction and 0 for x direction
            Octant o;
            if (dir == 1) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + sign * (j + 1);
            }
            if (dir == 0) {
                o.x = degree * face_cornerA.x + sign * (j + 1);
                o.y = degree * face_cornerA.y;
            }
            Integer sfc = mars::get_sfc_from_octant<Type>(o);

            return sfc;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION Integer enum_face_node(const ViewVectorType<Integer> &sfc_to_local,
                                                           const Octant &face_cornerA,
                                                           const int j,
                                                           const int dir) {
            return sfc_to_local(sfc_face_node<part, Type>(face_cornerA, j, dir));
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, void>
        build_starting_face_corner(Octant &face_cornerA, const int side, const Octant &oc) {
            // find the starting corner "face_cornerA" of the face and use the direction
            const int val = side ^ 1;
            if (dir == 0) {
                face_cornerA.x = oc.x + val;
                face_cornerA.y = oc.y;
                face_cornerA.z = oc.z;
            }
            if (dir == 1) {
                face_cornerA.x = oc.x;
                face_cornerA.y = oc.y + val;
                face_cornerA.z = oc.z;
            }
            if (dir == 2) {
                face_cornerA.x = oc.x;
                face_cornerA.y = oc.y;
                face_cornerA.z = oc.z + val;
            }
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, void>
        build_starting_face_corner(Octant &face_cornerA, const int side, const Octant &oc) {
            // find the starting corner "face_cornerA" of the face and use the direction
            const int val = side ^ 1;
            if (dir == 0) {
                face_cornerA.x = oc.x + val;
                face_cornerA.y = oc.y;
            }
            if (dir == 1) {
                face_cornerA.x = oc.x;
                face_cornerA.y = oc.y + val;
            }
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION Integer
        process_face_corner(Octant &face_cornerA, const Mesh *mesh, const int side, const Octant &oc) {
            Integer face_nr;
            if (side == 0)
                face_nr = 2 * dir + 1;
            else
                face_nr = 2 * dir;

            Octant nbh_oc = mesh->get_octant_face_nbh(oc, face_nr);
            Integer owner_proc = -1;

            // find owner proc method returns an invalid result with invalid octant.
            if (nbh_oc.is_valid()) {
                Integer enc_oc = mars::get_sfc_from_octant<Type>(nbh_oc);
                owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, mesh->get_proc());
            }

            build_starting_face_corner<Type, dir>(face_cornerA, side, oc);

            return owner_proc;
        }

        template <Integer Type, Integer dir, typename F>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Integer>
        build_face_node(F f, const Octant &face_cornerA, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
            }
            return f(o);
        }

        template <Integer Type, Integer dir, typename F>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer>
        build_face_node(F f, const Octant &face_cornerA, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
                o.z = degree * face_cornerA.z + j + 1;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
                o.z = degree * face_cornerA.z + j + 1;
            }
            if (dir == 2) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y + j + 1;
                o.z = degree * face_cornerA.z;
            }
            return f(o);
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION Integer process_face_node(const Octant &face_cornerA, const int j) {
            return build_face_node<Type, dir>(mars::get_sfc_from_octant<Type>, face_cornerA, j);
        }

        template <bool Ghost, Integer Label = DofLabel::lAll>
        struct RankBoundary {
            ViewMatrixType<bool> rank_boundary;
            ViewVectorType<Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            RankBoundary(ViewMatrixType<bool> rb, ViewVectorType<Integer> sl, ViewVectorType<Integer> m, Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            MARS_INLINE_FUNCTION
            void all_rank_boundary(const Mesh *mesh,
                                    const Integer i,
                                    const Integer sfc,
                                    const Integer owner_proc,
                                    std::false_type) const {
                const Integer ghost_proc = find_owner_processor(mesh->get_view_scan_boundary(), i, 1, proc);

                /* if (proc >= owner_proc && owner_proc >= 0) { */
                if (proc >= owner_proc) {
                    Integer index = sfc_to_locally_owned(sfc);
                    rank_boundary(index, map(ghost_proc)) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                all_rank_boundary(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct RankBoundary<Ghost, DofLabel::lVolume> {
            ViewMatrixType<bool> rank_boundary;
            ViewVectorType<Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            RankBoundary(ViewMatrixType<bool> rb,
                               ViewVectorType<Integer> sl,
                               ViewVectorType<Integer> m,
                               Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            MARS_INLINE_FUNCTION
            void volume_rank_boundary(const Mesh *mesh, const Integer i, const Integer sfc, std::false_type) const {
                const Integer ghost_proc = find_owner_processor(mesh->get_view_scan_boundary(), i, 1, proc);
                Integer index = sfc_to_locally_owned(sfc);
                rank_boundary(index, map(ghost_proc)) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh, const Integer i, const Integer dof_sfc) const {
                volume_rank_boundary(mesh, i, dof_sfc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <Integer Type, typename F>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer>
        process_edge_node(F f, const Octant &face_cornerA, const Integer dir, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
                o.z = degree * face_cornerA.z;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
                o.z = degree * face_cornerA.z;
            }
            if (dir == 2) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y;
                o.z = degree * face_cornerA.z + j + 1;
            }
            return f(o);
        }

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T == ElementType::Hex8, void> edge_iterate(const Integer sfc,
                                                                                                const Mesh *mesh,
                                                                                                const Integer i,
                                                                                                F f) {
            Octant oc = mesh->octant_from_sfc(sfc);
            for (int edge = 0; edge < 4 * ManifoldDim; ++edge) {
                // get the starting node of the edge to be used to compute the dofs contained in that edge.
                const Integer direction = oc.get_edge_direction(edge);
                const Octant start = mesh->get_octant_edge_start(oc, edge);

                Integer max_proc = -1;
                mesh->get_one_ring_edge_nbhs(
                    start, direction, MARS_LAMBDA_REF(const Octant &nbh_oc) {
                        // find owner proc method returns an invalid result with invalid octant.
                        if (nbh_oc.is_valid()) {
                            Integer enc_oc = mars::get_sfc_from_octant<ElemType>(nbh_oc);
                            Integer owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, mesh->get_proc());
                            if (owner_proc > max_proc) {
                                max_proc = owner_proc;
                            }
                        }
                    });

                for (int j = 0; j < edge_dofs; j++) {
                    Integer dof_sfc =
                        process_edge_node<ElemType>(mars::get_sfc_from_octant<ElemType>, start, direction, j);
                    f(mesh, i, dof_sfc, max_proc, direction);
                }
            }
        }

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T != ElementType::Hex8, void> edge_iterate(const Integer sfc,
                                                                                                const Mesh *mesh,
                                                                                                const Integer i,
                                                                                                F f) {}

        template <typename CornerF, typename FaceF, typename VolumeF, typename EdgeF>
        static MARS_INLINE_FUNCTION void elem_dof_iterate(const Integer sfc,
                                                          const Mesh *mesh,
                                                          const Integer i,
                                                          CornerF cf,
                                                          FaceF ff,
                                                          VolumeF vf,
                                                          EdgeF ef) {
            if (corner_dofs > 0) {
                corner_iterate(sfc, mesh, i, cf);
            }
            if (face_dofs > 0) {
                face_dir_iterate(sfc, mesh, i, ff);
            }
            if (volume_dofs > 0) {
                volume_iterate(sfc, mesh, i, vf);
            }
            if (edge_dofs > 0) {
                edge_iterate(sfc, mesh, i, ef);
            }
        }

        struct IdentifyBoundaryDofPerRank {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = mesh->get_boundary_sfc(i);
                const Integer proc = mesh->get_proc();

                constexpr bool BoundaryIter = false;
                // used for corners and edges that are based on one ring nbh logic.
                auto rb = RankBoundary<BoundaryIter>(rank_boundary, sfc_to_locally_owned, map, proc);
                auto vrb =
                    RankBoundary<BoundaryIter, DofLabel::lVolume>(rank_boundary, sfc_to_locally_owned, map, proc);
                // use frb for the edge logic as well. Same logic different input from the iterate.
                elem_dof_iterate(sfc, mesh, i, rb, rb, vrb, rb);
            }

            IdentifyBoundaryDofPerRank(Mesh *m,
                                       ViewMatrixType<bool> npb,
                                       ViewVectorType<Integer> mp,
                                       ViewVectorType<Integer> l)
                : mesh(m), rank_boundary(npb), map(mp), sfc_to_locally_owned(l) {}

            Mesh *mesh;

            ViewMatrixType<bool> rank_boundary;
            ViewVectorType<Integer> map;
            ViewVectorType<Integer> sfc_to_locally_owned;
        };

        void build_boundary_dof_sets(ViewVectorType<Integer> &scan_boundary,
                                     ViewVectorType<Integer> &boundary_lsfc_index,
                                     ViewVectorType<Integer> &sender_ranks_scan,
                                     const Integer nbh_rank_size) {
            using namespace Kokkos;

            const Integer size = get_mesh_manager().get_host_mesh()->get_view_boundary().extent(0);

            Integer xDim = get_mesh_manager().get_host_mesh()->get_XDim();
            Integer yDim = get_mesh_manager().get_host_mesh()->get_YDim();
            Integer zDim = get_mesh_manager().get_host_mesh()->get_ZDim();

            const Integer chunk_size_ = global_dof_enum.get_elem_size();
            ViewMatrixType<bool> rank_boundary("count_per_proc", chunk_size_, nbh_rank_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("identify_boundary_predicate",
                                 size,
                                 IdentifyBoundaryDofPerRank(get_mesh_manager().get_mesh(),
                                                            rank_boundary,
                                                            sender_ranks_scan,
                                                            global_dof_enum.get_view_sfc_to_local()));

            /* perform a scan for each row with the sum at the end for each rank */
            ViewMatrixType<Integer> rank_scan("rank_scan", chunk_size_ + 1, nbh_rank_size);
            for (int i = 0; i < nbh_rank_size; ++i) {
                auto subpredicate = subview(rank_boundary, ALL, i);
                auto subscan = subview(rank_scan, ALL, i);
                incl_excl_scan(0, chunk_size_, subpredicate, subscan);
            }

            scan_boundary = ViewVectorType<Integer>("scan_boundary_dof", nbh_rank_size + 1);
            // perform a scan on the last column to get the total sum.
            row_scan(nbh_rank_size, chunk_size_, rank_scan, scan_boundary);

            auto index_subview = subview(scan_boundary, nbh_rank_size);
            auto h_ic = create_mirror_view(index_subview);
            // Deep copy device view to host view.
            deep_copy(h_ic, index_subview);

            boundary_dofs_sfc = ViewVectorType<Integer>("boundary_sfc_dofs", h_ic());
            boundary_lsfc_index = ViewVectorType<Integer>("boundary_lsfc_index_dofs", h_ic());
            /*   parallel_for(
                "print scan", rank_size + 1, KOKKOS_LAMBDA(const int i) {
                    printf(" scan boundary: %i-%li\n", i, scan_boundary(i));
                }); */

            ViewVectorType<Integer> g_elems = global_dof_enum.get_view_elements();
            ViewVectorType<Integer> boundary_ds = boundary_dofs_sfc;
            /* We use this strategy so that the compacted elements from the local_sfc
            would still be sorted and unique. */
            parallel_for(
                MDRangePolicy<Rank<2>>({0, 0}, {chunk_size_, nbh_rank_size}),
                KOKKOS_LAMBDA(const Integer i, const Integer j) {
                    if (rank_boundary(i, j) == 1) {
                        Integer index = scan_boundary(j) + rank_scan(i, j);
                        boundary_ds(index) = g_elems(i);
                        boundary_lsfc_index(index) = i;
                    }
                });

            /* parallel_for(
                "print set", h_ic(), KOKKOS_LAMBDA(const Integer i) {
                    Integer proc = get_mesh_manager().get_host_mesh()->get_proc();
                    const Integer rank = find_owner_processor(scan_boundary, i, 1, proc);

                    printf("i:%li - boundary_ : %i - %li (%li) - proc: %li - rank: %li\n", i, boundary_dofs_sfc(i),
               boundary_lsfc_index(i), get_octant_from_sfc<simplex_type::ElemType>(boundary_dofs_sfc(i)).template
               get_global_index<simplex_type::ElemType>(xDim, yDim), rank, proc);
                }); */
        }

        /*
        template <bool Ghost>
        static typename std::enable_if<Ghost == false, void>::type MARS_INLINE_FUNCTION
        set_volume_predicate(ViewVectorType<bool> local_predicate,
                ViewVectorType<bool> global_predicate, const Integer enc_oc) {
            local_predicate(enc_oc) = 1;
            global_predicate(enc_oc) = 1;
        }


        template <bool Ghost>
        static typename std::enable_if<Ghost == true, void>::type MARS_INLINE_FUNCTION
        set_volume_predicate(ViewVectorType<bool> local_predicate,
                ViewVectorType<bool> global_predicate, const Integer enc_oc) {
            local_predicate(enc_oc) = 1;
        } */

        template <Integer Type>
        static MARS_INLINE_FUNCTION Integer process_corner(const Mesh *mesh, const Octant &o) {
            Integer xDim = mesh->get_XDim();
            Integer yDim = mesh->get_YDim();
            Integer zDim = mesh->get_ZDim();

            Integer max_proc = -1;

            Octant one_ring[simplex_type::ElemType];
            o.one_ring_corner_nbhs<Type>(one_ring, xDim, yDim, zDim, mesh->is_periodic());

            for (int k = 0; k < Type; k++) {
                if (one_ring[k].is_valid()) {
                    Integer enc_oc = mars::get_sfc_from_octant<Type>(one_ring[k]);
                    /* check the proc that owns the corner and then decide on the predicate.
                        This works because the corner defines an element and this element is
                        always on the largest proc number due to the z order partitioning*/
                    Integer owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, mesh->get_proc());
                    /* assert(owner_proc >= 0); */

                    if (owner_proc > max_proc) {
                        max_proc = owner_proc;
                    }
                }
            }

            return max_proc;
        }

        template <typename F>
        static MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            Octant oc = mesh->octant_from_sfc(sfc);
            corner_octant_transform(oc, mesh, index, f);
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        corner_octant_transform(const Octant &oc, const Mesh *mesh, const Integer index, F f) {
            for (int i = 0; i < corner_dofs + 1; i++) {
                for (int j = 0; j < corner_dofs + 1; j++) {
                    Octant o;
                    // get the next corner using the elem sfc
                    o.x = oc.x + i;
                    o.y = oc.y + j;

                    Integer max_proc = process_corner<ElemType>(mesh, o);

                    // convert the octant value into the new nodal sfc system
                    o.x *= degree;
                    o.y *= degree;
                    Integer dof_sfc = mars::get_sfc_from_octant<ElemType>(o);

                    f(mesh, index, dof_sfc, max_proc, 0);
                }
            }
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        corner_octant_transform(const Octant &oc, const Mesh *mesh, const Integer index, F f) {
            for (int i = 0; i < corner_dofs + 1; i++) {
                for (int j = 0; j < corner_dofs + 1; j++) {
                    for (int k = 0; k < corner_dofs + 1; k++) {
                        Octant o;
                        // get the next corner using the elem sfc
                        o.x = oc.x + i;
                        o.y = oc.y + j;
                        o.z = oc.z + k;

                        Integer max_proc = process_corner<ElemType>(mesh, o);

                        // convert the octant value into the new nodal sfc system
                        o.x *= degree;
                        o.y *= degree;
                        o.z *= degree;
                        Integer dof_sfc = mars::get_sfc_from_octant<ElemType>(o);

                        f(mesh, index, dof_sfc, max_proc, 0);
                    }
                }
            }
        }

        //generic local predicate functor to be used for edge, face and corners. The volume is specialized.
        template <bool Ghost, Integer Label>
        struct LocalGlobalPredicate {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<Integer> local_label;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<Integer> global_label;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;
            Integer proc;

            MARS_INLINE_FUNCTION
            LocalGlobalPredicate(ViewVectorType<bool> lp,
                                 ViewVectorType<Integer> llp,
                                 ViewVectorType<bool> gp,
                                 ViewVectorType<Integer> lgp,
                                 ViewVectorType<bool> npps,
                                 ViewVectorType<bool> nppr,
                                 Integer p)
                : local_predicate(lp),
                  local_label(llp),
                  global_predicate(gp),
                  global_label(lgp),
                  nbh_proc_predicate_send(npps),
                  nbh_proc_predicate_recv(nppr),
                  proc(p) {}

            MARS_INLINE_FUNCTION
            void lg_predicate(const Mesh *mesh,
                                const Integer i,
                                const Integer sfc,
                                const Integer owner_proc,
                                std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh->get_view_gp(), mesh->get_ghost_sfc(i), 2, mesh->get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    local_predicate(sfc) = 1;
                    local_label(sfc) = Label;
                }
            }

            MARS_INLINE_FUNCTION
            void lg_predicate(const Mesh *mesh,
                              const Integer i,
                              const Integer sfc,
                              const Integer owner_proc,
                              std::false_type) const {
                /* if the face neighbor element is ghost then check if the processor
                    is less than the owner. This is how the dofs are partitioned*/
                if (proc >= owner_proc) {
                    global_predicate(sfc) = 1;
                    global_label(sfc) = Label;
                }
                if (proc != owner_proc) {
                    nbh_proc_predicate_send(owner_proc) = 1;
                    nbh_proc_predicate_recv(owner_proc) = 1;
                }
                local_predicate(sfc) = 1;
                local_label(sfc) = Label;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                lg_predicate(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct LocalGlobalPredicate<Ghost, DofLabel::lVolume> {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<Integer> local_label;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<Integer> global_label;

            MARS_INLINE_FUNCTION
            LocalGlobalPredicate(ViewVectorType<bool> lp,
                                 ViewVectorType<Integer> llp,
                                 ViewVectorType<bool> gp,
                                 ViewVectorType<Integer> lgp)
                : local_predicate(lp), local_label(llp), global_predicate(gp), global_label(lgp) {}

            MARS_INLINE_FUNCTION
            void volume_predicate(const Integer sfc, std::true_type) const {
                local_predicate(sfc) = 1;
                local_label(sfc) = DofLabel::lVolume;
            }

            MARS_INLINE_FUNCTION
            void volume_predicate(const Integer sfc, std::false_type) const {
                local_predicate(sfc) = 1;
                local_label(sfc) = DofLabel::lVolume;
                global_predicate(sfc) = 1;
                global_label(sfc) = DofLabel::lVolume;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh, const Integer i, const Integer enc_oc) const {
                volume_predicate(enc_oc, std::integral_constant<bool, Ghost>{});
            }
        };

        // sfc | octant_from_sfc | volume_octant_transform. Composable functions.
        template <typename F, Integer ET = ElemType>
        static void volume_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            Octant oc = mesh->octant_from_sfc(sfc);
            volume_octant_transform(oc, mesh, index, f, mars::get_sfc_from_octant<ElemType>);
        }

        template <typename F, typename G, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        volume_octant_transform(const Octant &oc, const Mesh *mesh, const Integer index, F f, G g) {
            Octant o;
            // go through all the inside dofs for the current element
            for (int i = 0; i < degree - 1; i++) {
                for (int j = 0; j < degree - 1; j++) {
                    o.x = degree * oc.x + i + 1;
                    o.y = degree * oc.y + j + 1;
                    f(mesh, index, g(o));
                }
            }
        }

        template <typename F, typename G, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        volume_octant_transform(const Octant &oc, const Mesh *mesh, const Integer index, F f, G g) {
            Octant o;
            // go through all the inside dofs for the current element
            for (int i = 0; i < degree - 1; i++) {
                for (int j = 0; j < degree - 1; j++) {
                    for (int k = 0; k < degree - 1; k++) {
                        o.x = degree * oc.x + i + 1;
                        o.y = degree * oc.y + j + 1;
                        o.z = degree * oc.z + k + 1;
                        f(mesh, index, g(o));
                    }
                }
            }
        }

        template <typename F, Integer Type = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, void>
        face_dir_iterate(const Integer sfc, const Mesh *mesh, const Integer i, F f) {
            // dir 0 = x, 1 = y
            face_iterate<0>(sfc, mesh, i, f);
            face_iterate<1>(sfc, mesh, i, f);
        }

        template <typename F, Integer Type = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, void>
        face_dir_iterate(const Integer sfc, const Mesh *mesh, const Integer i, F f) {
            // dir 0 = x, 1 = y, 2 = z
            face_iterate<0>(sfc, mesh, i, f);
            face_iterate<1>(sfc, mesh, i, f);
            face_iterate<2>(sfc, mesh, i, f);
        }

        // careful: do not use this function for face iteration when expensive operations are involved. Useful mainly
        // for predicate builder. instead use the locally_owned_face_dof array to iterate over face dof sfcs.
        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh->octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_dofs; j++) {
                    Integer dof_sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(mesh, index, dof_sfc, owner_proc, dir);
                }
            }
        }

        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const Integer sfc, const Mesh *mesh, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh->octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_dofs; j++) {
                    Integer dof_sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(dof_sfc, owner_proc, dir);
                }
            }
        }

        static MARS_INLINE_FUNCTION Integer get_ghost_sfc(const Mesh *mesh, const Integer index) {
            return mesh->get_ghost_sfc(index);
        }
        static MARS_INLINE_FUNCTION Integer get_local_sfc(const Mesh *mesh, const Integer index) {
            return mesh->get_sfc(index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == true, Integer>::type get_sfc_ghost_or_local(
            const Mesh *mesh,
            const Integer index) {
            return get_ghost_sfc(mesh, index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == false, Integer>::type get_sfc_ghost_or_local(
            const Mesh *mesh,
            const Integer index) {
            return get_local_sfc(mesh, index);
        }

        template <bool Ghost>
        struct BuildLocalGlobalPredicate {
            BuildLocalGlobalPredicate(Mesh *m,
                                      ViewVectorType<bool> lp,
                                      ViewVectorType<Integer> llp,
                                      ViewVectorType<bool> gp,
                                      ViewVectorType<Integer> lgp,
                                      ViewVectorType<bool> npbs,
                                      ViewVectorType<bool> npbr)
                : mesh(m),
                  local_predicate(lp),
                  local_label(llp),
                  global_predicate(gp),
                  global_label(lgp),
                  nbh_proc_predicate_send(npbs),
                  nbh_proc_predicate_recv(npbr) {}

            Mesh *mesh;
            ViewVectorType<bool> local_predicate;
            ViewVectorType<Integer> local_label;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<Integer> global_label;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh->get_proc();
                auto cp = LocalGlobalPredicate<Ghost, DofLabel::lCorner>(local_predicate,
                                                                         local_label,
                                                                         global_predicate,
                                                                         global_label,
                                                                         nbh_proc_predicate_send,
                                                                         nbh_proc_predicate_recv,
                                                                         proc);

                auto fp = LocalGlobalPredicate<Ghost, DofLabel::lFace>(local_predicate,
                                                                       local_label,
                                                                       global_predicate,
                                                                       global_label,
                                                                       nbh_proc_predicate_send,
                                                                       nbh_proc_predicate_recv,
                                                                       proc);

                // fp logic is the same as the edge one. Just different label.
                auto ep = LocalGlobalPredicate<Ghost, DofLabel::lEdge>(local_predicate,
                                                                       local_label,
                                                                       global_predicate,
                                                                       global_label,
                                                                       nbh_proc_predicate_send,
                                                                       nbh_proc_predicate_recv,
                                                                       proc);

                auto vp = LocalGlobalPredicate<Ghost, DofLabel::lVolume>(
                    local_predicate, local_label, global_predicate, global_label);

                elem_dof_iterate(sfc, mesh, i, cp, fp, vp, ep);
            }
        };

        void print_dofs(const int rank) {
            auto handler = *this;

            SFC<ElemType> gdof = get_global_dof_enum();
            Kokkos::parallel_for(
                "for", gdof.get_elem_size(), MARS_LAMBDA(const int i) {
                    const Integer sfc_elem = handler.owned_to_sfc(i);
                    const Dof global_dof = handler.sfc_to_global_dof(sfc_elem);

                    /* Integer gid = -1, proc = -1;

                    if (locally_owned_dof(sfc_elem)) {
                        proc = get_mesh_manager().get_mesh()->get_proc();
                        const Integer sfc_lid = global_dof_enum.get_view_sfc_to_local()(sfc_elem);
                        gid = sfc_lid + global_dof_offset(proc);
                    } */

                    Octant o = get_octant_from_sfc(sfc_elem);

                    printf("i: %i global sfc: %li gdof: %li octant: [%li, %li, %li] -  rank: %i\n",
                           i,
                           sfc_elem,
                           global_dof.get_gid(),
                           o.x,
                           o.y,
                           o.z,
                           global_dof.get_proc());
                });

            /* SFC<ElemType> dof = get_local_dof_enum();
            Kokkos::parallel_for(
                "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
                    const Integer sfc_elem = handler.local_to_sfc(i);
                    const Dof global_dof = handler.local_to_global_dof(i);

                    Octant o = get_octant_from_sfc(sfc_elem);

                    printf("i: %i, local sfc: %li gdof: %li --- octant: [%li, %li, %li] -  rank: %i\n",
                           i,
                           sfc_elem,
                           global_dof.get_gid(),
                           o.x,
                           o.y,
                           o.z,
                           global_dof.get_proc());
                }); */
        }

        void build_lg_predicate(const context &context, ViewVectorType<bool> &npbs, ViewVectorType<bool> &npbr) {
            Integer xDim = get_mesh_manager().get_host_mesh()->get_XDim();
            Integer yDim = get_mesh_manager().get_host_mesh()->get_YDim();
            Integer zDim = get_mesh_manager().get_host_mesh()->get_ZDim();

            const Integer size = get_mesh_manager().get_host_mesh()->get_chunk_size();
            const Integer ghost_size = get_mesh_manager().get_host_mesh()->get_view_ghost().extent(0);

            /* TODO: xdim and ydim should be changed to max xdim and ydim
             * of the local partition to reduce memory footprint */
            local_dof_enum = SFC<simplex_type::ElemType>(degree * xDim, degree * yDim, degree * zDim);
            global_dof_enum = SFC<simplex_type::ElemType>(degree * xDim, degree * yDim, degree * zDim);

            const Integer lall_range = local_dof_enum.get_all_range();
            const Integer gall_range = global_dof_enum.get_all_range();

            ViewVectorType<bool> local_predicate("lpred", lall_range);
            ViewVectorType<Integer> local_label("llabel", lall_range);
            ViewVectorType<bool> global_predicate("gpred", gall_range);
            ViewVectorType<Integer> global_label("glabel", gall_range);

            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("lg_predicate",
                                 size,
                                 BuildLocalGlobalPredicate<false>(get_mesh_manager().get_mesh(),
                                                                  local_predicate,
                                                                  local_label,
                                                                  global_predicate,
                                                                  global_label,
                                                                  npbs,
                                                                  npbr));

            // Iterate through ghost sfc and enumerate
            Kokkos::parallel_for("lg_predicate_from_ghost",
                                 ghost_size,
                                 BuildLocalGlobalPredicate<true>(get_mesh_manager().get_mesh(),
                                                                 local_predicate,
                                                                 local_label,
                                                                 global_predicate,
                                                                 global_label,
                                                                 npbs,
                                                                 npbr));

            local_dof_enum.compact_element_and_labels(local_predicate, local_label);
            global_dof_enum.compact_element_and_labels(global_predicate, global_label);

            const int rank_size = num_ranks(context);

            ViewVectorType<Integer> global_dof_size_per_rank("global_dof_size_per_rank", rank_size);
            global_dof_offset = ViewVectorType<Integer>("global_dof_offset", rank_size + 1);

            const Integer global_size = global_dof_enum.get_elem_size();
            context->distributed->gather_all_view(global_size, global_dof_size_per_rank);

            incl_excl_scan(0, rank_size, global_dof_size_per_rank, global_dof_offset);
        }

        void exchange_ghost_counts(const context &context,
                                   ViewVectorType<Integer> &scan_boundary,
                                   ViewVectorType<Integer> &sender_ranks,
                                   ViewVectorType<Integer> &recver_ranks,
                                   std::vector<Integer> &send_count,
                                   std::vector<Integer> &receive_count)

        {
            using namespace Kokkos;

            Kokkos::Timer timer;

            int proc_num = rank(context);
            int size = num_ranks(context);

            auto scan_snd = create_mirror_view(scan_boundary);
            Kokkos::deep_copy(scan_snd, scan_boundary);

            auto s_rank_mirror = create_mirror_view(sender_ranks);
            Kokkos::deep_copy(s_rank_mirror, sender_ranks);
            const Integer nbh_rank_size = sender_ranks.extent(0);

            /* Integer proc_count = 0; */
            for (int i = 0; i < nbh_rank_size; ++i) {
                Integer count = scan_snd(i + 1) - scan_snd(i);
                if (count > 0) {
                    send_count[s_rank_mirror(i)] = count;
                    /* ++proc_count;
                    std::cout << "DM:****ToProc: " << s_rank_mirror(i) << " count:" << count
                              << " Proc: " << proc_num << std::endl; */
                }
            }

            auto r_rank_mirror = create_mirror_view(recver_ranks);
            Kokkos::deep_copy(r_rank_mirror, recver_ranks);
            const Integer nbh_rank_recv_size = recver_ranks.extent(0);

            /* Integer proc_count_r = 0; */
            for (int i = 0; i < nbh_rank_recv_size; ++i) {
                receive_count[r_rank_mirror(i)] = 1;
            }

            context->distributed->i_send_recv_vec(send_count, receive_count);

            /* for (int i = 0; i < size; ++i)
            {
                if (receive_count[i] > 0)
                {
                    std::cout << "DM:-----FromProc: " << i << " count:" << receive_count[i]
                              << " Proc: " << proc_num << std::endl;
                }
            } */
        }

        void exchange_ghost_dofs(const context &context,
                                 ViewVectorType<Integer> boundary_lsfc_index,
                                 ViewVectorType<Integer> &ghost_dofs_index) {
            int rank_size = num_ranks(context);
            Integer ghost_size = scan_recv_mirror(rank_size);

            ghost_dofs_sfc = ViewVectorType<Integer>("ghost_dofs", ghost_size);
            ghost_dofs_index = ViewVectorType<Integer>("ghost_dofs_index", ghost_size);
            // do it again to have all the process range to make it fit the i_send_recv_view

            context->distributed->i_send_recv_view(
                ghost_dofs_sfc, scan_recv_mirror.data(), boundary_dofs_sfc, scan_send_mirror.data());
            /* std::cout << "DM:Ending mpi send receive for the ghost sfc dofs " << std::endl; */

            context->distributed->i_send_recv_view(
                ghost_dofs_index, scan_recv_mirror.data(), boundary_lsfc_index, scan_send_mirror.data());
            std::cout << "DofHandler:MPI send receive for the ghost dofs layer done." << std::endl;
        }

        template <bool Ghost>
        struct FaceOrientDof {
            ViewVectorType<Integer> dir;
            ViewVectorType<Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            FaceOrientDof(ViewVectorType<Integer> d, ViewVectorType<Integer> l, Integer p)
                : dir(d), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            void face_orient_dof(const Mesh *mesh,
                                 const Integer i,
                                 const Integer sfc,
                                 const Integer owner_proc,
                                 const Integer d,
                                 std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh->get_view_gp(), mesh->get_ghost_sfc(i), 2, mesh->get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    Integer index = sfc_to_local(sfc);
                    dir(index) = d;
                }
            }

            MARS_INLINE_FUNCTION
            void face_orient_dof(const Mesh *mesh,
                                 const Integer i,
                                 const Integer sfc,
                                 const Integer owner_proc,
                                 const Integer d,
                                 std::false_type) const {
                Integer index = sfc_to_local(sfc);
                dir(index) = d;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                face_orient_dof(mesh, i, dof_sfc, owner_proc, dir, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct OrientDofs {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh->get_proc();

                if (face_dofs > 0) {
                    FaceOrientDof<Ghost> fp = FaceOrientDof<Ghost>(face_dir, sfc_to_local, proc);
                    face_dir_iterate(sfc, mesh, i, fp);
                }
            }

            OrientDofs(Mesh *m, ViewVectorType<Integer> sd, ViewVectorType<Integer> sl)
                : mesh(m), face_dir(sd), sfc_to_local(sl) {}

            Mesh *mesh;
            ViewVectorType<Integer> face_dir;
            ViewVectorType<Integer> sfc_to_local;
        };

        void build_local_orientation() {
            const Integer size = get_mesh_manager().get_host_mesh()->get_chunk_size();
            const Integer ghost_size = get_mesh_manager().get_host_mesh()->get_view_ghost().extent(0);

            const Integer local_size = get_local_dof_enum().get_elem_size();

            local_dof_enum.init_element_orientations(local_size);

            /* generate the sfc for the local and global dofs containing the generation locally
                     for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("orient_dofs",
                                 size,
                                 OrientDofs<false>(get_mesh_manager().get_mesh(),
                                                   local_dof_enum.get_view_element_orientations(),
                                                   get_local_dof_enum().get_view_sfc_to_local()));

            Kokkos::parallel_for("orient_ghost_dofs",
                                 ghost_size,
                                 OrientDofs<true>(get_mesh_manager().get_mesh(),
                                                  local_dof_enum.get_view_element_orientations(),
                                                  get_local_dof_enum().get_view_sfc_to_local()));
        }

        virtual void enumerate_dofs() {
            const auto &context = get_context();

            const Integer rank_size = num_ranks(context);
            const int proc_num = rank(context);

            ViewVectorType<bool> nbh_proc_predicate_send("send_to", rank_size);
            ViewVectorType<bool> nbh_proc_predicate_recv("receive_from", rank_size);
            ViewVectorType<Integer> proc_scan_send("nbh_scan_send", rank_size + 1);
            ViewVectorType<Integer> proc_scan_recv("nbh_scan_recv", rank_size + 1);

            build_lg_predicate(context, nbh_proc_predicate_send, nbh_proc_predicate_recv);
            build_local_orientation();

            incl_excl_scan(0, rank_size, nbh_proc_predicate_send, proc_scan_send);
            incl_excl_scan(0, rank_size, nbh_proc_predicate_recv, proc_scan_recv);

            auto ps = subview(proc_scan_send, rank_size);
            auto pr = subview(proc_scan_recv, rank_size);
            auto h_ps = create_mirror_view(ps);
            auto h_pr = create_mirror_view(pr);
            // Deep copy device view to host view.
            deep_copy(h_ps, ps);
            deep_copy(h_pr, pr);

            ViewVectorType<Integer> sender_ranks("sender_ranks", h_ps());
            ViewVectorType<Integer> recver_ranks("receiver_ranks", h_pr());
            compact_scan(nbh_proc_predicate_send, proc_scan_send, sender_ranks);
            compact_scan(nbh_proc_predicate_recv, proc_scan_recv, recver_ranks);

            ViewVectorType<Integer> scan_boundary, boundary_dofs_index;
            build_boundary_dof_sets(scan_boundary, boundary_dofs_index, proc_scan_send, h_ps());

            std::vector<Integer> s_count(rank_size, 0);
            std::vector<Integer> r_count(rank_size, 0);

            exchange_ghost_counts(context, scan_boundary, sender_ranks, recver_ranks, s_count, r_count);

            /* create the scan recv mirror view from the receive count */
            scan_recv = ViewVectorType<Integer>("scan_recv_", rank_size + 1);
            scan_recv_mirror = create_mirror_view(scan_recv);
            make_scan_index_mirror(scan_recv_mirror, r_count);
            Kokkos::deep_copy(scan_recv, scan_recv_mirror);

            /*create the scan send mirror view from the send count*/
            scan_send = ViewVectorType<Integer>("scan_send_", rank_size + 1);
            scan_send_mirror = create_mirror_view(scan_send);
            make_scan_index_mirror(scan_send_mirror, s_count);
            Kokkos::deep_copy(scan_send, scan_send_mirror);

            ViewVectorType<Integer> ghost_dofs_index;
            exchange_ghost_dofs(context, boundary_dofs_index, ghost_dofs_index);

            build_ghost_local_global_map(context, ghost_dofs_index);
        }

        /* build a map only for the ghost dofs as for the local ones the global dof enum can be used */
        struct BuildLocalToGlobalGhostMap {
            Mesh *mesh;
            UnorderedMap<Integer, Dof> glgm;
            ViewVectorType<Integer> global_dof_offset;
            ViewVectorType<Integer> scan_recv_proc;

            ViewVectorType<Integer> ghost_dofs;
            ViewVectorType<Integer> ghost_dofs_index;

            BuildLocalToGlobalGhostMap(Mesh *m,
                                       UnorderedMap<Integer, Dof> g,
                                       ViewVectorType<Integer> gdo,
                                       ViewVectorType<Integer> srm,
                                       ViewVectorType<Integer> gd,
                                       ViewVectorType<Integer> gdi)
                : mesh(m),
                  glgm(g),
                  global_dof_offset(gdo),
                  scan_recv_proc(srm),
                  ghost_dofs(gd),
                  ghost_dofs_index(gdi) {}

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                /* read the sfc and the local dof id from the ghost */
                const Integer ghost_sfc = ghost_dofs(i);
                const Integer ghost_sfc_lid = ghost_dofs_index(i);

                /* find the process by binary search in the scan_recv_proc view of size rank_size
                 * and calculate the global id by adding the ghost local id to the global offset for ghost process*/
                const Integer owner_proc = find_owner_processor(scan_recv_proc, i, 1, mesh->get_proc());
                const Integer gid = ghost_sfc_lid + global_dof_offset(owner_proc);

                // build the ghost dof object and insert into the map
                const auto result = glgm.insert(ghost_sfc, Dof(gid, owner_proc));
                assert(result.success());
            }
        };

        void build_ghost_local_global_map(const context &context, const ViewVectorType<Integer> &ghost_dofs_index) {
            int rank_size = num_ranks(context);
            Integer size = scan_recv_mirror(rank_size);

            ghost_local_to_global_map = UnorderedMap<Integer, Dof>(size);
            /* iterate through the unique ghost dofs and build the map */
            Kokkos::parallel_for("BuildLocalGlobalmap",
                                 size,
                                 BuildLocalToGlobalGhostMap(get_mesh_manager().get_mesh(),
                                                            ghost_local_to_global_map,
                                                            global_dof_offset,
                                                            scan_recv,
                                                            ghost_dofs_sfc,
                                                            ghost_dofs_index));
            /* In the end the size of the map should be as the size of the ghost_dofs.
             * Careful map size  is not capacity */
            assert(size == ghost_local_to_global_map.size());
        }

        MARS_INLINE_FUNCTION
        bool is_local(const Integer sfc) const {
            auto stl = local_dof_enum.get_view_sfc_to_local();
            if ((sfc + 1) >= stl.extent(0)) return false;
            /* use the sfc to local which is the scan of the predicate.
             * To get the predicate value the difference with the successive index is needed. */
            const Integer pred_value = stl(sfc + 1) - stl(sfc);
            return (pred_value > 0);
        }

        MARS_INLINE_FUNCTION
        bool locally_owned_dof(const Integer sfc) const {
            // use the sfc to local which the scan of the predicate. To get the predicate
            // value the difference with the successive index is needed.
            const Integer pred_value =
                global_dof_enum.get_view_sfc_to_local()(sfc + 1) - global_dof_enum.get_view_sfc_to_local()(sfc);
            return (pred_value > 0);
        }

        MARS_INLINE_FUNCTION
        Integer local_to_owned_dof(const Integer local) const {
            const Integer sfc = local_to_sfc(local);
            if (locally_owned_dof(sfc))
                return get_global_dof_enum().get_view_sfc_to_local()(sfc);
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Dof sfc_to_global_dof(const Integer sfc) const {
            Dof dof;
            if (locally_owned_dof(sfc)) {
                const Integer proc = get_mesh_manager().get_mesh()->get_proc();
                const Integer sfc_lid = global_dof_enum.get_view_sfc_to_local()(sfc);
                const Integer gid = sfc_lid + global_dof_offset(proc);
                dof.set_gid(gid);
                dof.set_proc(proc);
            } else {
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
        Dof local_to_global_dof(const Integer local) const {
            // use the local to sfc view (elements) to get the sfc of the local numbering.
            const Integer local_sfc = local_dof_enum.get_view_elements()(local);
            return sfc_to_global_dof(local_sfc);
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_global(const Integer sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        const Integer local_to_owned(const Integer local) const {
            Dof dof = local_to_global_dof(local);
            if (dof.is_valid())
                return dof.get_gid() - global_dof_offset(dof.get_proc());
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Integer local_to_global(const Integer local) const {
            Dof dof = local_to_global_dof(local);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_global_proc(const Integer sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Integer local_to_global_proc(const Integer local) const {
            Dof dof = local_to_global_dof(local);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Integer get_local_label(const Integer local) const { return local_dof_enum.get_view_element_labels()(local); }

        MARS_INLINE_FUNCTION
        Integer local_to_sfc(const Integer local) const { return local_dof_enum.get_view_elements()(local); }

        MARS_INLINE_FUNCTION
        Integer owned_to_sfc(const Integer owned) const { return global_dof_enum.get_view_elements()(owned); }

        MARS_INLINE_FUNCTION
        Integer sfc_to_owned(const Integer sfc) const { return global_dof_enum.get_view_sfc_to_local()(sfc); }

        MARS_INLINE_FUNCTION
        Integer sfc_to_local(const Integer sfc) const { return local_dof_enum.get_view_sfc_to_local()(sfc); }

        template <Integer Type>
        MARS_INLINE_FUNCTION void get_dof_coordinates_from_sfc(const Integer sfc, double *point) const {
            get_vertex_coordinates_from_sfc<Type>(sfc,
                                                  point,
                                                  get_local_dof_enum().get_XDim(),
                                                  get_local_dof_enum().get_YDim(),
                                                  get_local_dof_enum().get_ZDim());
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION void get_dof_coordinates_from_local(const Integer local, double *point) const {
            const Integer sfc = local_to_sfc(local);
            return get_dof_coordinates_from_sfc<Type>(sfc, point);
        }

        MARS_INLINE_FUNCTION void get_local_dof_coordinates(const Integer local, double *point) const {
            get_dof_coordinates_from_local<ElemType>(local, point);
        }

        template <Integer Type, Integer FaceNr = -1>
        MARS_INLINE_FUNCTION bool is_boundary(const Integer local) const {
            const Integer sfc = local_to_sfc(local);
            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            return is_boundary_sfc<Type, FaceNr>(sfc, xdim, ydim, zdim);
        }

        template <Integer FaceNr = -1>
        MARS_INLINE_FUNCTION bool is_boundary_dof(const Integer local) const {
            return is_boundary<ElemType, FaceNr>(local);
        }

        template <Integer face_nr = -1, typename F>
        void boundary_owned_dof_iterate(F f) {
            using namespace Kokkos;
            constexpr Integer Type = simplex_type::ElemType;

            const Integer size = get_global_dof_enum().get_elem_size();
            ViewVectorType<Integer> global_to_sfc = get_global_dof_enum().get_view_elements();

            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            Kokkos::parallel_for(
                "boundary_owned_dof_iterate", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    if (is_boundary_sfc<Type, face_nr>(sfc, xdim, ydim, zdim)) {
                        f(i, sfc);
                    }
                });
        }

        template <Integer face_nr = -1, typename F>
        void boundary_dof_iterate(F f) {
            using namespace Kokkos;
            constexpr Integer Type = simplex_type::ElemType;

            const Integer size = get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_local_dof_enum().get_view_sfc_to_local();

            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            Kokkos::parallel_for(
                "boundary_iterate", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);

                    if (is_boundary_sfc<Type, face_nr>(sfc, xdim, ydim, zdim)) {
                        f(local);
                    }
                });
        }

        MARS_INLINE_FUNCTION
        const Integer get_dof_size() const { return local_dof_enum.get_elem_size(); }

        MARS_INLINE_FUNCTION
        const Integer get_local_dof(const Integer i) const {
            assert(i < get_dof_size());
            return i;
        }

        // get the local dof of the owned index.
        MARS_INLINE_FUNCTION
        const Integer get_owned_dof(const Integer i) const {
            const Integer sfc = get_global_dof_enum().get_view_elements()(i);
            return sfc_to_local(sfc);
        }

        MARS_INLINE_FUNCTION
        const Integer get_owned_dof_size() const { return global_dof_enum.get_elem_size(); }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_local_dof_enum() const { return local_dof_enum; }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_global_dof_enum() const { return global_dof_enum; }

        MARS_INLINE_FUNCTION
        const Integer get_global_dof_offset(const Integer proc) const { return global_dof_offset(proc); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_global_dof_offset() const { return global_dof_offset; }

        MARS_INLINE_FUNCTION
        const Integer get_global_dof_size() const {
            const Integer rank_size = num_ranks(get_context());

            auto ss = subview(global_dof_offset, rank_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            return h_ss();
        }

        MARS_INLINE_FUNCTION
        const UnorderedMap<Integer, Dof> &get_ghost_lg_map() const { return ghost_local_to_global_map; }

        MARS_INLINE_FUNCTION
        Integer get_boundary_dof(const Integer i) const { return sfc_to_local(get_boundary_dofs()(i)); }

        MARS_INLINE_FUNCTION
        Integer get_boundary_dof_size() const { return get_boundary_dofs().extent(0); }

        MARS_INLINE_FUNCTION
        Integer get_ghost_dof(const Integer i) const { return sfc_to_local(get_ghost_dofs()(i)); }

        MARS_INLINE_FUNCTION
        Integer get_ghost_dof_size() const { return get_ghost_dofs().extent(0); }

        template <typename F>
        void ghost_iterate(F f) const {
            Kokkos::parallel_for("ghost_dof_iter", get_ghost_dof_size(), f);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_local_dof_map() const { return local_dof_enum.get_view_sfc_to_local(); }

        MARS_INLINE_FUNCTION
        const Integer get_local_dof_map(const Integer local_dof) const {
            return local_dof_enum.get_view_sfc_to_local()(local_dof);
        }

        MM get_mesh_manager() const { return mesh_manager; }

        MARS_INLINE_FUNCTION
        const Integer get_proc() const { return get_mesh_manager().get_mesh()->get_proc(); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_recv() const { return scan_recv; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_send() const { return scan_send; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_recv_mirror() const { return scan_recv_mirror; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_send_mirror() const { return scan_send_mirror; }

        MARS_INLINE_FUNCTION
        const Integer get_orientation(const Integer local_dof) const {
            return get_local_dof_enum().get_orientation(local_dof);
        }

        MARS_INLINE_FUNCTION
        const Integer get_label(const Integer local_dof) const { return get_local_dof_enum().get_label(local_dof); }

        MARS_INLINE_FUNCTION
        const Integer get_owned_label(const Integer owned_dof) const {
            return get_global_dof_enum().get_label(owned_dof);
        }

        const context &get_context() const { return ctx; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_boundary_dofs() const { return boundary_dofs_sfc; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_ghost_dofs() const { return ghost_dofs_sfc; }

        MARS_INLINE_FUNCTION Octant get_octant_from_local(const Integer local) const {
            const Integer sfc = local_to_sfc(local);
            return get_octant_from_sfc(sfc);
        }

        MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const Integer sfc) const {
            return mars::get_octant_from_sfc<ElemType>(sfc);
        }

        MARS_INLINE_FUNCTION Integer get_sfc_from_octant(const Octant &o) const {
            return mars::get_sfc_from_octant<ElemType>(o);
        }

        MARS_INLINE_FUNCTION Integer get_global_from_octant(const Octant &o) const {
            const Integer sfc = get_sfc_from_octant(o);
            const Integer local = is_local(sfc) ? sfc_to_local(sfc) : INVALID_INDEX;
            return local_to_global(local);
        }

        MARS_INLINE_FUNCTION Integer get_local_sfc_from_octant(const Octant &o) const {
            const Integer sfc = get_sfc_from_octant(o);
            return is_local(sfc) ? sfc : INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_local_from_octant(const Octant &o) const {
            const Integer sfc = get_sfc_from_octant(o);
            return is_local(sfc) ? sfc_to_local(sfc) : INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        const Integer get_XMax() const { return Degree * get_mesh_manager().get_mesh()->get_XDim(); }

        MARS_INLINE_FUNCTION
        const Integer get_YMax() const { return Degree * get_mesh_manager().get_mesh()->get_YDim(); }

        MARS_INLINE_FUNCTION
        const Integer get_ZMax() const { return Degree * get_mesh_manager().get_mesh()->get_ZDim(); }

    private:
        // data associated to the mesh elements (sfc) within the context.
        MM mesh_manager;
        const context &ctx;

        // local dof enumeration containing the local to sfc and sfc to local views.
        SFC<simplex_type::ElemType> local_dof_enum;
        // locally owned dof enumeration containing the global to sfc and sfc to global views.
        SFC<simplex_type::ElemType> global_dof_enum;
        // global offset used to calc the global numbering of the dofs.
        ViewVectorType<Integer> global_dof_offset;

        // the local to global mapp for the locally owned is mapped from the sfc_to global view
        // for the ghost dofs is used the following kokkos unordered map.
        UnorderedMap<Integer, Dof> ghost_local_to_global_map;

        // dofs sfc for the boundary dofs.
        ViewVectorType<Integer> boundary_dofs_sfc;
        // mirror view on the dof scan boundary view used as offset for the mpi send dofs.
        ViewVectorType<Integer> scan_send;
        ViewVectorType<Integer>::HostMirror scan_send_mirror;

        // ghost dofs received from other procs's boundary sfcs.
        ViewVectorType<Integer> ghost_dofs_sfc;
        // mirror view on the scan_ghost view used as an offset to receive ghost dofs.
        ViewVectorType<Integer> scan_recv;
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;
    };

}  // namespace mars

#endif
#endif

#endif
