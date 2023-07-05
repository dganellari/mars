#ifndef GENERATION_MARS_DISTRIBUTED_DofHandler_HPP_
#define GENERATION_MARS_DISTRIBUTED_DofHandler_HPP_

#include <cstdlib>
#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_dof.hpp"
#include "mars_distributed_mesh_management.hpp"

namespace mars {

    template <Integer degree, Integer Label, Integer Type>
    struct NumDofs {};

    template <Integer degree, Integer Label>
    struct NumDofs<degree, Label, ElementType::Quad4> {
        static constexpr Integer Dim = 2;

        static constexpr Integer num_corners = power_of_2(Dim);
        static constexpr Integer num_faces = 2 * Dim;
        static constexpr Integer num_edges = 0;

        static constexpr Integer corner_dofs = 1;
        static constexpr Integer face_dofs = (degree - 1);
        static constexpr Integer volume_dofs = (degree - 1) * (degree - 1);
        static constexpr Integer edge_dofs = 0;
        static constexpr Integer total_elem_dofs = (degree + 1) * (degree + 1);

        static constexpr Integer elem_dofs() noexcept {
            auto result = 0;

            if (Label & DofLabel::lCorner) {
                result += corner_dofs * num_corners;
            }
            if (Label & DofLabel::lVolume) {
                result += volume_dofs;
            }
            if (Label & DofLabel::lFace) {
                result += face_dofs * num_faces;
            }
            return result;
        }
    };

    template <Integer degree, Integer Label>
    struct NumDofs<degree, Label, ElementType::Hex8> {
        static constexpr Integer Dim = 3;

        static constexpr Integer num_corners = power_of_2(Dim);
        static constexpr Integer num_faces = 2 * Dim;
        static constexpr Integer num_edges = 4 * Dim;

        static constexpr Integer corner_dofs = 1;
        static constexpr Integer face_dofs = (degree - 1) * (degree - 1);
        static constexpr Integer volume_dofs = (degree - 1) * (degree - 1) * (degree - 1);
        static constexpr Integer edge_dofs = (degree - 1);
        static constexpr Integer total_elem_dofs = (degree + 1) * (degree + 1) * (degree + 1);

        static constexpr Integer elem_dofs() noexcept {
            auto result = 0;

            if (Label & DofLabel::lCorner) {
                result += corner_dofs * num_corners;
            }
            if (Label & DofLabel::lEdge) {
                result += edge_dofs * num_edges;
            }
            if (Label & DofLabel::lVolume) {
                result += volume_dofs;
            }
            if (Label & DofLabel::lFace) {
                result += face_dofs * num_faces;
            }
            return result;
        }
    };

    class IDofHandler {
    public:
        virtual MARS_INLINE_FUNCTION ~IDofHandler() {}
        // virtual MARS_FUNCTION const Integer get_ZMax() const;
    };

    template <class Mesh_, Integer degree, Integer Block_ = 1>
    class DofHandler : public IDofHandler {
    public:
        using Mesh = Mesh_;

        /* using MM = MeshManager<Mesh>; */
        using simplex_type = typename Mesh::Elem;

        static constexpr Integer ElemType = simplex_type::ElemType;

        static constexpr Integer dofLabel = DofLabel::lAll;

        static constexpr Integer Block = Block_;

        static constexpr Integer Degree = degree;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        using NDofs = NumDofs<degree, dofLabel, ElemType>;

        static constexpr Integer corner_dofs = NDofs::corner_dofs;
        static constexpr Integer num_corners = NDofs::num_corners;

        static constexpr Integer face_dofs = NDofs::face_dofs;
        static constexpr Integer num_faces = NDofs::num_faces;

        static constexpr Integer volume_dofs = NDofs::volume_dofs;

        static constexpr Integer edge_dofs = NDofs::edge_dofs;
        static constexpr Integer num_edges = NDofs::num_edges;

        static constexpr Integer elem_dofs = NDofs::elem_dofs();

        MARS_INLINE_FUNCTION
        constexpr Integer get_elem_type() const { return simplex_type::ElemType; }

        MARS_INLINE_FUNCTION
        DofHandler(Mesh m) : mesh(m) {}
        /* DofHandler(Mesh mesh) : mesh_manager(MM(mesh)) {} */

        // computes the component (block) from the block local
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Integer>::type compute_component(
            const Integer local) const {
            return 0;
        }

        // B==0 goes for runtime block input. Useful for interfacing with UtopiaFE.
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Integer>::type compute_component(
            const Integer local) const {
            return local % get_block<B>();
        }

        // computes the base (scalar) local id from the block local
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Integer>::type compute_base(const Integer local) const {
            return local;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Integer>::type compute_base(const Integer local) const {
            return local / get_block<B>();
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer compute_block_index(const Integer base_local, const Integer component) const {
            return base_local * get_block<B>() + component;
        }

        template <typename H>
        void owned_iterate(H f) const {
            owned_dof_iterate(f);
        }

        template <typename H, Integer B = Block_>
        void owned_dof_iterate(H f) const {
            Kokkos::parallel_for("init_initial_cond", get_block<B>() * global_dof_enum.get_elem_size(), f);
        }

        template <typename H, Integer B = Block_>
        void dof_iterate(H f) const {
            Kokkos::parallel_for("init_initial_cond", get_block<B>() * local_dof_enum.get_elem_size(), f);
        }

        template <typename H>
        void base_dof_iterate(H f) const {
            Kokkos::parallel_for("init_initial_cond", local_dof_enum.get_elem_size(), f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void elem_iterate(H f) const {
            get_mesh().elem_iterate(f);
        }

        MARS_INLINE_FUNCTION Integer get_elem_size() const {
            return get_mesh().get_chunk_size();
        }

        template <typename H>
        MARS_INLINE_FUNCTION void face_iterate(H f) const {
            get_mesh().face_iterate(f);
        }

        MARS_INLINE_FUNCTION
        Integer get_sfc_face_nbh(const Octant &oc, const Integer face_nr) const {
            Octant o = oc.sfc_face_nbh<simplex_type::ElemType>(face_nr);
            return get_sfc_from_octant(o);
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Integer> enum_corner(const Octant &oc,
                                                                                               const int i,
                                                                                               const int j) const {
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

        template <Integer Type>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer> enum_corner(const Octant &oc,
                                                                                              const int i,
                                                                                              const int j,
                                                                                              const int k) const {
            Octant o;
            // get the next corner using the elem sfc
            o.x = oc.x + i;
            o.y = oc.y + j;
            o.z = oc.z + k;
            // convert the octant value into the new nodal sfc system
            o.x *= degree;
            o.y *= degree;
            o.z *= degree;
            Integer sfc = mars::get_sfc_from_octant<Type>(o);

            return sfc_to_local(sfc);
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Octant> enum_face_corner(
            const Octant &oc,
            const int dir) {
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
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Octant> enum_face_corner(
            const Octant &oc,
            const int dir) {
            Octant face_cornerA;
            build_starting_face_corner<Type, part>(face_cornerA, dir, oc);
            return face_cornerA;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Integer>
        sfc_face_node(const Octant &face_cornerA, const int j, const int dir) {
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
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer>
        sfc_face_node(const Octant &face_cornerA, const int j, const int dir) {
            return build_face_node<Type, part>(mars::get_sfc_from_octant<Type>, face_cornerA, j);
        }

        template <Integer part, Integer Type>
        MARS_INLINE_FUNCTION Integer enum_face_node(const Octant &face_cornerA, const int j, const int dir) const {
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
        process_face_corner(Octant &face_cornerA, const Mesh mesh, const int side, const Octant &oc) {
            Integer face_nr;
            if (side == 0)
                face_nr = 2 * dir + 1;
            else
                face_nr = 2 * dir;

            Octant nbh_oc = mesh.get_octant_face_nbh(oc, face_nr);
            Integer owner_proc = -1;

            // find owner proc method returns an invalid result with invalid octant.
            if (nbh_oc.is_valid()) {
                Integer enc_oc = mars::get_sfc_from_octant<Type>(nbh_oc);
                owner_proc = find_owner_processor(mesh.get_view_gp(), enc_oc, 2, mesh.get_proc());
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
            ViewMatrixTypeLeft<bool> rank_boundary;
            UnorderedMap<Integer, Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            RankBoundary(ViewMatrixTypeLeft<bool> rb,
                         UnorderedMap<Integer, Integer> sl,
                         ViewVectorType<Integer> m,
                         Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            MARS_INLINE_FUNCTION
            void all_rank_boundary(const Mesh mesh,
                                   const Integer i,
                                   const Integer sfc,
                                   const Integer owner_proc,
                                   std::false_type) const {
                const Integer ghost_proc = find_owner_processor(mesh.get_view_scan_boundary(), i, 1, proc);

                /* if (proc >= owner_proc && owner_proc >= 0) { */
                if (proc >= owner_proc) {
                    Integer index = get_value_in_map(sfc_to_locally_owned, sfc);
                    rank_boundary(index, map(ghost_proc)) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer dof_sfc, const Integer owner_proc) const {
                all_rank_boundary(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct RankBoundary<Ghost, DofLabel::lVolume> {
            ViewMatrixTypeLeft<bool> rank_boundary;
            UnorderedMap<Integer, Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            RankBoundary(ViewMatrixTypeLeft<bool> rb,
                         UnorderedMap<Integer, Integer> sl,
                         ViewVectorType<Integer> m,
                         Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            MARS_INLINE_FUNCTION
            void volume_rank_boundary(const Mesh mesh, const Integer i, const Integer sfc, std::false_type) const {
                const Integer ghost_proc = find_owner_processor(mesh.get_view_scan_boundary(), i, 1, proc);
                Integer index = get_value_in_map(sfc_to_locally_owned, sfc);
                rank_boundary(index, map(ghost_proc)) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer dof_sfc) const {
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

        // TODO: Use this instead of the following functor when switching to c++17 and Cuda 11.
        /* auto max_lambda = MARS_LAMBDA_REF(const Octant &nbh_oc) {
            // find owner proc method returns an invalid result with invalid octant.
            if (nbh_oc.is_valid()) {
                Integer enc_oc = mars::get_sfc_from_octant<ElemType>(nbh_oc);
                Integer owner_proc = find_owner_processor(mesh.get_view_gp(), enc_oc, 2, mesh.get_proc());
                if (owner_proc > max_proc) {
                    max_proc = owner_proc;
                }
            }
        }; */

        struct MaxOneRingProc {
            Integer &max_proc;
            ViewVectorType<Integer> gp;
            Integer proc;

            MARS_INLINE_FUNCTION
            MaxOneRingProc(ViewVectorType<Integer> g, Integer p, Integer &mx) : gp(g), proc(p), max_proc(mx) {}

            MARS_INLINE_FUNCTION
            void operator()(const Octant &nbh_oc) const {
                // find owner proc method returns an invalid result with invalid octant.
                if (nbh_oc.is_valid()) {
                    Integer enc_oc = mars::get_sfc_from_octant<ElemType>(nbh_oc);
                    Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);
                    if (owner_proc > max_proc) {
                        max_proc = owner_proc;
                    }
                }
            }
        };

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T == ElementType::Hex8, void> edge_iterate(const Integer sfc,
                                                                                                const Mesh mesh,
                                                                                                const Integer i,
                                                                                                F f) {
            Octant oc = mesh.octant_from_sfc(sfc);
            for (int edge = 0; edge < 4 * ManifoldDim; ++edge) {
                // get the starting node of the edge to be used to compute the dofs contained in that edge.
                const Integer direction = oc.get_edge_direction(edge);
                const Octant start = mesh.get_octant_edge_start(oc, edge);

                Integer max_proc = -1;
                mesh.get_one_ring_edge_nbhs(
                    start, direction, MaxOneRingProc(mesh.get_view_gp(), mesh.get_proc(), max_proc));

                for (int j = 0; j < edge_dofs; j++) {
                    Integer dof_sfc =
                        process_edge_node<ElemType>(mars::get_sfc_from_octant<ElemType>, start, direction, j);
                    f(mesh, i, dof_sfc, max_proc, direction);
                }
            }
        }

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T != ElementType::Hex8, void> edge_iterate(const Integer sfc,
                                                                                                const Mesh mesh,
                                                                                                const Integer i,
                                                                                                F f) {}

        template <typename CornerF, typename FaceF, typename VolumeF, typename EdgeF>
        static MARS_INLINE_FUNCTION void elem_dof_iterate(const Integer sfc,
                                                          const Mesh mesh,
                                                          const Integer i,
                                                          CornerF cf,
                                                          FaceF ff,
                                                          VolumeF vf,
                                                          EdgeF ef) {
            if (corner_dofs > 0) {
                corner_iterate(sfc, mesh, i, cf);
            }
            if (face_dofs > 0) {
                // templetized with 1 value means filtered on the direction (No dir is forwarded to the f call).
                face_dir_iterate<1>(sfc, mesh, i, ff);
            }
            if (volume_dofs > 0) {
                volume_iterate(sfc, mesh, i, vf);
            }
            if (edge_dofs > 0) {
                edge_iterate(sfc, mesh, i, IterateFilter<EdgeF, 1>(ef));
            }
        }

        struct IdentifyBoundaryDofPerRank {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = mesh.get_boundary_sfc(i);
                const Integer proc = mesh.get_proc();

                constexpr bool BoundaryIter = false;
                // used for corners and edges that are based on one ring nbh logic.
                auto rb = RankBoundary<BoundaryIter>(rank_boundary, sfc_to_locally_owned, map, proc);
                auto vrb =
                    RankBoundary<BoundaryIter, DofLabel::lVolume>(rank_boundary, sfc_to_locally_owned, map, proc);
                // use frb for the edge logic as well. Same logic different input from the iterate.
                elem_dof_iterate(sfc, mesh, i, rb, rb, vrb, rb);
            }

            MARS_INLINE_FUNCTION
            IdentifyBoundaryDofPerRank(Mesh m,
                                       ViewMatrixTypeLeft<bool> npb,
                                       ViewVectorType<Integer> mp,
                                       UnorderedMap<Integer, Integer> l)
                : mesh(m), rank_boundary(npb), map(mp), sfc_to_locally_owned(l) {}

            Mesh mesh;

            ViewMatrixTypeLeft<bool> rank_boundary;
            ViewVectorType<Integer> map;
            UnorderedMap<Integer, Integer> sfc_to_locally_owned;
        };

        void build_boundary_dof_sets(ViewVectorType<Integer> &scan_boundary,
                                     ViewVectorType<Integer> &boundary_lsfc_index,
                                     ViewVectorType<Integer> &sender_ranks_scan,
                                     const Integer nbh_rank_size) {
            using namespace Kokkos;

            printf("Bilding the Boundary DOF Sets!\n");
            const Integer size = get_mesh().get_view_boundary().extent(0);

            Integer xDim = get_mesh().get_XDim();
            Integer yDim = get_mesh().get_YDim();
            Integer zDim = get_mesh().get_ZDim();

            const Integer chunk_size_ = global_dof_enum.get_elem_size();
            ViewMatrixTypeLeft<bool> rank_boundary("count_per_proc", chunk_size_, nbh_rank_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("identify_boundary_predicate",
                                 size,
                                 IdentifyBoundaryDofPerRank(get_mesh(),
                                                            rank_boundary,
                                                            sender_ranks_scan,
                                                            get_global_dof_enum().get_sfc_to_local_map()));

            /* perform a scan for each row with the sum at the end for each rank */
            ViewMatrixTypeLeft<Integer> rank_scan("rank_scan", chunk_size_ + 1, nbh_rank_size);
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
        static MARS_INLINE_FUNCTION Integer process_corner(const Mesh mesh, const Octant &o) {
            Integer xDim = mesh.get_XDim();
            Integer yDim = mesh.get_YDim();
            Integer zDim = mesh.get_ZDim();

            Integer max_proc = -1;

            Octant one_ring[simplex_type::ElemType];
            o.one_ring_corner_nbhs<Type>(one_ring, xDim, yDim, zDim, mesh.is_periodic());

            for (int k = 0; k < Type; k++) {
                if (one_ring[k].is_valid()) {
                    Integer enc_oc = mars::get_sfc_from_octant<Type>(one_ring[k]);
                    /* check the proc that owns the corner and then decide on the predicate.
                        This works because the corner defines an element and this element is
                        always on the largest proc number due to the z order partitioning*/
                    Integer owner_proc = find_owner_processor(mesh.get_view_gp(), enc_oc, 2, mesh.get_proc());
                    /* assert(owner_proc >= 0); */

                    if (owner_proc > max_proc) {
                        max_proc = owner_proc;
                    }
                }
            }

            return max_proc;
        }

        template <typename F>
        static MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc, const Mesh mesh, const Integer index, F f) {
            Octant oc = mesh.octant_from_sfc(sfc);
            corner_octant_transform(oc, mesh, index, f);
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        corner_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f) {
            for (int i = 0; i < corner_dofs + 1; i++) {
                for (int j = 0; j < corner_dofs + 1; j++) {
                    Octant o;
                    // get the next corner using the elem sfc
                    o.x = oc.x + i;
                    o.y = oc.y + j;

                    /* Integer mp = process_corner<ElemType>(mesh, o); */
                    Integer max_proc = -1;
                    mesh.get_one_ring_corner_nbhs(o, MaxOneRingProc(mesh.get_view_gp(), mesh.get_proc(), max_proc));

                    /* assert(mp == max_proc); */
                    // convert the octant value into the new nodal sfc system
                    o.x *= degree;
                    o.y *= degree;
                    Integer dof_sfc = mars::get_sfc_from_octant<ElemType>(o);

                    f(mesh, index, dof_sfc, max_proc);
                }
            }
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        corner_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f) {
            for (int i = 0; i < corner_dofs + 1; i++) {
                for (int j = 0; j < corner_dofs + 1; j++) {
                    for (int k = 0; k < corner_dofs + 1; k++) {
                        Octant o;
                        // get the next corner using the elem sfc
                        o.x = oc.x + i;
                        o.y = oc.y + j;
                        o.z = oc.z + k;

                        /* Integer mp = process_corner<ElemType>(mesh, o); */
                        Integer max_proc = -1;
                        mesh.get_one_ring_corner_nbhs(o,
                                                       MaxOneRingProc(mesh.get_view_gp(), mesh.get_proc(), max_proc));
                        /* assert(mp == max_proc); */

                        // convert the octant value into the new nodal sfc system
                        o.x *= degree;
                        o.y *= degree;
                        o.z *= degree;
                        Integer dof_sfc = mars::get_sfc_from_octant<ElemType>(o);

                        f(mesh, index, dof_sfc, max_proc);
                    }
                }
            }
        }

        // generic local predicate functor to be used for edge, face and corners. The volume is specialized.
        template <bool Ghost, Integer Label>
        struct LocalGlobalPredicate {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;
            Integer proc;

            MARS_INLINE_FUNCTION
            LocalGlobalPredicate(ViewVectorType<bool> lp,
                                 ViewVectorType<bool> gp,
                                 ViewVectorType<bool> npps,
                                 ViewVectorType<bool> nppr,
                                 Integer p)
                : local_predicate(lp),
                  global_predicate(gp),
                  nbh_proc_predicate_send(npps),
                  nbh_proc_predicate_recv(nppr),
                  proc(p) {}

            MARS_INLINE_FUNCTION
            void lg_predicate(const Mesh mesh,
                              const Integer i,
                              const Integer sfc,
                              const Integer owner_proc,
                              std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    local_predicate(sfc) = 1;
                    nbh_proc_predicate_send(elem_sfc_proc) = 1;
                    nbh_proc_predicate_recv(elem_sfc_proc) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void lg_predicate(const Mesh mesh,
                              const Integer i,
                              const Integer sfc,
                              const Integer owner_proc,
                              std::false_type) const {
                /* if the neighbor element is ghost then check if the processor
                    is less than the owner. This is how the dofs are partitioned*/
                if (proc >= owner_proc) {
                    global_predicate(sfc) = 1;
                }
                local_predicate(sfc) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer dof_sfc, const Integer owner_proc) const {
                lg_predicate(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct LocalGlobalPredicate<Ghost, DofLabel::lVolume> {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;

            MARS_INLINE_FUNCTION
            LocalGlobalPredicate(ViewVectorType<bool> lp, ViewVectorType<bool> gp)
                : local_predicate(lp), global_predicate(gp) {}

            MARS_INLINE_FUNCTION
            void volume_predicate(const Integer sfc, std::true_type) const { local_predicate(sfc) = 1; }

            MARS_INLINE_FUNCTION
            void volume_predicate(const Integer sfc, std::false_type) const {
                local_predicate(sfc) = 1;
                global_predicate(sfc) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer enc_oc) const {
                volume_predicate(enc_oc, std::integral_constant<bool, Ghost>{});
            }
        };

        // sfc | octant_from_sfc | volume_octant_transform. Composable functions.
        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc, const Mesh mesh, const Integer index, F f) {
            Octant oc = mesh.octant_from_sfc(sfc);
            volume_octant_transform(oc, mesh, index, f, mars::get_sfc_from_octant<ElemType>);
        }

        template <typename F, typename G, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        volume_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f, G g) {
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
        volume_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f, G g) {
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

        // Iterate filter redirecting the call with filtered number of parameters to different functors.
        template <typename F, Integer Params = 0>
        struct IterateFilter {
            F f;

            MARS_INLINE_FUNCTION
            IterateFilter(F fun) : f(fun) {}

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                f(mesh, i, dof_sfc, owner_proc, dir);
            }
        };

        // Could have specialized on a bool for dir or not dir but left it open for further filters.
        // For example one might need a face_iterate with only mesh, i, dofs_sfc parameteres.
        template <typename F>
        struct IterateFilter<F, 1> {
            F f;

            MARS_INLINE_FUNCTION
            IterateFilter(F fun) : f(fun) {}

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                f(mesh, i, dof_sfc, owner_proc);
            }
        };

        template <Integer Params, typename F, Integer Type = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, void>
        face_dir_iterate(const Integer sfc, const Mesh mesh, const Integer i, F f) {
            auto filter = IterateFilter<F, Params>(f);
            // dir 0 = x, 1 = y
            face_iterate<0>(sfc, mesh, i, filter);
            face_iterate<1>(sfc, mesh, i, filter);
        }

        template <Integer Params, typename F, Integer Type = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, void>
        face_dir_iterate(const Integer sfc, const Mesh mesh, const Integer i, F f) {
            auto filter = IterateFilter<F, Params>(f);
            // dir 0 = x, 1 = y
            // dir 0 = x, 1 = y, 2 = z
            face_iterate<0>(sfc, mesh, i, filter);
            face_iterate<1>(sfc, mesh, i, filter);
            face_iterate<2>(sfc, mesh, i, filter);
        }

        // careful: do not use this function for face iteration when expensive operations are involved. Useful mainly
        // for predicate builder. instead use the locally_owned_face_dof array to iterate over face dof sfcs.
        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const Integer sfc, const Mesh mesh, const Integer index, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh.octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_dofs; j++) {
                    Integer dof_sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(mesh, index, dof_sfc, owner_proc, dir);
                }
            }
        }

        /* template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const Integer sfc, const Mesh mesh, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh.octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_dofs; j++) {
                    Integer dof_sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(dof_sfc, owner_proc, dir);
                }
            }
        } */

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == true, Octant>::type get_octant_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return mesh.get_ghost_octant(index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == false, Octant>::type get_octant_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return mesh.get_octant(index);
        }

        static MARS_INLINE_FUNCTION Integer get_ghost_sfc(const Mesh mesh, const Integer index) {
            return mesh.get_ghost_sfc(index);
        }
        static MARS_INLINE_FUNCTION Integer get_local_sfc(const Mesh mesh, const Integer index) {
            return mesh.get_sfc(index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == true, Integer>::type get_sfc_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return get_ghost_sfc(mesh, index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == false, Integer>::type get_sfc_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return get_local_sfc(mesh, index);
        }

        template <bool Ghost>
        struct BuildLocalGlobalPredicate {

            MARS_INLINE_FUNCTION
            BuildLocalGlobalPredicate(Mesh m,
                                      ViewVectorType<bool> lp,
                                      ViewVectorType<bool> gp,
                                      ViewVectorType<bool> npbs,
                                      ViewVectorType<bool> npbr)
                : mesh(m),
                  local_predicate(lp),
                  global_predicate(gp),
                  nbh_proc_predicate_send(npbs),
                  nbh_proc_predicate_recv(npbr) {}

            Mesh mesh;
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh.get_proc();
                auto cp = LocalGlobalPredicate<Ghost, DofLabel::lCorner>(
                    local_predicate, global_predicate, nbh_proc_predicate_send, nbh_proc_predicate_recv, proc);

                auto fp = LocalGlobalPredicate<Ghost, DofLabel::lFace>(
                    local_predicate, global_predicate, nbh_proc_predicate_send, nbh_proc_predicate_recv, proc);

                // fp logic is the same as the edge one. Just different label.
                auto ep = LocalGlobalPredicate<Ghost, DofLabel::lEdge>(
                    local_predicate, global_predicate, nbh_proc_predicate_send, nbh_proc_predicate_recv, proc);

                auto vp = LocalGlobalPredicate<Ghost, DofLabel::lVolume>(local_predicate, global_predicate);

                elem_dof_iterate(sfc, mesh, i, cp, fp, vp, ep);
            }
        };

        void print_dofs(const int rank = 0) {
            auto handler = *this;

            SFC<ElemType> gdof = get_global_dof_enum();
            Kokkos::parallel_for(
                "for", gdof.get_elem_size(), MARS_LAMBDA(const int i) {
                    const Integer sfc_elem = handler.owned_to_sfc(i);
                    const Dof global_dof = handler.sfc_to_global_dof(sfc_elem);

                    /*Integer gid = -1, proc = -1;

                    if (is_owned_dof_sfc(sfc_elem)) {
                        proc = get_mesh_manager().get_mesh()->get_proc();
                        const Integer sfc_lid = global_dof_enum.get_view_sfc_to_local()(sfc_elem);
                        gid = sfc_lid + global_dof_offset(proc);
                    }*/

                    Octant o = handler.get_octant_from_sfc(sfc_elem);

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
                "for", dof.get_elem_size(), MARS_CLASS_LAMBDA(const int i) {
                    const Integer sfc_elem = local_to_sfc(i);
                    const Dof global_dof = local_to_global_dof(i);

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
            Integer xDim = get_mesh().get_XDim();
            Integer yDim = get_mesh().get_YDim();
            Integer zDim = get_mesh().get_ZDim();

            const Integer size = get_mesh().get_chunk_size();
            const Integer ghost_size = get_mesh().get_view_ghost().extent(0);

            /* TODO: xdim and ydim should be changed to max xdim and ydim
             * of the local partition to reduce memory footprint */
            local_dof_enum = SFC<simplex_type::ElemType>(degree * xDim, degree * yDim, degree * zDim);
            global_dof_enum = SFC<simplex_type::ElemType>(degree * xDim, degree * yDim, degree * zDim);

            ViewVectorType<bool> local_predicate("lpred", local_dof_enum.get_all_range());
            ViewVectorType<bool> global_predicate("gpred", global_dof_enum.get_all_range());

            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("lg_predicate",
                                 size,
                                 BuildLocalGlobalPredicate<false>(
                                     get_mesh(), local_predicate, global_predicate, npbs, npbr));

            // Iterate through ghost sfc and enumerate
            Kokkos::parallel_for("lg_predicate_from_ghost",
                                 ghost_size,
                                 BuildLocalGlobalPredicate<true>(
                                     get_mesh(), local_predicate, global_predicate, npbs, npbr));

            build_sfc_from_predicate(local_dof_enum, local_predicate);
            build_sfc_from_predicate(global_dof_enum, global_predicate);

            const int rank_size = num_ranks(context);

            ViewVectorType<Integer> global_dof_size_per_rank("global_dof_size_per_rank", rank_size);
            global_dof_offset = ViewVectorType<Integer>("global_dof_offset", rank_size + 1);

            const Integer global_size = global_dof_enum.get_elem_size();
            context->distributed->gather_all_view(global_size, global_dof_size_per_rank);

            incl_excl_scan(0, rank_size, global_dof_size_per_rank, global_dof_offset);
        }

        void build_sfc_from_predicate(SFC<ElemType> &sfc_enum, const ViewVectorType<bool> &predicate) {
            ViewVectorType<Integer> sfc_to_local("sfc_to_local dof handler", sfc_enum.get_all_range());
            sfc_enum.compact_elements(sfc_to_local, predicate);
            sfc_enum.generate_sfc_to_local_map();
        }

        // Memory efficient implementation but poor performance due to the segmented scan.
        /* void build_sfc_from_predicate(SFC<ElemType> &sfc_enum,
                                      const ViewVectorType<bool> &predicate,
                                      const ViewVectorType<Integer> &label) {
            auto all_range = sfc_enum.get_all_range();

            Integer nr_elements = 0;
            Kokkos::parallel_reduce(
                "Number of SFC codes",
                all_range,
                MARS_LAMBDA(const Integer i, Integer &update) { update += predicate(i); },
                nr_elements);

            sfc_enum.reserve_elements(nr_elements);
            sfc_enum.reserve_element_labels(nr_elements);

            //compact using the segmented scan to avoid storing all_range vector.
            auto elms = sfc_enum.get_view_elements();
            auto lbl = sfc_enum.get_view_element_labels();
            segmented_scan(
                all_range, predicate, MARS_LAMBDA(const Integer count, const Integer index) {
                    if (predicate(index) == 1) {
                        elms(count) = index;
                        lbl(count) = label(index);
                    }
                });

            sfc_enum.generate_sfc_to_local_map();
        } */

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
            std::cout << "DofHandler:MPI send receive for the ghost dofs layer started..." << std::endl;

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

        template <class UMAP>
        static MARS_INLINE_FUNCTION Integer get_value_in_map(const UMAP map, const Integer sfc) {
            const auto it = map.find(sfc);
            return map.value_at(it);
        }

        template <bool Ghost>
        struct FaceOrientDof {
            ViewVectorType<Integer> dir;
            UnorderedMap<Integer, Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            FaceOrientDof(ViewVectorType<Integer> d, UnorderedMap<Integer, Integer> l, Integer p)
                : dir(d), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            Integer get_local_from_sfc(const Integer sfc) const { return get_value_in_map(sfc_to_local, sfc); }

            MARS_INLINE_FUNCTION
            void face_orient_dof(const Mesh mesh,
                                 const Integer i,
                                 const Integer sfc,
                                 const Integer owner_proc,
                                 const Integer d,
                                 std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    Integer index = get_local_from_sfc(sfc);
                    dir(index) = d;
                }
            }

            MARS_INLINE_FUNCTION
            void face_orient_dof(const Mesh mesh,
                                 const Integer i,
                                 const Integer sfc,
                                 const Integer owner_proc,
                                 const Integer d,
                                 std::false_type) const {
                Integer index = get_local_from_sfc(sfc);
                dir(index) = d;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
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
                const Integer proc = mesh.get_proc();

                if (face_dofs > 0) {
                    FaceOrientDof<Ghost> fp = FaceOrientDof<Ghost>(face_dir, stl_map, proc);
                    // templatized with 0 values means no filter
                    face_dir_iterate<0>(sfc, mesh, i, fp);
                }
            }

            MARS_INLINE_FUNCTION
            OrientDofs(Mesh m, ViewVectorType<Integer> sd, UnorderedMap<Integer, Integer> sl)
                : mesh(m), face_dir(sd), stl_map(sl) {}

            Mesh mesh;
            ViewVectorType<Integer> face_dir;
            UnorderedMap<Integer, Integer> stl_map;
        };

        void build_local_orientation() {
            const Integer size = get_mesh().get_chunk_size();
            const Integer ghost_size = get_mesh().get_view_ghost().extent(0);

            const Integer local_size = get_local_dof_enum().get_elem_size();

            local_dof_enum.init_element_orientations(local_size);

            /* generate the sfc for the local and global dofs containing the generation locally
                     for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("orient_dofs",
                                 size,
                                 OrientDofs<false>(get_mesh(),
                                                   get_local_dof_enum().get_view_element_orientations(),
                                                   get_local_dof_enum().get_sfc_to_local_map()));

            Kokkos::parallel_for("orient_ghost_dofs",
                                 ghost_size,
                                 OrientDofs<true>(get_mesh(),
                                                  get_local_dof_enum().get_view_element_orientations(),
                                                  get_local_dof_enum().get_sfc_to_local_map()));
        }

        // generic local predicate functor to be used for edge, face and corners. The volume is specialized.
        template <bool Ghost, Integer Label>
        struct LocalGlobalLabel {
            SFC<ElemType> local_enum;
            SFC<ElemType> global_enum;
            Integer proc;

            MARS_INLINE_FUNCTION
            LocalGlobalLabel(SFC<ElemType> le, SFC<ElemType> ge, Integer p)
                : local_enum(le), global_enum(ge), proc(p) {}

            MARS_INLINE_FUNCTION
            void set_label_from_sfc(const Mesh mesh,
                                    const Integer i,
                                    const Integer sfc,
                                    const Integer owner_proc,
                                    std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                    local_enum.set_label(index, Label);
                }
            }

            MARS_INLINE_FUNCTION
            void set_label_from_sfc(const Mesh mesh,
                                    const Integer i,
                                    const Integer sfc,
                                    const Integer owner_proc,
                                    std::false_type) const {
                /* if the neighbor element is ghost then check if the processor
                    is less than the owner. This is how the dofs are partitioned*/
                if (proc >= owner_proc) {
                    auto gindex = get_value_in_map(global_enum.get_sfc_to_local_map(), sfc);
                    global_enum.set_label(gindex, Label);
                }
                auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                local_enum.set_label(index, Label);
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer dof_sfc, const Integer owner_proc) const {
                set_label_from_sfc(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct LocalGlobalLabel<Ghost, DofLabel::lVolume> {
            SFC<ElemType> local_enum;
            SFC<ElemType> global_enum;

            MARS_INLINE_FUNCTION
            LocalGlobalLabel(SFC<ElemType> le, SFC<ElemType> ge) : local_enum(le), global_enum(ge) {}

            MARS_INLINE_FUNCTION
            void volume_label(const Integer sfc, std::true_type) const {
                auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                local_enum.set_label(index, DofLabel::lVolume);
            }

            MARS_INLINE_FUNCTION
            void volume_label(const Integer sfc, std::false_type) const {
                auto gindex = get_value_in_map(global_enum.get_sfc_to_local_map(), sfc);
                global_enum.set_label(gindex, DofLabel::lVolume);
                auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                local_enum.set_label(index, DofLabel::lVolume);
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer enc_oc) const {
                volume_label(enc_oc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct BuildLocalGlobalLabel {

            MARS_INLINE_FUNCTION
            BuildLocalGlobalLabel(Mesh m, SFC<ElemType> le, SFC<ElemType> ge)
                : mesh(m), local_enum(le), global_enum(ge) {}

            Mesh mesh;
            SFC<ElemType> local_enum;
            SFC<ElemType> global_enum;

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh.get_proc();
                auto cp = LocalGlobalLabel<Ghost, DofLabel::lCorner>(local_enum, global_enum, proc);

                auto fp = LocalGlobalLabel<Ghost, DofLabel::lFace>(local_enum, global_enum, proc);

                // fp logic is the same as the edge one. Just different label.
                auto ep = LocalGlobalLabel<Ghost, DofLabel::lEdge>(local_enum, global_enum, proc);

                auto vp = LocalGlobalLabel<Ghost, DofLabel::lVolume>(local_enum, global_enum);

                elem_dof_iterate(sfc, mesh, i, cp, fp, vp, ep);
            }
        };

        void label_local_global_dofs() {
            const Integer size = get_mesh().get_chunk_size();
            const Integer ghost_size = get_mesh().get_view_ghost().extent(0);

            local_dof_enum.reserve_element_labels(local_dof_enum.get_elem_size());
            global_dof_enum.reserve_element_labels(global_dof_enum.get_elem_size());
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for("label_local",
                                 size,
                                 BuildLocalGlobalLabel<false>(
                                     get_mesh(), get_local_dof_enum(), get_global_dof_enum()));

            // Iterate through ghost sfc and enumerate
            Kokkos::parallel_for("label_ghost",
                                 ghost_size,
                                 BuildLocalGlobalLabel<true>(
                                     get_mesh(), get_local_dof_enum(), get_global_dof_enum()));
        }

        virtual void enumerate_dofs() {
            printf("Enumerating the Degrees of Freedom on the Mesh with Degree: %li\n", degree);

            const auto &context = get_context();

            const Integer rank_size = num_ranks(context);
            const int proc_num = rank(context);

            ViewVectorType<bool> nbh_proc_predicate_send("send_to", rank_size);
            ViewVectorType<bool> nbh_proc_predicate_recv("receive_from", rank_size);
            ViewVectorType<Integer> proc_scan_send("nbh_scan_send", rank_size + 1);
            ViewVectorType<Integer> proc_scan_recv("nbh_scan_recv", rank_size + 1);

            build_lg_predicate(context, nbh_proc_predicate_send, nbh_proc_predicate_recv);
            label_local_global_dofs();
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

            Kokkos::fence();

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
            Mesh mesh;
            UnorderedMap<Integer, Dof> glgm;
            ViewVectorType<Integer> global_dof_offset;
            ViewVectorType<Integer> scan_recv_proc;

            ViewVectorType<Integer> ghost_dofs;
            ViewVectorType<Integer> ghost_dofs_index;

            MARS_INLINE_FUNCTION
            BuildLocalToGlobalGhostMap(Mesh m,
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
                const Integer owner_proc = find_owner_processor(scan_recv_proc, i, 1, mesh.get_proc());
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
                                 BuildLocalToGlobalGhostMap(get_mesh(),
                                                            ghost_local_to_global_map,
                                                            global_dof_offset,
                                                            scan_recv,
                                                            ghost_dofs_sfc,
                                                            ghost_dofs_index));
            /* In the end the size of the map should be as the size of the ghost_dofs.
             * Careful map size  is not capacity */
            assert(size == ghost_local_to_global_map.size());
            printf("Built the Ghost Local to Global Map.\n");
        }

        /* MARS_INLINE_FUNCTION
        bool is_local(const Integer sfc) const {
            auto stl = local_dof_enum.get_view_sfc_to_local();
            if ((sfc + 1) >= stl.extent(0)) return false;
            [>use the sfc to local which is the scan of the predicate.
             * To get the predicate value the difference with the successive index is needed.<]
            const Integer pred_value = stl(sfc + 1) - stl(sfc);
            return (pred_value > 0);
        } */

        template <typename MT>
        static MARS_INLINE_FUNCTION bool exists_in_map(const UnorderedMap<MT, MT> map, const MT value) {
            const auto it = map.find(value);
            return map.valid_at(it);
        }

        MARS_INLINE_FUNCTION
        bool is_local(const Integer sfc) const {
            auto map = get_local_dof_enum().get_sfc_to_local_map();
            return exists_in_map(map, sfc);
        }

        MARS_INLINE_FUNCTION
        bool is_owned_dof_sfc(const Integer sfc) const {
            auto map = get_global_dof_enum().get_sfc_to_local_map();
            return exists_in_map(map, sfc);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION bool is_owned(const Integer local) const {
            const Integer sfc = local_to_sfc<B>(local);
            return is_owned_dof_sfc(sfc);
        }

        MARS_INLINE_FUNCTION Integer get_eval_value_in_local_map(const Integer sfc) const {
            auto map = get_local_dof_enum().get_sfc_to_local_map();
            const auto it = map.find(sfc);
            if (map.valid_at(it)) {
                return map.value_at(it);
            }
            return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_eval_value_in_global_map(const Integer sfc) const {
            auto map = get_global_dof_enum().get_sfc_to_local_map();
            const auto it = map.find(sfc);
            if (map.valid_at(it)) {
                return map.value_at(it);
            }
            return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Integer>::type local_to_owned_dof(
            const Integer local) const {
            const Integer comp_local = compute_component<B>(local);
            const Integer sfc = local_to_sfc<B>(local);
            const Integer base_owned = get_eval_value_in_global_map(sfc);
            return compute_block_index<B>(base_owned, comp_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Integer>::type local_to_owned_dof(
            const Integer local) const {
            const Integer sfc = local_to_sfc<B>(local);
            return get_eval_value_in_global_map(sfc);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_owned_index(const Integer local) const {
            return local_to_owned_dof<B>(local);
        }

        MARS_INLINE_FUNCTION
        Dof sfc_to_global_dof(const Integer sfc) const {
            Dof dof;
            const Integer sfc_lid = get_eval_value_in_global_map(sfc);
            if (sfc_lid != INVALID_INDEX) {
                const Integer proc = get_mesh().get_proc();
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

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Dof>::type local_to_global_dof(const Integer local) const {
            const Integer local_sfc = local_to_sfc<B>(local);
            return sfc_to_global_dof(local_sfc);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Dof>::type local_to_global_dof(const Integer local) const {

            // use the local to sfc view (elements) to get the sfc of the local numbering.
            const Integer component = compute_component<B>(local);

            const Integer local_sfc = local_to_sfc<B>(local);
            auto dof = sfc_to_global_dof(local_sfc);
            dof.set_gid(compute_block_index<B>(dof.get_gid(), component));
            return dof;
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_global(const Integer sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_global(const Integer local) const {
            Dof dof = local_to_global_dof<B>(local);
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

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_global_proc(const Integer local) const {
            Dof dof = local_to_global_dof<B>(local);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_local_label(const Integer local) const {
            const Integer base_local = compute_base<B>(local);
            return get_local_dof_enum().get_view_element_labels()(base_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_sfc(const Integer local) const {
            auto base = compute_base<B>(local);
            return get_local_dof_enum().get_view_elements()(base);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer owned_to_sfc(const Integer owned) const {
            auto base = compute_base<B>(owned);
            return get_global_dof_enum().get_view_elements()(base);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_owned(const Integer local) const {
            return local_to_owned_dof<B>(local);
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_owned(const Integer sfc) const { return get_eval_value_in_global_map(sfc); }

        MARS_INLINE_FUNCTION
        Integer sfc_to_owned(const Integer sfc, const Integer component) const {
            return compute_block_index<Block>(get_eval_value_in_global_map(sfc), component);
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_local(const Integer sfc) const { return get_eval_value_in_local_map(sfc); }

        MARS_INLINE_FUNCTION
        Integer sfc_to_local(const Integer sfc, const Integer component) const {
            return compute_block_index<Block>(get_eval_value_in_local_map(sfc), component);
        }

        MARS_INLINE_FUNCTION
        Integer get_local_index(const Integer sfc) const { return sfc_to_local(sfc); }

        MARS_INLINE_FUNCTION
        Integer get_local_index(const Integer sfc, const Integer component) const {
            return sfc_to_local(sfc, component);
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION void get_dof_coordinates_from_sfc(const Integer sfc, double *point) const {
            get_vertex_coordinates_from_sfc<Type>(sfc,
                                                  point,
                                                  get_local_dof_enum().get_XDim(),
                                                  get_local_dof_enum().get_YDim(),
                                                  get_local_dof_enum().get_ZDim());
        }

        template <Integer Type, Integer B = Block_>
        MARS_INLINE_FUNCTION void get_dof_coordinates_from_local(const Integer local, double *point) const {
            const Integer sfc = local_to_sfc<B>(local);
            return get_dof_coordinates_from_sfc<Type>(sfc, point);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION void get_local_dof_coordinates(const Integer local, double *point) const {
            get_dof_coordinates_from_local<ElemType, B>(local, point);
        }

        MARS_INLINE_FUNCTION void get_base_local_dof_coordinates(const Integer local, double *point) const {
            get_dof_coordinates_from_local<ElemType, 1>(local, point);
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION bool boundary_sfc(const Integer sfc, const Integer FaceNr = -1) const {
            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            return is_boundary_sfc<Type>(sfc, xdim, ydim, zdim, FaceNr);
        }

        template <Integer Type, Integer B = Block_>
        MARS_INLINE_FUNCTION bool is_boundary(const Integer local, const Integer FaceNr = -1) const {
            const Integer sfc = local_to_sfc<B>(local);
            return boundary_sfc<Type>(sfc, FaceNr);
        }

        MARS_INLINE_FUNCTION bool is_boundary_dof(const Integer local, const Integer FaceNr = -1) const {
            return is_boundary<ElemType>(local, FaceNr);
        }

        template <typename F>
        void boundary_owned_dof_iterate(F f, const std::string side = "all") {
            const Integer side_value = map_side_to_value<ElemType>(side);

            auto handler = *this;
            owned_dof_iterate(MARS_LAMBDA(const Integer i) {
                const Integer sfc = handler.owned_to_sfc(i);
                if (handler.template boundary_sfc<ElemType>(sfc, side_value)) {
                    f(i, sfc);
                }
            });
        }

        MARS_INLINE_FUNCTION
        static bool is_separated_block(const Integer block, const Integer component) {
            return (block == -1 || block == component);
        }

        template <typename F>
        void boundary_dof_iterate(F f, const Integer bl) {
            boundary_dof_iterate(f, "all", bl);
        }

        template <typename F>
        void boundary_dof_iterate(F f, const std::string side = "all", const Integer bl = -1) {
            const Integer side_value = map_side_to_value<ElemType>(side);

            auto handler = *this;
            owned_dof_iterate(MARS_LAMBDA(const Integer i) {
                const Integer sfc = handler.owned_to_sfc(i);
                auto comp = handler.compute_component(i);
                if (is_separated_block(bl, comp) && handler.template boundary_sfc<ElemType>(sfc, side_value)) {
                    const Integer local = handler.sfc_to_local(sfc);
                    f(handler.compute_block_index(local, comp));
                }
            });
        }

        template <Integer B = Block_>
        Integer get_dof_size() const {
            return get_block<B>() * local_dof_enum.get_elem_size();
        }

        MARS_INLINE_FUNCTION Integer get_base_dof_size() const { return local_dof_enum.get_elem_size(); }

        // no need to be vector valued since it will give the same result.
        MARS_INLINE_FUNCTION
        Integer get_local_dof(const Integer i) const {
            assert(i < get_dof_size());
            return i;
        }

        // get the local dof of the owned index.
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_owned_dof(const Integer owned_dof) const {
            auto base = compute_base<B>(owned_dof);
            auto component = compute_component<B>(owned_dof);
            const Integer sfc = get_global_dof_enum().get_view_elements()(base);
            return compute_block_index<B>(sfc_to_local(sfc), component);
        }

        MARS_INLINE_FUNCTION
        Integer get_base_owned_dof_size() const { return global_dof_enum.get_elem_size(); }

        template <Integer B = Block_>
        Integer get_owned_dof_size() const {
            return get_block<B>() * global_dof_enum.get_elem_size();
        }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_local_dof_enum() const { return local_dof_enum; }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_global_dof_enum() const { return global_dof_enum; }

        MARS_INLINE_FUNCTION
        Integer get_global_dof_offset(const Integer proc) const { return global_dof_offset(proc); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_global_dof_offset() const { return global_dof_offset; }

        Integer get_global_base_dof_size() const {
            const Integer rank_size = num_ranks(get_context());

            auto ss = subview(global_dof_offset, rank_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            return h_ss();
        }

        template <Integer B = Block_>
        Integer get_global_dof_size() const {
            return get_block<B>() * get_global_base_dof_size();
        }

        MARS_INLINE_FUNCTION
        const UnorderedMap<Integer, Dof> &get_ghost_lg_map() const { return ghost_local_to_global_map; }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_boundary_dof(const Integer i) const {
            const auto base = compute_base<B>(i);
            const auto component = compute_component<B>(i);
            const auto local = sfc_to_local(get_boundary_dofs()(base));
            assert(INVALID_INDEX != local);
            return compute_block_index<B>(local, component);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_boundary_local_index(const Integer i) const {
            return get_boundary_dof<B>(i);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_boundary_dof_size() const {
            return get_block<B>() * get_boundary_dofs().extent(0);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_ghost_dof(const Integer i) const {
            const auto base = compute_base<B>(i);
            const auto component = compute_component<B>(i);
            const auto local = sfc_to_local(get_ghost_dofs()(base));
            assert(INVALID_INDEX != local);
            return compute_block_index<B>(local, component);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_ghost_local_index(const Integer i) const {
            return get_ghost_dof<B>(i);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_ghost_dof_size() const {
            return get_block<B>() * get_ghost_dofs().extent(0);
        }

        template <typename F>
        void ghost_iterate(F f) const {
            Kokkos::parallel_for("ghost_dof_iter", get_ghost_dof_size(), f);
        }

        MARS_INLINE_FUNCTION
        const UnorderedMap<Integer, Integer> get_local_dof_map() const {
            return get_local_dof_enum().get_sfc_to_local_map();
        }

        /* MARS_INLINE_FUNCTION
        const Integer get_local_dof_map(const Integer sfc) const { return sfc_to_local(sfc); } */
        /* MM get_mesh_manager() const { return mesh_manager; } */
        MARS_INLINE_FUNCTION
        Mesh get_mesh() const { return mesh; }

        MARS_INLINE_FUNCTION
        const Integer get_proc() const { return get_mesh().get_proc(); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_recv() const { return scan_recv; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_send() const { return scan_send; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_recv_mirror() const { return scan_recv_mirror; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_send_mirror() const { return scan_send_mirror; }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_orientation(const Integer local_dof) const {
            const Integer base_local = compute_base<B>(local_dof);
            return get_local_dof_enum().get_orientation(base_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_label(const Integer local_dof) const {
            const Integer base_local = compute_base<B>(local_dof);
            return get_local_dof_enum().get_label(base_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_owned_label(const Integer owned_dof) const {
            const Integer base_owned = compute_base<B>(owned_dof);
            return get_global_dof_enum().get_label(base_owned);
        }

        const context &get_context() const { return get_mesh().get_context(); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_boundary_dofs() const { return boundary_dofs_sfc; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_ghost_dofs() const { return ghost_dofs_sfc; }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Octant get_octant_from_local(const Integer local) const {
            const Integer sfc = local_to_sfc<B>(local);
            return get_octant_from_sfc(sfc);
        }

        MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const Integer sfc) const {
            return mars::get_octant_from_sfc<ElemType>(sfc);
        }

        MARS_INLINE_FUNCTION Integer get_sfc_from_octant(const Octant &o) const {
            return mars::get_sfc_from_octant<ElemType>(o);
        }

        /* MARS_INLINE_FUNCTION Integer get_global_from_octant(const Octant &o) const {
            const Integer sfc = get_sfc_from_octant(o);
            const Integer local = is_local(sfc) ? sfc_to_local(sfc) : INVALID_INDEX;
            return local_to_global(local);
        } */

        MARS_INLINE_FUNCTION Integer get_local_sfc_from_octant(const Octant &o) const {
            const Integer sfc = get_sfc_from_octant(o);
            return is_local(sfc) ? sfc : INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_local_from_octant(const Octant &o) const {
            const Integer sfc = get_sfc_from_octant(o);
            return is_local(sfc) ? sfc_to_local(sfc) : INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_local_from_octant(const Octant &o, const Integer component) const {
            auto local = get_local_from_octant(o);
            return compute_block_index<Block>(local, component);
        }

        MARS_INLINE_FUNCTION
        const Integer get_XMax() const { return Degree * get_mesh().get_XDim(); }

        MARS_INLINE_FUNCTION
        const Integer get_YMax() const { return Degree * get_mesh().get_YDim(); }

        MARS_INLINE_FUNCTION
        const Integer get_ZMax() const { return Degree * get_mesh().get_ZDim(); }

        MARS_INLINE_FUNCTION
        void set_block(const Integer b) {
            assert(b > 0);
            block = b;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 0, Integer>::type get_block() const {
            assert(block > 0);
            return block;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION constexpr typename std::enable_if<(B > 0), Integer>::type get_block() const {
            return B;
        }

        template <class VIEW, class V>
        MARS_INLINE_FUNCTION void vector_apply_constraints(const Integer row, VIEW v, const V value) const {
            auto diag_row = local_to_owned_index(row);
            v(diag_row, 0) = value;
        }

        template <class VIEW>
        MARS_INLINE_FUNCTION void apply_zero_constraints(const Integer row, VIEW v) const {
            vector_apply_constraints(row, v, 0);
        }

    private:
        // manages host and device mesh objects.
        /* MM mesh_manager; */
        Mesh mesh;

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

        Integer block;
    };

}  // namespace mars

#endif
#endif

#endif
