#ifndef GENERATION_MARS_DISTRIBUTED_DM_HPP_
#define GENERATION_MARS_DISTRIBUTED_DM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_dof.hpp"
#include "mars_distributed_user_data.hpp"

namespace mars {

    template <class Mesh, Integer degree, typename... T>
    class DM {
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
        DM(Mesh *mesh, const context &c) : data(UD(mesh)) { create_ghost_layer<UD, simplex_type::ElemType>(c, data); }

        template <Integer... dataidx>
        MARS_INLINE_FUNCTION void reserve_user_data(user_tuple &tuple, std::string view_desc, const Integer size) {
            expand_tuple<resize_view_functor, dataidx...>(resize_view_functor(view_desc, size), tuple);
        }

        MARS_INLINE_FUNCTION
        const user_tuple &get_user_data() const { return user_data; }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, user_tuple>::type>
        MARS_INLINE_FUNCTION const H get_dof_data() const {
            return std::get<idx>(user_data);
        }

        template <std::size_t idx, typename H = typename std::tuple_element<idx, tuple>::type>
        MARS_INLINE_FUNCTION H &get_dof_data(const Integer i) const {
            return std::get<idx>(user_data)(i);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void parallel_for_data(const Integer size, H f) {
            Kokkos::parallel_for("init_initial_cond", size, f);
        }
        /*
            template <typename H, typename S>
            MARS_INLINE_FUNCTION void elem_iterate_reduce(H f, S s)
            {
                const Integer size = host_mesh->get_chunk_size();
                Kokkos::parallel_reduce("elem_reduce", size, f, s);
            } */

        template <typename H>
        MARS_INLINE_FUNCTION void owned_dof_iterate(H f) {
            const Integer size = global_dof_enum.get_elem_size();
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void dof_iterate(H f) {
            const Integer size = local_dof_enum.get_elem_size();
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void elem_iterate(H f) const {
            get_data().elem_iterate(f);
        }

        MARS_INLINE_FUNCTION Integer get_elem_size() const { return get_data().get_host_mesh()->get_chunk_size(); }

        template <typename H>
        MARS_INLINE_FUNCTION void face_iterate(H f) const {
            get_data().face_iterate(f);
        }

        MARS_INLINE_FUNCTION
        Integer get_sfc_face_nbh(const Octant &oc, const Integer face_nr) const {
            Octant o = oc.sfc_face_nbh<simplex_type::ElemType>(face_nr);
            return get_sfc_from_octant<simplex_type::ElemType>(o);
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
            Integer sfc = get_sfc_from_octant<Type>(o);

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
            Integer sfc = get_sfc_from_octant<Type>(o);

            return sfc;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION Integer enum_face_node(const ViewVectorType<Integer> &sfc_to_local,
                                                           const Octant &face_cornerA,
                                                           const int j,
                                                           const int dir) {
            return sfc_to_local(sfc_face_node<part, Type>(face_cornerA, j, dir));
        }

        template <Integer Type, Integer ManifoldDim>
        static MARS_INLINE_FUNCTION Integer process_corner(Integer &max_proc,
                                                           Integer one_ring_owners[Type],
                                                           const Mesh *mesh,
                                                           const Octant &oc,
                                                           const int i,
                                                           const int j) {
            Integer xDim = mesh->get_XDim();
            Integer yDim = mesh->get_YDim();
            Integer zDim = mesh->get_ZDim();

            Octant o;
            // get the next corner using the elem sfc
            o.x = oc.x + i;
            o.y = oc.y + j;

            Octant one_ring[simplex_type::ElemType];
            o.one_ring_nbh<Type, ManifoldDim>(one_ring, xDim, yDim, zDim, mesh->is_periodic());

            for (int k = 0; k < Type; k++) {
                one_ring_owners[k] = -1;

                if (one_ring[k].is_valid()) {
                    Integer enc_oc = get_sfc_from_octant<Type>(one_ring[k]);
                    /* check the proc that owns the corner and then decide on the predicate.
                        This works because the corner defines an element and this element is
                        always on the largest proc number due to the z order partitioning*/
                    Integer owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, mesh->get_proc());
                    assert(owner_proc >= 0);

                    one_ring_owners[k] = owner_proc;
                    if (owner_proc > max_proc) {
                        max_proc = owner_proc;
                    }
                }
            }

            // convert the octant value into the new nodal sfc system
            o.x *= degree;
            o.y *= degree;
            Integer sfc = get_sfc_from_octant<Type>(o);

            return sfc;
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION Integer
        process_face_corner(Octant &face_cornerA, const Mesh *mesh, const int side, const Octant& oc) {
            Integer face_nr;
            if (side == 0)
                face_nr = 2 * dir + 1;
            else
                face_nr = 2 * dir;

            Octant nbh_oc = mesh->get_octant_face_nbh(oc, face_nr);
            /* if (nbh_oc.is_valid())
                    { */
            Integer enc_oc = get_sfc_from_octant<Type>(nbh_oc);
            Integer owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, mesh->get_proc());
            /* assert(owner_proc >= 0); */

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

            return owner_proc;
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION Integer process_face_node(const Octant &face_cornerA, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
            }

            Integer sfc = get_sfc_from_octant<Type>(o);

            return sfc;
        }


        template <bool Ghost>
        struct CornerRankBoundary {
            ViewMatrixType<bool> rank_boundary;
            ViewVectorType<Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            CornerRankBoundary(ViewMatrixType<bool> rb,
                               ViewVectorType<Integer> sl,
                               ViewVectorType<Integer> m,
                               Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            /* void corner_rank_boundary(const Mesh *mesh,
                                  const Integer i,
                                  const Integer sfc,
                                  const Integer max_proc,
                                  std::true_type) const {} */

            MARS_INLINE_FUNCTION
            void corner_rank_boundary(const Mesh *mesh,
                                  const Integer i,
                                  const Integer sfc,
                                  const Integer max_proc,
                                  std::false_type) const {
                const Integer ghost_proc =
                    find_owner_processor(mesh->get_view_scan_boundary(), i, 1, proc);

                if (proc >= max_proc)  // Send it only if it is local to the proc
                {
                    Integer index = sfc_to_locally_owned(sfc);
                    rank_boundary(index, map(ghost_proc)) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh* mesh, const Integer i, const Integer dof_sfc, const Integer max_proc,
                    Integer oro[simplex_type::ElemType]) const {
                corner_rank_boundary(mesh, i, dof_sfc, max_proc, std::integral_constant<bool, Ghost>{});
            }
        };


        template <bool Ghost>
        struct FaceRankBoundary {
            ViewMatrixType<bool> rank_boundary;
            ViewVectorType<Integer> sfc_to_locally_owned;
            ViewVectorType<Integer> map;
            Integer proc;

            MARS_INLINE_FUNCTION
            FaceRankBoundary(ViewMatrixType<bool> rb,
                               ViewVectorType<Integer> sl,
                               ViewVectorType<Integer> m,
                               Integer p)
                : rank_boundary(rb), sfc_to_locally_owned(sl), map(m), proc(p) {}

            MARS_INLINE_FUNCTION
            void face_rank_boundary(const Mesh *mesh,
                                const Integer i,
                                const Integer sfc,
                                const Integer owner_proc,
                                std::false_type) const {
                const Integer ghost_proc =
                    find_owner_processor(mesh->get_view_scan_boundary(), i, 1, proc);

                /* if (proc >= owner_proc && owner_proc >= 0) { */
                if (proc >= owner_proc) {
                    Integer index = sfc_to_locally_owned(sfc);
                    rank_boundary(index, map(ghost_proc)) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh* mesh, const Integer i, const Integer dof_sfc, const Integer owner_proc) const {
                face_rank_boundary(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
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

        struct IdentifyBoundaryDofPerRank {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = mesh->get_boundary_sfc(i);
                /* const Integer sfc = mesh->get_sfc(i); */
                const Integer proc = mesh->get_proc();

                constexpr bool BoundaryIter = false;
                corner_iterate(
                    sfc, mesh, i, CornerRankBoundary<BoundaryIter>(rank_boundary, sfc_to_locally_owned, map, proc));

                if (face_nodes > 0) {
                    FaceRankBoundary<BoundaryIter> fp =
                        FaceRankBoundary<BoundaryIter>(rank_boundary, sfc_to_locally_owned, map, proc);
                    face_iterate<0>(sfc, mesh, i, fp);
                    face_iterate<1>(sfc, mesh, i, fp);
                }

                if (volume_nodes > 0) {
                    volume_iterate(
                        sfc, mesh, i, VolumeRankBoundary<BoundaryIter>(rank_boundary, sfc_to_locally_owned, map, proc));
                }
                // TODO: 3D part
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

            /* const Integer size = data.get_host_mesh()->get_chunk_size(); */
            const Integer size = data.get_view_boundary().extent(0);

            Integer xDim = data.get_host_mesh()->get_XDim();
            Integer yDim = data.get_host_mesh()->get_YDim();
            Integer zDim = data.get_host_mesh()->get_ZDim();

            const Integer chunk_size_ = global_dof_enum.get_elem_size();
            ViewMatrixType<bool> rank_boundary("count_per_proc", chunk_size_, nbh_rank_size);
            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for(
                "identify_boundary_predicate",
                size,
                IdentifyBoundaryDofPerRank(
                    data.get_mesh(), rank_boundary, sender_ranks_scan, global_dof_enum.get_view_sfc_to_local()));

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
                    Integer proc = data.get_host_mesh()->get_proc();
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

        template <bool Ghost>
        struct CornerPredicate {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;
            Integer proc;

            MARS_INLINE_FUNCTION
            CornerPredicate(ViewVectorType<bool> lp,
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
            void corner_predicate(const Mesh *mesh,
                                  const Integer i,
                                  const Integer sfc,
                                  const Integer max_proc,
                                  Integer one_ring_owners[simplex_type::ElemType],
                                  std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh->get_view_gp(), mesh->get_ghost_sfc(i), 2, mesh->get_proc());

                if (elem_sfc_proc >= max_proc) {
                    local_predicate(sfc) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void corner_predicate(const Mesh *mesh,
                                  const Integer i,
                                  const Integer sfc,
                                  const Integer max_proc,
                                  Integer one_ring_owners[simplex_type::ElemType],
                                  std::false_type) const {
                local_predicate(sfc) = 1;


                if (proc >= max_proc) {
                    global_predicate(sfc) = 1;
                }

                //is the one ring owners needed here? Maybe using max proc is enough?
                //It will end up to be a cleaner code as well! FIXME
                for (int k = 0; k < simplex_type::ElemType; k++) {
                    /* if the owner is strictly smaller than the larger should send some data to it */
                    const Integer oro = one_ring_owners[k];
                    if (proc != oro && oro >= 0) {
                        nbh_proc_predicate_send(oro) = 1;
                        nbh_proc_predicate_recv(oro) = 1;
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh* mesh, const Integer i, const Integer dof_sfc, const Integer max_proc,
                    Integer oro[simplex_type::ElemType]) const {
                corner_predicate(mesh, i, dof_sfc, max_proc, oro, std::integral_constant<bool, Ghost>{});
            }
        };

        template <typename F>
        static MARS_INLINE_FUNCTION void
        corner_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            Octant oc = mesh->octant_from_sfc(sfc);
            // go through all the corner dofs for the current element
            for (int i = 0; i < Mesh::ManifoldDim; i++) {
                for (int j = 0; j < Mesh::ManifoldDim; j++)
                /* for (int z = 0; z < Mesh::ManifoldDim; z++) // 3D case*/
                {
                    Integer max_proc = -1;
                    Integer one_ring_owners[simplex_type::ElemType];
                    Integer dof_sfc = process_corner<simplex_type::ElemType, Mesh::ManifoldDim>(
                        max_proc, one_ring_owners, mesh, oc, i, j);

                    f(mesh, index, dof_sfc, max_proc, one_ring_owners);
                }
            }
        }

        template <bool Ghost>
        struct VolumePredicate {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;

            MARS_INLINE_FUNCTION
            VolumePredicate(ViewVectorType<bool> lp, ViewVectorType<bool> gp)
                : local_predicate(lp), global_predicate(gp) {}

            MARS_INLINE_FUNCTION
            void volume_predicate(const Integer sfc, std::true_type) const {
                local_predicate(sfc) = 1; }

            MARS_INLINE_FUNCTION
            void volume_predicate(const Integer sfc, std::false_type) const {
                local_predicate(sfc) = 1;
                global_predicate(sfc) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh* mesh, const Integer i, const Integer enc_oc) const {
                volume_predicate(enc_oc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <typename F>
        static MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            Octant o;
            Octant oc = mesh->octant_from_sfc(sfc);
            // go through all the inside dofs for the current element
            for (int i = 0; i < degree - 1; i++) {
                for (int j = 0; j < degree - 1; j++) {
                    o.x = degree * oc.x + i + 1;
                    o.y = degree * oc.y + j + 1;
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);
                    f(mesh, index, enc_oc);
                }
            }
        }

        template <bool Ghost>
        struct FacePredicate {
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;
            Integer proc;

            MARS_INLINE_FUNCTION
            FacePredicate(ViewVectorType<bool> lp,
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
            void face_predicate(const Mesh *mesh,
                                const Integer i,
                                const Integer sfc,
                                const Integer owner_proc,
                                std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh->get_view_gp(), mesh->get_ghost_sfc(i), 2, mesh->get_proc());

                //if the ghost elem owns the dof then he is able to send it.
                 if (elem_sfc_proc >= owner_proc) {
                    local_predicate(sfc) = 1;
                }
            }

            MARS_INLINE_FUNCTION
            void face_predicate(const Mesh *mesh,
                                const Integer i,
                                const Integer sfc,
                                const Integer owner_proc,
                                std::false_type) const {
                /* if the face neighbor element is ghost then check if the processor
                    is less than the owner. This is how the dofs are partitioned*/
                if (proc >= owner_proc) {
                    global_predicate(sfc) = 1;
                }
                if (proc != owner_proc) {
                    nbh_proc_predicate_send(owner_proc) = 1;
                    nbh_proc_predicate_recv(owner_proc) = 1;
                }
                local_predicate(sfc) = 1;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh *mesh,
                            const Integer i,
                            const Integer dof_sfc,
                            const Integer owner_proc) const {
                face_predicate(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        //careful: do not use this function for face iteration when expensive operations are involved. Useful mainly for predicate builder.
        //instead use the locally_owned_face_dof array to iterate over face dof sfcs.
        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const Integer sfc, const Mesh *mesh, const Integer index, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh->octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>
                    (face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_nodes; j++) {
                    Integer dof_sfc = process_face_node<simplex_type::ElemType, dir>
                        (face_cornerA, j);
                    f(mesh, index, dof_sfc, owner_proc);
                }
            }
        }

        template <bool Ghost>
        struct BuildLocalGlobalPredicate {

            BuildLocalGlobalPredicate(Mesh *m,
                                      ViewVectorType<bool> lp,
                                      ViewVectorType<bool> gp,
                                      ViewVectorType<bool> npbs,
                                      ViewVectorType<bool> npbr)
                : mesh(m),
                  local_predicate(lp),
                  global_predicate(gp),
                  nbh_proc_predicate_send(npbs),
                  nbh_proc_predicate_recv(npbr) {}

            Mesh *mesh;
            ViewVectorType<bool> local_predicate;
            ViewVectorType<bool> global_predicate;
            ViewVectorType<bool> nbh_proc_predicate_send;
            ViewVectorType<bool> nbh_proc_predicate_recv;

            MARS_INLINE_FUNCTION Integer get_ghost_sfc(const Integer index) const { return mesh->get_ghost_sfc(index); }
            MARS_INLINE_FUNCTION Integer get_local_sfc(const Integer index) const { return mesh->get_sfc(index); }

            template <bool G>
            MARS_INLINE_FUNCTION typename std::enable_if<G == true, Integer>::type get_sfc_ghost_or_local(
                const Integer index) const {
                return get_ghost_sfc(index);
            }

            template <bool G>
            MARS_INLINE_FUNCTION typename std::enable_if<G == false, Integer>::type get_sfc_ghost_or_local(
                const Integer index) const {
                return get_local_sfc(index);
            }

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const Integer sfc = get_sfc_ghost_or_local<Ghost>(i);
                const Integer proc = mesh->get_proc();
                corner_iterate(
                    sfc,
                    mesh,
                    i,
                    CornerPredicate<Ghost>(
                        local_predicate, global_predicate, nbh_proc_predicate_send, nbh_proc_predicate_recv, proc));

                if (face_nodes > 0) {
                    FacePredicate<Ghost> fp = FacePredicate<Ghost>(
                        local_predicate, global_predicate, nbh_proc_predicate_send, nbh_proc_predicate_recv, proc);
                    face_iterate<0>(sfc, mesh, i, fp);
                    face_iterate<1>(sfc, mesh, i, fp);
                }

                if (volume_nodes > 0) {
                    volume_iterate(sfc, mesh, i, VolumePredicate<Ghost>(local_predicate, global_predicate));
                }
                // TODO: 3D part
            }
        };

        template <Integer Type>
        void print_dofs(const int rank) {
            SFC<Type> dof = get_local_dof_enum();
            Kokkos::parallel_for(
                "for", dof.get_elem_size(), MARS_LAMBDA(const int i) {
                    const Integer sfc_elem = local_to_sfc(i);
                    const Integer global_dof = local_to_global(i);

                    double point[3];
                    get_vertex_coordinates_from_sfc<Type>(
                        sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

                    printf("dof: %li - gdof: %li --- (%lf, %lf) - rank: %i\n", i, global_dof, point[0], point[1], rank);
                });
        }

        void build_lg_predicate(const context &context, ViewVectorType<bool> &npbs, ViewVectorType<bool> &npbr) {
            Integer xDim = data.get_host_mesh()->get_XDim();
            Integer yDim = data.get_host_mesh()->get_YDim();
            Integer zDim = data.get_host_mesh()->get_ZDim();

            const Integer size = data.get_host_mesh()->get_chunk_size();
            const Integer ghost_size = data.get_view_ghost().extent(0);

            /* TODO: xdim and ydim should be changed to max xdim and ydim
             * of the local partition to reduce memory footprint */
            local_dof_enum = SFC<simplex_type::ElemType>(degree * xDim, degree * yDim, degree * zDim);
            global_dof_enum = SFC<simplex_type::ElemType>(degree * xDim, degree * yDim, degree * zDim);

            const Integer lall_range = local_dof_enum.get_all_range();
            const Integer gall_range = global_dof_enum.get_all_range();

            ViewVectorType<bool> local_predicate("lpred", lall_range);
            ViewVectorType<bool> global_predicate("gpred", gall_range);

            /* generate the sfc for the local and global dofs containing the generation locally
            for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
            Kokkos::parallel_for(
                "lg_predicate",
                size,
                BuildLocalGlobalPredicate<false>(
                    data.get_mesh(), local_predicate, global_predicate, npbs, npbr));

            //Iterate through ghost sfc and enumerate
            Kokkos::parallel_for(
                "lg_predicate_from_ghost",
                ghost_size,
                BuildLocalGlobalPredicate<true>(
                    data.get_mesh(), local_predicate, global_predicate, npbs, npbr));

            local_dof_enum.compact_elements(local_predicate);
            global_dof_enum.compact_elements(global_predicate);

            // reserve the user data tuple with the size as the local dof size.
            reserve_user_data(user_data, "user_data tuple", local_dof_enum.get_elem_size());

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
            std::cout << "DM:MPI send receive for the ghost dofs layer done." << std::endl;
        }

        template <typename F, Integer... dataidx>
        MARS_INLINE_FUNCTION void expand_tuple(const F &f, user_tuple &t) {
            if (sizeof...(dataidx) == 0) {
                apply_impl(f, t);
            } else {
                for_each_arg<F, 0, dataidx...>(f, t);
            }
        }

        template <typename F, Integer... dataidx>
        MARS_INLINE_FUNCTION void expand_tuple(const F &f, user_tuple &t, user_tuple &v) {
            if (sizeof...(dataidx) == 0) {
                apply_impl(f, t, v);
            } else {
                for_each_arg<F, 0, dataidx...>(f, t, v);
            }
        }

        template <typename ElementType, Integer Op = 0>
        struct FillBufferData {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                buffer_data(i) = user_data(sfc_index);
            }
        };

        template <typename ElementType>
        struct FillBufferData<ElementType, 1> {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                Kokkos::atomic_add(&user_data(sfc_index), buffer_data(i));
            }
        };

        template <typename ElementType>
        struct FillBufferData<ElementType, 2> {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                Kokkos::atomic_fetch_max(&user_data(sfc_index), buffer_data(i));
            }
        };

        template <typename ElementType>
        struct FillBufferData<ElementType, 3> {
            ElementType buffer_data;
            ElementType user_data;
            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillBufferData(ElementType bf, ElementType ud, ViewVectorType<Integer> bd, ViewVectorType<Integer> stl)
                : buffer_data(bf), user_data(ud), boundary_sfc(bd), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer sfc_index = sfc_to_local(boundary_sfc(i));
                Kokkos::atomic_fetch_min(&user_data(sfc_index), buffer_data(i));
            }
        };

        template <Integer Op>
        struct FillBufferDataFunctor {
            FillBufferDataFunctor(std::string d, size_t s, ViewVectorType<Integer> b, ViewVectorType<Integer> stl)
                : desc(d), size(s), boundary_sfc(b), sfc_to_local(stl) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                Kokkos::parallel_for(
                    desc, size, FillBufferData<ElementType, Op>(el_1, el_2, boundary_sfc, sfc_to_local));
            }

            std::string desc;
            size_t size;

            ViewVectorType<Integer> boundary_sfc;
            ViewVectorType<Integer> sfc_to_local;
        };

        template <Integer Op, Integer... dataidx>
        MARS_INLINE_FUNCTION void fill_buffer_data(user_tuple &buffer_data, const ViewVectorType<Integer> &boundary) {
            const Integer size = boundary.extent(0);
            expand_tuple<FillBufferDataFunctor<Op>, dataidx...>(
                FillBufferDataFunctor<Op>("fill_buffer_data", size, boundary, local_dof_enum.get_view_sfc_to_local()),
                buffer_data,
                user_data);
        }

        struct ExchangeGhostDofsData {
            ExchangeGhostDofsData(const context &c, Integer *sr, Integer *ss)
                : con(c), sc_rcv_mirror(sr), sc_snd_mirror(ss) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                con->distributed->i_send_recv_view(el_1, sc_rcv_mirror, el_2, sc_snd_mirror);
            }

            Integer *sc_rcv_mirror;
            Integer *sc_snd_mirror;
            const context &con;
        };

        // gather operation: fill the data from the received ghost data
        template <typename ElementType, bool Op = 0>
        struct FillUserData {
            ElementType ghost_data;
            ElementType user_data;
            ViewVectorType<Integer> ghost_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillUserData(ElementType gd, ElementType ud, ViewVectorType<Integer> gs, ViewVectorType<Integer> stl)
                : ghost_data(gd), user_data(ud), ghost_sfc(gs), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer local_sfc_index = sfc_to_local(ghost_sfc(i));
                user_data(local_sfc_index) = ghost_data(i);
            }
        };

        // scatter operation: fill the to be sent ghost buffer from data
        template <typename ElementType>
        struct FillUserData<ElementType, 1> {
            ElementType ghost_data;
            ElementType user_data;
            ViewVectorType<Integer> ghost_sfc;
            ViewVectorType<Integer> sfc_to_local;

            FillUserData(ElementType gd, ElementType ud, ViewVectorType<Integer> gs, ViewVectorType<Integer> stl)
                : ghost_data(gd), user_data(ud), ghost_sfc(gs), sfc_to_local(stl) {}

            MARS_INLINE_FUNCTION
            void operator()(Integer i) const {
                const Integer local_sfc_index = sfc_to_local(ghost_sfc(i));
                ghost_data(i) = user_data(local_sfc_index);
            }
        };

        template <bool Op>
        struct FillUserDataFunctor {
            FillUserDataFunctor(std::string d, size_t s, ViewVectorType<Integer> b, ViewVectorType<Integer> stl)
                : desc(d), size(s), ghost_sfc(b), sfc_to_local(stl) {}

            template <typename ElementType>
            void operator()(ElementType &el_1, ElementType &el_2) const {
                Kokkos::parallel_for(desc, size, FillUserData<ElementType, Op>(el_1, el_2, ghost_sfc, sfc_to_local));
            }

            std::string desc;
            size_t size;

            ViewVectorType<Integer> ghost_sfc;
            ViewVectorType<Integer> sfc_to_local;
        };

        /* fill the user data views from the ghost data that were filled in FillBufferData and received through mpi. */
        template <bool Op, Integer... dataidx>
        MARS_INLINE_FUNCTION void fill_user_data(user_tuple &ghost_user_data,
                                                 const ViewVectorType<Integer> &ghost_sfc) {
            const Integer size = ghost_sfc.extent(0);
            expand_tuple<FillUserDataFunctor<Op>, dataidx...>(
                FillUserDataFunctor<Op>("fill_user_data", size, ghost_sfc, local_dof_enum.get_view_sfc_to_local()),
                ghost_user_data,
                user_data);
        }

        template <Integer... dataidx>
        void gather_ghost_data(const context &context) {
            using namespace Kokkos;

            Kokkos::Timer timer;

            // exchange the ghost dofs first since it will be used to find the address
            // of the userdata based on the sfc code.

            int proc_num = rank(context);
            int size = num_ranks(context);

            Integer ghost_size = scan_recv_mirror(size);
            user_tuple ghost_user_data;
            reserve_user_data<dataidx...>(ghost_user_data, "ghost_user_data", ghost_size);

            // prepare the buffer to send the boundary data
            const Integer buffer_size = get_boundary_dofs().extent(0);
            user_tuple buffer_data;
            reserve_user_data<dataidx...>(buffer_data, "buffer_data", buffer_size);

            fill_buffer_data<0, dataidx...>(buffer_data, get_boundary_dofs());

            expand_tuple<ExchangeGhostDofsData, dataidx...>(
                ExchangeGhostDofsData(context, scan_recv_mirror.data(), scan_send_mirror.data()),
                ghost_user_data,
                buffer_data);

            // use the received ghost data and the sfc to put them to the unified local data
            fill_user_data<0, dataidx...>(ghost_user_data, get_ghost_dofs());

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

            Integer ghost_size = scan_recv_mirror(size);
            user_tuple ghost_buffer_data;
            reserve_user_data<dataidx...>(ghost_buffer_data, "ghost_user_data", ghost_size);

            fill_user_data<1, dataidx...>(ghost_buffer_data, get_ghost_dofs());

            const Integer boundary_size = get_boundary_dofs().extent(0);
            user_tuple boundary_user_data;
            reserve_user_data<dataidx...>(boundary_user_data, "boundary_user_data", boundary_size);

            // prepare the buffer to send the boundary data
            expand_tuple<ExchangeGhostDofsData, dataidx...>(
                ExchangeGhostDofsData(context, scan_send_mirror.data(), scan_recv_mirror.data()),
                boundary_user_data,
                ghost_buffer_data);
            /* print_nth_tuple<1>(proc_num); */

            return boundary_user_data;
        }

        template <Integer... dataidx>
        void scatter_add(user_tuple &boundary_user_data) {
            fill_buffer_data<1, dataidx...>(boundary_user_data, get_boundary_dofs());
        }

        template <Integer... dataidx>
        void scatter_max(user_tuple &boundary_user_data) {
            fill_buffer_data<2, dataidx...>(boundary_user_data, get_boundary_dofs());
        }

        template <Integer... dataidx>
        void scatter_min(user_tuple &boundary_user_data) {
            fill_buffer_data<3, dataidx...>(boundary_user_data, get_boundary_dofs());
        }

        struct EnumLocalDofs {
            MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc_index, Integer &index) const {
                Octant o;
                Octant oc = mesh->get_octant(sfc_index);
                // go through all the inside dofs for the current element
                // TODO: maybe a topolgical order would be neccessary
                for (int j = 0; j < degree - 1; j++) {
                    for (int i = 0; i < degree - 1; i++) {
                        o.x = degree * oc.x + i + 1;
                        o.y = degree * oc.y + j + 1;
                        Integer sfc = get_sfc_from_octant<simplex_type::ElemType>(o);
                        Integer localid = sfc_to_local(sfc);
                        elem_dof_enum(sfc_index, index++) = localid;
                    }
                }
            }

            MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc_index, Integer &index) const {
                Octant oc = mesh->get_octant(sfc_index);

                // go through all the corner dofs for the current element counterclockwise
                /* for (const auto &x : {{0,0}, {1, 0}, {1, 1}, {0, 1}}) */
                /* for (i,j from 0 to 2) first = (i + j) % 2; second = i;*/

                Integer localid = enum_corner<simplex_type::ElemType>(sfc_to_local, oc, 0, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = enum_corner<simplex_type::ElemType>(sfc_to_local, oc, 1, 0);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = enum_corner<simplex_type::ElemType>(sfc_to_local, oc, 1, 1);
                elem_dof_enum(sfc_index, index++) = localid;

                localid = enum_corner<simplex_type::ElemType>(sfc_to_local, oc, 0, 1);
                elem_dof_enum(sfc_index, index++) = localid;
            }

            template <Integer part>
            MARS_INLINE_FUNCTION void face_iterate(const Integer sfc_index, Integer &index) const {
                Octant oc = mesh->get_octant(sfc_index);
                // side  0 means origin side and 1 destination side.
                for (int dir = 0; dir < 2; ++dir) {
                    Octant face_cornerA = enum_face_corner<part>(oc, dir);

                    for (int j = 0; j < face_nodes; j++) {
                        Integer localid =
                            enum_face_node<part, simplex_type::ElemType>(sfc_to_local, face_cornerA, j, dir);
                        elem_dof_enum(sfc_index, index++) = localid;
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                Integer index = 0;
                // topological order within the element
                corner_iterate(i, index);
                face_iterate<0>(i, index);
                face_iterate<1>(i, index);
                volume_iterate(i, index);
                // TODO: 3D part
            }

            EnumLocalDofs(Mesh *m, ViewMatrixType<Integer> ede, ViewVectorType<Integer> stl)
                : mesh(m), elem_dof_enum(ede), sfc_to_local(stl) {}

            Mesh *mesh;
            ViewMatrixType<Integer> elem_dof_enum;
            ViewVectorType<Integer> sfc_to_local;
        };

        void enumerate_local_dofs() {
            const Integer size = data.get_host_mesh()->get_chunk_size();
            elem_dof_enum = ViewMatrixType<Integer>("elem_dof_enum", size, elem_nodes);

            /* enumerates the dofs within each element topologically */
            Kokkos::parallel_for("enum_local_dofs",
                                 size,
                                 EnumLocalDofs(data.get_mesh(), elem_dof_enum, local_dof_enum.get_view_sfc_to_local()));
        }

        virtual void enumerate_dofs(const context &context) {
            const Integer rank_size = num_ranks(context);
            const int proc_num = rank(context);

            ViewVectorType<bool> nbh_proc_predicate_send("send_to", rank_size);
            ViewVectorType<bool> nbh_proc_predicate_recv("receive_from", rank_size);
            ViewVectorType<Integer> proc_scan_send("nbh_scan_send", rank_size + 1);
            ViewVectorType<Integer> proc_scan_recv("nbh_scan_recv", rank_size + 1);

            build_lg_predicate(context, nbh_proc_predicate_send, nbh_proc_predicate_recv);

            /* print_dofs<simplex_type::ElemType>(proc_num); */

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

            enumerate_local_dofs();
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
                                 BuildLocalToGlobalGhostMap(data.get_mesh(),
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
            // use the sfc to local which the scan of the predicate. To get the predicate
            // value the difference with the successive index is needed.
            const Integer pred_value =
                local_dof_enum.get_view_sfc_to_local()(sfc + 1) - local_dof_enum.get_view_sfc_to_local()(sfc);
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
        Dof sfc_to_global_dof(const Integer sfc) const {
            Dof dof;
            if (locally_owned_dof(sfc)) {
                const Integer proc = data.get_mesh()->get_proc();
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
        Integer local_to_owned(const Integer local) const {
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
        Integer local_to_sfc(const Integer local) const { return local_dof_enum.get_view_elements()(local); }

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

        template <Integer Type, Integer FaceNr = -1>
        MARS_INLINE_FUNCTION bool is_boundary(const Integer local) const {
            const Integer sfc = local_to_sfc(local);
            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            return is_boundary_sfc<Type, FaceNr>(sfc, xdim, ydim, zdim);
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void get_locally_owned_data(const ViewVectorType<H> &x) {
            using namespace Kokkos;

            assert(get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);
                    x(i) = dof_data(local);
                });
        }

        template <Integer idx, typename H = typename std::tuple_element<idx, tuple>::type>
        void set_locally_owned_data(const ViewVectorType<H> &x) {
            using namespace Kokkos;

            assert(get_global_dof_enum().get_elem_size() == x.extent(0));
            const Integer size = get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);
                    dof_data(local) = x(i);
                });
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
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    if (is_boundary_sfc<Type, face_nr>(sfc, xdim, ydim, zdim)) {
                        f(i, sfc);
                    }
                });
        }

        template <Integer idx,
                  Integer face_nr = -1,
                  typename F,
                  typename H = typename std::tuple_element<idx, tuple>::type>
        void boundary_dof_iterate(F f) {
            using namespace Kokkos;
            constexpr Integer Type = simplex_type::ElemType;

            const Integer size = get_global_dof_enum().get_elem_size();

            ViewVectorType<Integer> global_to_sfc = get_global_dof_enum().get_view_elements();
            ViewVectorType<Integer> sfc_to_local = get_local_dof_enum().get_view_sfc_to_local();
            ViewVectorType<H> dof_data = get_dof_data<idx>();

            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            Kokkos::parallel_for(
                "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                    const Integer sfc = global_to_sfc(i);
                    const Integer local = sfc_to_local(sfc);

                    if (is_boundary_sfc<Type, face_nr>(sfc, xdim, ydim, zdim)) {
                        f(local, dof_data(local));
                    }
                });
        }

        /* :TODO stencil use the face iterate on the dof sfc. */
        /* The face nbh will give an sfc code. To get the local code from it
         * sfc_to_local can be used */

        // Two way of iterations: face and element. You can also build your stencil yourself.

        MARS_INLINE_FUNCTION
        const Integer get_dof_size() const { return local_dof_enum.get_elem_size(); }

        MARS_INLINE_FUNCTION
        const Integer get_locally_owned_dof_size() const { return global_dof_enum.get_elem_size(); }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_local_dof_enum() const { return local_dof_enum; }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType> &get_global_dof_enum() const { return global_dof_enum; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_global_dof_offset() const { return global_dof_offset; }

        MARS_INLINE_FUNCTION
        const UnorderedMap<Integer, Dof> &get_ghost_lg_map() const { return ghost_local_to_global_map; }

        MARS_INLINE_FUNCTION
        const ViewMatrixType<Integer> get_elem_dof_enum() const { return elem_dof_enum; }

        MARS_INLINE_FUNCTION
        const Integer get_elem_local_dof(const Integer elem_index, const Integer i) const {
            return elem_dof_enum(elem_index, i);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_boundary_dofs() const { return boundary_dofs_sfc; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_ghost_dofs() const { return ghost_dofs_sfc; }

        UD get_data() const { return data; }

        const Integer get_proc() { return data.get_host_mesh->get_proc(); }

    private:
        // data associated to the mesh elements (sfc).
        UD data;

        // local dof enumeration containing the local to sfc and sfc to local views.
        SFC<simplex_type::ElemType> local_dof_enum;
        // locally owned dof enumeration containing the global to sfc and sfc to global views.
        SFC<simplex_type::ElemType> global_dof_enum;
        // global offset used to calc the global numbering of the dofs.
        ViewVectorType<Integer> global_dof_offset;

        // the local to global mapp for the locally owned is mapped from the sfc_to global view
        // for the ghost dofs is used the following kokkos unordered map.
        UnorderedMap<Integer, Dof> ghost_local_to_global_map;

        // local enumeration of the dofs topologically foreach element
        ViewMatrixType<Integer> elem_dof_enum;

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

        // data associated to the dof data.
        user_tuple user_data;
    };

}  // namespace mars

#endif
#endif

#endif
