#ifndef GENERATION_MARS_DISTRIBUTED_DM_HPP_
#define GENERATION_MARS_DISTRIBUTED_DM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_user_data.hpp"

namespace mars
{

template <class Mesh, Integer degree>
class DM
{
    /* using UD = UserData<Mesh, double>; */
    using UD = UserData<Mesh>;
    using simplex_type = typename Mesh::Elem;
    static constexpr Integer dofs_per_elem = (degree + 1) ^ 2;

public:
    MARS_INLINE_FUNCTION
    DM(Mesh *mesh, const context &c) : data(UD(mesh))
    {
        create_ghost_layer<UD, simplex_type::ElemType>(c, data);
    }

    /* template <typename H>
    MARS_INLINE_FUNCTION void parallel_for_data(const Integer size, H f)
    {
        Kokkos::parallel_for("init_initial_cond", size, f);
    }

    template <typename H, typename S>
    MARS_INLINE_FUNCTION void elem_iterate_reduce(H f, S s)
    {
        const Integer size = host_mesh->get_chunk_size();
        Kokkos::parallel_reduce("elem_reduce", size, f, s);
    }

    template <typename H>
    MARS_INLINE_FUNCTION void elem_iterate(H f)
    {
        const Integer size = host_mesh->get_chunk_size();
        Kokkos::parallel_for("elem_iterate", size, f);
    }
 */

    template <Integer Type, Integer ManifoldDim>
    static MARS_INLINE_FUNCTION Integer process_corner(Integer &max_proc, Integer one_ring_owners[Type], const Mesh* mesh, const Octant& oc, const int i, const int j)
    {
        Integer xDim = mesh->get_XDim();
        Integer yDim = mesh->get_YDim();
        Integer zDim = mesh->get_ZDim();

        Octant o;
        //get the next corner using the elem sfc
        o.x = oc.x + i;
        o.y = oc.y + j;

        Octant one_ring[simplex_type::ElemType];
        o.one_ring_nbh<Type, ManifoldDim>(one_ring, xDim, yDim,
                zDim, mesh->is_periodic());

        for (int k = 0; k < Type; k++)
        {
            one_ring_owners[k] = -1;

            if (one_ring[k].is_valid())
            {
                Integer enc_oc = get_sfc_from_octant<Type>(one_ring[k]);
                /* check the proc that owns the corner and then decide on the predicate.
                    This works because the corner defines an element and this element is
                    always on the largest proc number due to the z order partitioning*/
                Integer owner_proc = find_owner_processor(mesh->get_view_gp(),
                        enc_oc, 2, mesh->get_proc());
                assert(owner_proc >= 0);

                one_ring_owners[k] = owner_proc;
                if (owner_proc > max_proc)
                {
                    max_proc = owner_proc;
                }
            }
        }


        //convert the octant value into the new nodal sfc system
        o.x *= degree;
        o.y *= degree;
        Integer sfc = get_sfc_from_octant<Type>(o);

        return sfc;
    }

    template<Integer Type, Integer dir>
    static MARS_INLINE_FUNCTION Integer process_face_corner(Octant& face_cornerA, const Mesh * mesh, const int side, const Integer i)
    {
        Octant oc = mesh->get_octant(i);

        Integer face_nr;
        if (side == 0)
            face_nr = 2 * dir + 1;
        else
            face_nr = 2 * dir;

        Octant nbh_oc = mesh->get_octant_face_nbh(i, face_nr);
        /* if (nbh_oc.is_valid())
                { */
        Integer enc_oc = get_sfc_from_octant<Type>(nbh_oc);
        Integer owner_proc =
            find_owner_processor(mesh->get_view_gp(), enc_oc, 2, mesh->get_proc());
        /* assert(owner_proc >= 0); */

        //find the starting corner "face_cornerA" of the face and use the direction
        const int val = side ^ 1;
        if (dir == 0)
        {
            face_cornerA.x = oc.x + val;
            face_cornerA.y = oc.y;
        }
        if (dir == 1)
        {
            face_cornerA.x = oc.x;
            face_cornerA.y = oc.y + val;
        }

        return owner_proc;
    }

    template<Integer Type, Integer dir>
    static MARS_INLINE_FUNCTION Integer process_face_node(const Octant& face_cornerA, const int j)
    {
        Octant o;
        if (dir == 0)
        {
            o.x = degree * face_cornerA.x;
            o.y = degree * face_cornerA.y + j + 1;
        }
        if (dir == 1)
        {
            o.x = degree * face_cornerA.x + j + 1;
            o.y = degree * face_cornerA.y;
        }

        Integer sfc = get_sfc_from_octant<Type>(o);

        return sfc;
    }

    struct IdentifyBoundaryDofPerRank
    {
        static constexpr Integer volume_nodes = (degree - 1) ^ 2;
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;

        MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc_index) const
        {
            Octant oc = mesh->get_octant(sfc_index);
            //go through all the corner dofs for the current element
            for (int i = 0; i < Mesh::ManifoldDim; i++)
            {
                for (int j = 0; j < Mesh::ManifoldDim; j++)
                /* for (int z = 0; z < Mesh::ManifoldDim; z++) // 3D case*/
                {
                    Integer max_proc = -1;
                    Integer one_ring_owners[simplex_type::ElemType];
                    Integer sfc = process_corner<simplex_type::ElemType, Mesh::ManifoldDim>(max_proc,
                            one_ring_owners, mesh, oc, i, j);

                    if (proc >= max_proc)
                    {
                        for (int k = 0; k < simplex_type::ElemType; k++)
                        {
                            //if the owner is strictly smaller than the larger should send some data to it
                            if (proc > one_ring_owners[k] && one_ring_owners[k] >= 0)
                            {
                                Integer index = sfc_to_locally_owned(sfc);
                                rank_boundary(map(one_ring_owners[k]), index) = 1;
                            }
                        }
                    }
                }
            }
        }

        template <Integer dir>
        MARS_INLINE_FUNCTION void face_iterate(const Integer i) const
        {
            //side  0 means origin side and 1 destination side.
            for (int side = 0; side < 2; ++side)
            {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA,
                        mesh, side, i);

                for (int j = 0; j < face_nodes; j++)
                {
                    Integer sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    //if the owner is strictly smaller than the larger should send some data to it
                    if (proc > owner_proc && owner_proc >= 0)
                    {
                        Integer index = sfc_to_locally_owned(sfc);
                        rank_boundary(map(owner_proc), index) = 1;
                    }
                }
            }
        }

        MARS_INLINE_FUNCTION
        void operator()(const Integer i) const
        {
            //topological order within the element
            corner_iterate(i);
            face_iterate<0>(i);
            face_iterate<1>(i);
            // TODO: 3D part
        }

        IdentifyBoundaryDofPerRank(Mesh *m, ViewMatrixType<bool> npb, ViewVectorType<Integer> mp,
                                   ViewVectorType<Integer> l, Integer p, Integer x, Integer y, Integer z)
            : mesh(m), rank_boundary(npb), map(mp), sfc_to_locally_owned(l), proc(p), xDim(x), yDim(y), zDim(z) {}

        Mesh *mesh;

        ViewMatrixType<bool> rank_boundary;
        ViewVectorType<Integer> map;
        ViewVectorType<Integer> sfc_to_locally_owned;

        Integer proc;
        Integer xDim;
        Integer yDim;
        Integer zDim;
    };

    void build_boundary_dof_sets(ViewVectorType<Integer> &scan_boundary, ViewVectorType<Integer> &boundary_sfc,
            ViewVectorType<Integer> &boundary_lsfc_index, ViewVectorType<Integer> &sender_ranks_scan)
    {
        using namespace Kokkos;

        const Integer size = data.get_host_mesh()->get_chunk_size();

        Integer xDim = data.get_host_mesh()->get_XDim();
        Integer yDim = data.get_host_mesh()->get_YDim();
        Integer zDim = data.get_host_mesh()->get_ZDim();

        //scan is always rank_scan +1 so we need rank scan size
        const Integer nbh_rank_size = sender_ranks_scan(sender_ranks_scan.extent(0) - 1);
        const Integer chunk_size_ = global_dof_enum.get_elem_size();

        printf("nbh rank size: %li - %li - %li\n", nbh_rank_size, chunk_size_, size);
        ViewMatrixType<bool> rank_boundary("count_per_proc", nbh_rank_size, chunk_size_);
        /* generate the sfc for the local and global dofs containing the generation locally
        for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
        Kokkos::parallel_for(
            "identify_boundary_predicate", size,
            IdentifyBoundaryDofPerRank(data.get_mesh(), rank_boundary, sender_ranks_scan,
                                       global_dof_enum.get_view_sfc_to_local(), data.get_host_mesh()->get_proc(),
                                       xDim, yDim, zDim));

        /* perform a scan for each row with the sum at the end for each rank */
        ViewMatrixType<Integer> rank_scan("rank_scan", nbh_rank_size, chunk_size_ + 1);
        for (int i = 0; i < nbh_rank_size; ++i)
        {
            auto row_predicate = subview(rank_boundary, i, ALL);
            auto row_scan = subview(rank_scan, i, ALL);
            incl_excl_scan(0, chunk_size_, row_predicate, row_scan);
        }

        scan_boundary = ViewVectorType<Integer>("scan_boundary_dof", nbh_rank_size + 1);
        //perform a scan on the last column to get the total sum.
        column_scan(nbh_rank_size, chunk_size_, rank_scan, scan_boundary);

        auto index_subview = subview(scan_boundary, nbh_rank_size);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);
        /* std::cout << "DM: boundary_ count result: " << h_ic() << std::endl; */

        boundary_sfc = ViewVectorType<Integer>("boundary_sfc_dofs", h_ic());
        boundary_lsfc_index = ViewVectorType<Integer>("boundary_lsfc_index_dofs", h_ic());

        /*   parallel_for(
            "print scan", rank_size + 1, KOKKOS_LAMBDA(const int i) {
                printf(" scan boundary: %i-%li\n", i, scan_boundary(i));
            }); */

        /* We use this strategy so that the compacted elements from the local_sfc
        would still be sorted and unique. */
        parallel_for(
            MDRangePolicy<Rank<2>>({0, 0}, {nbh_rank_size, chunk_size_}),
            KOKKOS_LAMBDA(const Integer i, const Integer j) {
                if (rank_boundary(i, j) == 1)
                {
                    Integer index = scan_boundary(i) + rank_scan(i, j);
                    boundary_sfc(index) = global_dof_enum.get_view_elements()(j);
                    boundary_lsfc_index(index) = j;
                }
            });

        parallel_for(
            "print set", h_ic(), KOKKOS_LAMBDA(const Integer i) {
                Integer proc = data.get_host_mesh()->get_proc();
                const Integer rank = find_owner_processor(scan_boundary, i, 1, proc);

                printf("i:%li - boundary_ : %i - %li (%li) - proc: %li - rank: %li\n", i, boundary_sfc(i), boundary_lsfc_index(i),
                       get_octant_from_sfc<simplex_type::ElemType>(boundary_sfc(i)).template get_global_index<simplex_type::ElemType>(xDim, yDim),
                       rank, proc);
            });
    }

    struct BuildLocalGlobalPredicate
    {
        static constexpr Integer volume_nodes = (degree - 1) ^ 2;
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;

        MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc_index) const
        {
            Octant o;
            Octant oc = mesh->get_octant(sfc_index);
            //go through all the inside dofs for the current element
            for (int i = 0; i < degree - 1; i++)
            {
                for (int j = 0; j < degree - 1; j++)
                {
                    o.x = degree * oc.x + i + 1;
                    o.y = degree * oc.y + j + 1;
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);

                    local_predicate(enc_oc) = 1;
                    global_predicate(enc_oc) = 1;
                }
            }
        }

        MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc_index) const
        {
            Octant oc = mesh->get_octant(sfc_index);
            //go through all the corner dofs for the current element
            for (int i = 0; i < Mesh::ManifoldDim; i++)
            {
                for (int j = 0; j < Mesh::ManifoldDim; j++)
                /* for (int z = 0; z < Mesh::ManifoldDim; z++) // 3D case*/
                {
                    Integer max_proc = -1;
                    Integer one_ring_owners[simplex_type::ElemType];
                    Integer sfc = process_corner<simplex_type::ElemType, Mesh::ManifoldDim>(max_proc,
                            one_ring_owners, mesh, oc, i, j);

                    /* if the face neighbor element is ghost then check if the processor
                        is less than the owner. This is how the dofs are partitioned */
                    /* if (!o.shares_boundary() && proc >= owner_proc) */
                    if (proc >= max_proc)
                    {
                        global_predicate(sfc) = 1;
                        for (int k = 0; k < simplex_type::ElemType; k++)
                        {
                            //if the owner is strictly smaller than the larger should send some data to it
                            if (proc > one_ring_owners[k] && one_ring_owners[k] >= 0)
                            {
                                nbh_proc_predicate_send(one_ring_owners[k]) = 1;
                            }
                        }
                    }
                    else if (proc < max_proc)
                    {
                        nbh_proc_predicate_recv(max_proc) = 1;
                    }
                    local_predicate(sfc) = 1;
                }
            }
        }

        template <Integer dir>
        MARS_INLINE_FUNCTION void face_iterate(const Integer i) const
        {
            //side  0 means origin side and 1 destination side.
            for (int side = 0; side < 2; ++side)
            {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA,
                        mesh, side, i);

               for (int j = 0; j < face_nodes; j++)
                {
                    Integer sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                   /* if the face neighbor element is ghost then check if the processor
                         is less than the owner. This is how the dofs are partitioned*/
                    if (proc >= owner_proc)
                    {
                        global_predicate(sfc) = 1;
                    }
                    /* if the owner is strictly smaller than the larger should send some
                    data to it otherwise it should receive from it */
                    if (proc > owner_proc)
                    {
                        nbh_proc_predicate_send(owner_proc) = 1;
                    }
                    else if (proc < owner_proc)
                    {
                        nbh_proc_predicate_recv(owner_proc) = 1;
                    }
                    local_predicate(sfc) = 1;
                }
            }
        }

        MARS_INLINE_FUNCTION
        void operator()(const Integer i) const
        {
            //topological order within the element
            corner_iterate(i);
            face_iterate<0>(i);
            face_iterate<1>(i);
            volume_iterate(i);
            // TODO: 3D part
        }

        BuildLocalGlobalPredicate(Mesh *m, ViewVectorType<bool> lp, ViewVectorType<bool> gp,
                                  ViewVectorType<bool> npbs, ViewVectorType<bool> npbr, Integer p, Integer x, Integer y, Integer z)
            : mesh(m), local_predicate(lp), global_predicate(gp), nbh_proc_predicate_send(npbs), nbh_proc_predicate_recv(npbr), proc(p), xDim(x), yDim(y), zDim(z) {}

        Mesh *mesh;

        ViewVectorType<bool> local_predicate;
        ViewVectorType<bool> global_predicate;
        ViewVectorType<bool> nbh_proc_predicate_send;
        ViewVectorType<bool> nbh_proc_predicate_recv;

        Integer proc;
        Integer xDim;
        Integer yDim;
        Integer zDim;
    };

    void build_lg_predicate(ViewVectorType<bool> &npbs, ViewVectorType<bool> &npbr)
    {
        Integer xDim = data.get_host_mesh()->get_XDim();
        Integer yDim = data.get_host_mesh()->get_YDim();
        Integer zDim = data.get_host_mesh()->get_ZDim();

        const Integer size = data.get_host_mesh()->get_chunk_size();

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
            "lg_predicate", size,
            BuildLocalGlobalPredicate(data.get_mesh(), local_predicate, global_predicate,
                                      npbs, npbr, data.get_host_mesh()->get_proc(),
                                      xDim, yDim, zDim));

        local_dof_enum.compact_elements(local_predicate);
        global_dof_enum.compact_elements(global_predicate);
    }

    void exchange_ghost_counts(const context &context, ViewVectorType<Integer> &scan_boundary, ViewVectorType<Integer> &sender_ranks, ViewVectorType<Integer> &recver_ranks, std::vector<Integer> &send_count, std::vector<Integer> &receive_count)

    {
        using namespace Kokkos;

        Kokkos::Timer timer;

        int proc_num = rank(context);
        int size = num_ranks(context);

        auto scan_send_mirror = create_mirror_view(scan_boundary);
        Kokkos::deep_copy(scan_send_mirror, scan_boundary);

        auto s_rank_mirror = create_mirror_view(sender_ranks);
        Kokkos::deep_copy(s_rank_mirror, sender_ranks);
        const Integer nbh_rank_size = sender_ranks.extent(0);

        /* Integer proc_count = 0; */
        for (int i = 0; i < nbh_rank_size; ++i)
        {
            Integer count = scan_send_mirror(i + 1) - scan_send_mirror(i);
            if (count > 0)
            {
                send_count[s_rank_mirror(i)] = count;
                /* ++proc_count;
                std::cout << "DM:****ToProc: " << i << " count:" << count
                          << " Proc: " << proc_num << std::endl; */
            }
        }

        auto r_rank_mirror = create_mirror_view(recver_ranks);
        Kokkos::deep_copy(r_rank_mirror, recver_ranks);
        const Integer nbh_rank_recv_size = recver_ranks.extent(0);

        /* Integer proc_count_r = 0; */
        for (int i = 0; i < nbh_rank_recv_size; ++i)
        {
            receive_count[r_rank_mirror(i)] = 1;
            /* ++proc_count_r;
            std::cout << "DM:****ToProc: " << i << " count:" << 1
                      << " Proc: " << proc_num << std::endl; */
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

    void exchange_ghost_dofs(const context &context, ViewVectorType<Integer> boundary, ViewVectorType<Integer> boundary_lsfc_index, ViewVectorType<Integer> &ghost_dofs, ViewVectorType<Integer> &ghost_dofs_index, const Integer *recv_scan, const Integer *send_scan)
    {
        int rank_size = num_ranks(context);
        Integer ghost_size = recv_scan[rank_size];

        ghost_dofs = ViewVectorType<Integer>("ghost_dofs", ghost_size);
        ghost_dofs_index = ViewVectorType<Integer>("ghost_dofs_index", ghost_size);
        //do it again to have all the process range to make it fit the i_send_recv_view

        context->distributed->i_send_recv_view(
            ghost_dofs, recv_scan, boundary, send_scan);
        std::cout << "DM:Ending mpi send receive for the ghost sfc dofs " << std::endl;

        context->distributed->i_send_recv_view(
            ghost_dofs_index, recv_scan, boundary_lsfc_index, send_scan);
        std::cout << "DM:Ending mpi send receive for the ghost dofs local index" << std::endl;
    }

    void enumerate_dofs(const context &context)
    {
        const Integer rank_size = num_ranks(context);
        const int proc_num = rank(context);

        ViewVectorType<bool> nbh_proc_predicate_send("send_to", rank_size);
        ViewVectorType<bool> nbh_proc_predicate_recv("receive_from", rank_size);
        ViewVectorType<Integer> proc_scan_send("nbh_scan_send", rank_size + 1);
        ViewVectorType<Integer> proc_scan_recv("nbh_scan_recv", rank_size + 1);

        build_lg_predicate(nbh_proc_predicate_send, nbh_proc_predicate_recv);

        incl_excl_scan(0, rank_size, nbh_proc_predicate_send, proc_scan_send);
        incl_excl_scan(0, rank_size, nbh_proc_predicate_recv, proc_scan_recv);

        ViewVectorType<Integer> sender_ranks("sender_ranks", proc_scan_send(rank_size));
        ViewVectorType<Integer> recver_ranks("receiver_ranks", proc_scan_recv(rank_size));
        compact_scan(nbh_proc_predicate_send, proc_scan_send, sender_ranks);
        compact_scan(nbh_proc_predicate_recv, proc_scan_recv, recver_ranks);

        ViewVectorType<Integer> scan_boundary, boundary_dofs, boundary_dofs_index;
        build_boundary_dof_sets(scan_boundary, boundary_dofs, boundary_dofs_index, proc_scan_send);

        std::vector<Integer> s_count(rank_size, 0);
        std::vector<Integer> r_count(rank_size, 0);

        exchange_ghost_counts(context, scan_boundary, sender_ranks, recver_ranks, s_count, r_count);

        /* create the scan recv mirror view from the receive count */
        ViewVectorType<Integer> scan_ghost("scan_ghost_", rank_size + 1);
        auto scan_recv_mirror = create_mirror_view(scan_ghost);

        make_scan_index_mirror(scan_recv_mirror, r_count);
        Kokkos::deep_copy(scan_ghost, scan_recv_mirror);

        /*create the scan recv mirror view from the receive count*/
        auto scan_send_mirror = make_scan_index(s_count);

        ViewVectorType<Integer> ghost_dofs;
        ViewVectorType<Integer> ghost_dofs_index;
        exchange_ghost_dofs(context, boundary_dofs, boundary_dofs_index, ghost_dofs, ghost_dofs_index,
                            scan_recv_mirror.data(), scan_send_mirror.data());

        using namespace Kokkos;
        Integer ghost_size = scan_recv_mirror(rank_size);
        parallel_for(
            "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer rank =
                    find_owner_processor(scan_ghost, i, 1, data.get_host_mesh()->get_proc());

                printf(" i: %li ghost: %i - %li - owner rank: %li - rank: %li\n", i,
                       ghost_dofs(i), ghost_dofs_index(i), rank, data.get_host_mesh()->get_proc());
            });
        /* build_map();
        enum_elem_dof(); */
    }

    const Integer get_dof_size() const
    {
        return dof_size;
    }

    const SFC<simplex_type::ElemType> &get_local_dof_enum() const
    {
        return local_dof_enum;
    }

    const SFC<simplex_type::ElemType> &get_global_dof_enum() const
    {
        return global_dof_enum;
    }

private:
    UD data;

    SFC<simplex_type::ElemType> local_dof_enum;
    SFC<simplex_type::ElemType> global_dof_enum;

    UnorderedMap<Integer, Integer> local_to_global_enum;
    ViewMatrixType<Integer> elem_dof_enum;

    Integer dof_size;
    /* Integer proc_count; */
};

} // namespace mars

#endif
#endif

#endif
