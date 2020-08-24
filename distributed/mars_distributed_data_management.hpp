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

    struct IdentifyBoundaryDofPerRank
    {
        static constexpr Integer volume_nodes = (degree - 1) ^ 2;
        static constexpr Integer face_nodes = (degree - 1);
        static constexpr Integer corner_nodes = 1;

        MARS_INLINE_FUNCTION void corner_iterate(const Integer sfc_index) const
        {
            Octant o;
            Octant oc = mesh->get_octant(sfc_index);
            //go through all the corner dofs for the current element
            for (int i = 0; i < Mesh::ManifoldDim; i++)
            {
                for (int j = 0; j < Mesh::ManifoldDim; j++)
                /* for (int z = 0; z < Mesh::ManifoldDim; z++) // 3D case*/
                {
                    //get the next corner using the elem sfc
                    o.x = oc.x + i;
                    o.y = oc.y + j;

                    Integer max_proc = -1;
                    Octant one_ring[simplex_type::ElemType];
                    o.one_ring_nbh<simplex_type::ElemType, Mesh::ManifoldDim>
                        (one_ring, xDim, yDim, zDim, mesh->is_periodic());

                    for (int k = 0; k < simplex_type::ElemType; k++)
                    {
                        if (one_ring[k].is_valid())
                        {
                            Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(one_ring[k]);

                            /* check the proc that owns the corner and then decide on the predicate.
                    This works because the corner defines an element and this element is
                    always on the largest proc number due to the z order partitioning*/
                            Integer owner_proc =
                                find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                            assert(owner_proc >= 0);

                            if(owner_proc > max_proc)
                                max_proc = owner_proc;
                        }
                    }

                    //convert the octant value into the new nodal sfc system
                    o.x *= degree;
                    o.y *= degree;
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);

                    //if the owner is strictly smaller than the larger should send some data to it
                    if (proc > max_proc)
                    {
                        Integer index = locally_owned(enc_oc);
                        nbh_proc_predicate(map(max_proc), index) = 1;
                    }
                }
            }
        }

        template <Integer dir>
        MARS_INLINE_FUNCTION void face_iterate(const Integer i) const
        {
            Octant oc = mesh->get_octant(i);
            //side  0 means origin side and 1 destination side.
            for (int side = 0; side < 2; ++side)
            {
                Integer face_nr;
                if (side == 0)
                    face_nr = 2 * dir + 1;
                else
                    face_nr = 2 * dir;

                Octant nbh_oc = mesh->get_octant_face_nbh(i, face_nr);
                /* if (nbh_oc.is_valid())
                { */
                Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(nbh_oc);
                Integer owner_proc =
                    find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                /* assert(owner_proc >= 0); */

                //find the starting corner of the face and use the direction
                Octant face_cornerA;
                Octant face_cornerB;
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

                Octant o;
                for (int j = 0; j < face_nodes; j++)
                {
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
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);

                   //if the owner is strictly smaller than the larger should send some data to it
                    if(proc >  owner_proc)
                    {
                        Integer index = sfc_to_locally_owned(enc_oc);
                        nbh_proc_predicate(map(max_proc), index) = 1;
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

        BuildLocalGlobalPredicate(Mesh *m, ViewMatrixType<bool> npb, ViewVectorType<Integer> m,
                                  ViewVectorTypeC<Integer> l, Integer p, Integer x, Integer y, Integer z)
            : mesh(m), nbh_proc_predicate(npb), map(m), sfc_to_locally_owned(l), proc(p), xDim(x), yDim(y), zDim(z) {}

        Mesh *mesh;

        ViewMatrixType<bool> nbh_proc_predicate;
        ViewVectorType<Integer> map;
        ViewVectorType<Integer> sfc_to_locally_owned;

        Integer proc;
        Integer xDim;
        Integer yDim;
        Integer zDim;
    };

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
            Octant o;
            Octant oc = mesh->get_octant(sfc_index);
            //go through all the corner dofs for the current element
            for (int i = 0; i < Mesh::ManifoldDim; i++)
            {
                for (int j = 0; j < Mesh::ManifoldDim; j++)
                /* for (int z = 0; z < Mesh::ManifoldDim; z++) // 3D case*/
                {
                    //get the next corner using the elem sfc
                    o.x = oc.x + i;
                    o.y = oc.y + j;

                    Integer max_proc = -1;
                    Octant one_ring[simplex_type::ElemType];
                    o.one_ring_nbh<simplex_type::ElemType, Mesh::ManifoldDim>
                        (one_ring, xDim, yDim, zDim, mesh->is_periodic());

                    for (int k = 0; k < simplex_type::ElemType; k++)
                    {
                        if (one_ring[k].is_valid())
                        {
                            Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(one_ring[k]);

                            /* check the proc that owns the corner and then decide on the predicate.
                    This works because the corner defines an element and this element is
                    always on the largest proc number due to the z order partitioning*/
                            Integer owner_proc =
                                find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                            assert(owner_proc >= 0);

                            if(owner_proc > max_proc)
                                max_proc = owner_proc;
                        }
                    }

                    //convert the octant value into the new nodal sfc system
                    o.x *= degree;
                    o.y *= degree;
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);

                    /* if the face neighbor element is ghost then check if the processor
                        is less than the owner. This is how the dofs are partitioned */
                    /* if (!o.shares_boundary() && proc >= owner_proc) */
                    if (proc >= max_proc)
                    {
                        global_predicate(enc_oc) = 1;
                    }
                    //if the owner is strictly smaller than the larger should send some data to it
                    if (proc > max_proc)
                    {
                        nbh_proc_predicate(max_proc) = 1;
                    }
                    local_predicate(enc_oc) = 1;
                }
            }
        }

        template <Integer dir>
        MARS_INLINE_FUNCTION void face_iterate(const Integer i) const
        {
            Octant oc = mesh->get_octant(i);
            //side  0 means origin side and 1 destination side.
            for (int side = 0; side < 2; ++side)
            {
                Integer face_nr;
                if (side == 0)
                    face_nr = 2 * dir + 1;
                else
                    face_nr = 2 * dir;

                Octant nbh_oc = mesh->get_octant_face_nbh(i, face_nr);
                /* if (nbh_oc.is_valid())
                { */
                Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(nbh_oc);
                Integer owner_proc =
                    find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                /* assert(owner_proc >= 0); */

                //find the starting corner of the face and use the direction
                Octant face_cornerA;
                Octant face_cornerB;
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

                Octant o;
                for (int j = 0; j < face_nodes; j++)
                {
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
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);
                    /* if the face neighbor element is ghost then check if the processor
                         is less than the owner. This is how the dofs are partitioned*/
                    if (proc >= owner_proc)
                    {
                        global_predicate(enc_oc) = 1;
                    }
                    //if the owner is strictly smaller than the larger should send some data to it
                    if(proc >  owner_proc)
                    {
                        nbh_proc_predicate(owner_proc) =1;
                    }
                    local_predicate(enc_oc) = 1;
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
                                  ViewVectorType<bool> npb, Integer p, Integer x, Integer y, Integer z)
            : mesh(m), local_predicate(lp), global_predicate(gp), nbh_proc_predicate(npb), proc(p), xDim(x), yDim(y), zDim(z) {}

        Mesh *mesh;

        ViewVectorType<bool> local_predicate;
        ViewVectorType<bool> global_predicate;
        ViewVectorType<bool> nbh_proc_predicate;

        Integer proc;
        Integer xDim;
        Integer yDim;
        Integer zDim;
    };

    void build_lg_predicate(ViewVectorType<Integer>& nbh_proc_predicate)
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
                                      nbh_proc_predicate, data.get_host_mesh()->get_proc(),
                                      xDim, yDim, zDim));

        local_dof_enum.compact_elements(local_predicate);
        global_dof_enum.compact_elements(global_predicate);
   }

    void build_boundary_dof_sets(ViewVectorType<Integer>& nbh_proc_predicate, const Integer rank_size)
    {
        const Integer nbh_rank_size = nbh_proc_scan(rank_size);

        ViewVectorType<Integer> nbh_proc_scan("nbh_scan_pred", rank_size + 1);
        incl_excl_scan(0, rank_size, nbh_proc_predicate, nbh_proc_scan);

        const Integer chunk_size_ = global_dof_enum.get_elem_size();
        ViewMatrixType<bool> rank_boundary("count_per_proc", nbh_rank_size, chunk_size_);
        /* generate the sfc for the local and global dofs containing the generation locally
        for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
        Kokkos::parallel_for(
            "identify_boundary_predicate", size,
            IdentifyBoundaryDofPerRank(data.get_mesh(), rank_boundary, nbh_proc_scan,
                                       global_dof_enum.get_view_sfc_to_local(), data.get_host_mesh()->get_proc(),
                                       xDim, yDim, zDim));


        /* perform a scan for each row with the sum at the end for each rank */
        ViewMatrixType<Integer> rank_scan("rank_scan", nbh_rank_size, chunk_size_ + 1);
        for (int i = 0; i < nbh_rank_size; ++i)
        {
            const Integer remote_proc = nbh_proc_scan(i);
            if (remote_proc != proc)
            {
                auto row_predicate = subview(rank_boundary, i, ALL);
                auto row_scan = subview(rank_scan, i, ALL);
                incl_excl_scan(0, chunk_size_, row_predicate, row_scan);
/*
                parallel_for(
                    "print scan", chunk_size_, KOKKOS_LAMBDA(const int i) {
                        printf(" boundary -inside: %i-%i", i, row_predicate(i));
                    });

                printf("\n"); */

            }
        }


        ViewVectorType<Integer> scan_boundary_("scan_boundary_", nbh_rank_size + 1);
        //perform a scan on the last column to get the total sum.
        column_scan(nbh_rank_size, chunk_size_, rank_scan, scan_boundary_);

        auto index_subview = subview(scan_boundary_, nbh_rank_size);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);
        std::cout << "boundary_ count result: " << h_ic() << std::endl;

        ViewVectorType<Integer> boundary_sfc("boundary_", h_ic());
        ViewVectorType<Integer> boundary_lsfc_index("boundary_lsfc_index_", h_ic());

        /*   parallel_for(
            "print scan", rank_size + 1, KOKKOS_LAMBDA(const int i) {
                printf(" scan boundary: %i-%li\n", i, scan_boundary_(i));
            }); */

        /* We use this strategy so that the compacted elements from the local_sfc
        would still be sorted and unique. */
        /* compact_boundary_elements(scan_boundary_, rank_boundary, rank_scan, nbh_rank_size); */
        parallel_for(
            MDRangePolicy<Rank<2>>({0, 0}, {nbh_rank_size, chunk_size_}),
            KOKKOS_LAMBDA(const Integer i, const Integer j) {
                if (rank_boundary(i, j) == 1)
                {
                    Integer index = scan_boundary_(i) + rank_scan(i, j);
                    boundary(index) = global_dof_enum.get_view_elements(j);
                    boundary_lsfc_index(index) = j;
                }
            });

        /* parallel_for(
            "print set", h_ic(), KOKKOS_LAMBDA(const Integer i) {

                const Integer rank = find_owner_processor(scan_boundary_, i, 1, proc);

                printf(" boundary_ : %i - %li (%li) - proc: %li - rank: %li\n", i, boundary_(i),
                            get_octant_from_sfc<Type>(boundary_(i)).template get_global_index<Type>(xDim, yDim),
                            rank , proc);
            }); */
    }

    void enumerate_dofs()
    {

        const Integer rank_size = mesh->get_view_gp.extent(0) / 2 - 1;
        ViewVectorType<bool> nbh_proc_predicate("npred", rank_size);


        build_lg_predicate(nbh_proc_predicate);
        build_boundary_dof_sets(nbh_proc_predicate, rank_size);
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

    const SFC<simplex_type::ElemType>& get_global_dof_enum() const
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
    Integer proc_count;
};

} // namespace mars

#endif
#endif

#endif
