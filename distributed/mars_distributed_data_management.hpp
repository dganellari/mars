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
    using simplex_type = typename Mesh::Elem;
    static constexpr Integer dofs_per_elem = (degree + 1) ^ 2;

public:
    MARS_INLINE_FUNCTION
    DM(Mesh *mesh, const context &c)
    {
        context = c;
        data = UserData<Mesh>(mesh);
        create_ghost_layer(context, data);
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

    struct BuildLocalGlobalPredicate
    {
        constexpr Integer volume_nodes = (degree - 1) ^ 2;
        constexpr Integer face_nodes = (degree - 1);
        constexpr Integer corner_nodes = 1;

        MARS_INLINE_FUNCTION void volume_iterate(const Integer sfc_index) const
        {
            Octant o;
            Octant oc = get_octant(sfc_index);
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
            Octant oc = get_octant(sfc_index);
            //go through all the corner dofs for the current element
            for (int i = 0; i < Mesh::ManifoldDim; i++)
            {
                for (int j = 0; j < Mesh::ManifoldDim; j++)
                /* for (int z = 0; z < Mesh::ManifoldDim; z++) // 3D case*/
                {
                    //get the next corner using the elem sfc
                    o.x = oc.x + i;
                    o.y = oc.y + j;
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);

                    /* check the proc that owns the corner and then decide on the predicate.
                    This works because the corner defines an element and this element is
                    always on the largest proc number due to the z order partitioning*/
                    Integer owner_proc =
                        find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                    assert(owner_proc >= 0);
                    //convert the octant value into the new nodal sfc system
                    o.x *= degree;
                    o.y *= degree;
                    enc_oc = get_sfc_from_octant<simplex_type::ElemType>(o);

                    /* if the face neighbor element is ghost then check if the processor
                      is less than the owner. This is how the dofs are partitioned*/
                    if (proc >= owner_proc)
                    {
                        global_predicate(enc_oc) = 1;
                    }
                    local_predicate(enc_oc) = 1;
                }
            }
        }

        template <Integer dir>
        MARS_INLINE_FUNCTION void face_iterate(const Integer i) const
        {
            Octant oc = get_octant(i);
            //side  0 means origin side and 1 destination side.
            for (int side = 0; side < 2; ++side)
            {
                Integer face_nr;
                if (side == 0)
                    face_nr = 2 * dir + 1;
                else
                    face_nr = 2 * dir;

                Octant nbh_oc = mesh->get_octant_face_nbh(i, face_nr);
                if (nbh_oc.is_valid())
                {
                    Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(nbh_oc);
                    Integer owner_proc =
                        find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                    assert(owner_proc >= 0);

                    //find the starting corner of the face and use the direction
                    Octant face_cornerA;
                    /* Octant face_cornerB; */
                    const int val = origin_side ^ 1;
                    if (dir == 0)
                    {
                        face_cornerA.x = oc.x + val;
                        face_cornerA.y = oc.y;
                        /* face_cornerB.x = oc.x + val;
                        face_cornerB.y = ocy + 1; */
                    }
                    if (dir == 1)
                    {
                        face_cornerA.x = oc.x;
                        face_cornerA.y = oc.y + val
                        /* face_cornerB.x = oc.x + 1;
                        face_cornerB.y = oc.y + val; */
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
                        local_predicate(enc_oc) = 1;
                    }
                    /* printf("Index: %li, o.x: %li, y: %li, elem-index: %li, owner_proc:
           * %li, proc: %li , o.x: %li, y: %li, index: %li, ghost: %i\n", index,
           * ref_octant.x, ref_octant.y, elem_index(ref_octant.x, ref_octant.y,
           * ref_octant.z, xDim, yDim), owner_proc, proc, o.x, o.y,
           * elem_index(o.x, o.y, o.z, xDim, yDim),
           * face.get_second_side().is_ghost()); */
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

        BuildLocalGlobalPredicate(Mesh *m,  ViewVectorType<bool> lp, ViewVectorType<bool> gp,
                Integer p, Integer x, Integer y, Integer z)
            : mesh(m), local_predicate(lp), global_predicate(gp), proc(p), xDim(x), yDim(y), zDim(z) {}

        Mesh *mesh;

        ViewVectorType<bool> local_predicate;
        ViewVectorType<bool> global_predicate;

        Integer proc;
        Integer xDim;
        Integer yDim;
        Integer zDim;
    };

    void build_lg_predicate()
    {
        Integer xDim = host_mesh->get_XDim();
        Integer yDim = host_mesh->get_YDim();
        Integer zDim = host_mesh->get_ZDim();

        const Integer size = host_mesh->get_chunk_size();

        /* TODO: xdim and ydim should be changed to max xdim and ydim
         * of the local partition to reduce memory footprint */
        local_dof_enum = SFC<simplex_type::ElemType>(xDim, yDim, zDim);
        global_dof_enum = SFC<simplex_type::ElemType>(xDim, yDim, zDim);

        ViewVectorType<bool> local_predicate("lpred", local_dof_enum.get_all_range());
        ViewVectorType<bool> global_predicate("gpred", global_dof_enum.get_all_range());

        /* generate the sfc for the local and global dofs containing the generation locally
        for each partition of the mesh using the existing elem sfc to build this nodal sfc. */
        Kokkos::parallel_for(
            "lg_predicate", size,
            BuildLocalGlobalPredicate(mesh, local_predicate, global_predicate,
                                      host_mesh->get_proc(), xDim, yDim, zDim));

        local_dof_enum.compact_elements(local_predicate);
        global_dof_enum.compact_elements(global_predicate);
    }

    const Integer get_dof_size() const
    {
        return dof_size;
    }

private:
    context &context;

    UserData<Mesh> data;

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
