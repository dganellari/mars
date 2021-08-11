#ifndef GENERATION_MARS_DISTRIBUTED_MM_HPP_
#define GENERATION_MARS_DISTRIBUTED_MM_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    /* template <Integer Dim, Integer ManifoldDim, Integer Type, typename ...T> */
    template <class Mesh>
    class MeshManager {
        using simplex_type = typename Mesh::Elem;

    public:
        MARS_INLINE_FUNCTION MeshManager(Mesh *mesh) : mesh(nullptr), host_mesh(mesh) { copy_mesh_to_device(); }

        template <typename H>
        MARS_INLINE_FUNCTION void elem_iterate(H f) {
            const Integer size = host_mesh->get_chunk_size();
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        struct FaceIterate {
            FaceIterate(Mesh *m,
                        H f,
                        ViewVectorType<Integer> gl,
                        ViewVectorType<Integer> sg,
                        Integer p,
                        Integer x,
                        Integer y,
                        Integer z)
                : mesh(m), func(f), ghost_layer(gl), scan_ghost(sg), proc(p), xDim(x), yDim(y), zDim(z) {}

            template <Integer dir>
            MARS_INLINE_FUNCTION void iterate(const Integer i) const {
                // side  0 means origin side and 1 destination side.
                for (int side = 0; side < 2; ++side) {
                    Integer face_nr;

                    if (side == 0)
                        face_nr = 2 * dir + 1;
                    else
                        face_nr = 2 * dir;

                    /* Octant nbh_oc = face_nbh<simplex_type::ElemType>(ref_octant, face_nr,
                     * mesh); */
                    Octant nbh_oc = mesh->get_octant_face_nbh(i, face_nr);

                    bool ghost = false;
                    Integer index;

                    if (nbh_oc.is_valid()) {
                        Integer enc_oc = get_sfc_from_octant<simplex_type::ElemType>(nbh_oc);

                        Integer owner_proc = find_owner_processor(mesh->get_view_gp(), enc_oc, 2, proc);
                        assert(owner_proc >= 0);

                        /* if the face neighbor element is ghost then do a binary search
                         * on the ghost layer to find the index */
                        if (proc != owner_proc) {
                            ghost = true;

                            /* to narrow down the range of search we use the scan ghost
                and the owner proc of the ghost. */
                            const int start_index = scan_ghost(owner_proc);
                            const int last_index = scan_ghost(owner_proc + 1) - 1;

                            /* as opposed to the whole range: */
                            /* const int start_index = 0;
                const int last_index = ghost_layer.extent(0) -1; */

                            index = binary_search(ghost_layer.data(), start_index, last_index, enc_oc);
                            assert(index >= 0);
                        } else {
                            // using the sfc (global) to local mapping of the mesh.
                            index = mesh->get_index_of_sfc_elem(enc_oc);
                            assert(index >= 0);
                        }

                        /* printf("Index: %li, o.x: %li, y: %li, elem-index: %li, owner_proc:
                         * %li, proc: %li , o.x: %li, y: %li, index: %li, ghost: %i\n", index,
                         * ref_octant.x, ref_octant.y, elem_index(ref_octant.x, ref_octant.y,
                         * ref_octant.z, xDim, yDim), owner_proc, proc, o.x, o.y,
                         * elem_index(o.x, o.y, o.z, xDim, yDim),
                         * face.get_second_side().is_ghost()); */
                    }

                    bool boundary = nbh_oc.shares_boundary();

                    /* constructed valid for period and non-periodic. */
                    Face<simplex_type::ElemType, dir> face;

                    /* build only faces from face nr 1 and 3 (right and up) face sides
            meaning: side 0 only if the nbc_oc is not ghost to avoid a boundary face
            been called twice. Check the validate_nbh func.*/
                    if (side == 1 && mesh->is_periodic() && boundary && !ghost) {
                        nbh_oc.set_invalid();
                        face.invalidate();
                    }

                    if (face.is_valid() && ((side == 0 && nbh_oc.is_valid()) || ghost || boundary)) {
                        int origin_side = side;

                        if (boundary && !mesh->is_periodic()) origin_side = 0;

                        face.get_side(origin_side).set_elem_id(i);
                        face.get_side(origin_side).set_boundary(boundary);

                        /* if it is the side element of the ref octant. */
                        face.get_side(origin_side).set_origin();

                        if (!boundary || mesh->is_periodic()) {
                            int otherside = origin_side ^ 1;

                            face.get_side(otherside).set_elem_id(index);
                            face.get_side(otherside).set_ghost(ghost);
                        }

                        if (boundary && mesh->is_periodic()) face.swap_sides();

                        func(face);
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                /* const Integer oc = mesh->get_view_sfc()(i);
          Octant ref_octant = get_octant_from_sfc<simplex_type::ElemType>(oc); */
                iterate<0>(i);
                iterate<1>(i);
                // TODO: 3D part
            }

            Mesh *mesh;
            H func;

            ViewVectorType<Integer> ghost_layer;
            ViewVectorType<Integer> scan_ghost;

            Integer proc;
            Integer xDim;
            Integer yDim;
            Integer zDim;
        };

        template <typename H>
        void face_iterate(H f) {
            Integer xDim = host_mesh->get_XDim();
            Integer yDim = host_mesh->get_YDim();
            Integer zDim = host_mesh->get_ZDim();

            const Integer size = host_mesh->get_chunk_size();

            Kokkos::parallel_for("face_iterate",
                                 size,
                                 FaceIterate<H>(mesh,
                                                f,
                                                host_mesh->get_view_ghost(),
                                                host_mesh->get_view_scan_ghost(),
                                                host_mesh->get_proc(),
                                                xDim,
                                                yDim,
                                                zDim));
        }

        Mesh *get_mesh() const { return mesh; }

        Mesh *get_host_mesh() const { return host_mesh; }

        /* Does not perform a deep copy of the view containted in the mesh. Just the mesh object. */
        void copy_mesh_to_device() {
            Mesh *tmp = (Mesh *)Kokkos::kokkos_malloc(sizeof(Mesh));
            Mesh mCopy = *host_mesh;
            Mesh *oldDeviceMesh = mesh;
            Kokkos::parallel_for(
                "CreateDistributedMeshObject", 1, KOKKOS_LAMBDA(const int &) {
                    /* two local copies for m and tmp since this->m, this->mesh host pointers */
                    new ((Mesh *)tmp) Mesh(mCopy);
                    if (oldDeviceMesh) oldDeviceMesh->~Mesh();
                    /* it only works on a copy since **m is still a host pointer that fails on the device. */
                });

            /* make the mesh pointer a device one so that this init func is not neccessary anymore. */
            mesh = tmp;
        }

    private:
        Mesh *mesh;       // device mesh
        Mesh *host_mesh;  // host mesh that is copied to device.
    };
}  // namespace mars

#endif
#endif

#endif
