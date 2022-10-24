#ifndef GENERATION_MARS_DISTRIBUTED_MM_HPP_
#define GENERATION_MARS_DISTRIBUTED_MM_HPP_

#ifdef MARS_ENABLE_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    /* template <Integer Dim, Integer ManifoldDim, Integer Type, typename ...T> */
    template <class Mesh>
    class MeshManager {

    public:
        MARS_INLINE_FUNCTION MeshManager(Mesh *mesh) : host_mesh(mesh), mesh(nullptr) { copy_mesh_to_device(); }

        MARS_INLINE_FUNCTION
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
