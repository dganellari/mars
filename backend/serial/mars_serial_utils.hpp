#ifndef MARS_SERIAL_UTILS_HPP
#define MARS_SERIAL_UTILS_HPP

#include "mars_vector.hpp"

namespace mars {

    template <class SerialMesh, class ParallelMesh>
    void convert_parallel_mesh_to_serial(SerialMesh& mesh, const ParallelMesh& pMesh) {
        ViewMatrixType<Integer>::HostMirror h_el = Kokkos::create_mirror_view(pMesh.get_view_elements());
        ViewMatrixType<Real>::HostMirror h_pt = Kokkos::create_mirror_view(pMesh.get_view_points());
        ViewVectorType<bool>::HostMirror h_ac = Kokkos::create_mirror_view(pMesh.get_view_active());

        // Deep copy device views to host views.
        Kokkos::deep_copy(h_el, pMesh.get_view_elements());
        Kokkos::deep_copy(h_pt, pMesh.get_view_points());
        Kokkos::deep_copy(h_ac, pMesh.get_view_active());

        mesh.reserve(pMesh.n_elements(), pMesh.n_nodes());

        Vector<Real, SerialMesh::Dim> p;

        for (Integer i = 0; i < pMesh.n_nodes(); ++i) {
            for (Integer j = 0; j < SerialMesh::Dim; ++j) {
                p[j] = h_pt(i, j);
            }

            mesh.add_point(p);
        }

        for (Integer i = 0; i < pMesh.n_elements(); ++i) {
            auto& e = mesh.add_elem();

            for (Integer j = 0; j < e.nodes.size(); ++j) {
                e.nodes[j] = h_el(i, j);
            }
        }

        // add_elem sets all elements to active. In case there is any non-active...
        for (Integer i = 0; i < pMesh.n_elements(); ++i) {
            if (!h_ac(i)) mesh.set_active(i, false);
        }
    }

    template <class SerialMesh, class ParallelMesh>
    void convert_serial_mesh_to_parallel(ParallelMesh& pMesh, const SerialMesh& mesh) {
        pMesh.reserve(mesh.n_elements(), mesh.n_nodes());

        ViewMatrixType<Integer>::HostMirror h_el = Kokkos::create_mirror_view(pMesh.get_view_elements());
        ViewMatrixType<Real>::HostMirror h_pt = Kokkos::create_mirror_view(pMesh.get_view_points());
        ViewVectorType<bool>::HostMirror h_ac = Kokkos::create_mirror_view(pMesh.get_view_active());

        for (Integer i = 0; i < mesh.n_nodes(); ++i) {
            for (Integer j = 0; j < SerialMesh::Dim; ++j) {
                h_pt(i, j) = mesh.point(i)[j];
            }
        }

        for (Integer i = 0; i < mesh.n_elements(); ++i) {
            for (Integer j = 0; j < mesh.elem(i).nodes.size(); ++j) {
                h_el(i, j) = mesh.elem(i).nodes[j];
            }

            h_ac(i) = mesh.is_active(i);
        }

        // Deep copy host views to device views.
        Kokkos::deep_copy(pMesh.get_view_elements(), h_el);
        Kokkos::deep_copy(pMesh.get_view_points(), h_pt);
        Kokkos::deep_copy(pMesh.get_view_active(), h_ac);
    }

}  // namespace mars

#endif  // MARS_SERIAL_UTILS_HPP