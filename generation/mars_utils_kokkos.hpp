#ifndef GENERATION_MARS_UTILS_KOKKOS_HPP_
#define GENERATION_MARS_UTILS_KOKKOS_HPP_

#include "mars_mesh_kokkos.hpp"
#include "mars_mesh.hpp"

#include <Kokkos_Core.hpp>

namespace mars {
namespace generation {

using namespace kokkos;

template<Integer Dim, Integer ManifoldDim>
void convertParallelMeshToSerial(mars::Mesh<Dim, ManifoldDim>& mesh,
		mars::generation::kokkos::Parallel_Mesh<Dim, ManifoldDim>& pMesh) {

	ViewMatrixType<Integer>::HostMirror h_el = Kokkos::create_mirror_view(
			pMesh.get_view_elems());
	ViewMatrixType<Real>::HostMirror h_pt = Kokkos::create_mirror_view(
			pMesh.get_view_points());
	ViewVectorType<bool>::HostMirror h_ac = Kokkos::create_mirror_view(
			pMesh.get_view_active());

	// Deep copy device views to host views.
	Kokkos::deep_copy(h_el,pMesh.get_view_elems());
	Kokkos::deep_copy(h_pt,pMesh.get_view_points());
	Kokkos::deep_copy(h_ac,pMesh.get_view_active());

	mesh.reserve(pMesh.n_elements(), pMesh.n_nodes());

	Vector<Real, Dim> p;

	for (Integer i = 0; i < pMesh.n_nodes(); ++i) {
		for (Integer j = 0; j < Dim; ++j) {
			p[j] = h_pt(i, j);
		}

		mesh.add_point(p);
	}

	for (Integer i = 0; i < pMesh.n_elements(); ++i) {
		auto& e = mesh.add_elem();

		for (Integer j = 0; j < ManifoldDim + 1; ++j) {
			e.nodes[j] = h_el(i, j);
		}

	}

	//add_elem sets all elements to active. In case there is any non-active...
	for (Integer i = 0; i < pMesh.n_elements(); ++i) {
		if (!h_ac(i))
			mesh.set_active(i, false);
	}

}

}
}
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
