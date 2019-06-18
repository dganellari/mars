#ifndef GENERATION_MARS_UTILS_KOKKOS_HPP_
#define GENERATION_MARS_UTILS_KOKKOS_HPP_

#include "mars_mesh_kokkos.hpp"
#include "mars_mesh.hpp"

#include <Kokkos_Core.hpp>

namespace mars {

constexpr int hex_n_sides = 6; // 6 faces in total for the hex27.
constexpr int hex_n_nodes = 27; // 27 nodes for the hex27.
constexpr int hex_side_n_nodes = 9; // 9 nodes per face for the hex27.

//libmesh method to map the sides to nodes.
const std::vector<std::vector<unsigned int>> hex_side_nodes{ { 0, 3, 2,
		1, 11, 10, 9, 8, 20 }, // Side 0
		{ 0, 1, 5, 4, 8, 13, 16, 12, 21 }, // Side 1
		{ 1, 2, 6, 5, 9, 14, 17, 13, 22 }, // Side 2
		{ 2, 3, 7, 6, 10, 15, 18, 14, 23 }, // Side 3
		{ 3, 0, 4, 7, 11, 12, 19, 15, 24 }, // Side 4
		{ 4, 5, 6, 7, 16, 17, 18, 19, 25 }  // Side 5
};

//return size of an array as a compile-time constant.
template<typename T, std::size_t N, std::size_t M>
constexpr std::size_t arraySize(T (&)[N][M]) noexcept
{
	return N;
}

#ifdef MARS_USE_CUDA
#define KokkosSpace Kokkos::CudaSpace
#define KokkosLayout Kokkos::LayoutLeft
#else
#define KokkosSpace Kokkos::OpenMP
#define KokkosLayout Kokkos::LayoutRight
#endif


template<typename T>
using ViewVectorType = Kokkos::View<T*,KokkosLayout,KokkosSpace>;

template<typename T>
using ViewMatrixType = Kokkos::View<T**,KokkosLayout,KokkosSpace>;

template<typename T>
using ViewMatrixTexture = Kokkos::View<T**,KokkosLayout,KokkosSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

template<typename T>
using ViewVectorTexture = Kokkos::View<T*,KokkosLayout,KokkosSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

class KokkosImplementation {};

template<Integer Dim, Integer ManifoldDim>
void convert_parallel_mesh_to_serial(mars::Mesh<Dim, ManifoldDim>& mesh,
		mars::Mesh<Dim, ManifoldDim,KokkosImplementation>& pMesh) {

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

template<typename T>
void copy_matrix_from_host(std::vector<std::vector<T>> hostData,
		ViewMatrixTexture<T> map_side_to_nodes, const int xDim,
		const int yDim) {

	using namespace Kokkos;

	typename ViewMatrixTexture<T>::HostMirror h_view = create_mirror_view(
			map_side_to_nodes);

	parallel_for(MDRangePolicy<Rank<2>, OpenMP>( {0, 0}, {xDim, yDim}),
			KOKKOS_LAMBDA (int i, int j) {

				h_view(i,j) = hostData[i][j];

			});

	Kokkos::deep_copy(map_side_to_nodes, h_view);
}

}
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
