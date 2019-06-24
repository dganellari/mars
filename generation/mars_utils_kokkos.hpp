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

template<typename T, Integer YDim_>
using ViewMatrixTexture = Kokkos::View<T*[YDim_],KokkosLayout,KokkosSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

template<typename T, Integer XDim_, Integer YDim_>
using ViewMatrixTextureC = Kokkos::View<T[XDim_][YDim_],KokkosLayout,KokkosSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

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

//copy matrix from host data to the host mirror view and then deep copy to the device texture view.
template<typename T, Integer xDim_, Integer yDim_>
void copy_matrix_from_host(std::vector<std::vector<T>> hostData,
		ViewMatrixTextureC<T, xDim_, yDim_> map_side_to_nodes, const int xDim,
		const int yDim) {

	using namespace Kokkos;

	typename ViewMatrixTextureC<T, xDim_, yDim_>::HostMirror h_view = create_mirror_view(
			map_side_to_nodes);

	parallel_for(MDRangePolicy<Rank<2>, OpenMP>( {0, 0}, {xDim, yDim}),
			KOKKOS_LAMBDA (int i, int j) {

				h_view(i,j) = hostData[i][j];

			});

	Kokkos::deep_copy(map_side_to_nodes, h_view);
}


KOKKOS_INLINE_FUNCTION
Integer index(const Integer xDim, const Integer yDim, const Integer i,
		const Integer j, const Integer k) {
	//return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
	return i + (2 * xDim + 1) * (j + k * (2 * yDim + 1));
}

void add_side(std::vector<Integer>& side, const Integer a, const Integer b,
		const Integer index) {

	if (a != 1 && b != 1) { //add only nodes which are not mid faces or mid edges
		side.push_back(index);
	} else if (a == 1 && b == 1) { // then add only mid faces
		side.push_back(index);
	}
}

/*template<Integer Dim, Integer ManifoldDim, class Point_>
void remove_extra_nodes(Mesh<Dim, ManifoldDim, Point_>& mesh,
		std::vector<Vector<Real, Dim> >& np, const std::vector<bool>& active) {

	int count = 0;
	for (unsigned int i = 0; i < active.size(); ++i) {
		if (active[i]) {
			np[count] = mesh.point(i);
			++count;
		}

	}

	mesh.setPoints(move(np));

}*/


//TODO: Repeated code. Try to optimize.
//Maybe using traits and a wrapper for the operator () in kokkos views to have a unified operator []
KOKKOS_INLINE_FUNCTION
void build_hex27(ViewMatrixTexture<Integer,hex_n_nodes> nodes, const int cube_index, const Integer xDim,
		const Integer yDim, const int i, const int j, const int k) {

	nodes(cube_index,0) = index(xDim, yDim, i, j, k);
	nodes(cube_index,1) = index(xDim, yDim, i + 2, j, k);
	nodes(cube_index,2) = index(xDim, yDim, i + 2, j + 2, k);
	nodes(cube_index,3) = index(xDim, yDim, i, j + 2, k);
	nodes(cube_index,4) = index(xDim, yDim, i, j, k + 2);
	nodes(cube_index,5) = index(xDim, yDim, i + 2, j, k + 2);
	nodes(cube_index,6) = index(xDim, yDim, i + 2, j + 2, k + 2);
	nodes(cube_index,7) = index(xDim, yDim, i, j + 2, k + 2);
	nodes(cube_index,8) = index(xDim, yDim, i + 1, j, k);
	nodes(cube_index,9) = index(xDim, yDim, i + 2, j + 1, k);
	nodes(cube_index,10) = index(xDim, yDim, i + 1, j + 2, k);
	nodes(cube_index,11) = index(xDim, yDim, i, j + 1, k);
	nodes(cube_index,12) = index(xDim, yDim, i, j, k + 1);
	nodes(cube_index,13) = index(xDim, yDim, i + 2, j, k + 1);
	nodes(cube_index,14) = index(xDim, yDim, i + 2, j + 2, k + 1);
	nodes(cube_index,15) = index(xDim, yDim, i, j + 2, k + 1);
	nodes(cube_index,16) = index(xDim, yDim, i + 1, j, k + 2);
	nodes(cube_index,17) = index(xDim, yDim, i + 2, j + 1, k + 2);
	nodes(cube_index,18) = index(xDim, yDim, i + 1, j + 2, k + 2);
	nodes(cube_index,19) = index(xDim, yDim, i, j + 1, k + 2);
	nodes(cube_index,20) = index(xDim, yDim, i + 1, j + 1, k);
	nodes(cube_index,21) = index(xDim, yDim, i + 1, j, k + 1);
	nodes(cube_index,22) = index(xDim, yDim, i + 2, j + 1, k + 1);
	nodes(cube_index,23) = index(xDim, yDim, i + 1, j + 2, k + 1);
	nodes(cube_index,24) = index(xDim, yDim, i, j + 1, k + 1);
	nodes(cube_index,25) = index(xDim, yDim, i + 1, j + 1, k + 2);
	nodes(cube_index,26) = index(xDim, yDim, i + 1, j + 1, k + 1);
}

void build_hex27(std::array<Integer, hex_n_nodes>& nodes, const Integer xDim,
		const Integer yDim, const int i, const int j, const int k) {

	nodes[0] = index(xDim, yDim, i, j, k);
	nodes[1] = index(xDim, yDim, i + 2, j, k);
	nodes[2] = index(xDim, yDim, i + 2, j + 2, k);
	nodes[3] = index(xDim, yDim, i, j + 2, k);
	nodes[4] = index(xDim, yDim, i, j, k + 2);
	nodes[5] = index(xDim, yDim, i + 2, j, k + 2);
	nodes[6] = index(xDim, yDim, i + 2, j + 2, k + 2);
	nodes[7] = index(xDim, yDim, i, j + 2, k + 2);
	nodes[8] = index(xDim, yDim, i + 1, j, k);
	nodes[9] = index(xDim, yDim, i + 2, j + 1, k);
	nodes[10] = index(xDim, yDim, i + 1, j + 2, k);
	nodes[11] = index(xDim, yDim, i, j + 1, k);
	nodes[12] = index(xDim, yDim, i, j, k + 1);
	nodes[13] = index(xDim, yDim, i + 2, j, k + 1);
	nodes[14] = index(xDim, yDim, i + 2, j + 2, k + 1);
	nodes[15] = index(xDim, yDim, i, j + 2, k + 1);
	nodes[16] = index(xDim, yDim, i + 1, j, k + 2);
	nodes[17] = index(xDim, yDim, i + 2, j + 1, k + 2);
	nodes[18] = index(xDim, yDim, i + 1, j + 2, k + 2);
	nodes[19] = index(xDim, yDim, i, j + 1, k + 2);
	nodes[20] = index(xDim, yDim, i + 1, j + 1, k);
	nodes[21] = index(xDim, yDim, i + 1, j, k + 1);
	nodes[22] = index(xDim, yDim, i + 2, j + 1, k + 1);
	nodes[23] = index(xDim, yDim, i + 1, j + 2, k + 1);
	nodes[24] = index(xDim, yDim, i, j + 1, k + 1);
	nodes[25] = index(xDim, yDim, i + 1, j + 1, k + 2);
	nodes[26] = index(xDim, yDim, i + 1, j + 1, k + 1);
}

}
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
