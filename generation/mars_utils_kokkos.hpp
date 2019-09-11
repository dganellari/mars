#ifndef GENERATION_MARS_UTILS_KOKKOS_HPP_
#define GENERATION_MARS_UTILS_KOKKOS_HPP_

#include "mars_mesh.hpp"
#include "mars_globals.hpp"

namespace mars {


//return size of an array as a compile-time constant.
template<typename T, std::size_t N, std::size_t M>
constexpr std::size_t arraySize(T (&)[N][M]) noexcept
{
	return N;
}

#ifdef KOKKOS_ENABLE_OPENMP // for the case when cuda is used but for the host copy the openmp could still be used.
	#define KokkosHostSpace Kokkos::OpenMP
#else
	// #ifdef KOKKOS_ENABLE_THREADS
	// 	#define KokkosHostSpace Kokkos::Threads
	// #else //KOKKOS_ENABLE_THREADS
	#define KokkosHostSpace Kokkos::Serial
	// #endif //KOKKOS_ENABLE_THREADS
#endif

#ifdef MARS_USE_CUDA
	#define KokkosSpace Kokkos::CudaSpace
	#define KokkosLayout Kokkos::LayoutLeft
#else //MARS_USE_CUDA
	#ifdef KOKKOS_ENABLE_OPENMP
		#define KokkosSpace Kokkos::OpenMP
		#define KokkosLayout Kokkos::LayoutRight
	#else //KOKKOS_ENABLE_OPENMP
		// #ifdef KOKKOS_ENABLE_THREADS
		// 	#define KokkosSpace Kokkos::Threads
		// 	#define KokkosLayout Kokkos::LayoutRight
		// #else //KOKKOS_ENABLE_THREADS
		#define KokkosSpace Kokkos::Serial
		#define KokkosLayout Kokkos::LayoutRight
		// #endif //KOKKOS_ENABLE_THREADS
	#endif //KOKKOS_ENABLE_OPENMP
#endif //MARS_USE_CUDA


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

template<typename T, Integer XDim_>
using ViewVectorTypeC = Kokkos::View<T[XDim_],KokkosLayout,KokkosSpace>;

template<typename T, Integer YDim_>
struct IndexView {
	ViewMatrixTexture<T, YDim_> view;
	int index;

	IndexView(ViewMatrixTexture<T, YDim_> v, int idx) :
			view(v), index(idx) {
	}

	KOKKOS_INLINE_FUNCTION
	T& operator[](int i) {
		return view(index, i);
	}

};

class KokkosImplementation {
	std::string name = "kokkos";
};

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

	parallel_for(MDRangePolicy<Rank<2>, KokkosHostSpace>( {0, 0}, {xDim, yDim}),
			KOKKOS_LAMBDA (int i, int j) {

				h_view(i,j) = hostData[i][j];

			});

	Kokkos::deep_copy(map_side_to_nodes, h_view);
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



}
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
