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

template<typename T>
using ViewVectorTypeU = Kokkos::View<T*, KokkosSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template<typename T>
using ViewObject = Kokkos::View<T[1], KokkosSpace>;

template<typename T, class space>
using ViewObjectU = Kokkos::View<T[1], space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
/*

template<typename T>
using ViewObjectUH = Kokkos::View<T[1], KokkosHostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
*/

template<typename Key, typename Value>
using UnorderedMap = Kokkos::UnorderedMap<Key,Value,KokkosSpace>;

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
		const mars::Mesh<Dim, ManifoldDim,KokkosImplementation>& pMesh) {

	ViewMatrixType<Integer>::HostMirror h_el = Kokkos::create_mirror_view(
			pMesh.get_view_elements());
	ViewMatrixType<Real>::HostMirror h_pt = Kokkos::create_mirror_view(
			pMesh.get_view_points());
	ViewVectorType<bool>::HostMirror h_ac = Kokkos::create_mirror_view(
			pMesh.get_view_active());

	// Deep copy device views to host views.
	Kokkos::deep_copy(h_el,pMesh.get_view_elements());
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

template<Integer Dim, Integer ManifoldDim>
void convert_serial_mesh_to_parallel(
		mars::Mesh<Dim, ManifoldDim, KokkosImplementation>& pMesh,
		const mars::Mesh<Dim, ManifoldDim>& mesh)
{

	pMesh.reserve(mesh.n_elements(), mesh.n_nodes());

	ViewMatrixType<Integer>::HostMirror h_el = Kokkos::create_mirror_view(
			pMesh.get_view_elements());
	ViewMatrixType<Real>::HostMirror h_pt = Kokkos::create_mirror_view(
			pMesh.get_view_points());
	ViewVectorType<bool>::HostMirror h_ac = Kokkos::create_mirror_view(
			pMesh.get_view_active());

	for (Integer i = 0; i < mesh.n_nodes(); ++i)
	{
		for (Integer j = 0; j < Dim; ++j)
		{
			h_pt(i, j) = mesh.point(i)[j];
		}
	}

	for (Integer i = 0; i < mesh.n_elements(); ++i)
	{

		for (Integer j = 0; j < ManifoldDim + 1; ++j)
		{
			h_el(i, j) = mesh.elem(i).nodes[j];
		}

		h_ac(i) = mesh.is_active(i);
	}

	// Deep copy host views to device views.
	Kokkos::deep_copy(pMesh.get_view_elements(), h_el);
	Kokkos::deep_copy(pMesh.get_view_points(), h_pt);
	Kokkos::deep_copy(pMesh.get_view_active(), h_ac);
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

template<class Mesh_>
struct Hypersphere
{
	static constexpr Integer Dim = Mesh_::Dim;
	static constexpr Integer ManifoldDim = Mesh_::ManifoldDim;

	Mesh_* mesh;
	ViewVectorTypeC<Integer, Dim> center;
	Real radius;

	ViewVectorType<Integer> count;

	Hypersphere(Mesh_* ms, ViewVectorTypeC<Integer, Dim> c, Real rd,
			ViewVectorType<Integer> pc) :
			mesh(ms), center(c), radius(rd), count(pc)
	{
	}

	Hypersphere(Mesh_* ms, ViewVectorTypeC<Integer, Dim> c, Real rd) :
				mesh(ms), center(c), radius(rd)
	{
	}


	MARS_INLINE_FUNCTION
	void operator()(const int element_id) const
	{
		Integer hyper_count=0;

		bool inside = false;
		bool outside = false;

		if (mesh->is_active(element_id)) //TODO: remove if by compacting on active elems.
		{
			auto e = mesh->elem(element_id);

			for(Integer i=0; i<ManifoldDim+1; ++i){

				auto dir = mesh->point(e.nodes[i]);
				auto d = dir.norm();

				if(d < radius) {
					inside = true;
				} else if(d > radius) {
					outside = true;
				} else if(std::abs(d) < 1e-16) {
					inside = true;
					outside = true;
					break;
				}
			}

			if(inside && outside) {
				++hyper_count;
			}
		}

		count(element_id+1) = hyper_count;
		//+1 for leaving the first cell 0 and performing an inclusive scan on the rest
		//to have both exclusive and inclusive and the total on the last cell.
	}
};


template<class Mesh_>
struct ScatterHypersphere : Hypersphere<Mesh_>
{
	static constexpr Integer Dim = Mesh_::Dim;
	static constexpr Integer ManifoldDim = Mesh_::ManifoldDim;

	ViewVectorType<Integer> elements;
	ViewVectorType<Integer> scan;

	ScatterHypersphere(Mesh_* ms, ViewVectorTypeC<Integer, Dim> c, Real rd,
			ViewVectorType<Integer> sc, ViewVectorType<Integer> elems) :
			Hypersphere<Mesh_>(ms, c, rd), scan(sc), elements(elems)
	{
	}

	MARS_INLINE_FUNCTION
	void operator()(const int element_id) const
	{
		bool inside = false;
		bool outside = false;

		if (mesh->is_active(element_id)) //TODO: remove if by compacting on active elems.
		{
			auto e = mesh->elem(element_id);

			for(Integer i=0; i<ManifoldDim+1; ++i){

				auto dir = mesh->point(e.nodes[i]);
				auto d = dir.norm();

				if(d < radius) {
					inside = true;
				} else if(d > radius) {
					outside = true;
				} else if(std::abs(d) < 1e-16) {
					inside = true;
					outside = true;
					break;
				}
			}

			if(inside && outside) {
				elements(scan(element_id)) = element_id;
			}
		}
	}
};


void inclusive_scan(const Integer start, const Integer end,
			ViewVectorType<Integer> in_)
{
	using namespace Kokkos;

	parallel_scan (RangePolicy<>(start , end ),	KOKKOS_LAMBDA (const int& i,
				Integer& upd, const bool& final)
	{
		// Load old value in case we update it before accumulating
		const float val_i = in_(i);

		upd += val_i;

		if (final)
		{
			in_(i) = upd; // only update array on final pass
		}
	});
}

void inclusive_scan(const Integer start, const Integer end,
			const ViewVectorType<Integer> in_, ViewVectorType<Integer> out_))
{
	using namespace Kokkos;

	parallel_scan (RangePolicy<>(start , end ),	KOKKOS_LAMBDA (const int& i,
				Integer& upd, const bool& final)
	{
		// Load old value in case we update it before accumulating
		const float val_i = in_(i);

		upd += val_i;

		if (final)
		{
			out_(i) = upd; // only update array on final pass
		}
	});
}

void incl_excl_scan(const Integer start, const Integer end,
			const ViewVectorType<Integer> in_, ViewVectorType<Integer> out_))
{
	using namespace Kokkos;

	parallel_scan (RangePolicy<>(start , end ),	KOKKOS_LAMBDA (const int& i,
				Integer& upd, const bool& final)
	{
		// Load old value in case we update it before accumulating
		const float val_i = in_(i);

		upd += val_i;

		if (final)
		{
			out_(i+1) = upd; // To have both ex and inclusive in the same output.
		}
	});
}

void exclusive_scan(const Integer start, const Integer end,
			ViewVectorType<Integer> in_)
{
	using namespace Kokkos;
	parallel_scan (RangePolicy<>(start , end ),	KOKKOS_LAMBDA (const int& i,
				Integer& upd, const bool& final)
	{
		// Load old value in case we update it before accumulating
		const float val_i = in_(i);

		if (final)
		{
			in_(i) = upd; // only update array on final pass
		}

		upd += val_i;

	});
}

template<class Mesh_>
ViewVectorType<Integer> mark_hypersphere_for_refinement(const Mesh_* mesh,
ViewVectorTypeC<Integer, Mesh_::Dim> center, Real radius, const Integer nr_elements)
{
	using namespace Kokkos;

	Timer timer1;

	ViewVectorType<Integer> count_ = ViewVectorType<Integer>("count_hypers",
		nr_elements + 1);

	parallel_for(nr_elements,
			Hypersphere(mesh, center, radius, count_));

	double time1 = timer1.seconds();
	std::cout << "Hypersphere Count took: " << time1 << " seconds." << std::endl;


	Timer timer2;

	inclusive_scan(1, nr_elements + 1, count_);

	double time2 = timer2.seconds();
	std::cout << "Inclusive_scan took: " << time2 << " seconds." << std::endl;

	auto index_subview = subview(count_, nr_elements);
	auto h_ic = create_mirror_view(index_subview);

	// Deep copy device view to host view.
	deep_copy(h_ic, index_subview);

	ViewVectorType<Integer> elements = ViewVectorType<Integer>("hyper_elements",
		h_ic(0));

	Timer timer3;

	parallel_for(nr_elements,
		ScatterHypersphere(mesh, center, radius, count_, elements));

	double time3 = timer3.seconds();
	std::cout << "ScatterHypersphere took: " << time3 << " seconds." << std::endl;

	fence();

	return elements;
}

//mesh is a device pointer and can not be used here for getting the nr_elements.
template<class Mesh_>
ViewVectorType<Integer> mark_active(const Mesh_* mesh, const Integer nr_elements)
{
	using namespace Kokkos;


	Timer timer2;

	ViewVectorType<Integer> scan = ViewVectorType<Integer>("count_hypers",
		nr_elements + 1);

	incl_excl_scan(0, nr_elements, mesh->get_view_active(), scan);

	double time2 = timer2.seconds();
	std::cout << "Inclusive_scan took: " << time2 << " seconds." << std::endl;

	auto index_subview = subview(scan, nr_elements);
	auto h_ic = create_mirror_view(index_subview);

	// Deep copy device view to host view.
	deep_copy(h_ic, index_subview);

	Timer timer3;

	ViewVectorType<Integer> elements = ViewVectorType<Integer>("hyper_elements",
		h_ic(0));

	parallel_for(nr_elements, KOKKOS_LAMBDA(const Integer element_id){

		if (mesh->is_active(element_id))
		{
			elements(scan(element_id)) = element_id;
		}
	});

	double time3 = timer3.seconds();
	std::cout << "ScatterHypersphere took: " << time3 << " seconds." << std::endl;

	fence();

	return elements;
}

}
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
