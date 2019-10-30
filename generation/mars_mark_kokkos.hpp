#ifndef MARS_MARK_KOKKOS_HPP_
#define MARS_MARK_KOKKOS_HPP_

#include "mars_mesh_kokkos.hpp"

namespace mars {

template<class Mesh_>
struct Hypersphere
{
	static constexpr Integer Dim = Mesh_::Dim;
	static constexpr Integer ManifoldDim = Mesh_::ManifoldDim;

	Mesh_* mesh;
	ViewVectorTypeC<Real, Dim> center;
	Real radius;

	ViewVectorType<Integer> count;

	Hypersphere(Mesh_* ms, ViewVectorTypeC<Real, Dim> c, Real rd,
			ViewVectorType<Integer> pc) :
			mesh(ms), center(c), radius(rd), count(pc)
	{
	}

	Hypersphere(Mesh_* ms, ViewVectorTypeC<Real, Dim> c, Real rd) :
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

				Point<Real, Dim> pt_center(center);
				auto dir = mesh->point(e.nodes[i]) - pt_center;
				auto d = dir.norm();

				if(d < radius) {
					inside = true;
				} else if(d > radius) {
					outside = true;
				} else if(::abs(d) < 1e-16) {
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

	ScatterHypersphere(Mesh_* ms, ViewVectorTypeC<Real, Dim> c, Real rd,
			ViewVectorType<Integer> sc, ViewVectorType<Integer> elems) :
			Hypersphere<Mesh_>(ms, c, rd), scan(sc), elements(elems)
	{
	}

	MARS_INLINE_FUNCTION
	void operator()(const int element_id) const
	{
		bool inside = false;
		bool outside = false;

		if (this->mesh->is_active(element_id)) //TODO: remove if by compacting on active elems.
		{
			auto e = this->mesh->elem(element_id);

			for(Integer i=0; i<ManifoldDim+1; ++i){

				Point<Real, Dim> pt_center(this->center);

				auto dir = this->mesh->point(e.nodes[i]) - pt_center;
				auto d = dir.norm();

				if(d < this->radius) {
					inside = true;
				} else if(d > this->radius) {
					outside = true;
				} else if(::abs(d) < 1e-16) {
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



template<class Mesh_>
ViewVectorType<Integer> mark_hypersphere_for_refinement(Mesh_* mesh,
const ViewVectorTypeC<Real, Mesh_::Dim> center, const Real radius, const Integer nr_elements)
{
	using namespace Kokkos;

	Timer timer1;

	ViewVectorType<Integer> count_ = ViewVectorType<Integer>("count_hypers",
		nr_elements + 1);

	parallel_for(nr_elements,
			Hypersphere<Mesh_>(mesh, center, radius, count_));

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
		ScatterHypersphere<Mesh_>(mesh, center, radius, count_, elements));

	double time3 = timer3.seconds();
	std::cout << "ScatterHypersphere took: " << time3 << " seconds." << std::endl;

	fence();

	return elements;
}

//mesh is a device pointer and can not be used here for getting the nr_elements.
template<class Mesh_>
ViewVectorType<Integer> mark_active(const Mesh_* mesh, const ViewVectorType<bool> active, const Integer nr_elements)
{
	using namespace Kokkos;


	Timer timer2;

	ViewVectorType<Integer> scan = ViewVectorType<Integer>("count_hypers",
		nr_elements + 1);

	incl_excl_scan(0, nr_elements, active, scan);

	double time2 = timer2.seconds();
	std::cout << "Inclusive_scan took: " << time2 << " seconds." << std::endl;

	auto index_subview = subview(scan, nr_elements);
	auto h_ic = create_mirror_view(index_subview);

	// Deep copy device view to host view.
	deep_copy(h_ic, index_subview);
	std::cout << "Hyper count result: " << h_ic(0)<< std::endl;

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

template<typename T, Integer Dim>
void fill_view(ViewVectorTypeC<T, Dim> point, const T value)
{
	using namespace Kokkos;

	typename ViewVectorTypeC<T, Dim>::HostMirror h_pt = create_mirror_view(point);

	for (Integer i=0; i<Dim; ++i)
		h_pt(i) = value;

	deep_copy(point, h_pt);
}

}
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
