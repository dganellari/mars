#ifndef MARS_MESH_KOKKOS_HPP
#define MARS_MESH_KOKKOS_HPP

#include "mars_simplex.hpp"
/*
 #include "mars_edge_element_map.hpp"
 #include "mars_edge_node_map.hpp"
 #include "mars_dual_graph.hpp"
 #include "mars_red_green_refinement.hpp"
 */

#include "mars_visualization.hpp"
#include "mars_point.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include "mars_imesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

template<Integer Dim_, Integer ManifoldDim_>
class Mesh<Dim_,ManifoldDim_,KokkosImplementation>: public ParallelIMesh<Dim_> {
public:
	static const Integer Dim = Dim_;
	static const Integer ManifoldDim = ManifoldDim_;

	using Elem = mars::Simplex<Dim, ManifoldDim>;
	using SideElem = mars::Simplex<Dim, ManifoldDim-1>;
	using Point = mars::Vector<Real, Dim>;

	void reserve(const std::size_t n_elements, const std::size_t n_points) override
	{
		elements_size_ = n_elements;
		points_size_ = n_points;

		elements_ = ViewMatrixType<Integer>("elems", n_elements,
				ManifoldDim + 1);
		active_ = ViewVectorType<bool>("active_", n_elements);
		points_ = ViewMatrixType<Real>("pts", n_points, Dim);
	}

	void reserve_points(const std::size_t n_points)
	{
		points_size_ = n_points;
		points_ = ViewMatrixType<Real>("pts", n_points, Dim);
	}

	void reserve_elements(const std::size_t n_elements)
	{
		elements_size_ = n_elements;
		elements_ = ViewMatrixType<Integer>("elems", n_elements,
				ManifoldDim + 1);
		active_ = ViewVectorType<bool>("active_", n_elements);
	}

	/*inline Elem &elem(const Integer id) override
	 {
	 assert(id >= 0);
	 assert(id < n_elements());
	 return elements_k(id);
	 }

	 inline const Elem &elem(const Integer id) const override
	 {
	 assert(id >= 0);
	 assert(id < n_elements());
	 return elements_k(id);
	 }*/

	inline bool is_active(const Integer id) const override
	{
		assert(id >= 0);
		assert(id < n_elements());
		return active_(id);
	}

	/*inline bool is_valid(const Integer id) const
	 {
	 return id >= 0 && id < n_elements();
	 }

	 inline bool is_node_valid(const Integer id) const
	 {
	 return id >= 0 && id < n_nodes();
	 }

	 inline bool is_child(
	 const Integer parent_id,
	 const Integer child_id) const
	 {
	 return std::find(
	 elem(parent_id).children.begin(),
	 elem(parent_id).children.end(),
	 child_id) != elem(parent_id).children.end();
	 }

	 inline void set_active(const Integer id, const bool val)
	 {
	 assert(id >= 0);
	 assert(id < active_.size());
	 active_[id] = val;
	 }*/

	inline Integer add_point(const Point &point, const int row) override
	{
		for (unsigned int i = 0; i < Dim; ++i) {
			points_(row, i) = point[i];
		}

		return row;

	}

	//add point functor
	struct AddCenterHexPoint {

		ViewMatrixType<Real> points;
		Integer xDim;
		Integer yDim;
		Integer zDim;

		AddCenterHexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm,
				Integer zdm) :
				points(pts), xDim(xdm), yDim(ydm), zDim(zdm) {
		}

		KOKKOS_INLINE_FUNCTION
		void operator()(int z, int y, int x) const {

			Integer cube_index = xDim * yDim * z + xDim * y + x;

			const int n_cube_nodes = (2 * xDim + 1) * (2 * yDim + 1)
					* (2 * zDim + 1);

			//add center of the hex to the new points.
			int centerHex = n_cube_nodes / 2 + cube_index + 1;

			points(centerHex, 0) = static_cast<Real>(2 * x + 1)
					/ static_cast<Real>(2 * xDim);
			points(centerHex, 1) = static_cast<Real>(2 * y + 1)
					/ static_cast<Real>(2 * yDim);
			points(centerHex, 2) = static_cast<Real>(2 * z + 1)
					/ static_cast<Real>(2 * zDim);

		}
	};

	//add point functor
	struct AddPoint {

		ViewMatrixType<Real> points;
		Integer xDim;
		Integer yDim;
		Integer zDim;

		AddPoint(ViewMatrixType<Real> pts, Integer xdm) :
				points(pts), xDim(xdm) {
		}

		AddPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm) :
				points(pts), xDim(xdm), yDim(ydm) {
		}

		AddPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm, Integer zdm) :
				points(pts), xDim(xdm), yDim(ydm), zDim(zdm) {
		}


		KOKKOS_INLINE_FUNCTION
		void operator()(int i) const {

			for (int l = 0; l < Dim; ++l) {
				points(i, l) = static_cast<Real>(i) / static_cast<Real>(xDim);
			}

		}

		KOKKOS_INLINE_FUNCTION
		void operator()(int i, int j) const {

			int index = i * (yDim + 1) + j;
			points(index, 0) = static_cast<Real>(i) / static_cast<Real>(xDim);
			points(index, 1) = static_cast<Real>(j) / static_cast<Real>(yDim);

		}

		KOKKOS_INLINE_FUNCTION
		void operator()(int z, int y, int x) const {

			Integer in = index(xDim, yDim, x, y, z);
/*
			 If not center hex point existed then the rule of dividing by two holds.
			 The multiple of 2 indexes are the ones to be added to the mesh. The others discarded.
*/
			if (in % 2 == 0) {

				Integer i = in / 2;

				points(i, 0) = static_cast<Real>(x)
						/ static_cast<Real>(2 * xDim);
				points(i, 1) = static_cast<Real>(y)
						/ static_cast<Real>(2 * yDim);
				points(i, 2) = static_cast<Real>(z)
						/ static_cast<Real>(2 * zDim);
			}
		}
	};

	inline bool generate_points(const int xDim, const int yDim, const int zDim) {

		using namespace Kokkos;

		switch (ManifoldDim_) {

		case 1: {

			assert(xDim != 0);
			assert(yDim == 0);
//			assert(zDim == 0);

			const int n_nodes = xDim + 1;
			reserve_points(n_nodes);

			parallel_for(n_nodes, AddPoint(points_, xDim));

			return true;
		}
		case 2: {

			assert(xDim != 0);
			assert(yDim != 0);
//			assert(zDim == 0);

			const int n_nodes = (xDim + 1) * (yDim + 1);
			reserve_points(n_nodes);

			parallel_for(
					MDRangePolicy<Rank<2> >( { 0, 0 }, { xDim + 1, yDim + 1 }),
					AddPoint(points_, xDim, yDim));

			return true;
		}
		case 3: {
			assert(xDim != 0);
			assert(yDim != 0);
			assert(zDim != 0);

			const int n_tetra_nodes = 5 * (xDim * yDim * zDim)
					+ 2 * (xDim * yDim + xDim * zDim + yDim * zDim)
					+ (xDim + yDim + zDim) + 1;

			reserve_points(n_tetra_nodes);

			parallel_for(
					MDRangePolicy<Rank<3> >( {0, 0, 0}, {zDim, yDim, xDim}),
					AddCenterHexPoint(points_, xDim, yDim,zDim));

			fence();

			parallel_for(
					MDRangePolicy<Rank<3> >( { 0, 0, 0 },
							{ 2 * zDim + 1, 2 * yDim + 1, 2 * xDim + 1 }),
					AddPoint(points_, xDim, yDim, zDim));

			return true;
		}
		default: {
			return false;
		}
		}

	}

	//add elem functor
	struct AddElem {

		ViewMatrixType<Integer> elem;
		ViewVectorType<bool> active;

		ViewMatrixTexture<Integer,hex_n_nodes> hexes;
		ViewMatrixTextureC<unsigned int,hex_n_sides, hex_side_n_nodes> map_side_nodes;

		Integer xDim;
		Integer yDim;
		Integer zDim;

		AddElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
				Integer xdm) :
				elem(el), active(ac), xDim(xdm) {
		}

		AddElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
				Integer xdm, Integer ydm) :
				elem(el), active(ac), xDim(xdm), yDim(ydm) {
		}

		AddElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
				ViewMatrixTexture<Integer,hex_n_nodes> hxs, ViewMatrixTextureC<unsigned int,hex_n_sides, hex_side_n_nodes> map,
				Integer xdm, Integer ydm, Integer zdm) :
				elem(el), active(ac), hexes(hxs), map_side_nodes(map), xDim(xdm), yDim(ydm), zDim(zdm) {
		}

		KOKKOS_INLINE_FUNCTION
		void operator()(int index) const {

			switch (ManifoldDim_) {

			case 1: {

				for (int i = 0; i < ManifoldDim + 1; ++i) {
					elem(index, i) = index + i;
				}
				active(index) = true;

				break;
			}
			case 2: {

				const int offset = yDim + 1;

				//extracting i and j from the global index from the parallel for to make it coalesced.
				int i = (index / 2) / yDim;
				int j = (index / 2) % yDim;

				int add_to_i = (index +1) % 2;
				int add_to_j = index % 2;

				elem(index, 0) = i * offset + j;
				elem(index, 1) = (i + 1) * offset + (j + 1);
				elem(index, 2) = (i + add_to_i) * offset + (j + add_to_j);
				active(index) = true;

				break;
			}
			case 3: {

				const int n_cube_nodes = (2 * xDim + 1) * (2 * yDim + 1)
						* (2 * zDim + 1);

				//const Integer cube_index = cube_indexes[index];
				const Integer cube_index = index/24;

				//add center of the hex to the new points.
				int centerHex = n_cube_nodes / 2 + cube_index + 1;

				//const int index = 24 * cube_index + 4 * i + k;
				const int k = (index - 24 * cube_index) % 4;
				const int i = (index - 24 * cube_index) / 4;

				elem(index, 0) = hexes(cube_index, map_side_nodes(i, k)) / 2;
				elem(index, 1) = hexes(cube_index, map_side_nodes(i, 8)) / 2; // mid-face point always the last element.
				elem(index, 2) = ( // rotation to catch all combinations.
						k == 3 ?
								hexes(cube_index, map_side_nodes(i, 0)) :
								hexes(cube_index, map_side_nodes(i, k + 1))) / 2;
				elem(index, 3) = centerHex; // the center of the cube.

				active(index) = true;

				break;
			}
			default:
				break;
			}
		}

		//alternative approach for the 2D generation with a stride access of 2 as opposed to the first totally coalesced case.
		/*
		 KOKKOS_INLINE_FUNCTION
		 void operator()(int i, int j) const {

		 const int offset = yDim + 1;

		 int index = 2 * (i * xDim + j);

		 elem(index, 0) = i * offset + j;
		 elem(index, 1) = (i + 1) * offset + (j + 1);
		 elem(index, 2) = i * offset + (j + 1);

		 active(index) = true;

		 index = index + 1; // stride access of 2.

		 elem(index, 0) = i * offset + j;
		 elem(index, 1) = (i + 1) * offset + (j + 1); //just to write it clearer
		 elem(index, 2) = (i + 1) * offset + j;

		 active(index) = true;
		 }
		 */

		KOKKOS_INLINE_FUNCTION
		void operator()(int z, int y, int x) const {

			Integer cube_index = xDim * yDim * z + xDim * y + x;

			//build the hex27 element which serves to generalize the idea to many hex27.
			//without generating the hex27 element first there is no easy way to create the sides.
			build_hex27(IndexView<Integer,hex_n_nodes>(hexes, cube_index), xDim, yDim, 2 * x, 2 * y, 2 * z);

		}
	};


	inline bool generate_elements(const int xDim, const int yDim,
			const int zDim) {

		using namespace Kokkos;

		switch (ManifoldDim_) {

		case 1: {
			const int n_elements = xDim;
			reserve_elements(n_elements);

			parallel_for(n_elements, AddElem(elements_, active_, xDim));

			return true;
		}
		case 2: {
			const int n_elements = 2 * xDim * yDim;
			reserve_elements(n_elements);

			parallel_for(n_elements, AddElem(elements_, active_, xDim, yDim));

			/*parallel_for(MDRangePolicy<Rank<2> >( { 0, 0 }, { xDim, yDim }),
			 AddElem(elements_, active_, xDim, yDim));*/
			return true;
		}
		case 3: {
			const int n_elements = xDim * yDim * zDim * 24; //24 tetrahedrons on one hex27
			reserve_elements(n_elements);

			const int n_cube_nodes = (2 * xDim + 1) * (2 * yDim + 1)
					* (2 * zDim + 1);

			//texture view for random access
			ViewMatrixTexture<Integer, hex_n_nodes> hexes("hexes",	xDim * yDim * zDim);

			//compile time defined texture view
			ViewMatrixTextureC<unsigned int, hex_n_sides, hex_side_n_nodes> map_side_to_nodes("cube_sides");

			//copy the map from the host to the device, putting it to the compile time defined view.
			copy_matrix_from_host(hex_side_nodes, map_side_to_nodes,
					hex_n_sides, hex_side_n_nodes);

			AddElem el(elements_, active_, hexes, map_side_to_nodes, xDim, yDim, zDim);

			parallel_for(MDRangePolicy<Rank<3> >( { 0, 0, 0 }, { zDim, yDim,
					xDim }), el);

			fence();

			parallel_for(n_elements, el);

			return true;
		}
		default: {
			return false;
		}
		}

	}

	/*	inline __device__ __host__ void add_point1(const size_t row,const Integer xDim) {

	 for (int i = 0; i < Dim; ++i)
	 points_(row, i) = static_cast<Real>(row) / static_cast<Real>(xDim);

	 }*/

	/*inline Point &point(const Integer i) override
	 {
	 assert(i >= 0);
	 assert(i < points_size);
	 return points_k(i);
	 }

	 inline const Point &point(const Integer i) const override
	 {
	 assert(i >= 0);
	 assert(i < points_.size());
	 return points_(i);
	 }*/

	/*inline void points(const Integer id, std::vector<Point> &pts) const override
	 {
	 assert(id >= 0);
	 assert(id < n_elements());

	 auto &e = elements_[id];
	 pts.resize(ManifoldDim + 1);

	 for (Integer i = 0; i < ManifoldDim + 1; ++i) {
	 pts[i] = points_(e.nodes[i]);
	 }
	 }*/

	const ViewMatrixType<Real> &points() const //override
	{
		return points_;
	}

	ViewMatrixType<Real> get_view_points() const //override
	{
		return points_;
	}

	ViewMatrixType<Integer> get_view_elems() const //override
	{
		return elements_;
	}

	ViewVectorType<bool> get_view_active() const //override
	{
		return active_;
	}

	/*
	 void setPoints(std::vector<Point>&& points)
	 {
	 points_ = std::forward<std::vector<Point>>(points);
	 }

	 template<typename Iter>
	 void remove_point(const Iter pos){
	 points_.erase(pos);
	 }

	 inline Integer add_elem(const Elem &elem)
	 {
	 auto id = elements_.size();
	 elements_.push_back(elem);
	 elements_.back().id = id;
	 active_.push_back(true);
	 assert(elements_.back().id == id);
	 return elements_.back().id;
	 }*/

	template<std::size_t NNodes>
	inline Integer add_elem(const std::array<Integer, NNodes> &nodes,
			const int row) //override
			{
		/*Elem elem;
		 assert(nodes.size() == mars::n_nodes(elem));*/

		//for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
		for (int i = 0; i < ManifoldDim + 1; ++i) {
			elements_(row, i) = nodes[i];
		}

		return row;
	}

	/*inline Integer add_elem(const IElem &elem) override
	 {
	 assert(elem.type() == ManifoldDim);

	 const Elem * elem_ptr = dynamic_cast<const Elem *>(&elem);
	 if(elem_ptr) {
	 return add_elem(*elem_ptr);
	 }

	 // fallback for other types of elements
	 Elem elem_copy;

	 std::vector<Integer> e_nodes;
	 elem.get_nodes(e_nodes);

	 assert(e_nodes.size() == ManifoldDim + 1);

	 for(std::size_t i = 0; i < mars::n_nodes(elem_copy); ++i) {
	 elem_copy.nodes[i] = e_nodes[i];
	 }

	 return add_elem(elem_copy);
	 }

	 template<std::size_t NNodes>
	 Integer add_elem(const std::array<Integer, NNodes> &nodes)
	 {
	 static_assert(NNodes == std::size_t(ManifoldDim + 1), "does not have the correct number of nodes");
	 elements_.emplace_back();
	 auto &e = elements_.back();
	 e.id = elements_.size() - 1;
	 e.nodes = nodes;
	 active_.push_back(true);
	 assert(e.id == elements_.size() - 1);
	 return e.id;
	 }

	 Elem& add_elem()
	 {
	 elements_.emplace_back();
	 auto &e = elements_.back();
	 e.id = elements_.size() - 1;
	 //e.nodes = nodes;
	 active_.push_back(true);
	 return e;
	 }


	 inline void deactivate_children(const Integer id)
	 {
	 assert(id >= 0);
	 assert(id < n_elements());

	 for(auto c : elem(id).children) {
	 active_[c] = false;
	 }
	 }

	 void repair_element(const Integer element_id, const bool verbose = false)
	 {
	 assert(element_id >= 0);
	 assert(element_id < n_elements());

	 if(sorted_elements_) {
	 auto &e = elem(element_id);
	 std::sort(e.nodes.begin(), e.nodes.end());
	 // return;
	 }

	 auto &e = elem(element_id);
	 const Real vol = mars::volume(e, points_);

	 if(vol < 0.) {
	 if(verbose) {
	 std::cout << element_id << " has negative volume" << std::endl;
	 }

	 std::swap(e.nodes[ManifoldDim-1], e.nodes[ManifoldDim]);
	 assert(mars::volume(e, points_) > 0.);
	 }
	 }

	 template<typename Iter>
	 void find_elements_by_nodes(
	 const Iter nodes_begin,
	 const Iter nodes_end,
	 std::vector<Integer> &elements,
	 const bool match_all = true,
	 const bool skip_inactive = true
	 ) const
	 {
	 elements.clear();

	 const Integer nn = std::distance(nodes_begin, nodes_end);

	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(!is_active(i) && skip_inactive) continue;

	 const auto &e = elem(i);

	 if(match_all) {
	 Integer n_matching = 0;
	 for(auto e_n : e.nodes) {
	 for(auto n_it = nodes_begin; n_it != nodes_end; ++n_it) {
	 if(*n_it == e_n) {
	 ++n_matching;
	 break;
	 }
	 }
	 }

	 if(n_matching == nn) {
	 elements.push_back(i);
	 }

	 } else {
	 bool found = false;
	 for(auto e_n : e.nodes) {
	 for(auto n_it = nodes_begin; n_it != nodes_end; ++n_it) {
	 if(*n_it == e_n) {
	 elements.push_back(i);
	 found = true;
	 break;
	 }
	 }

	 if(found) break;
	 }
	 }
	 }
	 }

	 void repair(const bool verbose = false)
	 {
	 for(std::size_t i = 0; i < elements_.size(); ++i) {
	 repair_element(i, verbose);
	 }
	 }

	 bool is_boundary(const Integer id) const {
	 const auto &adj = dual_graph_.adj(id);

	 for(auto a : adj) {
	 if(a == INVALID_INDEX) return true;
	 }

	 return false;
	 }

	 bool is_boundary(const Integer id, const Integer side_num) const {
	 const auto &adj = dual_graph_.adj(id);
	 return adj[side_num] == INVALID_INDEX;
	 }

	 bool is_interface(const Integer id) const {
	 const auto &adj = dual_graph_.adj(id);

	 for(auto a : adj) {
	 if(a < INVALID_INDEX) return true;
	 }

	 return false;
	 }

	 void describe_boundary_elements(std::ostream &os)
	 {
	 std::cout << "-------------------------\n";
	 for(std::size_t i = 0; i < elements_.size(); ++i) {
	 if(active_[i] && is_boundary(i)) {
	 dual_graph().describe_adj(i, os);
	 }
	 }
	 std::cout << "-------------------------\n";
	 }

	 void describe_element(const Integer i, std::ostream &os, const bool print_sides = false) const
	 {
	 const auto &e = elem(i);
	 const Real vol = mars::volume(e, points_);
	 const auto b   = barycenter(e, points_);

	 os << "---------------------------------\n";
	 os << "[" << i << "]: vol: " << vol << ", ";
	 for(auto v : e.nodes) {
	 os << " " << v;
	 }

	 os << "\n";

	 if(print_sides) {
	 Simplex<Dim, ManifoldDim-1> side;
	 Matrix<Real, Dim, ManifoldDim-1> J;

	 os << "sides:\n";
	 for(Integer k = 0; k < n_sides(e); ++k) {
	 e.side(k, side);
	 os << "==============\n";
	 jacobian(side, points_, J);

	 const auto n = normal(side, points_);
	 const auto sign = dot(points_[side.nodes[0]] - b, n) > 0? 1 : -1;
	 const Real u_area = mars::unsigned_volume(side, points_);
	 const Real area   = sign * u_area;

	 // J.describe(os);
	 os << area << " == " << u_area << std::endl;
	 }
	 }

	 os << "---------------------------------\n";
	 os << "\n";
	 }

	 void describe(std::ostream &os, const bool print_sides = false) const
	 {
	 for(std::size_t i = 0; i < elements_.size(); ++i) {
	 if(!active_[i]) continue;

	 describe_element(i, os, print_sides);

	 // const auto &e = elements_[i];
	 // const Real vol = mars::volume(e, points_);
	 // const auto b   = barycenter(e, points_);

	 // os << "---------------------------------\n";
	 // os << "[" << i << "]: vol: " << vol << ", ";
	 // for(auto v : e.nodes) {
	 // 	os << " " << v;
	 // }

	 // os << "\n";

	 // if(print_sides) {
	 // 	Simplex<Dim, ManifoldDim-1> side;
	 // 	Matrix<Real, Dim, Dim-1> J;

	 // 	os << "sides:\n";
	 // 	for(Integer k = 0; k < n_sides(e); ++k) {
	 // 		e.side(k, side);
	 // 		os << "==============\n";
	 // 		jacobian(side, points_, J);

	 // 		const auto n = normal(side, points_);
	 // 		const auto sign = dot(points_[side.nodes[0]] - b, n) > 0? 1 : -1;
	 // 		const Real u_area = mars::unsigned_volume(side, points_);
	 // 		const Real area   = sign * u_area;

	 // 		J.describe(os);
	 // 		os << area << " == " << u_area << std::endl;
	 // 	}
	 // }

	 // os << "---------------------------------\n";
	 // os << "\n";
	 }

	 for(std::size_t i = 0; i < points_.size(); ++i) {
	 os << i << ") ";
	 points_[i].describe(os);
	 }
	 }*/

	inline Integer n_nodes() const override
	{
		return points_size_;
	}

	inline Integer n_elements() const override
	{
		return elements_size_;
	}

	inline Integer n_active_elements() const override
	{
		Integer ret = 0;
		for (Integer i = 0; i < elements_size_; ++i) {
			ret += active_(i);
		}

		return ret;
	}

	/*bool have_common_sub_surface(
	 const Integer e_index_1,
	 const Integer e_index_2,
	 const Integer required_common_nodes = ManifoldDim) const
	 {
	 const auto &e1 = elem(e_index_1);
	 const auto &e2 = elem(e_index_2);

	 Integer n_common_nodes = 0;
	 for(Integer i = 0; i < ManifoldDim + 1; ++i) {
	 for(Integer j = 0; j < ManifoldDim + 1; ++j) {
	 n_common_nodes += e1.nodes[i] == e2.nodes[j];
	 }
	 }

	 assert(n_common_nodes <= ManifoldDim);
	 return n_common_nodes == required_common_nodes;
	 }

	 bool have_common_side(
	 const Integer e_index_1,
	 const Integer e_index_2) const
	 {
	 return have_common_sub_surface(e_index_1, e_index_2, ManifoldDim);
	 }

	 Integer common_side_num(
	 const Integer e_index_1,
	 const Integer e_index_2) const
	 {
	 const auto &e1 = elem(e_index_1);
	 const auto &e2 = elem(e_index_2);

	 Simplex<Dim, ManifoldDim-1> side;
	 for(Integer k = 0; k < n_sides(e1); ++k) {
	 e1.side(k, side);

	 Integer nn = 0;

	 for(Integer i = 0; i < ManifoldDim; ++i) {
	 const auto side_node = side.nodes[i];

	 for(Integer j = 0; j < ManifoldDim + 1; ++j) {
	 if(side_node == e2.nodes[j]) {
	 nn++;
	 break;
	 }
	 }
	 }

	 if(nn == ManifoldDim) {
	 return k;
	 }
	 }

	 assert(false);
	 return INVALID_INDEX;
	 }

	 void describe_dual_graph(std::ostream &os) const
	 {
	 dual_graph_.describe(os);
	 }

	 Integer n_boundary_sides() const
	 {
	 assert( !dual_graph_.empty() && "requires that build_dual_graph is called first");

	 Integer ret = 0;
	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(!active_[i]) continue;

	 const auto &e = elem(i);
	 const auto &e_adj = dual_graph_.adj(i);
	 for(Integer k = 0; k < e_adj.size(); ++k) {
	 const Integer j = e_adj[k];
	 if(j == INVALID_INDEX) {
	 ret++;
	 }
	 }
	 }

	 return ret;
	 }

	 bool check_side_ordering() const
	 {
	 assert( !dual_graph_.empty() && "requires that build_dual_graph is called first");

	 if(ManifoldDim == 4) {
	 std::cerr << "not implemented for 4d yet" << std::endl;
	 return false;
	 }

	 Simplex<Dim, ManifoldDim-1> side, other_side;

	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(!active_[i]) continue;

	 const auto &e = elem(i);
	 const auto &e_adj = dual_graph_.adj(i);
	 for(Integer k = 0; k < e_adj.size(); ++k) {
	 const Integer j = e_adj[k];
	 if(j == INVALID_INDEX) continue;
	 e.side(k, side);

	 const auto &other = elem(j);
	 const auto &other_adj = dual_graph_.adj(j);


	 Integer other_side_index = 0;
	 {
	 auto it = std::find(other_adj.begin(), other_adj.end(), i);


	 if(it == other_adj.end()) {
	 std::cerr << "Bad dual graph for " <<  i << " <-> " << j << std::endl;
	 assert(it != other_adj.end());
	 return false;
	 }

	 other_side_index = std::distance(other_adj.begin(), it);
	 other.side(other_side_index, other_side);
	 }

	 auto it = std::find(other_side.nodes.begin(), other_side.nodes.end(), side.nodes[0]);
	 assert(it != other_side.nodes.end());

	 Integer other_offset = std::distance(other_side.nodes.begin(), it);

	 for(Integer q = 0; q < ManifoldDim; ++q) {
	 Integer other_q = other_offset - q;

	 if(other_q < 0) {
	 other_q += ManifoldDim;
	 }

	 if(side.nodes[q] != other_side.nodes[other_q]) {
	 std::cerr << "common face not matching for (" << i << ", " << k << ") and (" << j << ", " << other_side_index << ")" << std::endl;
	 std::cerr << "[ ";
	 for(auto s : side.nodes) {
	 std::cerr << s << " ";
	 }
	 std::cerr << " ]\n";

	 std::cerr << "[ ";
	 for(auto s : other_side.nodes) {
	 std::cerr << s << " ";
	 }
	 std::cerr << " ]\n";
	 break;
	 }
	 }
	 }
	 }

	 return true;
	 }

	 DualGraph<ManifoldDim> &dual_graph() { return dual_graph_; }
	 const DualGraph<ManifoldDim> &dual_graph() const { return dual_graph_; }

	 void update_dual_graph(const bool force = false)
	 {
	 dual_graph_.update(*this, force);
	 }

	 void build_dual_graph()
	 {
	 update_dual_graph();
	 }

	 Real volume() const
	 {
	 Real ret = 0.;
	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(is_active(i)) {
	 ret += mars::volume(elem(i), points());
	 }
	 }

	 return ret;
	 }

	 void renumber_nodes()
	 {
	 assert(n_elements() == n_active_elements());

	 Integer edge_index = 0;
	 EdgeNodeMap enm;
	 for(Integer i = 0; i < n_elements(); ++i) {
	 const auto &e = elem(i);

	 for(Integer k = 0; k < n_edges(e); ++k) {
	 Integer v1, v2;
	 e.edge(k, v1, v2);

	 Integer edge_id = enm.get(v1, v2);

	 if(edge_id == INVALID_INDEX) {
	 enm.update(v1, v2, edge_index++);
	 }
	 }
	 }

	 std::vector< std::pair<Real, Integer> > weight(
	 n_nodes(),
	 std::pair<Real, Integer>(0, INVALID_INDEX)
	 );

	 std::vector<Integer> hits(n_nodes(), 0);

	 for(auto en : enm) {
	 Integer v1 = en.first[0];
	 Integer v2 = en.first[1];

	 Real d = (point(v1) - point(v2)).norm();

	 weight[v1].first += d;
	 weight[v1].second = v1;

	 weight[v2].first += d;
	 weight[v2].second = v2;

	 ++hits[v1];
	 ++hits[v2];
	 }

	 for(std::size_t i = 0; i < weight.size(); ++i) {
	 weight[i].first /= hits[i];
	 }

	 std::sort(std::begin(weight), std::end(weight));
	 std::reverse(std::begin(weight), std::end(weight));

	 std::vector<Integer> new_index(n_nodes(), INVALID_INDEX);

	 for(std::size_t i = 0; i < weight.size(); ++i) {
	 new_index[weight[i].second] = i;
	 }

	 {
	 std::vector<Point> points(n_nodes());

	 for(std::size_t i = 0; i < weight.size(); ++i) {
	 points[i] = points_[weight[i].second];
	 }

	 points_ = std::move(points);
	 }

	 for(Integer i = 0; i < n_elements(); ++i) {
	 auto &e = elem(i);

	 for(Integer k = 0; k < e.nodes.size(); ++k) {
	 e.nodes[k] = new_index[e.nodes[k]];
	 }
	 }
	 }

	 void reorder_nodes(const bool descending_order = true)
	 {

	 if(descending_order) {
	 for(Integer i = 0; i < n_elements(); ++i) {
	 std::sort(elem(i).nodes.begin(), elem(i).nodes.end(), [](const Integer v1, const Integer v2) -> bool {
	 return v2 < v1;
	 });
	 }
	 } else {
	 for(Integer i = 0; i < n_elements(); ++i) {
	 std::sort(elem(i).nodes.begin(), elem(i).nodes.end(), [](const Integer v1, const Integer v2) -> bool {
	 return v1 < v2;
	 });
	 }
	 }
	 }


	 void clean_up()
	 {
	 std::vector< Elem > elements;
	 std::vector<Integer> tags;

	 elements.reserve(n_active_elements());
	 tags.reserve(elements.capacity());

	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(!is_active(i)) continue;

	 assert(elem(i).children.empty());

	 elements.push_back(elem(i));
	 elements.back().id        = elements.size() - 1;
	 elements.back().parent_id = INVALID_INDEX;

	 if(!tags_.empty()) {
	 tags.push_back(tags_[i]);
	 }
	 }

	 elements_ = std::move(elements);
	 tags_  	  = std::move(tags);

	 dual_graph_.clear();
	 active_.resize(n_elements());
	 std::fill(active_.begin(), active_.end(), true);

	 update_dual_graph();
	 }

	 void clear()
	 {
	 elements_.clear();
	 points_.clear();
	 tags_.clear();
	 dual_graph_.clear();
	 active_.clear();
	 }

	 inline Integer root(const Integer id) const
	 {
	 if(id == INVALID_INDEX) return INVALID_INDEX;

	 Integer current_id = id;
	 while(elem(current_id).parent_id != INVALID_INDEX) {
	 current_id = elem(current_id).parent_id;
	 }

	 return current_id;
	 }


	 void scale(const Real factor)
	 {
	 for(auto &p : points_) {
	 p *= factor;
	 }
	 }

	 std::vector<Integer> &tags() { return tags_; }
	 const std::vector<Integer> &tags() const { return tags_; }

	 Mesh(const bool sorted_elements = false)
	 : sorted_elements_(sorted_elements)
	 {}


	 void collect_sides(
	 const Integer tag,
	 std::vector<Simplex<Dim, ManifoldDim-1>> &sides,
	 const bool active_only = true)
	 {
	 sides.clear();
	 sides.reserve(n_elements());

	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(active_only && !is_active(i)) continue;

	 const auto &e = elem(i);

	 for(Integer k = 0; k < n_sides(e); ++k) {
	 if(e.side_tags[k] == tag) {
	 sides.emplace_back();
	 e.side(k, sides.back());
	 }
	 }
	 }
	 }

	 void remove_elements(const std::vector<Integer> &elems)
	 {
	 Integer n_elems = n_elements();
	 std::vector<bool> to_remove(n_elems, false);

	 Integer min_elem_index = n_elems;

	 for(auto e : elems) {
	 to_remove[e] = true;
	 min_elem_index = std::min(e, min_elem_index);

	 Integer parent_id = elem(e).parent_id;
	 if(parent_id != INVALID_INDEX) {
	 elem(parent_id).children.clear();
	 }
	 }

	 bool is_contiguous = true;
	 for(Integer i = min_elem_index + 1; i < n_elems; ++i) {
	 if(to_remove[i] != to_remove[i-1]) {
	 is_contiguous = false;
	 assert(false);
	 }

	 assert(to_remove[i]);
	 }

	 assert(is_contiguous);

	 std::vector<bool> is_node_referenced(n_nodes(), false);

	 Integer max_node_index = 0;
	 for(Integer i = 0; i < n_elems; ++i) {
	 if(to_remove[i]) continue;

	 for(auto n : elem(i).nodes) {
	 is_node_referenced[n] = true;
	 max_node_index = std::max(n, max_node_index);
	 }
	 }


	 elements_.resize(min_elem_index);
	 active_.resize(elements_.size());
	 tags_.resize(elements_.size());

	 //assume contiguous for runtime efficiency
	 points_.resize(max_node_index + 1);

	 dual_graph_.clear();
	 update_dual_graph();
	 }

	 bool is_conforming() const
	 {
	 for(Integer i = 0; i < n_elements(); ++i) {
	 if(!is_active(i)) continue;

	 auto &e   = elem(i);
	 auto &adj = dual_graph().adj(i);

	 for(Integer k = 0; k < n_sides(e); ++k) {
	 if(adj[k] == INVALID_INDEX) {
	 if(e.side_tags[k] == INVALID_INDEX) {
	 return false;
	 }
	 }
	 }
	 }

	 return true;
	 }*/

	inline Integer type() const override
	{
		return ManifoldDim;
	}

private:

	ViewMatrixType<Integer> elements_;
	ViewMatrixType<Real> points_;
	ViewVectorType<bool> active_;
	Integer elements_size_;
	Integer points_size_;

	std::vector<Integer> tags_;
	DualGraph<ManifoldDim> dual_graph_;
	bool sorted_elements_;
};

/*template<Integer Dim, Integer ManifoldDim>
 bool read_mesh(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false)
 {
 std::ifstream is(path);
 if(!is.good()) {
 return false;
 }

 int dim = -1;
 int n_elements = -1;
 int n_nodes = -1;
 int n_coords = -1;

 std::string line;
 while(is.good()) {
 std::getline(is, line);

 if(line == "dimension") {
 std::getline(is, line);
 dim = atoi(line.c_str());
 assert(dim == ManifoldDim);
 } else if(line == "elements") {
 std::getline(is, line);
 n_elements = atoi(line.c_str());

 for(Integer i = 0; i < n_elements; ++i) {
 assert(is.good());
 std::getline(is, line);
 std::stringstream ss(line);
 int attr, type;

 std::array<Integer, ManifoldDim+1> nodes;
 ss >> attr >> type;

 for(Integer k = 0; k < ManifoldDim+1; ++k) {
 ss >> nodes[k];
 }

 mesh.add_elem(nodes);
 }
 } else if(line == "vertices") {
 std::getline(is, line);
 n_nodes = atoi(line.c_str());
 std::getline(is, line);
 n_coords = atoi(line.c_str());
 assert(n_coords == Dim);

 Vector<Real, Dim> p;
 p.zero();
 for(Integer i = 0; i < n_nodes; ++i) {
 assert(is.good());

 for(Integer k = 0; k < n_coords; ++k) {
 is >> p(k);
 }

 mesh.add_point(p);
 }

 }
 }

 is.close();

 mesh.repair(verbose);
 return true;
 }

 template<Integer Dim, Integer ManifoldDim>
 bool write_mesh_MFEM(const std::string &path,const Mesh<Dim, ManifoldDim> &mesh)
 {
 std::ofstream os;
 os.open(path.c_str());
 if (!os.good()) {
 os.close();
 return false;
 }

 writeHeader(mesh, os);
 writeElements(mesh, os);
 writeVertices(mesh, os);

 os.close();
 //	clear();
 return true;

 }

 template<Integer Dim, Integer ManifoldDim>
 void writeHeader(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) {
 Integer dim = mesh.ManifoldDim;
 os << "MFEM mesh v1.0\n\ndimension\n";
 os << dim<<"\n\n";
 }

 template<Integer Dim, Integer ManifoldDim>
 void writeElements(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) {
 os << "elements\n";
 os << mesh.n_elements() << "\n";

 for (Integer k = 0; k < mesh.n_elements(); ++k) {
 if (!mesh.is_active(k))
 continue;

 if(mesh.tags().size() == mesh.n_elements())
 os<<mesh.tags()[k]<< " " << mesh.elem(k).type()<<" ";
 else
 os<<INVALID_INDEX<< " " << INVALID_INDEX<<" ";

 for (Integer i = 0; i < ManifoldDim + 1; ++i) {
 const Integer v = mesh.elem(k).nodes[i];
 os << v;
 if (i < ManifoldDim) {
 os << " ";
 }
 }
 os << "\n";
 }
 os << "\n";
 }

 template<Integer Dim, Integer ManifoldDim>
 void writeVertices(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) {
 Integer dim = mesh.Dim;
 os << "vertices\n";

 const std::vector<Vector<Real, Dim>> points = mesh.points();

 os << points.size() << "\n";
 os << dim << "\n";

 for (Integer i = 0; i < points.size(); ++i) {
 for (Integer d = 0; d < Dim; ++d) {
 os << points[i](d);
 if (d < Dim - 1) {
 os << " ";
 }
 }
 os << "\n";
 }

 }*/

/*inline bool mesh_hyper_cube(
 const std::array<Integer, 4> &dims,
 const Vector<Real, 4> &lobo,
 const Vector<Real, 4> &upbo,
 const Mesh<4, 4> &mesh)
 {

 return false;
 }*/

using ParallelMesh1 = mars::Mesh<1, 1, KokkosImplementation>;
using ParallelMesh2 = mars::Mesh<2, 2, KokkosImplementation>;
using ParallelMesh3 = mars::Mesh<3, 3, KokkosImplementation>;
using ParallelMesh4 = mars::Mesh<4, 4, KokkosImplementation>;
using ParallelMesh5 = mars::Mesh<5, 5, KokkosImplementation>;
using ParallelMesh6 = mars::Mesh<6, 6, KokkosImplementation>;

}
#endif //MARS_MESH_HPP
