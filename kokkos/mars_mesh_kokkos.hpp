#ifndef MARS_MESH_KOKKOS_HPP
#define MARS_MESH_KOKKOS_HPP

#include "mars_visualization.hpp"
#include "mars_point.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include "mars_imesh_kokkos.hpp"
#include "mars_simplex_kokkos.hpp"
#include "mars_non_simplex_kokkos.hpp"

namespace mars {

template<Integer Dim_, Integer ManifoldDim_, class Simplex_>
class Mesh<Dim_,ManifoldDim_,KokkosImplementation, Simplex_>: public ParallelIMesh<Dim_> {
public:
	static constexpr Integer Dim = Dim_;
	static constexpr Integer ManifoldDim = ManifoldDim_;

	using Elem = Simplex_;
	using SideElem = mars::Simplex<Dim, ManifoldDim-1,KokkosImplementation>;
	using Point = mars::Point<Real,Dim>;
	using SerialMesh = mars::Mesh<Dim_,ManifoldDim_>;
	using Edge = mars::Side<2, KokkosImplementation>;
	using Comb = Combinations<ManifoldDim + 1, 2, KokkosImplementation>;

	MARS_INLINE_FUNCTION Mesh()
		: ParallelIMesh<Dim_>()
		, elements_size_(0)
		, points_size_(0)
		, children_size_(0)
		, sorted_elements_(false)
		//, combinations(nullptr)
	{}

    void reserve(const std::size_t n_elements, const std::size_t n_points) override
    {
        elements_size_ = n_elements;
        points_size_ = n_points;
        children_size_ = 0;

        elements_ = ViewMatrixType<Integer>("elems", n_elements, Elem::ElemType);
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
		children_size_=0;

		elements_ = ViewMatrixType<Integer>("elems", n_elements,
				Elem::ElemType);
		active_ = ViewVectorType<bool>("active_", n_elements);
        reserve_parents(n_elements);
	}

    void reserve_parents(const Integer n_elements)
    {
        parent_ = ViewVectorType<Integer>("parent_", n_elements);
        Kokkos::parallel_for("init_parent", n_elements,
                MARS_LAMBDA(const Integer i){
                    parent_(i) = -1;
                });
    }

    void reserve_elem_children_map(const Integer size)
    {
        if(elem_children_map_.is_allocated())
        {
            elem_children_map_.rehash(size);
        }
        else
        {
            elem_children_map_ = UnorderedMap<Integer, Integer>(size);
        }
    }

	/*inline Elem &elem(const Integer id) override
	{
		assert(id >= 0);
		assert(id < n_elements());
		return simplices_(id);
	}

	inline const Elem &elem(const Integer id) const override
	{
		assert(id >= 0);
		assert(id < n_elements());
		return simplices_(id);
	}*/

	MARS_INLINE_FUNCTION Elem elem(const Integer id)// override
	{
		assert(id >= 0);
		assert(id < n_elements());
		Elem e = Elem(SubView<Integer,Elem::ElemType>(&elements_,id), combinations);
        e.set_view_parent(parent_);
		e.id = id;
		return e;
	}

	MARS_INLINE_FUNCTION const Elem elem(const Integer id) const //override
	{
		assert(id >= 0);
		assert(id < n_elements());
		Elem e = Elem(SubView<Integer,Elem::ElemType>(&elements_,id), combinations);
        e.set_view_parent(parent_);
		e.id = id;
		return e;
	}


    MARS_INLINE_FUNCTION const Elem elem_with_children(const Integer id) const //override
	{
		assert(id >= 0);
		assert(id < n_elements());

        auto children = get_children(id);
		Elem e = Elem(SubView<Integer,Elem::ElemType>(&elements_,id), children);
        e.set_view_parent(parent_);
		e.id = id;
		return e;
	}


    MARS_INLINE_FUNCTION
    const SubView<Integer, 2> get_children(const Integer element_id) const
    {
        SubView<Integer, 2> sv;

        const auto it = get_elem_children_map().find(element_id);
        if (get_elem_children_map().valid_at(it))
        {
            const Integer childrens_id = get_elem_children_map().value_at(it);
            sv = SubView<Integer, 2>(&children_, childrens_id);
        }
        else
        {
            sv.set_valid(false);
        }

        return sv;
    }

    MARS_INLINE_FUNCTION
    bool get_children(Integer children[2], const Integer element_id) const
    {
        SubView<Integer, 2> sv = get_children(element_id);
        if(sv.is_valid())
        {
            children[0] = sv(0);
            children[1] = sv(1);
            return true;
        }
        else
        {
            return false;
        }
    }
        /*MARS_INLINE_FUNCTION
	Elem elem(const Integer id,
			ViewMatrixTextureC<Integer, Comb::value, 2> combs)
	{
		assert(id >= 0);
		assert(id < n_elements());
		Elem e = Elem(SubView<Integer, ManifoldDim + 1>(elements_, id), combs);
		e.id = id;
		return e;
	}

	MARS_INLINE_FUNCTION
	const Elem elem(const Integer id,
			ViewMatrixTextureC<Integer, Comb::value, 2> combs) const
	{
		assert(id >= 0);
		assert(id < n_elements());
		Elem e = Elem(SubView<Integer, ManifoldDim + 1>(elements_, id), combs);
		e.id = id;
		return e;
	}*/

	MARS_INLINE_FUNCTION Elem elem(const Integer el_id, const Integer ch_id)// override
	{
		assert(el_id >= 0 && ch_id >=0);
		assert(el_id < n_elements()); // no assert for the children since children nr is not yet updated.
		Elem e = Elem(SubView<Integer,Elem::ElemType>(&elements_,el_id), SubView<Integer,2>(&children_,ch_id));
        e.set_view_parent(parent_);
		e.id = el_id;
		return e;
	}

	MARS_INLINE_FUNCTION const Elem elem(const Integer el_id, const Integer ch_id) const //override
	{
		assert(el_id >= 0 && ch_id >=0);
		assert(el_id < n_elements()); // no assert for the children since children nr is not yet updated.
		Elem e = Elem(SubView<Integer,Elem::ElemType>(&elements_,el_id), SubView<Integer,2>(&children_,ch_id));
        e.set_view_parent(parent_);
		e.id = el_id;
		return e;
	}

	MARS_INLINE_FUNCTION bool is_active(const Integer id) const override
	{
		assert(id >= 0);
		assert(id < n_elements());
		return active_(id);
	}

    MARS_INLINE_FUNCTION bool is_valid(const Integer id) const
    {
        return id >= 0 && id < n_elements();
    }

    MARS_INLINE_FUNCTION bool is_node_valid(const Integer id) const
    {
        return id >= 0 && id < n_nodes();
    }

    MARS_INLINE_FUNCTION
    const ViewMatrixTextureC<Integer, Comb::value, 2> & combs() const
	{
		return combinations;
	}

	MARS_INLINE_FUNCTION
	void set_combs(const ViewMatrixTextureC<Integer, Comb::value, 2>& c)
	{
		combinations = c;
	}

	 /* inline bool is_child(
	 const Integer parent_id,
	 const Integer child_id) const
	 {
	 return std::find(
	 elem(parent_id).children.begin(),
	 elem(parent_id).children.end(),
	 child_id) != elem(parent_id).children.end();
	 }
*/

	MARS_INLINE_FUNCTION
	void set_active(const Integer id, const bool val=true)
	{
		assert(id >= 0);
		assert(id < n_elements());
		active_(id) = val;
	}

    MARS_INLINE_FUNCTION
    void add_point(const Point &point, const Integer index)
    {
        for (Integer i = 0; i < Dim; ++i)
        {
            points_(index, i) = point[i];
        }
    }

    MARS_INLINE_FUNCTION
	void add_point(const TempArray<Real, Dim> &point, const Integer index)
	{
		for (Integer i = 0; i < Dim; ++i)
		{
			points_(index, i) = point[i];
		}
	}

	MARS_INLINE_FUNCTION
	Point point(const Integer i) override
	{
		assert(i >= 0);
		assert(i < points_size_);
		return Point(points_, i);
	}

	MARS_INLINE_FUNCTION
	const Point point(const Integer i) const override
	{
		assert(i >= 0);
		assert(i < points_size_);
		return Point(points_, i);
	}

	const ViewMatrixType<Real> &points() const //override
	{
		return points_;
	}

	ViewMatrixType<Real> get_view_points() const //override
	{
		return points_;
	}

	ViewMatrixType<Integer> get_view_elements() const //override
	{
		return elements_;
	}

	MARS_INLINE_FUNCTION
	ViewVectorType<bool> get_view_active() const //override
	{
		return active_;
	}

	ViewMatrixType<Integer> get_view_children() const //override
	{
		return children_;
	}

    MARS_INLINE_FUNCTION
    const UnorderedMap<Integer, Integer>& get_elem_children_map() const
    {
        return elem_children_map_;
    }

	void resize_points(const Integer size)
	{
		points_size_+= size;
		resize(points_, points_size_, Dim);
	}


	void resize_elements(const Integer size)
	{
		elements_size_+= size;
		resize(elements_, elements_size_, Elem::ElemType);
		resize(active_, elements_size_);
		resize(parent_, elements_size_);
	}

	void resize_active(const Integer size)
	{
		resize(active_, elements_size_+size);
	}

	void resize_children(const Integer size)
	{
		children_size_+= size;
		resize(children_, children_size_, 2);
	}

	/*
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


	MARS_INLINE_FUNCTION
	Integer n_nodes() const override
	{
		return points_size_;
	}

	MARS_INLINE_FUNCTION
	void set_n_nodes(Integer size)
	{
		points_size_ = size;
	}

	MARS_INLINE_FUNCTION
	void set_n_childrens(Integer size)
	{
		children_size_ = size;
	}

	MARS_INLINE_FUNCTION
	Integer n_elements() const override
	{
		return elements_size_;
	}

	MARS_INLINE_FUNCTION
	Integer n_childrens() const
	{
		return children_size_;
	}

	MARS_INLINE_FUNCTION
	void set_n_elements(Integer size)
	{
		elements_size_ = size;
	}

	/*Integer n_active_elements() const override
	{
		Integer ret = 0;
		for (Integer i = 0; i < elements_size_; ++i)
		{
			ret += active_(i);
		}

		return ret;
	}*/

	inline Integer n_active_elements(const Integer N)
	{
		using namespace Kokkos;

		Timer timer;

		ViewVectorType<bool> active_ = get_view_active();

		double result=0;
		parallel_reduce("Active_all", N, KOKKOS_LAMBDA (const int& i, double& lsum )
		{
			lsum += active_(i);
		},result);

		double time = timer.seconds();
		std::cout << "Active_all Reduction took: " << time << " seconds." << std::endl;

		printf("Result: %i %lf\n", N, result);

		return result;
	}


	inline Integer n_active_elements(const ViewVectorType<Integer> elements)
	{
		using namespace Kokkos;

		Timer timer;

		ViewVectorType<bool> active_ = get_view_active();

		Integer N = elements.extent(0);

		double result=0;
		parallel_reduce("Active_elem", N, KOKKOS_LAMBDA (const int& i, double& lsum )
		{
			lsum += active_(elements(i));
		},result);

		double time = timer.seconds();
		std::cout << "Active_elements Reduction took: " << time << " seconds." << std::endl;

		printf("Result: %i %lf\n", N, result);

		return result;
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

	 MARS_INLINE_FUNCTION Integer type() const override
	{
		return ManifoldDim;
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

//add point functor
	struct AddNonSimplexPoint {

		ViewMatrixType<Real> points;
		Integer xDim;
		Integer yDim;
		Integer zDim;


		AddNonSimplexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm) :
				points(pts), xDim(xdm), yDim(ydm) {
		}

		AddNonSimplexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm, Integer zdm) :
				points(pts), xDim(xdm), yDim(ydm), zDim(zdm) {
		}


		KOKKOS_INLINE_FUNCTION
		void operator()(int j, int i) const {

			int index = j * (xDim + 1) + i;

			points(index, 0) = static_cast<Real>(i) / static_cast<Real>(xDim);
			points(index, 1) = static_cast<Real>(j) / static_cast<Real>(yDim);
		}

		KOKKOS_INLINE_FUNCTION
		void operator()(int z, int y, int x) const {

			int index =  (xDim + 1) * (yDim +1) * z + (xDim +1) * y + x;

			points(index, 0) = static_cast<Real>(x)
					/ static_cast<Real>(xDim);
			points(index, 1) = static_cast<Real>(y)
					/ static_cast<Real>(yDim);
			points(index, 2) = static_cast<Real>(z)
					/ static_cast<Real>(zDim);
		}

	};

	inline bool generate_points(const int xDim, const int yDim, const int zDim,
				Integer type) {

		using namespace Kokkos;

		switch (ManifoldDim_) {

		case 2: {

			switch(type){

			case ElementType::Quad4: {
				assert(xDim != 0);
				assert(yDim != 0);
				assert(zDim == 0);

				const int n_nodes = (xDim + 1) * (yDim + 1);
				reserve_points(n_nodes);

				parallel_for(
						MDRangePolicy<Rank<2> >( { 0, 0 }, { yDim + 1, xDim + 1 }),
						AddNonSimplexPoint(points_, xDim, yDim));
				return true;
			}
			default: {
				return false;
			}
			}
		}
		case 3: {

			switch(type){

			case ElementType::Hex8: {
				assert(xDim != 0);
				assert(yDim != 0);
				assert(zDim != 0);

				const int n_nodes = (xDim + 1) * (yDim + 1)	* (zDim + 1);

				reserve_points(n_nodes);

				parallel_for(
						MDRangePolicy<Rank<3> >( { 0, 0, 0 },
								{ zDim + 1, yDim + 1, xDim + 1}),
						AddNonSimplexPoint(points_, xDim, yDim, zDim));
				return true;
			}
			default: {
				return false;
			}
			}
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
				elem(el), active(ac), xDim(xdm)
		{
		}

		AddElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
				Integer xdm, Integer ydm) :
				elem(el), active(ac), xDim(xdm), yDim(ydm)
		{
		}

		AddElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
				ViewMatrixTexture<Integer,hex_n_nodes> hxs, ViewMatrixTextureC<unsigned int,hex_n_sides, hex_side_n_nodes> map,
				Integer xdm, Integer ydm, Integer zdm) :
				elem(el), active(ac), hexes(hxs), map_side_nodes(map), xDim(xdm), yDim(ydm), zDim(zdm)
				{
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

				if (add_to_i == 0)
				{
					elem(index, 1) = (i + 1) * offset + (j + 1);
					elem(index, 2) = (i + add_to_i) * offset + (j + add_to_j);
				}
				else
				{
					elem(index, 1) = (i + add_to_i) * offset + (j + add_to_j);
					elem(index, 2) = (i + 1) * offset + (j + 1);
				}

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
		void operator()(int z, int y, int x) const
		{
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

	//add elem functor
	struct AddNonSimplexElem
	{

		ViewMatrixType<Integer> elem;
		ViewVectorType<bool> active;

		Integer xDim;
		Integer yDim;
		Integer zDim;

		AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
							Integer xdm) : elem(el), active(ac), xDim(xdm)
		{
		}

		AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
							Integer xdm, Integer ydm) : elem(el), active(ac), xDim(xdm), yDim(ydm)
		{
		}

		AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
							Integer xdm, Integer ydm, Integer zdm) : elem(el), active(ac), xDim(xdm), yDim(ydm), zDim(zdm)
		{
		}

        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const
        {
            const int offset = xDim + 1;
            int index = j * xDim + i;

            elem(index, 0) = i + offset * j;
            elem(index, 1) = (i + 1) + offset * j;
            elem(index, 2) = (i + 1) + offset * (j + 1);
            elem(index, 3) = i + offset * (j + 1);

            active(index) = true;
        }

        KOKKOS_INLINE_FUNCTION
		void operator()(int k, int j, int i) const
		{
			int index = k * xDim * yDim + j * xDim + i;

			elem(index, 0) = elem_index(i, j, k, xDim, yDim);
			elem(index, 1) = elem_index(i + 1, j, k, xDim, yDim);
			elem(index, 2) = elem_index(i + 1, j + 1, k, xDim, yDim);
			elem(index, 3) = elem_index(i, j + 1, k, xDim, yDim);
			elem(index, 4) = elem_index(i, j, k + 1, xDim, yDim);
			elem(index, 5) = elem_index(i + 1, j, k + 1, xDim, yDim);
			elem(index, 6) = elem_index(i + 1, j + 1, k + 1, xDim, yDim);
			elem(index, 7) = elem_index(i, j + 1, k + 1, xDim, yDim);

			active(index) = true;
		}
	};

    inline bool generate_elements(const int xDim, const int yDim,
			 const int zDim, Integer type) {

		using namespace Kokkos;

		switch (ManifoldDim_) {

		case 2: {

			switch(type){

			case ElementType::Quad4: {
				const int n_elements = xDim * yDim;
				reserve_elements(n_elements);

				parallel_for(MDRangePolicy<Rank<2> >( { 0, 0 }, {yDim, xDim}),
			 				AddNonSimplexElem(elements_, active_, xDim, yDim));
				return true;
			}
			default: {
				return false;
			}
			}
		}
		case 3: {

			switch(type){

			case ElementType::Hex8: {
				const int n_elements = xDim * yDim * zDim;
				reserve_elements(n_elements);

				parallel_for(MDRangePolicy<Rank<3> >( { 0, 0, 0 }, {zDim, yDim, xDim}),
							AddNonSimplexElem(elements_, active_, xDim, yDim, zDim));
				return true;
			}
			default: {
				return false;
			}
			}
		}
		default: {
			return false;
		}
		}
	}

    struct RefineMesh
    {
        Mesh *mesh;
        Integer xDim;

        void init(const Mesh &m)
        {
            Mesh *tmp = (Mesh *)Kokkos::kokkos_malloc(sizeof(Mesh));

            Kokkos::parallel_for("CreateMeshObject", 1, KOKKOS_LAMBDA(const int &) {
                new ((Mesh *)tmp) Mesh(m);
            });

            mesh = tmp; //otherwise "this" pointer is still a host pointer within the parallel_for.
        }

        RefineMesh(const Mesh &m, Integer xdm) : xDim(xdm)
        {
            init(m);
        }
        /*~RefineMesh() {
			Kokkos::kokkos_free(m);
		}*/
        KOKKOS_INLINE_FUNCTION
        void operator()(int x) const
        {
            LongestEdgeSelect<Mesh> es;
            es.select(*mesh, x);
        }
    };

    inline bool refine_mesh(const int xDim)
	{
		using namespace Kokkos;
		const int n_elements = xDim;

		RefineMesh ref= RefineMesh(*this, xDim);
		parallel_for(n_elements, ref);

		const Mesh* tmp = ref.mesh;

		parallel_for("DestroyMeshObject",1, KOKKOS_LAMBDA (const int&) {
			tmp->~Mesh();
		});

		return true;
	}

private:
	ViewMatrixTextureC<Integer, Comb::value, 2> combinations;

	ViewMatrixType<Integer> elements_;
	ViewMatrixType<Real> points_;
	ViewVectorType<bool> active_;
	Integer elements_size_;
	Integer points_size_;

	ViewMatrixType<Integer> side_tags_;
	ViewMatrixType<Integer> children_;
	Integer children_size_;
    UnorderedMap<Integer, Integer> elem_children_map_;
    ViewVectorType<Integer> parent_;

	ViewVectorType<Integer> tags_;
//	DualGraph<ManifoldDim> dual_graph_;
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

 std::array<Integer, Elem::ElemType> nodes;
 ss >> attr >> type;

 for(Integer k = 0; k < Elem::ElemType; ++k) {
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

 }

/*inline bool mesh_hyper_cube(
 const std::array<Integer, 4> &dims,
 const Vector<Real, 4> &lobo,
 const Vector<Real, 4> &upbo,
 const Mesh<4, 4> &mesh)
 {

 return false;
 }*/

	using ParallelMesh1 = mars::Mesh<1, 1, KokkosImplementation, Simplex<1,1, KokkosImplementation>>;
	using ParallelMesh2 = mars::Mesh<2, 2, KokkosImplementation, Simplex<2,2, KokkosImplementation>>;
	using ParallelMesh3 = mars::Mesh<3, 3, KokkosImplementation, Simplex<3,3, KokkosImplementation>>;
	using ParallelMesh4 = mars::Mesh<4, 4, KokkosImplementation, Simplex<4,4, KokkosImplementation>>;


	using ParallelQuad4Mesh = mars::Mesh<2, 2, KokkosImplementation, Quad4PElem>;
	using ParallelHex8Mesh = mars::Mesh<3, 3, KokkosImplementation, Hex8PElem>;

	template<Integer Type>
	using ParallelNSMesh2 = mars::Mesh<2, 2, KokkosImplementation, NonSimplex<Type, KokkosImplementation>>;
}
#endif //MARS_MESH_HPP
