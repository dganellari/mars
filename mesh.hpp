#ifndef MARS_MESH_HPP
#define MARS_MESH_HPP

#include "simplex.hpp"
#include "edge_element_map.hpp"
#include "edge_node_map.hpp"
#include "dual_graph.hpp"
#include "red_green_refinement.hpp"

#include "visualization.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>  

namespace mars {

	template<Integer Dim_, Integer ManifoldDim_ = Dim_>
	class Mesh {
	public:
		static const Integer Dim = Dim_;
		static const Integer ManifoldDim = ManifoldDim_;

		void reserve(
			const std::size_t n_elements,
			const std::size_t n_points)
		{
			elements_.reserve(n_elements);
			active_.reserve(n_elements);
			points_.reserve(n_points);
		}

		inline Simplex<Dim, ManifoldDim> &elem(const Integer id)
		{
			assert(id >= 0);
			assert(id < n_elements());
			return elements_[id];
		}

		inline const Simplex<Dim, ManifoldDim> &elem(const Integer id) const
		{
			assert(id >= 0);
			assert(id < n_elements());
			return elements_[id];
		}

		inline bool is_active(const Integer id) const
		{
			assert(id >= 0);
			assert(id < n_elements());
			return active_[id];
		}

		inline bool is_valid(const Integer id) const
		{
			return id >= 0 && id < n_elements();
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
		}

		inline Integer add_point(const Vector<Real, Dim> &point)
		{
			points_.push_back(point);
			return points_.size() - 1;
		}

		inline Vector<Real, Dim> &point(const Integer i)
		{
			assert(i >= 0);
			assert(i < points_.size());
			return points_[i];
		}
		
		inline const Vector<Real, Dim> &point(const Integer i) const
		{
			assert(i >= 0);
			assert(i < points_.size());
			return points_[i];
		}

		const std::vector<Vector<Real, Dim>> &points() const
		{
			return points_;
		}

		inline Integer add_elem(const Simplex<Dim, ManifoldDim> &elem)
		{
			auto id = elements_.size();
			elements_.push_back(elem);
			elements_.back().id = id;
			active_.push_back(true);
			assert(elements_.back().id == id);
			return elements_.back().id;
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

		inline void points(const Integer id, std::vector<Vector<Real, Dim>> &pts) const
		{
			assert(id >= 0);
			assert(id < n_elements());

			auto &e = elements_[id];
			pts.resize(ManifoldDim + 1);
			
			for(Integer i = 0; i < ManifoldDim + 1; ++i) {
				pts[i] = points_[e.nodes[i]];
			}
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
				Matrix<Real, Dim, Dim-1> J;

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
		}

		inline Integer n_nodes() const
		{
			return points_.size();
		}

		inline Integer n_elements() const
		{
			return elements_.size();
		}

		inline Integer n_active_elements() const
		{
			Integer ret = 0;
			for(auto a : active_) {
				ret += a;
			}

			return ret;
		}

		bool have_common_sub_surface(
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
				std::vector<Vector<Real, Dim>> points(n_nodes());

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
			std::vector< Simplex<Dim, ManifoldDim> > elements;
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
		}

	private:
		std::vector< Simplex<Dim, ManifoldDim> > elements_;
		std::vector< Vector<Real, Dim> > points_;
		std::vector<Integer> tags_;
		DualGraph<ManifoldDim> dual_graph_;
		std::vector<bool> active_;
		bool sorted_elements_;
	};


	template<Integer Dim, Integer ManifoldDim>
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

	inline bool mesh_hyper_cube(
		const std::array<Integer, 4> &dims,
		const Vector<Real, 4> &lobo,
		const Vector<Real, 4> &upbo,
		const Mesh<4, 4> &mesh)
	{

		return false;
	}
}

#endif //MARS_MESH_HPP
