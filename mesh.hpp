#ifndef MARS_MESH_HPP
#define MARS_MESH_HPP

#include "simplex.hpp"
#include "edge_element_map.hpp"
#include "edge_node_map.hpp"

namespace moonolith {
	using Integer = mars::Integer;
}

#include "moonolith_config.hpp"
#include "moonolith_mesh.hpp"
// #include "moonolith_svg_canvas.hpp"
#include "moonolith_eps_canvas.hpp"
#include "moonolith_plotter.hpp"
#include "moonolith_func_to_color.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>  

namespace mars {
	enum RefinementFlag {
		NONE = 0,
		RED = 1,
		GREEN_1 = 2,
		GREEN_2 = 3,
		GREEN_3 = 4,
		CHILD_OF_GREEN = 5,
		PARENT_PROMOTED_TO_RED = 6
	};

	template<Integer Dim, Integer ManifoldDim = Dim>
	class Mesh {
	public:
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
			return elements_[id];
		}

		inline const Simplex<Dim, ManifoldDim> &elem(const Integer id) const
		{
			return elements_[id];
		}

		inline bool is_active(const Integer id) const
		{
			return active_[id];
		}

		inline Integer add_point(const Vector<Real, Dim> &point)
		{
			points_.push_back(point);
			return points_.size() - 1;
		}

		inline const Vector<Real, Dim> &point(const Integer i) const
		{
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
			refinement_flag_.push_back(NONE);
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
			refinement_flag_.push_back(NONE);
			assert(e.id == elements_.size() - 1);
			return e.id;
		}

		inline void set_refinement_flag(const Integer &element_id, const Integer flag)
		{
			refinement_flag_[element_id] = flag;
		}

		inline Integer refinement_flag(const Integer &element_id) const
		{
			return refinement_flag_[element_id];
		}

		inline void points(const Integer element_id, std::vector<Vector<Real, Dim>> &pts)
		{
			auto &e = elements_[element_id];
			pts.resize(ManifoldDim + 1);
			
			for(Integer i = 0; i < ManifoldDim + 1; ++i) {
				pts[i] = points_[e.nodes[i]];
			}
		}

		inline bool is_green(const Integer id) const
		{
			if(INVALID_INDEX == id) return false;

			switch(refinement_flag(id)) 
			{
				case GREEN_1: return true;
				case GREEN_2: return true;
				case GREEN_3: return true;
				default: return false;
			}
		}	

		inline Integer n_adjacients(const Integer id) const
		{
			Integer ret = 0;
			for(auto a : dual_graph_[id]) {
				ret += a != INVALID_INDEX;
			}

			return ret;
		}


		inline void deactivate_children(const Integer id)
		{
			for(auto c : elem(id).children) {
				active_[c] = false;
			}
		}

		void propagate_flags(
			std::vector<Integer> &red_elements,
			std::vector<Integer> &green_elements)
		{
			std::vector<Integer> potential_green_elements;
			for(auto &e : red_elements) {
				const auto &adj = dual_graph_[e];
				for(auto a : adj) {
					if(a == INVALID_INDEX || !is_active(a)) continue;
					auto flag = refinement_flag(a);

					if(!elem(a).children.empty()) {
						std::cout << "promoting " << a << " to red" << std::endl;
						set_refinement_flag(a, RED);
						//deactivate children
						deactivate_children(a);
						red_elements.push_back(a); 
						continue;
					}

					const auto parent = elem(a).parent_id;
					if(is_green(parent)) {
						std::cout << "promoting " << a << " to red" << std::endl;
						set_refinement_flag(parent, RED);
						//deparentctivate children
						deactivate_children(parent);
						red_elements.push_back(parent);
					}

					switch(flag)
					{
						case RED:
						{
							continue;
						}
						case NONE:
						{
							if(n_adjacients(a) == 1) {
								set_refinement_flag(a, RED);
								red_elements.push_back(a);
							} else {
								set_refinement_flag(a, GREEN_1);
								potential_green_elements.push_back(a);
							}

							break;
						}
						case GREEN_1:
						{
							if(n_adjacients(a) == 2) {
								set_refinement_flag(a, RED);
								red_elements.push_back(a);
							} else {
								set_refinement_flag(a, GREEN_2);
								potential_green_elements.push_back(a);
							}

							break;
						}
						case GREEN_2:
						{
							if(n_adjacients(a) == 3) {
								set_refinement_flag(a, RED);
								red_elements.push_back(a);
							} else {
								//GREEN_3?
								set_refinement_flag(a, RED);
								red_elements.push_back(a);
							}

							break;
						}
						case GREEN_3:
						{
							assert(false && "implement me");
							break;
						}
					}
				}
			}

			green_elements.clear();
			green_elements.reserve(potential_green_elements.size());
			for(auto e : potential_green_elements) {
				if(refinement_flag(e) != RED && is_active(e)) {
					green_elements.push_back(e);
				}
			}
		}

		void red_green_refinement(const std::vector<Integer> &elements_to_refine)
		{
			std::vector<Integer> promoted_elements;
			std::vector<Integer> red_elements, green_elements;
			for(auto &e : elements_to_refine) {
				if(!is_active(e)) {
					std::cerr << e << " cannot refine inactive" << std::endl;
					continue;
				}

				auto parent = elem(e).parent_id;
				
				if(is_green(parent)) {
					std::cout << "promoting " << parent << " to red" << std::endl;
					//promote parent to red
					set_refinement_flag(parent, RED);
					//deactivate children
					deactivate_children(parent);
					active_[parent] = true;
					promoted_elements.push_back(parent);
					continue;
				}

				set_refinement_flag(e, RED);
				red_elements.push_back(e);
			}

			update_dual_graph();

			propagate_flags(promoted_elements, green_elements);

			for(auto e : promoted_elements) {
				red_refine_element(e);
			}

			for(auto e : green_elements) {
				green_refine_element(e);
			}

			update_dual_graph();
			edge_element_map_.build(*this);

			propagate_flags(red_elements, green_elements);
			for(auto e : red_elements) {
				red_refine_element(e);
			}

			for(auto e : green_elements) {
				green_refine_element(e);
			}
		}	

		inline void green_refine_element(const Integer element_id)
		{
			auto &e = elem(element_id);
			active_[element_id] = false;

			const auto &adj = dual_graph_[element_id];
			std::array<Integer, ManifoldDim> side_flags;
			std::vector<Integer> red_side_index;
			std::vector<Integer> green_side_index;

			Integer n_red_neighs = 0;
			Integer n_green_neighs = 0;
			
			Integer k = 0;
			for(auto a : adj) {
				if(a == INVALID_INDEX) {
					side_flags[k++] = INVALID_INDEX;
				 	continue;
				}

				if(refinement_flag(a) == RED) {
					red_side_index.push_back(k);
					++n_red_neighs;
				} else if(refinement_flag(a) == GREEN_1) {
					++n_green_neighs;
					green_side_index.push_back(k);
				}
				
				side_flags[k++] = refinement_flag(a);
			}

			Simplex<Dim, ManifoldDim-1> side_1, side_2;
			Simplex<Dim, ManifoldDim> child;
			child.parent_id = element_id;
			e.children.clear();

			switch(n_red_neighs) {
				case 1: 
				{	
					e.side(red_side_index[0], side_1);
					
					Integer n0 = side_1.nodes[0];
					Integer n1 = side_1.nodes[1];
					const Integer midpoint = edge_node_map_.get(n0, n1);
					const Integer opposite = e.vertex_opposite_to_side(red_side_index[0]);

					child.nodes[0] = n0;
					child.nodes[1] = opposite;
					child.nodes[2] = midpoint;

					e.children.push_back( add_elem(child) );

					child.nodes[0] = n1;
					child.nodes[1] = midpoint;
					child.nodes[2] = opposite;

					e.children.push_back( add_elem(child) );
					break;
				}

				case 2:
				{
					e.side(red_side_index[0], side_1);
					
					Integer n1_0 = side_1.nodes[0];
					Integer n1_1 = side_1.nodes[1];
					const Integer midpoint_1 = edge_node_map_.get(n1_0, n1_1);
					const Integer opposite_1 = e.vertex_opposite_to_side(red_side_index[0]);

					e.side(red_side_index[1], side_2);
					
					Integer n2_0 = side_2.nodes[0];
					Integer n2_1 = side_2.nodes[1];
					const Integer midpoint_2 = edge_node_map_.get(n2_0, n2_1);
					const Integer opposite_2 = e.vertex_opposite_to_side(red_side_index[1]);


					child.nodes[0] = midpoint_1;
					child.nodes[1] = opposite_2;
					child.nodes[2] = opposite_1;

					e.children.push_back( add_elem(child) ); 

					child.nodes[0] = midpoint_1;
					child.nodes[1] = midpoint_2;
					child.nodes[2] = n2_0;

					e.children.push_back( add_elem(child) ); 

					child.nodes[0] = midpoint_1;
					child.nodes[1] = n2_1;
					child.nodes[2] = midpoint_2;

					e.children.push_back( add_elem(child) ); 
					break;
				}

				case 3:
				{
					std::cerr << "[" << element_id << "] should have been promoted to RED" << std::endl;
					break;
				}
				default:
				{
					assert(false);
				}
			}
		}

		inline void red_refine_element(const Integer element_id)
		{
			static const Integer NSubs = NSubSimplices<ManifoldDim>::value;
			static_assert(NSubSimplices<ManifoldDim>::value > 0, "!");

			auto &parent_e = elements_[element_id];
			std::vector<Vector<Real, Dim>> parent_points;
			points(element_id, parent_points);

			std::array<Simplex<Dim, ManifoldDim>, NSubs> children;
			std::vector<Vector<Real, Dim>> children_points;
			auto interp = std::make_shared< SimplexInterpolator<ManifoldDim> >();

			Simplex<Dim, ManifoldDim> modified_e = parent_e;
			
			if(ManifoldDim == 4) {
				//4D hack
				std::sort(modified_e.nodes.begin(), modified_e.nodes.end());
			}

			red_refinement<Dim, ManifoldDim, NSubs>(
			    modified_e,
			    parent_points,
			    children, 
			    children_points,
			    *interp
			);

			if(interp_.size() <= element_id) {
				interp_.resize(element_id + 1);
			}

			std::vector<Integer> point_ids(interp->rows(), INVALID_INDEX);

			for(Integer i = 0; i < ManifoldDim + 1; ++i) {
				point_ids[i] = modified_e.nodes[i];

            	for(Integer j = i + 1; j < ManifoldDim + 1; ++j) {
                	Integer offset = midpoint_index<ManifoldDim>(i, j); 
                	point_ids[offset] = edge_node_map_.get(modified_e.nodes[i], modified_e.nodes[j]);

                	if(point_ids[offset] == INVALID_INDEX) {
                		const auto new_id = add_point(children_points[offset]);
                		edge_node_map_.update(
                			modified_e.nodes[i],
                			modified_e.nodes[j],
							new_id
						);

						point_ids[offset] = new_id;
						assert(new_id < this->n_nodes());
                	}
            	}
            }

			interp_[element_id] = interp;

			for(auto &c : children) {
				for(Integer i = 0; i < ManifoldDim + 1; ++i) {
					c.nodes[i] = point_ids[c.nodes[i]];
				}

				// std::sort(c.nodes.begin(), c.nodes.end()); add_elem(c);

				//4D hack
				repair_element(add_elem(c), ManifoldDim != 4);
			}

			active_[element_id] = false;
		}

		void repair_element(const Integer element_id, const bool verbose = false)
		{
			assert(element_id >= 0);

			auto &e = elements_[element_id];
			const Real vol = volume(e, points_);

			if(vol < 0.) {
				if(verbose) {
					std::cout << element_id << " has negative volume" << std::endl;
				}

				switch(Dim) {
					case 1:
					{
						std::swap(e.nodes[0], e.nodes[1]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					case 2:
					{
						std::swap(e.nodes[1], e.nodes[2]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					case 3:
					{
						std::swap(e.nodes[2], e.nodes[3]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					case 4:
					{
						std::swap(e.nodes[3], e.nodes[4]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					default: 
					{
						assert(false && "implement me");
						break;
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

		void describe(std::ostream &os, const bool print_sides = false) const
		{
			for(std::size_t i = 0; i < elements_.size(); ++i) {
				if(!active_[i]) continue;

				const auto &e = elements_[i];
				const Real vol = volume(e, points_);
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
						const Real u_area = unsigned_volume(side, points_);
						const Real area   = sign * u_area;

						J.describe(os);
						os << area << " == " << u_area << std::endl;
					}
				}

				os << "---------------------------------\n";
				os << "\n";
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

		bool have_common_side(const Integer e_index_1, const Integer e_index_2) const
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
			return n_common_nodes == ManifoldDim;
		}

		Integer common_side_num(const Integer e_index_1, const Integer e_index_2) const
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
			for(std::size_t i = 0; i < dual_graph_.size(); ++i) {
				os << "[" << i << "]:";
				for(std::size_t j = 0; j < dual_graph_[i].size(); ++j) {
					os << " " << dual_graph_[i][j];
				}
				os << "\n";
			}
		}

		

		Integer n_boundary_sides() const
		{
			assert( !dual_graph_.empty() && "requires that build_dual_graph is called first");

			Integer ret = 0;
			for(Integer i = 0; i < n_elements(); ++i) {
				if(!active_[i]) continue;

				const auto &e = elem(i);
				const auto &e_adj = dual_graph_[i];
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
				const auto &e_adj = dual_graph_[i];
				for(Integer k = 0; k < e_adj.size(); ++k) {
					const Integer j = e_adj[k];
					if(j == INVALID_INDEX) continue;
					e.side(k, side);

					const auto &other = elem(j);
					const auto &other_adj = dual_graph_[j];


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
							

							// assert(side.nodes[q] == other_side.nodes[other_q]);
							// return false;
							break;
						}
					}
				}
			}

			return true;
		}

		void uniform_refinement(const Integer n_levels = 1)
		{
			for(Integer l = 0; l < n_levels; ++l) {
				auto ne = n_elements();

				for(Integer i = 0; i < ne; ++i) {
					if(active_[i]) {
						red_refine_element(i);
					}
				}
			}

			// update_dual_graph();
		}

		void update_dual_graph(const bool force = false)
		{
			const Integer n_nodes    = this->n_nodes();
			const Integer n_elements = this->n_elements();
			Integer el_index_size = 0;

			if(dual_graph_.size() == n_elements && !force) {
				return;
			}

			std::vector< std::vector< Integer> > node_2_element(n_nodes);
			dual_graph_.resize(n_elements);

			for(Integer i = 0; i < n_elements; ++i) {
				if(!active_[i]) continue;

				const auto &e    = elem(i);
				const Integer nn = ManifoldDim + 1;

				std::fill(std::begin(dual_graph_[i]), std::end(dual_graph_[i]), INVALID_INDEX);

				for(Integer k = 0; k < nn; ++k) {
					node_2_element[e.nodes[k]].push_back(i);
				}
			}

			std::vector<Integer> el_offset(n_elements, 0);

			for(Integer i = 0; i < n_nodes; ++i) {
				const auto &elements = node_2_element[i];

				for(std::size_t e_i = 0; e_i < elements.size(); ++e_i) {
					const Integer e = elements[e_i];

					for(std::size_t e_i_adj = 0; e_i_adj < elements.size(); ++e_i_adj) {
						if(e_i == e_i_adj) continue;

						const Integer e_adj = elements[e_i_adj];

						bool must_add = true;
						for(Integer k = 0; k < el_offset[e]; ++k) {
							if(e_adj == dual_graph_[e][k]) {
								must_add = false;
							}
						}

						if(must_add && have_common_side(e, e_adj)) {
							assert(el_offset[e] < ManifoldDim + 1);
							dual_graph_[e][el_offset[e]] = e_adj;
							++el_offset[e];
						}
					}
				}
			}

			for(Integer i = 0; i < n_elements; ++i) {
				if(!active_[i]) continue;

				const std::array<Integer, ManifoldDim+1> neighs = dual_graph_[i];

				std::fill(std::begin(dual_graph_[i]), std::end(dual_graph_[i]), INVALID_INDEX);

				for(Integer j = 0; j < neighs.size(); ++j) {
					if(neighs[j] == INVALID_INDEX) break;

					const auto s = common_side_num(i, neighs[j]);
					assert(s != INVALID_INDEX);

					if(dual_graph_[i][s] != INVALID_INDEX) {
						std::cerr << "bad side numbering or creation" << std::endl;
						assert(dual_graph_[i][s] == INVALID_INDEX);
					}

					dual_graph_[i][s] = neighs[j];
				}
			}
		}

		void build_dual_graph()
		{
			// dual_graph_.clear();
			update_dual_graph();
		}

		inline Integer root(const Integer id) const
		{	
			Integer current_id = id;
			while(elem(current_id).parent_id != INVALID_INDEX) {
				current_id = elem(current_id).parent_id;
			}

			return current_id;
		}

	private:
		std::vector< Simplex<Dim, ManifoldDim> > elements_;
		std::vector< Vector<Real, Dim> > points_;
		
		//refinement
		std::vector<Integer> refinement_flag_;
		std::vector<std::array<Integer, ManifoldDim+1>> dual_graph_;
		std::vector<bool> active_;
		std::vector< std::shared_ptr<SimplexInterpolator<ManifoldDim>> > interp_;
		EdgeNodeMap edge_node_map_;
		EdgeElementMap edge_element_map_;

	};

	enum PlotFun
	{
		PLOT_ROOT = 0,
		PLOT_FLAG = 1,
		PLOT_ID = 2,
		PLOT_PARENT = 3,
		PLOT_PARENT_FLAG = 4
	};

	inline void flag_color(const Integer flag, Real &r, Real &g, Real &b)
	{
		r = 0.; g = 0.; b = 0;
		switch(flag)
		{
			case RED: {
				r = 1.;
				break;
			}
			case GREEN_1:
			case GREEN_2:
			case GREEN_3:
			{
				g = 1.;
				break;
			}
			default: {
				r = 1.; 
				g = 1.;
				b = 1.;
				break;
			}
		}
	}


	template<Integer Dim>
	bool write_mesh(
		const std::string &path,
		const Mesh<Dim, 2> &mesh,
		const Real scale_factor = 1.,
		const PlotFun plot_fun = PLOT_ROOT)
	{

		moonolith::Mesh m;
		m.dim = Dim;

		m.points.resize(mesh.n_nodes() * Dim);
		m.el_index.resize(mesh.n_active_elements() * 3);

		for(Integer i = 0; i < mesh.n_nodes(); ++i) {
			for(Integer d = 0; d < Dim; ++d) {
				m.points[i * Dim + d] = mesh.point(i)(d) * scale_factor;
			}
		}

		m.elem_type.resize(mesh.n_active_elements());
		std::fill(m.elem_type.begin(), m.elem_type.end(), moonolith::ElemType::TRI3);
		m.uniform_elem_type = moonolith::ElemType::TRI3;
		m.has_uniform_elem_type = true;

		m.el_ptr.resize(m.elem_type.size() + 1);

		m.el_ptr[0] = 0;
		for(std::size_t i = 1; i < m.el_ptr.size(); ++i) {
			m.el_ptr[i] = m.el_ptr[i - 1] + 3;
		}

		Integer k = 0;
		for(std::size_t i = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				for(Integer j = 0; j < 3; ++j) {
					m.el_index[k * 3 + j] = mesh.elem(i).nodes[j];
				}

				k++;
			}
		}

		// moonolith::SVGCanvas canvas;
		moonolith::EPSCanvas canvas;
		canvas.set_line_width(0.1/mesh.n_active_elements());

		std::vector<Real> f, hsv;

		if(plot_fun == PLOT_FLAG || plot_fun == PLOT_PARENT_FLAG) {
			f.resize(mesh.n_active_elements() * 3);
		} else {
			f.resize(mesh.n_active_elements());
		}
		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				switch(plot_fun) {
					case PLOT_ROOT: {
						f[k++] = mesh.root(i);
						break;
					}

					case PLOT_FLAG: {
						Real r = 0., g = 0., b = 0;
						flag_color(mesh.refinement_flag(i), r, g, b);
						f[k * 3]     = r;
						f[k * 3 + 1] = g;
						f[k * 3 + 2] = b; 

						k++;
						break;
					}

					case PLOT_PARENT_FLAG: {
						Real r = 0., g = 0., b = 0;
						if(mesh.elem(i).parent_id != INVALID_INDEX) {
							flag_color(mesh.refinement_flag(mesh.elem(i).parent_id), r, g, b);
						} else {
							flag_color(NONE, r, g, b);
						}

						f[k * 3]     = r;
						f[k * 3 + 1] = g;
						f[k * 3 + 2] = b; 

						k++;
						break;
					}

					case PLOT_PARENT: {
						if(mesh.elem(i).parent_id != INVALID_INDEX) {
							f[k++] = mesh.elem(i).parent_id;
						} else {
							f[k++] = mesh.elem(i).id;
						}
						break;
					}

					default:
					{
						f[k++] = mesh.elem(i).id;
						break;
					}
				}
			}
		}

		if(plot_fun == PLOT_FLAG || plot_fun == PLOT_PARENT_FLAG) {
			hsv = f;
		} else {
			moonolith::func_to_hsv(f, hsv);
		}

		canvas.fill_mesh(m, hsv);
		canvas.set_color(0,0,0);
		canvas.stroke_mesh(m);

		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				auto b = barycenter(mesh.elem(i), mesh.points());
				canvas.draw_text(b(0)*scale_factor, b(1)*scale_factor, 4./mesh.n_active_elements(), "Arial", std::to_string(i), true);
			}
		}

		return canvas.write(path);
	}

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

	bool mesh_hyper_cube(
		const std::array<Integer, 4> &dims,
		const Vector<Real, 4> &lobo,
		const Vector<Real, 4> &upbo,
		const Mesh<4, 4> &mesh)
	{

		return false;
	}
}

#endif //MARS_MESH_HPP
