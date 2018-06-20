#ifndef MARS_MESH_HPP
#define MARS_MESH_HPP

#include "simplex.hpp"
#include "edge_element_map.hpp"
#include "edge_node_map.hpp"

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
		CHILD_OF_GREEN = 5
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

		inline Integer refinement_flag(const Integer &element_id)
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

		void red_green_refinement(const std::vector<Integer> &elements_to_refine)
		{
			build_dual_graph();
			edge_element_map_.build(*this);

			for(auto &e : elements_to_refine) {
				set_refinement_flag(e, RED);
			}

			std::vector<Integer> red_green_elements = elements_to_refine;
			for(auto &e : elements_to_refine) {
				const auto &adj = dual_graph_[e];
				for(auto a : adj) {
					if(a == INVALID_INDEX) continue;

					switch(refinement_flag(a))
					{
						case RED:
						{
							continue;
						}
						case NONE:
						{
							set_refinement_flag(a, GREEN_1);
							red_green_elements.push_back(a);
							break;
						}
						case GREEN_1:
						case GREEN_2:
						case GREEN_3:
						{
							assert(false && "implement me");
							break;
						}
					}
				}
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

		void update_dual_graph()
		{
			const Integer n_nodes    = this->n_nodes();
			const Integer n_elements = this->n_elements();
			Integer el_index_size = 0;

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
			dual_graph_.clear();
			update_dual_graph();
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
