#ifndef MARS_RED_GREEN_REFINEMENT_HPP
#define MARS_RED_GREEN_REFINEMENT_HPP

#include <vector>
#include <array>
#include <algorithm>  

namespace mars {
	
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	enum RefinementFlag {
		NONE = 0,
		RED = 1,
		GREEN_1 = 2,
		GREEN_2 = 3,
		GREEN_3 = 4,
		CHILD_OF_GREEN = 5,
		PARENT_PROMOTED_TO_RED = 6
	};

	template<Integer Dim, Integer ManifoldDim>
	class RedGreenRefinement {
	public:
		class Refinement {
		public:
			class Trigger {
			public:
				//red|green|...
				Integer flag;

				//element id
				Integer element;

				//self=Manifold|side=Manifold-1|sub-side=Manifold-2|....|edge=1
				Integer type;

				//face index|edge index
				Integer index;
			};

			std::vector<Trigger> triggers;
		};

		RedGreenRefinement(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh)
		{}

		void update_children(const Integer id)
		{
			for(auto c : mesh.elem(id).children) {
				update_element(c);
			}
		}

		void update_element(const Integer id)
		{
			const auto &e = mesh.elem(id);
			const auto p_id = e.parent_id;
			const auto &parent_e = mesh.elem(p_id);
			const auto &adj = mesh.dual_graph().adj(p_id);

			std::map<Integer, std::vector<Integer> > local_node_2_element;

			for(auto c : mesh.elem(p_id).children) {
				for(auto n : mesh.elem(c).nodes) {
					local_node_2_element[n].push_back(c);
				}
			}

			for(auto a : adj) {
				if(a == INVALID_INDEX) continue;

				for(auto c : mesh.elem(a).children) {
					for(auto n : mesh.elem(c).nodes) {
						local_node_2_element[n].push_back(c);
					}
				}
			}

			bool updated = false;

			for(Integer i = 0; i < ManifoldDim + 1; ++i) {
				auto it = local_node_2_element.find(e.nodes[i]);
				if(it == local_node_2_element.end()) continue;

				for(auto other : it->second) {
					if(id == other) continue;

					if(mesh.have_common_side(id, other)) {
						updated = true;
						auto &e_adj = mesh.dual_graph().safe_adj(id);
						e_adj[mesh.common_side_num(id, other)] = other;

						auto &other_adj = mesh.dual_graph().safe_adj(other);
						other_adj[mesh.common_side_num(other, id)] = id;
					}
				}
			}

			if(!updated) {
				std::cerr << "element " << id << " with parent " << p_id << " not updated " << std::endl;
				assert(updated);
			}
		}

		void promote_to_red(const Integer green_parent)
		{
			set_refinement_flag(green_parent, RED);
			mesh.deactivate_children(green_parent);
			mesh.set_active(green_parent, true);
		}
		
		void red_refine(const std::vector<Integer> &elements_to_refine)
		{
			update_dual_graph();
			assert(mesh.n_elements() == mesh.dual_graph().size());

			if(elements_to_refine.empty()) return;

			if(refinement_flag_.empty()) {
				refinement_flag_.resize(mesh.n_elements(), NONE);
			}

			std::vector<Integer> promoted;
			std::vector<Integer> red_elements;
			for(auto &e : elements_to_refine) {
				if(!mesh.is_active(e)) {
					std::cerr << e << " cannot refine inactive" << std::endl;
					continue;
				}

				auto parent = mesh.elem(e).parent_id;
				
				if(is_green(parent)) {
					promote_to_red(parent);
					promoted.push_back(parent);
				} else {
					set_refinement_flag(e, RED);
					red_elements.push_back(e);
				}
			}

			for(auto &e : red_elements) {
				const auto &adj = mesh.dual_graph().adj(e);
				for(auto a : adj) {
					if(a == INVALID_INDEX || !mesh.is_active(a) || mesh.elem(a).parent_id == INVALID_INDEX) continue;
					auto parent = mesh.elem(a).parent_id;
					
					if(is_green(parent)) {
						promote_to_red(parent);
						promoted.push_back(parent);
					}
				}
			}

			for(auto &e : promoted) {
				red_refine_element(e);
				update_element(e);
			}

			for(auto &e : red_elements) {
				red_refine_element(e);
			}

			//update dual_graph for children
			for(auto e : promoted) {
				update_children(e);
			}

			for(auto e : red_elements) {
				update_children(e);
			}

			// assert(mesh.n_elements() == mesh.dual_graph().size());
			mesh.tags() = refinement_flag_;
		}

		void green_refine()
		{
			std::vector<Integer> green_elements, promoted_elements;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(mesh.is_active(i) && !is_red(i)) {
					const auto &adj = mesh.dual_graph().adj(i);

					Integer n_red_neighs = 0;
					for(auto a : adj) {
						if(a == INVALID_INDEX || mesh.is_active(a)) continue;
						if(is_red(a)) {
							n_red_neighs++;
						}
					}

					if(n_red_neighs == 0) continue;


					switch(n_red_neighs) {
						case 1: 
						{
							set_refinement_flag(i, GREEN_1); 
							green_elements.push_back(i);
							break;
						}
						case 2: 
						{
							set_refinement_flag(i, GREEN_2); 
							green_elements.push_back(i);
							break;
						}
						case 3:
						{
							set_refinement_flag(i, RED);
							promoted_elements.push_back(i);
							break;
						}
						default: 
						{
							assert(false);
							break;
						}
					}
				}
			}

			if(!promoted_elements.empty()) {
				std::cout << "handling promoted elements" << std::endl;
				red_refine(promoted_elements);
				green_refine();
				return;
			}

			for(auto e : green_elements) {
				green_refine_element(e);
			}

			update_dual_graph();

			mesh.tags() = refinement_flag_;
		}
		
		void propagate_flags(
			std::vector<Integer> &red_elements,
			std::vector<Integer> &green_elements)
		{
			std::vector<Integer> potential_green_elements;
			for(auto &e : red_elements) {
				const auto &adj = mesh.dual_graph().adj(e);
				for(auto a : adj) {
					if(a == INVALID_INDEX || !mesh.is_active(a)) continue;
					auto flag = refinement_flag(a);

					if(!mesh.elem(a).children.empty()) {
						std::cout << "promoting " << a << " to red" << std::endl;
						set_refinement_flag(a, RED);
						//deactivate children
						mesh.deactivate_children(a);
						red_elements.push_back(a); 
						continue;
					}

					const auto parent = mesh.elem(a).parent_id;
					if(is_green(parent)) {
						std::cout << "promoting " << a << " to red" << std::endl;
						set_refinement_flag(parent, RED);
						//deparentctivate children
						mesh.deactivate_children(parent);
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
							// if(mesh.dual_graph().n_adjacients(a) == 1) {
							// 	set_refinement_flag(a, RED);
							// 	red_elements.push_back(a);
							// } else {
								set_refinement_flag(a, GREEN_1);
								potential_green_elements.push_back(a);
							// }

							break;
						}
						case GREEN_1:
						{
							// if(mesh.dual_graph().n_adjacients(a) == 2) {
							// 	set_refinement_flag(a, RED);
							// 	red_elements.push_back(a);
							// } else {
								set_refinement_flag(a, GREEN_2);
								potential_green_elements.push_back(a);
							// }

							break;
						}
						case GREEN_2:
						{
							// if(mesh.dual_graph().n_adjacients(a) == 3) {
							// 	set_refinement_flag(a, RED);
							// 	red_elements.push_back(a);
							// } else {
								//GREEN_3?
								set_refinement_flag(a, RED);
								red_elements.push_back(a);
							// }

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
				if(refinement_flag(e) != RED && mesh.is_active(e)) {
					green_elements.push_back(e);
				}
			}
		}

		void refine(const std::vector<Integer> &elements_to_refine)
		{
			if(refinement_flag_.empty()) {
				refinement_flag_.resize(mesh.n_elements(), NONE);
			}

			std::vector<Integer> promoted_elements;
			std::vector<Integer> red_elements, green_elements;
			for(auto &e : elements_to_refine) {
				if(!mesh.is_active(e)) {
					std::cerr << e << " cannot refine inactive" << std::endl;
					continue;
				}

				auto parent = mesh.elem(e).parent_id;
				
				if(is_green(parent)) {
					std::cout << "promoting " << parent << " to red" << std::endl;
					//promote parent to red
					set_refinement_flag(parent, RED);
					//deactivate children
					mesh.deactivate_children(parent);
					mesh.set_active(parent, true);
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
			edge_element_map_.build(mesh);

			propagate_flags(red_elements, green_elements);
			for(auto e : red_elements) {
				red_refine_element(e);
			}

			for(auto e : green_elements) {
				green_refine_element(e);
			}

			mesh.tags() = refinement_flag_;
		}	

		inline void green_refine_element(const Integer element_id)
		{
			mesh.set_active(element_id, false);

			const auto &adj = mesh.dual_graph().adj(element_id);
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
			mesh.elem(element_id).children.clear();

			switch(n_red_neighs) {
				case 1: 
				{	
					mesh.elem(element_id).side(red_side_index[0], side_1);
					
					Integer n0 = side_1.nodes[0];
					Integer n1 = side_1.nodes[1];
					const Integer midpoint = edge_node_map_.get(n0, n1);
					const Integer opposite = mesh.elem(element_id).vertex_opposite_to_side(red_side_index[0]);

					child.nodes[0] = n0;
					child.nodes[1] = opposite;
					child.nodes[2] = midpoint;

					mesh.elem(element_id).children.push_back( add_elem(child) );

					child.nodes[0] = n1;
					child.nodes[1] = midpoint;
					child.nodes[2] = opposite;

					mesh.elem(element_id).children.push_back( add_elem(child) );
					break;
				}

				case 2:
				{
					mesh.elem(element_id).side(red_side_index[0], side_1);
					
					Integer n1_0 = side_1.nodes[0];
					Integer n1_1 = side_1.nodes[1];
					const Integer midpoint_1 = edge_node_map_.get(n1_0, n1_1);
					const Integer opposite_1 = mesh.elem(element_id).vertex_opposite_to_side(red_side_index[0]);

					mesh.elem(element_id).side(red_side_index[1], side_2);
					
					Integer n2_0 = side_2.nodes[0];
					Integer n2_1 = side_2.nodes[1];
					const Integer midpoint_2 = edge_node_map_.get(n2_0, n2_1);
					const Integer opposite_2 = mesh.elem(element_id).vertex_opposite_to_side(red_side_index[1]);


					child.nodes[0] = midpoint_1;
					child.nodes[1] = opposite_2;
					child.nodes[2] = opposite_1;

					mesh.elem(element_id).children.push_back( add_elem(child) ); 

					child.nodes[0] = midpoint_1;
					child.nodes[1] = midpoint_2;
					child.nodes[2] = n2_0;

					mesh.elem(element_id).children.push_back( add_elem(child) ); 

					child.nodes[0] = midpoint_1;
					child.nodes[1] = n2_1;
					child.nodes[2] = midpoint_2;

					mesh.elem(element_id).children.push_back( add_elem(child) ); 
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

		void update_dual_graph(const bool force = false)
		{
			mesh.update_dual_graph(force);
		}

		inline void red_refine_element(const Integer element_id)
		{
			static const Integer NSubs = NSubSimplices<ManifoldDim>::value;
			static_assert(NSubSimplices<ManifoldDim>::value > 0, "!");

			std::vector<Vector<Real, Dim>> parent_points;
			mesh.points(element_id, parent_points);

			std::array<Simplex<Dim, ManifoldDim>, NSubs> children;
			std::vector<Vector<Real, Dim>> children_points;
			auto interp = std::make_shared< SimplexInterpolator<ManifoldDim> >();

			Simplex<Dim, ManifoldDim> modified_e = mesh.elem(element_id);
			
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
                		const auto new_id = mesh.add_point(children_points[offset]);
                		edge_node_map_.update(
                			modified_e.nodes[i],
                			modified_e.nodes[j],
							new_id
						);

						point_ids[offset] = new_id;
						assert(new_id < mesh.n_nodes());
                	}
            	}
            }

			interp_[element_id] = interp;
			mesh.elem(element_id).children.clear();

			for(auto &c : children) {
				for(Integer i = 0; i < ManifoldDim + 1; ++i) {
					c.nodes[i] = point_ids[c.nodes[i]];
				}

				// std::sort(c.nodes.begin(), c.nodes.end()); add_elem(c);

				//4D hack
				const auto c_id = add_elem(c);
				mesh.repair_element(c_id, ManifoldDim != 4);
				mesh.elem(element_id).children.push_back(c_id);
			}

			assert(mesh.elem(element_id).children.size() == NSubs);
			mesh.set_active(element_id, false);
		}

		void uniformly_refine(const Integer n_levels = 1)
		{
			for(Integer l = 0; l < n_levels; ++l) {
				auto ne = mesh.n_elements();

				for(Integer i = 0; i < ne; ++i) {
					if(mesh.is_active(i)) {
						red_refine_element(i);
					}
				}
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

		inline bool is_red(const Integer id) const
		{
			return refinement_flag(id) == RED;
		}

		inline void set_refinement_flag(const Integer &element_id, const Integer flag)
		{
			refinement_flag_[element_id] = flag;
		}

		inline Integer refinement_flag(const Integer &element_id) const
		{
			return refinement_flag_[element_id];
		}

	private:

		inline Integer add_elem(const Simplex<Dim, ManifoldDim> &elem)
		{
			refinement_flag_.push_back(NONE);
			return mesh.add_elem(elem);
		}

		Mesh<Dim, ManifoldDim> &mesh;

		//refinement
		std::vector<Integer> refinement_flag_;
		std::vector< std::shared_ptr<SimplexInterpolator<ManifoldDim>> > interp_;
		EdgeNodeMap edge_node_map_;
		EdgeElementMap edge_element_map_;
	};
}

#endif //MARS_RED_GREEN_REFINEMENT_HPP
