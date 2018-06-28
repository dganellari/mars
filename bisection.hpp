#ifndef MARS_BISECTION_HPP
#define MARS_BISECTION_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;


	

	template<Integer Dim, Integer ManifoldDim>
	class Bisection {
	public:
		Bisection(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh), verbose(false), limit_to_1_level_diff(true)
		{}

		void set_limit_to_1_level_diff(const bool val)
		{
			limit_to_1_level_diff = val;
		}
		
		Integer add_elem(const Simplex<Dim, ManifoldDim> &e)
		{
			flags.push_back(NONE);
			Integer id = mesh.add_elem(e);
			mesh.repair_element(id);
			level.push_back((mesh.elem(id).parent_id == INVALID_INDEX)? 0 : (level[mesh.elem(id).parent_id] + 1)); 
			edge_element_map_.update(mesh.elem(id));
			return id;
		}

		void other_nodes(
			const std::array<Integer, ManifoldDim+1> &nodes,
			const Integer v1, 
			const Integer v2,
			std::array<Integer, ManifoldDim-1>  &opposite_nodes) const
		{
			Integer i = 0;
			for(auto n : nodes) {
				if(n != v1 && n != v2) {
					opposite_nodes[i++] = n;
				}
			}
		}

		void deactivated_element_dual_graph_update(const Integer id)
		{
			const auto &e = mesh.elem(id);

			if(e.children.empty()) {
				std::cerr << "calling element_deactivated on childless element " << id << std::endl;
				return;
			}

			std::fill(
				mesh.dual_graph().adj(id).begin(), 
				mesh.dual_graph().adj(id).end(), 
				INVALID_INDEX);

			std::map<Integer, std::vector<Integer> > local_node_2_element;

			for(auto c : mesh.elem(id).children) {
				for(auto n : mesh.elem(c).nodes) {
					local_node_2_element[n].push_back(c);
				}
			}

			for(auto a : mesh.dual_graph().adj(id)) {
				if(a == INVALID_INDEX) continue;

				for(auto n : mesh.elem(a).nodes) {
					local_node_2_element[n].push_back(a);
				}

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
				std::cerr << "element " << id  << " not updated " << std::endl;
				assert(updated);
			}
		}

		void bisect_element(
			const Integer element_id,
			const Edge &edge)
		{
			mesh.elem(element_id).children.clear();
			mesh.set_active(element_id, false);

			Simplex<Dim, ManifoldDim> s;
			s.parent_id = element_id;

			const Integer v1 = edge.nodes[0];
			const Integer v2 = edge.nodes[1];
			
			if(verbose) {
				std::cout << "bisect(" << v1 << ", " << v2 << ") for " << element_id << std::endl;
			}

			auto midpoint = edge_node_map_.get(v1, v2);

			if(midpoint == INVALID_INDEX) {
				midpoint = mesh.add_point(0.5 * (mesh.point(v1) + mesh.point(v2)));
				edge_node_map_.update(v1, v2, midpoint);
			}

			std::array<Integer, ManifoldDim-1> opposite_nodes;
			other_nodes(mesh.elem(element_id).nodes, v1, v2, opposite_nodes);

			for(Integer i = 0; i < ManifoldDim-1; ++i) {
				s.nodes[2+i] = opposite_nodes[i];
			}

			s.nodes[0] = v1;
			s.nodes[1] = midpoint;

			Integer new_id = add_elem(s);
			mesh.elem(element_id).children.push_back(new_id);

			s.nodes[0] = v2;
			s.nodes[1] = midpoint;

			new_id = add_elem(s);
			mesh.elem(element_id).children.push_back(new_id);
			return;
		}

		//one possible strategy for refining
		void longest_edge_ordering(const Simplex<Dim, ManifoldDim> &e, std::vector<Integer> &ordering) const
		{
			ordering.clear();
			
			std::vector< std::pair<Real, Integer> > len2edge;
			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);
				len2edge.emplace_back((mesh.point(v1) - mesh.point(v2)).norm(), i);
			}

			std::sort(len2edge.begin(), len2edge.end());

			ordering.reserve(n_edges(e));
			for(auto it = len2edge.rbegin(); it != len2edge.rend(); ++it) {
				ordering.push_back(it->second);
			}
		}

		void longest_edge_ordering_with_tol(const Simplex<Dim, ManifoldDim> &e, const Edge &edge, std::vector<Integer> &ordering)
		{
			ordering.clear();
			
			Real best_dist = 0.;
			Real edge_dist = 0.;
			Integer edge_index = INVALID_INDEX;

			std::vector< std::pair<Real, Integer> > len2edge;
			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);
				Real dist = (mesh.point(v1) - mesh.point(v2)).norm();

				best_dist = std::max(dist, best_dist);
				
				if(Edge(v1, v2) == edge) {
					edge_dist = dist;
					edge_index = i;
				}

				len2edge.emplace_back(dist, i);
			}

			std::sort(len2edge.begin(), len2edge.end());
			ordering.reserve(n_edges(e));

			bool use_edge = false;
			if(edge_dist/best_dist >= (0.99)) {
				use_edge = true;
				ordering.push_back(edge_index);
			}

			for(auto it = len2edge.rbegin(); it != len2edge.rend(); ++it) {
				if(use_edge && it->second == edge_index) {
					continue;
				} else {
					ordering.push_back(it->second);
				}
			}
		}

		bool refine_element_recursive(
			const Integer element_id,
			const Edge &edge,
			const Integer max_level)
		{
			assert(has_edge(mesh.elem(element_id), edge.nodes[0], edge.nodes[1]));

			std::vector<Integer> edges;
			// longest_edge_ordering(mesh.elem(element_id), edges);
			longest_edge_ordering_with_tol(mesh.elem(element_id), edge, edges);

			Edge new_edge;
			mesh.elem(element_id).edge(edges.front(), new_edge.nodes[0], new_edge.nodes[1]);
			new_edge.fix_ordering();

			if(edge == new_edge) {
				bisect_element(element_id, edge);
			} else if(limit_to_1_level_diff && level[element_id] >= max_level) {
				bisect_element(element_id, edge);
			} else {
				refine_edge(new_edge);
				assert(!mesh.is_active(element_id));
			}

			return false;
		}

		void refine_element(const Integer element_id)
		{
			std::vector<Integer> edges;
			longest_edge_ordering(mesh.elem(element_id), edges);

			Edge edge;
			mesh.elem(element_id).edge(edges.front(), edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();
			refine_edge(edge);
		}

		void refine_edge(const Edge &edge)
		{
			bool complete = false;

			auto incidents = edge_element_map_.elements(edge);

			while(!complete) {
				complete = true;
				
				for(auto i : incidents) {
					if(!mesh.is_active(i)) continue;

					assert(has_edge(mesh.elem(i), edge.nodes[0], edge.nodes[1]));
					refine_element_recursive(i, edge, level[i]);
					complete = false;
				}

				const auto &next_incidents = edge_element_map_.elements(edge);
				
				if(next_incidents.size() > incidents.size()) {
					incidents = next_incidents;
					complete = false;
				}
			}
		}

		void split_neighs(
			const std::vector<Integer> &neighs,
			const Integer split_level,
			std::vector<Integer> &before,
			std::vector<Integer> &after) const
		{
			before.clear();
			after.clear();
			for(auto n : neighs) {
				if(!mesh.is_active(n)) continue;

				if(level[n] >= split_level) {
					after.push_back(n);
				} else {
					before.push_back(n);
				}
			}
		}

		void uniform_refine(const Integer n_levels)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			for(Integer l = 0; l < n_levels; ++l) {
				const Integer n_elements = mesh.n_elements();
				for(Integer i = 0; i < n_elements; ++i) {
					if(!mesh.is_valid(i)) {
						std::cerr << "tried to refine non-existing element " << i << std::endl;
						continue;
					}

					if(!mesh.is_active(i)) continue;
					refine_element(i);
				}

			}

			mesh.update_dual_graph();
			mesh.tags() = level;
		}

		void refine(const std::vector<Integer> &elements)
		{
			std::vector<Integer> sorted_elements = elements;
			std::sort(sorted_elements.begin(), sorted_elements.end());

			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			for(auto i : sorted_elements) {
				if(!mesh.is_valid(i)) {
					std::cerr << "tried to refine non-existing element " << i << std::endl;
					continue;
				}


				if(!mesh.is_active(i)) {
					std::cerr << "tried to refine inactive element " << i << std::endl;
					continue;
				}

				refine_element(i);
			}

			mesh.update_dual_graph();
			mesh.tags() = level;
		}

		void refine_(const std::vector<Integer> &elements)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			std::vector<Integer> defferred;
			std::vector<std::pair<Integer, Edge>> to_refine;

			for(auto i : elements) {
				if(!mesh.is_valid(i)) {
					std::cerr << "tried to refine non-existing element " << i << std::endl;
					continue;
				}


				if(!mesh.is_active(i)) {
					std::cerr << "tried to refine inactive element " << i << std::endl;
					continue;
				}

				if(flags[i] == BISECTION) {
					continue;
				}


				Integer edge_to_split = INVALID_INDEX;
				

				std::vector<Integer> edges;
				longest_edge_ordering(mesh.elem(i), edges);
				
				for(auto edge_num : edges) {
				// for(Integer edge_num = 0; edge_num < n_edges(mesh.elem(i)); ++edge_num) {
					
					Edge edge;
					mesh.elem(i).edge(edge_num, edge.nodes[0], edge.nodes[1]);
					edge.fix_ordering();

					auto neighs = edge_element_map_.elements(edge.nodes[0], edge.nodes[1]);
					
					bool skip_edge = false;
					for(auto n : neighs) {
						if(!mesh.is_active(n)) continue;

						if(flags[n] == BISECTION) {
							skip_edge = true;
							break;
						}
					}

					if(!skip_edge)  {
						edge_to_split = edge_num;
						break;
					}
				}

				if(edge_to_split == INVALID_INDEX) {
					defferred.push_back(i);
				} else {
					Edge edge;
					mesh.elem(i).edge(edge_to_split, edge.nodes[0], edge.nodes[1]);
					edge.fix_ordering();

					auto neighs = edge_element_map_.elements(edge.nodes[0], edge.nodes[1]);

					for(auto n : neighs) {
						if(!mesh.is_active(n)) continue;
						flags[n] = BISECTION;
					}

					to_refine.emplace_back(i, edge);
				}

			}

			for(auto &tr : to_refine) {
				bisect_element(tr.first, tr.second);

				auto neighs = edge_element_map_.elements(tr.second.nodes[0], tr.second.nodes[1]);

				for(auto n : neighs) {
					if(!mesh.is_active(n)) continue;
					bisect_element(n, tr.second);
				}
			}

			if(!defferred.empty()) {
				return refine(defferred);
			}


			mesh.update_dual_graph();
			// mesh.tags() = flags;
			mesh.tags() = level;
		}

		const EdgeNodeMap &edge_node_map() const
		{
			return edge_node_map_;
		}

		const EdgeElementMap &edge_element_map() const
		{
			return edge_element_map_;
		}
		
		void set_verbose(const bool val)
		{
			verbose = val;
		}

	private:
		Mesh<Dim, ManifoldDim> &mesh;
		std::vector<Integer> flags;
		std::vector<Integer> level;
		std::vector<std::array<Integer, ManifoldDim+1> > side_flags;
		EdgeNodeMap edge_node_map_;
		EdgeElementMap edge_element_map_;
		bool verbose;
		bool limit_to_1_level_diff;
	};
}

#endif //MARS_BISECTION_HPP
