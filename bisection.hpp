#ifndef MARS_BISECTION_HPP
#define MARS_BISECTION_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;


	

	template<Integer Dim, Integer ManifoldDim>
	class Bisection {
	public:
		Bisection(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh)
		{}
		
		Integer add_elem(const Simplex<Dim, ManifoldDim> &e)
		{
			flags.push_back(NONE);
			Integer id = mesh.add_elem(e);
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
			mesh.repair_element(new_id);

			s.nodes[0] = v2;
			s.nodes[1] = midpoint;

			new_id = add_elem(s);
			mesh.elem(element_id).children.push_back(new_id);
			mesh.repair_element(new_id);
			return;
		}


		//one possible strategy for refining
		void longest_edge_ordering(const Simplex<Dim, ManifoldDim> &e, std::vector<Integer> &ordering) const
		{
			ordering.clear();

			Integer v1, v2;
			std::vector< std::pair<Real, Integer> > len2edge;
			for(Integer i = 0; i < n_edges(e); ++i) {
				e.edge(i, v1, v2);
				len2edge.emplace_back((mesh.point(v1) - mesh.point(v2)).norm(), i);
			}

			std::sort(len2edge.begin(), len2edge.end());

			ordering.reserve(n_edges(e));
			for(auto &l2e : len2edge) {
				ordering.push_back(l2e.second);
			}
		}

		void refine(const std::vector<Integer> &elements)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
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

					auto &neighs = edge_element_map_.elements(edge.nodes[0], edge.nodes[1]);
					
					bool skip_edge = false;
					for(auto n : neighs) {
						if(!mesh.is_active(n)) continue;

						if(flags[n] == BISECTION) {
							skip_edge = true;
							break;
						}
					}

					if(!skip_edge)  {
						//FIXME USE QUALITY CRITERION
						edge_to_split = edge_num;
					}
				}

				if(edge_to_split == INVALID_INDEX) {
					defferred.push_back(i);
				} else {
					Edge edge;
					mesh.elem(i).edge(edge_to_split, edge.nodes[0], edge.nodes[1]);
					auto &neighs = edge_element_map_.elements(edge.nodes[0], edge.nodes[1]);

					for(auto n : neighs) {
						if(!mesh.is_active(n)) continue;
						flags[n] = BISECTION;
					}

					to_refine.emplace_back(i, edge);
				}

			}

			for(auto &tr : to_refine) {
				bisect_element(tr.first, tr.second);

				auto &neighs = edge_element_map_.elements(tr.second.nodes[0], tr.second.nodes[1]);

				for(auto n : neighs) {
					if(!mesh.is_active(n)) continue;
					bisect_element(n, tr.second);
				}
			}

			if(!defferred.empty()) {
				return refine(defferred);
			}


			mesh.update_dual_graph();
			mesh.tags() = flags;
		}

		const EdgeNodeMap &edge_node_map() const
		{
			return edge_node_map_;
		}
		
	private:
		Mesh<Dim, ManifoldDim> &mesh;
		std::vector<Integer> flags;
		std::vector<std::array<Integer, ManifoldDim+1> > side_flags;
		EdgeNodeMap edge_node_map_;
		EdgeElementMap edge_element_map_;
	};
}

#endif //MARS_BISECTION_HPP
