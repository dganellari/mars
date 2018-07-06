#ifndef MARS_BISECTION_HPP
#define MARS_BISECTION_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer Dim, Integer ManifoldDim>
	class EdgeSelect {
	public:
		virtual ~EdgeSelect() {}
		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const
		{
			//first edge selected
			return 0;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const
		{
			//first edge selected
			return 0;
		}

		virtual bool is_recursive() const
		{
			return false;
		}

		virtual std::string name() const = 0;
	};

	template<Integer Dim, Integer ManifoldDim>
	class LongestEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		LongestEdgeSelect(const bool recursive = true, const bool use_tollerance = true)
		: recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			const auto &e = mesh.elem(element_id);

			Integer edge_num = 0;
			Real len = 0;

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				Real len_i = (mesh.point(v1) - mesh.point(v2)).norm();

				if(len_i > len) {
					len = len_i;
					edge_num = i;
				}
			}

			return edge_num;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{
			if(!use_tollerance_) {
				return select(mesh, element_id);
			}

			const auto &e = mesh.elem(element_id);

			Real len = 0.;
			Real neigh_len = 0.;
			Integer edge_num = 0;

			Integer neigh_edge_num = INVALID_INDEX;

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				const Real dist = (mesh.point(v1) - mesh.point(v2)).norm();

				if(dist > len) {
					len = dist;
					edge_num = i;
				}
				
				if(Edge(v1, v2) == neighbor_edge) {
					neigh_len 	   = dist;
					neigh_edge_num = i;
				}
			}

			if(neigh_len/len >= (0.99)) {
				edge_num = neigh_edge_num;
			}

			return edge_num;
		}

		void set_recursive(const bool recursive)
		{
			recursive_ = recursive;
		}

		bool is_recursive() const override
		{
			return recursive_;
		}

		virtual std::string name() const override
		{
			return "LongestEdge";
		}

	private:
		bool recursive_;
		bool use_tollerance_;
	};


	template<Integer Dim, Integer ManifoldDim>
	class NewestVertexEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		NewestVertexEdgeSelect(const bool recursive = true)
		: recursive_(recursive)
		{}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			const auto &e = mesh.elem(element_id);
			Integer edge_num = 0;

			auto node_ids = e.nodes;
			std::sort(node_ids.begin(), node_ids.end());
			Edge edge(node_ids[0], node_ids[1]);

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);
				
				if(Edge(v1, v2) == edge) {
					return i;
				}
			}

			assert(false);
			return edge_num;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{	
			(void) neighbor_edge;
			return select(mesh, element_id);
		}

		void set_recursive(const bool recursive)
		{
			recursive_ = recursive;
		}

		bool is_recursive() const override
		{
			return recursive_;
		}

		virtual std::string name() const override
		{
			return "NewestVertex";
		}

	private:
		bool recursive_;
	};


	template<Integer Dim, Integer ManifoldDim>
	class NewestVertexAndLongestEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		NewestVertexAndLongestEdgeSelect(const bool recursive = true)
		: recursive_(recursive)
		{}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			const auto &e = mesh.elem(element_id);
			Integer edge_num = 0;
			Real len = 0;

			auto node_ids = e.nodes;
			std::sort(node_ids.begin(), node_ids.end());
			const Integer vertex_to_skip = node_ids.back();

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				if(v1 == vertex_to_skip || v2 == vertex_to_skip) continue;

				Real len_i = (mesh.point(v1) - mesh.point(v2)).norm();

				if(len_i > len) {
					len = len_i;
					edge_num = i;
				}
			}

			return edge_num;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{	
			const auto &e = mesh.elem(element_id);

			Real len = 0.;
			Real neigh_len = 0.;
			Integer edge_num = 0;

			Integer neigh_edge_num = INVALID_INDEX;

			auto node_ids = e.nodes;
			std::sort(node_ids.begin(), node_ids.end());
			const Integer vertex_to_skip = node_ids.back();

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				if(v1 == vertex_to_skip || v2 == vertex_to_skip) continue;

				const Real dist = (mesh.point(v1) - mesh.point(v2)).norm();

				if(dist > len) {
					len = dist;
					edge_num = i;
				}
				
				if(Edge(v1, v2) == neighbor_edge) {
					neigh_len 	   = dist;
					neigh_edge_num = i;
				}
			}

			if(neigh_len/len >= (0.99)) {
				edge_num = neigh_edge_num;
			}

			return edge_num;
		}

		void set_recursive(const bool recursive)
		{
			recursive_ = recursive;
		}

		bool is_recursive() const override
		{
			return recursive_;
		}

		virtual std::string name() const override
		{
			return "NewestVertexAndLongestEdge";
		}

	private:
		bool recursive_;
	};
	

	template<Integer Dim, Integer ManifoldDim>
	class Bisection {
	public:
		Bisection(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh), edge_select_(std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>()),
		  verbose(false)
		{}

		void set_edge_select(const std::shared_ptr<EdgeSelect<Dim, ManifoldDim>> &edge_select)
		{
			edge_select_ = edge_select;
		}

		Integer add_elem(const Simplex<Dim, ManifoldDim> &e)
		{
			flags.push_back(NONE);
			Integer id = mesh.add_elem(e);
			mesh.repair_element(id);
			mesh.update_side_flags_from_parent(id);
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

		// //one possible strategy for refining
		// void longest_edge_ordering(const Simplex<Dim, ManifoldDim> &e, std::vector<Integer> &ordering) const
		// {
		// 	ordering.clear();
			
		// 	std::vector< std::pair<Real, Integer> > len2edge;
		// 	for(Integer i = 0; i < n_edges(e); ++i) {
		// 		Integer v1, v2;
		// 		e.edge(i, v1, v2);
		// 		len2edge.emplace_back((mesh.point(v1) - mesh.point(v2)).norm(), i);
		// 	}

		// 	std::sort(len2edge.begin(), len2edge.end());

		// 	ordering.reserve(n_edges(e));
		// 	for(auto it = len2edge.rbegin(); it != len2edge.rend(); ++it) {
		// 		ordering.push_back(it->second);
		// 	}
		// }

		// void longest_edge_ordering_with_tol(const Simplex<Dim, ManifoldDim> &e, const Edge &edge, std::vector<Integer> &ordering)
		// {
		// 	ordering.clear();
			
		// 	Real best_dist = 0.;
		// 	Real edge_dist = 0.;
		// 	Integer edge_index = INVALID_INDEX;

		// 	std::vector< std::pair<Real, Integer> > len2edge;
		// 	for(Integer i = 0; i < n_edges(e); ++i) {
		// 		Integer v1, v2;
		// 		e.edge(i, v1, v2);
		// 		Real dist = (mesh.point(v1) - mesh.point(v2)).norm();

		// 		best_dist = std::max(dist, best_dist);
				
		// 		if(Edge(v1, v2) == edge) {
		// 			edge_dist = dist;
		// 			edge_index = i;
		// 		}

		// 		len2edge.emplace_back(dist, i);
		// 	}

		// 	std::sort(len2edge.begin(), len2edge.end());
		// 	ordering.reserve(n_edges(e));

		// 	bool use_edge = false;
		// 	if(edge_dist/best_dist >= (0.99)) {
		// 		use_edge = true;
		// 		ordering.push_back(edge_index);
		// 	}

		// 	for(auto it = len2edge.rbegin(); it != len2edge.rend(); ++it) {
		// 		if(use_edge && it->second == edge_index) {
		// 			continue;
		// 		} else {
		// 			ordering.push_back(it->second);
		// 		}
		// 	}
		// }

		bool refine_element_recursive(
			const Integer element_id,
			const Edge &edge,
			const Integer max_level)
		{
			assert(has_edge(mesh.elem(element_id), edge.nodes[0], edge.nodes[1]));

			// Edge mandatory_edge = mandatory_edge_splitting(element_id);
			// if(mandatory_edge.is_valid()) {
			// 	bisect_element(element_id, mandatory_edge);
			// 	return false;
			// }

			const Integer edge_num = edge_select_->select(mesh, edge, element_id);

			Edge new_edge;
			mesh.elem(element_id).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
			new_edge.fix_ordering();

			if(edge == new_edge) {
				bisect_element(element_id, edge);
			} else if(!edge_select_->is_recursive()) {
				bisect_element(element_id, edge);
			} else {
				refine_edge(new_edge);
				assert(!mesh.is_active(element_id));
			}

			return false;
		}

		void refine_element(const Integer element_id)
		{
			// Edge mandatory_edge = mandatory_edge_splitting(element_id);
			// if(mandatory_edge.is_valid()) {
			// 	refine_edge(mandatory_edge);
			// 	return false;
			// }

			const Integer edge_num = edge_select_->select(mesh, element_id);

			Edge edge;
			mesh.elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();
			refine_edge(edge);
		}

		void refine_edge(const Edge &edge)
		{
			bool complete = false;
			bool has_refined = false;

			auto incidents = edge_element_map_.elements(edge);

			while(!complete) {
				complete = true;
				
				for(auto i : incidents) {
					if(!mesh.is_active(i)) continue;

					has_refined = true;

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

			if(has_refined) {
				bisected_edges_.push_back(edge);
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
					// std::cerr << "tried to refine inactive element " << i << std::endl;
					continue;
				}

				refine_element(i);
			}

			mesh.update_dual_graph();
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

	 	EdgeElementMap &edge_element_map()
		{
			return edge_element_map_;
		}
		
		void set_verbose(const bool val)
		{
			verbose = val;
		}

		const Mesh<Dim, ManifoldDim> &get_mesh() const
		{
			return mesh;
		}

		Mesh<Dim, ManifoldDim> &get_mesh()
		{
			return mesh;
		}

		const std::vector<Edge> &bisected_edges() const
		{
			return bisected_edges_;
		}

		void clear_bisected_edges()
		{
			bisected_edges_.clear();
		}

		Edge mandatory_edge_splitting(const Integer element_id)
		{
			if(mandatory_edge_splitting_.empty()) {
				return Edge();
			}

			assert(element_id < mandatory_edge_splitting_.size());
			assert(element_id >= 0);

			return mandatory_edge_splitting_[element_id];
		}

		// void refine_sides(const std::vector<Side<ManifoldDim>> &sides)
		// {

		// }

		void if_exist_refine_edges(const std::vector<Edge> &edges)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			// mandatory_edge_splitting_.resize(mesh.n_elements());

			for(auto e : edges) {
				refine_edge(e);
			}
		}

	private:
		Mesh<Dim, ManifoldDim> &mesh;
		std::vector<Integer> flags;
		std::vector<Integer> level;
		std::vector<std::array<Integer, ManifoldDim+1> > side_flags;
		EdgeNodeMap edge_node_map_;
		EdgeElementMap edge_element_map_;
		std::shared_ptr<EdgeSelect<Dim, ManifoldDim>> edge_select_;
		bool verbose;
		//tracking the refinement 
		std::vector<Edge> bisected_edges_;
		std::vector<std::vector<Edge>> mandatory_edge_splitting_;
	};
}

#endif //MARS_BISECTION_HPP
