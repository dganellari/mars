#ifndef MARS_BISECTION_HPP
#define MARS_BISECTION_HPP

#include "dof_map.hpp"
#include <iostream>

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

		virtual bool can_refine(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const
		{
			return true;
		}

		virtual void reorder_edge(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const
		{
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
	class UniqueLongestEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		UniqueLongestEdgeSelect(const bool recursive = true, const bool use_tollerance = true)
		: recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		void reorder_edge(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const override
		{
			if(v2 < v1) {
				std::swap(v1, v2);
			}
		}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			const auto &e = mesh.elem(element_id);
			std::vector< std::pair<Edge, Integer> > edge_pairs;
			edge_pairs.reserve(n_edges(e));

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);
				edge_pairs.emplace_back(Edge(v1, v2), i);
			}

			std::sort(edge_pairs.begin(), edge_pairs.end());

			Integer edge_num = 0;
			Real len = 0;

			for(auto &ep : edge_pairs) {
				const Integer v1 = ep.first.nodes[0];
				const Integer v2 = ep.first.nodes[1];

				Real len_i = (mesh.point(v1) - mesh.point(v2)).norm();

				if(len_i > len) {
					len = len_i;
					edge_num = ep.second;
				}
			}

			return edge_num;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{
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
			return "UniqueLongestEdgeSelect";
		}

	private:
		bool recursive_;
		bool use_tollerance_;
	};


	template<Integer Dim, Integer ManifoldDim>
	class GloballyUniqueLongestEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		GloballyUniqueLongestEdgeSelect(
			const Map &map,
			const bool recursive = true,
			const bool use_tollerance = true)
		: map(map), recursive_(recursive), use_tollerance_(use_tollerance)//, use_trick_(true)
		{}

		void reorder_edge(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const override
		{
			if(map.global(v2) < map.global(v1)) {
				std::swap(v1, v2);
			}
		}

		bool can_refine(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			for(auto n : mesh.elem(element_id).nodes) {
				assert(n != INVALID_INDEX);

				if(map.global(n) == INVALID_INDEX) return false;

				assert(map.local(map.global(n)) != INVALID_INDEX);
			}

			return true;
		}

		class EdgeDesc {
		public:
			EdgeDesc(const Edge &edge, const Integer edge_num) 
			: edge(edge), edge_num(edge_num)
			{
				assert(edge.is_valid());
			}

			inline bool operator<(const EdgeDesc &other) const
			{
				return edge < other.edge;
			}

			Edge edge;
			Integer edge_num;
		};

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			assert(can_refine(mesh, element_id));

			const auto &e = mesh.elem(element_id);
			std::vector<EdgeDesc> edge_pairs;
			edge_pairs.reserve(n_edges(e));

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);
				edge_pairs.emplace_back(Edge(map.global(v1), map.global(v2)), i);
			}

			std::sort(edge_pairs.begin(), edge_pairs.end(), 
				[](const EdgeDesc &e1, const EdgeDesc &e2) -> bool {
					return e2 < e1;
			});

			Integer edge_num = 0;
			Real len = 0;

			for(auto &ep : edge_pairs) {
				const Integer v1 = map.local(ep.edge[0]);
				const Integer v2 = map.local(ep.edge[1]);

				Real len_i = (mesh.point(v1) - mesh.point(v2)).norm();
				
				// if(use_trick_) {
				// 	len_i += ( len_i / ( ep.first[0] + ep.first[1] ) );
				// }

				if(len_i > len) {
					len = len_i;
					edge_num = ep.edge_num;
				}
			}

			return edge_num;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{
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
			return "GloballyUniqueLongestEdgeSelect";
		}

	private:
		const Map &map;
		bool recursive_;
		bool use_tollerance_;
		// bool use_trick_;
	};


	template<Integer Dim, Integer ManifoldDim>
	class GlobalNewestVertexEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		GlobalNewestVertexEdgeSelect(
			const Map &map,
			const bool recursive = true,
			const bool use_tollerance = true)
		: map(map), recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		void reorder_edge(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const override
		{
			if(map.global(v2) < map.global(v1)) {
				std::swap(v1, v2);
			}
		}

		bool can_refine(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			for(auto n : mesh.elem(element_id).nodes) {
				assert(n != INVALID_INDEX);

				if(map.global(n) == INVALID_INDEX) return false;

				assert(map.local(map.global(n)) != INVALID_INDEX);
			}

			return true;
		}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			assert(can_refine(mesh, element_id));

			const auto &e = mesh.elem(element_id);
			std::vector< std::pair<Edge, Integer> > edge_pairs;
			edge_pairs.reserve(n_edges(e));

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);
				edge_pairs.emplace_back(Edge(map.global(v1), map.global(v2)), i);
			}

			std::sort(edge_pairs.begin(), edge_pairs.end());
			return edge_pairs[0].second;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{
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
			return "GlobalNewestVertexEdgeSelect";
		}

	private:
		const Map &map;
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
		NewestVertexAndLongestEdgeSelect(const bool recursive = true, const bool use_tollerance = true)
		: recursive_(recursive), use_tollerance_(use_tollerance)
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
			if(!use_tollerance_) {
				return select(mesh, element_id);
			}

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
		bool use_tollerance_;
	};
	

	template<Integer Dim, Integer ManifoldDim>
	class Bisection {
	public:
		Bisection(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh), edge_select_(std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>()),
		verbose(false), fail_if_not_refine(false)
		{}

		void set_fail_if_not_refine(const bool val) {
			fail_if_not_refine = val;
		}

		void set_edge_select(const std::shared_ptr<EdgeSelect<Dim, ManifoldDim>> &edge_select)
		{
			edge_select_ = edge_select;
		}

		Integer add_elem(const Simplex<Dim, ManifoldDim> &e)
		{
			flags.push_back(NONE);
			Integer id = mesh.add_elem(e);
			mesh.repair_element(id);
			// mesh.update_side_flags_from_parent(id);
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

		void bisect_element(
			const Integer element_id,
			const Edge &edge)
		{
			Integer v1 = edge.nodes[0];
			Integer v2 = edge.nodes[1];

			edge_select_->reorder_edge(mesh, element_id, v1, v2);
			bisect_element(element_id, v1, v2);
		}

		void bisect_element(
			const Integer element_id,
			const Integer v1,
			const Integer v2)
		{
			mesh.elem(element_id).children.clear();
			mesh.set_active(element_id, false);

			Simplex<Dim, ManifoldDim> s;
			s.parent_id = element_id;
			
			if(verbose) {
				std::cout << "bisect(" << v1 << ", " << v2 << ") for " << element_id << std::endl;
			}

			auto midpoint = edge_node_map_.get(v1, v2);

			if(midpoint == INVALID_INDEX) {
				// midpoint = mesh.add_point(0.5 * (mesh.point(v1) + mesh.point(v2)));
				midpoint = mesh.add_point((mesh.point(v1) + mesh.point(v2))/2.);
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

			bisect_side_tags(element_id, Edge(v1, v2), midpoint);
			return;
		}

		inline Integer side_num(
			const Integer element_id,
			const Simplex<Dim, ManifoldDim-1> &side) const
		{
			auto nodes = side.nodes;
			std::sort(nodes.begin(), nodes.end());

			const auto &e = mesh.elem(element_id);

			Simplex<Dim, ManifoldDim-1> e_side;
			
			for(Integer i = 0; i < n_sides(e); ++i) {
				e.side(i, e_side);
				std::sort(std::begin(e_side.nodes), std::end(e_side.nodes));

				bool same_side = true;
				for(Integer k = 0; k < nodes.size(); ++k) {
					assert(nodes[k] != INVALID_INDEX);
					assert(e_side.nodes[k] != INVALID_INDEX);

					if(nodes[k] != e_side.nodes[k]) {
						same_side = false;
						break;
					} 
				}

				if(same_side) {
					return i;
				}
			}

			return INVALID_INDEX;
		}


		void bisect_side_tags(
			const Integer element_id,
			const Edge &edge,
			const Integer midpoint_id)
		{
			const auto &e = mesh.elem(element_id);

			for(auto c : e.children) {
				std::fill(mesh.elem(c).side_tags.begin(),
					mesh.elem(c).side_tags.end(),
					INVALID_INDEX);
			}	

			Simplex<Dim, ManifoldDim-1> side;
			Simplex<Dim, ManifoldDim-1> child_side;

			for(Integer i = 0; i < n_sides(e); ++i) {
				e.side(i, side);
				const Integer tag = e.side_tags[i];
				if(tag == INVALID_INDEX) continue;

				if(has_edge(side, edge[0], edge[1])) {
					//tag of split sides
					child_side.nodes[0] = midpoint_id;
					
					for(Integer j = 0; j < 2; ++j) {
						const Integer vj = edge[j];
						child_side.nodes[1] = vj;

						Integer local_ind = 2;
						for(auto s : side.nodes) {
							if(s == edge[0] || s == edge[1]) {
								continue;
							}
							child_side.nodes[local_ind++] = s;
						}

						bool found_side = false;
						for(auto c : e.children) {
							auto sn = side_num(c, child_side);

							if(INVALID_INDEX != sn) {
								assert(mesh.elem(c).side_tags[sn] == INVALID_INDEX || 
									mesh.elem(c).side_tags[sn] == tag);

								mesh.elem(c).side_tags[sn] = tag;
								found_side = true;
							}
						}

						assert(found_side);
					}
				} else {
					bool found_side = false;

					for(auto c : e.children) {
						auto sn = side_num(c, side);

						if(INVALID_INDEX != sn) {
							assert(mesh.elem(c).side_tags[sn] == INVALID_INDEX || 
								mesh.elem(c).side_tags[sn] == tag);

							mesh.elem(c).side_tags[sn] = tag;
							found_side = true;
						}
					}

					assert(found_side);
				}
			}
		}

		bool refine_element_recursive(
			const Integer element_id,
			const Edge &edge,
			const Integer max_level)
		{
			assert(has_edge(mesh.elem(element_id), edge.nodes[0], edge.nodes[1]));

			const Integer edge_num = edge_select_->select(mesh, edge, element_id);

			Edge new_edge;
			mesh.elem(element_id).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
			new_edge.fix_ordering();

			bool success = true;
			if(edge == new_edge) {
				bisect_element(element_id, edge);
			} else if(!edge_select_->is_recursive()) {
				bisect_element(element_id, edge);
			} else {
				success = refine_edge(new_edge);
				assert(!mesh.is_active(element_id));
			}

			return success;
		}

		void refine_element(const Integer element_id)
		{
			if(!edge_select_->can_refine(mesh, element_id)) {
				incomplete_elements_.push_back(element_id);
				assert(!fail_if_not_refine);
				return;
			}

			const Integer edge_num = edge_select_->select(mesh, element_id);
			Edge edge;
			mesh.elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();
			refine_edge(edge);
		}

		bool refine_edge(const Edge &edge)
		{
			bool visited_all_incidents = false;
			bool has_refined = false;
			bool success = true;
			bool all_done = true;
			auto incidents = edge_element_map_.elements(edge);

			while(!visited_all_incidents) {
				visited_all_incidents = true;

				for(auto i : incidents) {
					if(!mesh.is_active(i)) continue;

					assert(has_edge(mesh.elem(i), edge.nodes[0], edge.nodes[1]));
					all_done = false;

					const bool can_refine = edge_select_->can_refine(mesh, i);
					if(can_refine) {
						if(refine_element_recursive(i, edge, level[i])) {
							visited_all_incidents = false;
							has_refined = true;
							continue;
						} 

						assert(!fail_if_not_refine);
					}

					assert(!fail_if_not_refine);
					success = false;
				}

				const auto &next_incidents = edge_element_map_.elements(edge);

				if(next_incidents.size() > incidents.size()) {
					incidents = next_incidents;
					visited_all_incidents = false;
				}
			}

			if(!success) {
				incomplete_edges_.push_back(edge);
			}

			if(has_refined) {
				bisected_edges_.push_back(edge);
			}

			return success || all_done;
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
			// std::sort(sorted_elements.begin(), sorted_elements.end());

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

		void refine_edges(const std::vector<Edge> &edges)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			for(auto e : edges) {
				refine_edge(e);
			}

			mesh.update_dual_graph();
			mesh.tags() = level;
		}

		bool refine_incomplete()
		{
			std::cout << "incomplete elems: " << incomplete_elements_.size();
			std::cout << " edges: " << incomplete_edges_.size() << std::endl;

			auto elements_to_refine = std::move(incomplete_elements_);

			incomplete_elements_.clear();
			refine(elements_to_refine);

			auto edges_to_refine = std::move(incomplete_edges_);
			incomplete_edges_.clear();

			refine_edges(edges_to_refine);
			return incomplete_elements_.empty() && incomplete_edges_.empty();
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

		std::vector<Edge>    incomplete_edges_;
		std::vector<Integer> incomplete_elements_;
		bool fail_if_not_refine;
	};
}

#endif //MARS_BISECTION_HPP
