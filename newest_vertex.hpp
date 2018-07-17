#ifndef MARS_NEWEST_VERTEX_HPP
#define MARS_NEWEST_VERTEX_HPP

#include "edge_select.hpp"

namespace mars {

	template<class Mesh>
	class ImplicitOrderEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		ImplicitOrderEdgeSelect()
		{}

		virtual bool repair_element() override
		{
			return false;
		}

		void reorder_edge(
			const Mesh &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const override
		{
			const auto &e = mesh.elem(element_id);

			Integer l1 = e.nodes[n_nodes(e) - 2];
			Integer l2 = e.nodes[n_nodes(e) - 1];

			if(l1 != v1) {
				std::swap(v1, v2);
			}

			assert(v1 == l1);
			assert(v2 == l2);
		}

		Integer select(
			const Mesh &mesh,
			const Integer element_id) const override
		{
			const auto &e = mesh.elem(element_id);
			Integer edge_num, v1, v2;

			edge_num = n_edges(mesh.elem(element_id)) - 1;
			e.edge(edge_num, v1, v2);

			assert(v1 == e.nodes[n_nodes(e) - 2]);
			assert(v2 == e.nodes[n_nodes(e) - 1]);
			return edge_num;
		}

		virtual Integer select(
			const Mesh &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const override
		{
			return select(mesh, element_id);
		}

		bool is_recursive() const override
		{
			return true;
		}


		virtual std::string name() const override
		{
			return "ImplicitOrderEdgeSelect";
		}
	};

	template<class Mesh>
	class GlobalNewestVertexEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		GlobalNewestVertexEdgeSelect(
			const Map &map,
			const bool recursive = true,
			const bool use_tollerance = true)
		: map(map), recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		void reorder_edge(
			const Mesh &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const override
		{
			if(map.global(v2) < map.global(v1)) {
				std::swap(v1, v2);
			}
		}

		bool can_refine(
			const Mesh &mesh,
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
			const Mesh &mesh,
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
			const Mesh &mesh,
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



	template<class Mesh>
	class NewestVertexEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		NewestVertexEdgeSelect(const bool recursive = true)
		: recursive_(recursive)
		{}

		Integer select(
			const Mesh &mesh,
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
			const Mesh &mesh,
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


	template<class Mesh>
	class NewestVertexAndLongestEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		NewestVertexAndLongestEdgeSelect(const bool recursive = true, const bool use_tollerance = true)
		: recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		Integer select(
			const Mesh &mesh,
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
			const Mesh &mesh,
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
}

#endif //MARS_NEWEST_VERTEX_HPP
