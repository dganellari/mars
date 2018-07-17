#ifndef MARS_OLDEST_EDGE_HPP
#define MARS_OLDEST_EDGE_HPP

#include "edge_select.hpp"
#include "node_rank.hpp"

namespace mars {

	template<Integer Dim, Integer ManifoldDim>
	class OldestEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		class EdgeDesc {
		public:
			EdgeDesc(const Edge &edge, const Integer rank, const Integer edge_num) 
			: edge(edge), rank(rank), edge_num(edge_num)
			{
				// assert(edge.is_valid());
			}

			inline bool operator<(const EdgeDesc &other) const
			{
				if(rank < other.rank) {
					return true;
				}

				if(other.rank < rank) {
					return false;
				}

				return edge < other.edge;
			}

			Edge edge;
			Integer rank;
			Integer edge_num;
		};

		OldestEdgeSelect(
			const Map &map,
			const std::shared_ptr<NodeRank> &node_rank = std::make_shared<NodeRank>())
		: map(map), node_rank_(node_rank)
		{}

		bool can_refine(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			assert(mesh.is_valid(element_id));

			if(is_oldest_edge_uniquelly_ranked(mesh, element_id)) {
				// std::cout << "unique" << std::endl;
				return true;
			}

			for(auto n : mesh.elem(element_id).nodes) {
				assert(n != INVALID_INDEX);

				assert(mesh.is_node_valid(n));

				if(map.global(n) == INVALID_INDEX) {
					return false;
				}

				assert(map.local(map.global(n)) != INVALID_INDEX);
			}

			return true;
		}

		bool is_oldest_edge_uniquelly_ranked(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const
		{
			const auto &e = mesh.elem(element_id);

			std::vector<Integer> ranks; 
			ranks.reserve(n_edges(e));
			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				Integer r = rank(Edge(v1, v2));
				ranks.push_back(r);
			}

			std::sort(ranks.begin(), ranks.end());
			return ranks[0] != ranks[1];
		}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			assert(can_refine(mesh, element_id));

			const auto &e = mesh.elem(element_id);
			std::vector<EdgeDesc> ranked_edges;
			ranked_edges.reserve(n_edges(e));

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				Integer r = rank(Edge(v1, v2)); assert(r != INVALID_INDEX);
				
				ranked_edges.emplace_back(
					Edge(map.global(v1), map.global(v2)),
					r,
					i
				);
			}

			std::sort(ranked_edges.begin(), ranked_edges.end());
			return ranked_edges[0].edge_num;
		}

		virtual Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
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
			return "OldestEdgeSelect";
		}

		inline Integer rank(const Edge &local_edge) const
		{
			return node_rank_->rank(local_edge);
		}

		void update(const Mesh<Dim, ManifoldDim> &mesh) override
		{
			node_rank_->init(mesh);
		}

		void element_refined(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id,
			const Edge &edge,
			const Integer local_midpoint_id) override
		{
			set_rank(local_midpoint_id, rank(edge));
		}

		void describe(std::ostream &os) const override
		{
			node_rank_->describe(os);
		}

	private:
		void set_rank(const Integer &node_id, const Integer rank)
		{
			node_rank_->set_rank(node_id, rank);
		}

		const Map &map;
		std::shared_ptr<NodeRank> node_rank_;
	};

}

#endif //MARS_OLDEST_EDGE_HPP