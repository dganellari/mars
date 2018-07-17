#ifndef MARS_NODE_RANK_HPP
#define MARS_NODE_RANK_HPP

namespace mars {

	class NodeRank {
	public:
		inline Integer rank(const Integer node_id) const
		{
			assert(node_id < ranks_.size());
			assert(ranks_[node_id] != INVALID_INDEX);
			return ranks_[node_id];
		}

		inline Integer rank(const Edge &local_edge) const
		{
			assert(local_edge.is_valid());
			assert(local_edge[0] < ranks_.size());
			assert(local_edge[1] < ranks_.size());
			assert(ranks_[local_edge[0]] != INVALID_INDEX);
			assert(ranks_[local_edge[1]] != INVALID_INDEX);

			return ranks_[local_edge[0]] + ranks_[local_edge[1]];
		}

		void set_rank_to_midpoint(
			const Edge &edge,
			const Integer local_midpoint_id)
		{
			set_rank(local_midpoint_id, rank(edge));
		}

		void describe(std::ostream &os) const
		{
			for(auto nr : ranks_) {
				os << nr << "\n";
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void init(const Mesh<Dim, ManifoldDim> &mesh)
		{
			ranks_.reserve(mesh.n_nodes() * 2);
			ranks_.resize(mesh.n_nodes(), 1);
		}

		void set_rank(const Integer &node_id, const Integer rank)
		{
			if(ranks_.size() <= node_id) {
				ranks_.resize(node_id + 1, INVALID_INDEX);
			}

			ranks_[node_id] = rank;
		}

	private:

		std::vector<Integer> ranks_;
	};
}

#endif //MARS_NODE_RANK_HPP
