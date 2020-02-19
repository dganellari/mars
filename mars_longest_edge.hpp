#ifndef MARS_LONGEST_EDGE_SELECT_HPP
#define MARS_LONGEST_EDGE_SELECT_HPP

#include "mars_edge_select.hpp"
// #include "mars_edge_select_kokkos.hpp"
#include "mars_node_rank.hpp"

namespace mars {
	template<class Mesh>
	class LongestEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		using Edge = typename Mesh::Edge;

		MARS_INLINE_FUNCTION LongestEdgeSelect(const bool recursive = true, const bool use_tollerance = true)
		: recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		MARS_INLINE_FUNCTION Integer select(
			const Mesh &mesh,
			const Integer element_id) const override
		{
			const auto &e = mesh.elem(element_id);

			Integer edge_num = 0;
			Real len = 0;

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				Real len_i = (mesh.point(v1) - mesh.point(v2)).squared_norm();
				//Real len_i = mesh.point(v1).squared_distance(mesh.point(v2));

				if(len_i > len) {
					len = len_i;
					edge_num = i;
				}
			}

			//printf ("edge num = %i\n", edge_num);

			return edge_num;
		}

		MARS_INLINE_FUNCTION
		Integer stable_select(const Mesh &mesh, const Integer element_id) const
				override
		{
			const auto &e = mesh.elem(element_id);

			Integer edge_num = 0;
			Real len = 0;

			Integer otherV = 0, otherV2 = 0;
			for (Integer i = 0; i < n_edges(e); ++i)
			{
				Integer v1, v2;
				e.edge(i, v1, v2);

//				Real len_i = (mesh.point(v1) - mesh.point(v2)).squared_norm();
				Real len_i = mesh.point(v1).squared_distance(mesh.point(v2));

				if (len_i == len)
				{
					assert(!(min(v1, v2) == otherV && v1 + v2 == otherV2));
					if ((min(v1, v2) == otherV && v1 + v2 < otherV2)
							|| min(v1, v2) < otherV)
					{
						len = len_i;
						edge_num = i;
						otherV = min(v1, v2);
						otherV2 = v1 + v2;
					}
				}

				if (len_i > len)
				{
					len = len_i;
					edge_num = i;
					otherV = min(v1, v2);
					otherV2 = v1 + v2;
				}
			}

			return edge_num;
		}

		/*MARS_INLINE_FUNCTION
		Integer stable_select(const Mesh &mesh, const Integer element_id) const
				override
		{

			using Point = typename Mesh::Point;

			const auto &e = mesh.elem(element_id);

			Integer edge_num = 0;
			Real len = 0;

			Point smallestPoint;
			for (Integer i = 0; i < n_edges(e); ++i)
			{
				Integer v1, v2;
				e.edge(i, v1, v2);

				Real len_i = (mesh.point(v1) - mesh.point(v2)).squared_norm();

				if (len_i == len)
				{
					Point edgeSmallPoint = min(mesh.point(v1), mesh.point(v2));
					if (edgeSmallPoint < smallestPoint) {
						len = len_i;
						edge_num = i;
						smallestPoint = edgeSmallPoint;
					}
				}

				if (len_i > len)
				{
					len = len_i;
					edge_num = i;
					smallestPoint = min(mesh.point(v1), mesh.point(v2));
				}
			}

			return edge_num;
		}*/

		MARS_INLINE_FUNCTION
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

			for(Integer i = 0; i < n_edges(e); ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				const Real dist = (mesh.point(v1) - mesh.point(v2)).squared_norm();

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

		//check always the edges which are equals to the longest L
		//it then selects always the same edge as the stable select with the addition that it does it
		//only for the path with the smallest number of edges which are larger than L.
	/*	MARS_INLINE_FUNCTION
		Integer select(
					const Mesh &mesh,
					const Integer element_id,
					const EdgeElementMap &edge_element_map) const override
		{
			const auto &e = mesh.elem(element_id);

			Integer edge_num = 0;
			Real len = 0;

			constexpr Integer size = Combinations<Mesh::ManifoldDim + 1, 2>::value;
			std::array<Real,size> lengths;

			Integer otherV, otherV2;
			for(Integer i = 0; i < size; ++i) {
				Integer v1, v2;
				e.edge(i, v1, v2);

				Real len_i = (mesh.point(v1) - mesh.point(v2)).squared_norm();
				lengths[i]=len_i;

				if(len_i > len) {
					len = len_i;
					otherV = std::min(v1, v2);
					otherV2 = v1 + v2;
					//edge_num = i;
				}
			}


			Integer min_nr_edges = std::numeric_limits<Integer>::max();
			//Integer min_nr_edges = 0;
			for(Integer i=0; i<size;++i){

				if(lengths[i] == len){

					Integer v1, v2;
					Edge edge;
					mesh.elem(element_id).edge(i, edge.nodes[0], edge.nodes[1]);
					edge.fix_ordering();

					auto incidents = edge_element_map.elements(edge);

					Integer nr_edges_gt_len=0;
					for (auto l : incidents) {

						if(mesh.is_active(l)){

							const auto &el = mesh.elem(l);
							for(Integer k = 0; k < size; ++k) {
								el.edge(k, v1, v2);

								Real len_in = (mesh.point(v1) - mesh.point(v2)).squared_norm();

								if(len_in > len)
									++nr_edges_gt_len;
							}
						}
					}


					if(nr_edges_gt_len<min_nr_edges){
						min_nr_edges=nr_edges_gt_len;
						edge_num = i;
						otherV = edge.nodes[0];
						otherV2 = edge.nodes[1];
					}
					if (nr_edges_gt_len == min_nr_edges) { //TODO: stable for uniqueness
						v1 = edge.nodes[0];
						v2 = edge.nodes[1];
						assert(!(std::min(v1, v2) == otherV && v1 + v2 == otherV2));
						if ((std::min(v1, v2) == otherV && v1 + v2 < otherV2)
								|| std::min(v1, v2) < otherV) {
							edge_num = i;
							otherV = std::min(v1, v2);
							otherV2 = v1 + v2;
						}
					}
				}
			}

			return edge_num;
		}*/

		void set_recursive(const bool recursive)
		{
			recursive_ = recursive;
		}

		MARS_INLINE_FUNCTION bool  is_recursive() const override
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

	template<class Mesh>
	class UniqueLongestEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		UniqueLongestEdgeSelect(const bool recursive = true, const bool use_tollerance = true)
		: recursive_(recursive), use_tollerance_(use_tollerance)
		{}

		void reorder_edge(
			const Mesh &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const override
		{
			if(v2 < v1) {
				std::swap(v1, v2);
			}
		}

		Integer select(
			const Mesh &mesh,
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

				Real len_i = (mesh.point(v1) - mesh.point(v2)).squared_norm();

				if(len_i > len) {
					len = len_i;
					edge_num = ep.second;
				}
			}

			return edge_num;
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
			return "UniqueLongestEdgeSelect";
		}

	private:
		bool recursive_;
		bool use_tollerance_;
	};


	template<class Mesh>
	class GloballyUniqueLongestEdgeSelect final : public EdgeSelect<Mesh> {
	public:
		GloballyUniqueLongestEdgeSelect(
			const Map &map,
			const bool recursive = true,
			const bool use_tollerance = true)
		: map(map), recursive_(recursive), use_tollerance_(use_tollerance)//, use_trick_(true)
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
			const Mesh &mesh,
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
					return e1 < e2;
			});

			Integer edge_num = 0;
			Real len = 0;

			for(auto &ep : edge_pairs) {
				const Integer v1 = map.local(ep.edge[0]);
				const Integer v2 = map.local(ep.edge[1]);

				Real len_i = (mesh.point(v1) - mesh.point(v2)).squared_norm();

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

		// class EdgeDesc {
		// public:
		// 	EdgeDesc(
		// 		const Edge &edge,
		// 		const Real len,
		// 		const Integer edge_num)
		// 	: edge(edge), len(len), edge_num(edge_num)
		// 	{
		// 		assert(edge.is_valid());
		// 	}

		// 	inline bool operator<(const EdgeDesc &other) const
		// 	{
		// 		if(len < other.len) {
		// 			return true;
		// 		}

		// 		if(other.len < len) {
		// 			return false;
		// 		}

		// 		return edge < other.edge;
		// 	}

		// 	Edge edge;
		// 	Real len;
		// 	Integer edge_num;
		// };

		// Integer select(
		// 	const Mesh &mesh,
		// 	const Integer element_id) const override
		// {
		// 	assert(can_refine(mesh, element_id));

		// 	const auto &e = mesh.elem(element_id);
		// 	std::vector<EdgeDesc> edge_pairs;
		// 	edge_pairs.reserve(n_edges(e));

		// 	for(Integer i = 0; i < n_edges(e); ++i) {
		// 		Integer v1, v2;
		// 		e.edge(i, v1, v2);
		// 		edge_pairs.emplace_back(
		// 			Edge(map.global(v1), map.global(v2)),
		// 			(mesh.point(v1) - mesh.point(v2)).squared_norm(),
		// 		    i
		// 		);
		// 	}

		// 	std::sort(edge_pairs.begin(), edge_pairs.end());
		// 	return edge_pairs.back().edge_num;
		// }

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
			return "GloballyUniqueLongestEdgeSelect";
		}

		void element_refined(
			const Mesh &mesh,
			const Integer element_id,
			const Edge &edge,
			const Integer local_midpoint_id) override
		{
			if(node_rank_) {
				node_rank_->set_rank_to_midpoint(edge, local_midpoint_id);
			}
		}

		void set_node_rank(const std::shared_ptr<NodeRank> &node_rank)
		{
			node_rank_ = node_rank;
		}

		void update(const Mesh &mesh) override
		{
			if(node_rank_) {
				node_rank_->init(mesh);
			}
		}

	private:
		const Map &map;
		bool recursive_;
		bool use_tollerance_;
		std::shared_ptr<NodeRank> node_rank_;
	};
}

#endif //MARS_LONGEST_EDGE_SELECT_HPP
