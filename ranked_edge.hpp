#ifndef MARS_RANKED_EDGE_HPP
#define MARS_RANKED_EDGE_HPP

#include "edge_select.hpp"

namespace mars {

	template<Integer Dim, Integer ManifoldDim>
	class RankedEdgeSelect final : public EdgeSelect<Dim, ManifoldDim> {
	public:
		class EdgeDesc {
		public:
			EdgeDesc(const Edge &edge, const Integer rank, const Integer edge_num) 
			: edge(edge), rank(rank), edge_num(edge_num)
			{
				assert(edge.is_valid());
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

		RankedEdgeSelect(const Map &map, const bool online_update)
		: map(map), online_update(online_update)
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

			assert(mesh.is_valid(element_id));

			for(auto n : mesh.elem(element_id).nodes) {
				assert(n != INVALID_INDEX);

				assert(mesh.is_node_valid(n));

				if(map.global(n) == INVALID_INDEX) {
					// std::cout << "invalid node: " << n << " for elem: " << element_id << std::endl;
					return false;
				}

				assert(map.local(map.global(n)) != INVALID_INDEX);
			}

			if(!online_update) {
				return element_is_ranked(mesh, element_id);
			} else {
				return true;
			}
		}

		Integer select(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const override
		{
			assert(can_refine(mesh, element_id));
			assert(element_is_ranked(mesh, element_id));

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
			return "RankedEdgeSelect";
		}

		//update ranking of children
		void element_refined(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id,
			const Edge &edge,
			const Integer local_midpoint_id) override
		{
			if(!online_update) { return; }

			const auto &e = mesh.elem(element_id);
			auto r = edge_rank_.get(edge);
			assert(r != INVALID_INDEX);

			assert(element_is_ranked(mesh, element_id));

			//new sub edges
			const Edge mid1(local_midpoint_id, edge[0]);
			const Edge mid2(local_midpoint_id, edge[1]);

			edge_rank_.update(mid1, r * 2);
			edge_rank_.update(mid2, r * 2);

			//generate rank for new edges
			Edge edge_k;
			Integer mid_face_rank = 0;
			for(Integer k = 0; k < n_edges(e); ++k) {
				e.edge(k, edge_k[0], edge_k[1]); edge_k.fix_ordering();
				if(edge_k == edge) { continue; } 
				mid_face_rank = std::max(mid_face_rank, rank(edge_k));
			}

			assert(!e.children.empty());

			for(auto c : e.children) {
				const auto &child = mesh.elem(c);
				
				for(Integer k = 0; k < n_edges(e); ++k) {
					child.edge(k, edge_k[0], edge_k[1]);
					edge_k.fix_ordering();

					assert(edge_k.is_valid());

					auto rk = rank(edge_k);
					if(rk == INVALID_INDEX) {
						edge_rank_.update(edge_k, mid_face_rank);
					}
				}
			}

#ifndef NDEBUG
			for(auto c : e.children) {
				assert(element_is_ranked(mesh, c));
			}
#endif // NDEBUG

		}

		bool element_is_ranked(
			const Mesh<Dim, ManifoldDim> &mesh,
			const Integer element_id) const
		{
			const auto &e = mesh.elem(element_id);

			Edge edge_k;
			for(Integer k = 0; k < n_edges(e); ++k) {
				e.edge(k, edge_k[0], edge_k[1]); edge_k.fix_ordering();
				auto r = rank(edge_k);
				// assert(r != INVALID_INDEX);
				if(r == INVALID_INDEX) return false;
			}

			return true;
		}

		inline Integer rank(const Edge &local_edge) const
		{
			auto r = edge_rank_.get(local_edge);
			return r;
		}

		void update(const Mesh<Dim, ManifoldDim> &mesh) override
		{
			if(online_update) {
				if(edge_rank_.empty()) {
					init_uniform_ranking(mesh);
					// init_local_ranking(mesh);
				}

			} else {
				init_ranking(mesh);
			}
		}


		bool is_valid(const Mesh<Dim, ManifoldDim> &mesh)
		{
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;

				if(!element_is_ranked(mesh, i)) {
					assert(false);
					return false;
				}
			}

			return true;
		}

		void describe(std::ostream &os) const override
		{
			for(auto it = edge_rank_.begin(); it != edge_rank_.end(); ++it) {
				it->first.describe(os);
				os << " -> " << it->second << "\n";
			}
		}

	private:

		void init_uniform_ranking(const Mesh<Dim, ManifoldDim> &mesh)
		{
			Integer v1, v2;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				const auto &e = mesh.elem(i);
				for(Integer k = 0; k < n_edges(e); ++k) {
					e.edge(k, v1, v2);
					
					edge_rank_.update(
						map.global(v1),
						map.global(v2),
						1
						);
				}
			}

			assert(is_valid(mesh));
		}

		void init_local_ranking(const Mesh<Dim, ManifoldDim> &mesh)
		{
			edge_rank_.clear();
			std::map<Edge, Integer> edge_rank_count;

			Integer v1, v2;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;

				const auto &e = mesh.elem(i);

				std::vector< std::pair<Real, Edge> > local_len;
				local_len.reserve(n_edges(e));
				
				for(Integer k = 0; k < n_edges(e); ++k) {
					e.edge(k, v1, v2);
					const Real e_len = (mesh.point(v1) - mesh.point(v2)).squared_norm();

					//global indexing for global consistency
					local_len.emplace_back(
						e_len, 
						Edge(map.global(v1), map.global(v2))
					);
				}

				std::sort(local_len.begin(), local_len.end());

				for(Integer k = 0; k < local_len.size(); ++k) {
					const auto &edge = local_len[k].second;

					Edge local_edge(
						map.local(edge[0]),
						map.local(edge[1])
					);

					edge_rank_.mapping()[local_edge] += k;
					++edge_rank_count[local_edge];
				}
			}

			for(auto &er : edge_rank_.mapping()) {
				auto it = edge_rank_count.find(er.first);
				assert(it != edge_rank_count.end());
				er.second /= it->second;
			}

			assert(is_valid(mesh));
		}

		void init_ranking(const Mesh<Dim, ManifoldDim> &mesh)
		{
			edge_rank_.clear();

			EdgeNodeMap global_edge_rank;

			Integer v1, v2;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;

				const auto &e = mesh.elem(i);
				for(Integer k = 0; k < n_edges(e); ++k) {
					e.edge(k, v1, v2);

					assert(map.global(v1) != INVALID_INDEX);
					assert(map.global(v2) != INVALID_INDEX);
					
					global_edge_rank.update(
						map.global(v1),
						map.global(v2),
						INVALID_INDEX
					);
				}
			}

			const Integer ne = global_edge_rank.size();
			std::vector<std::pair<Real, Edge>> len;
			len.reserve(ne);

			for(auto it = global_edge_rank.begin(); it != global_edge_rank.end(); ++it) {
				v1 = map.local(it->first[0]);
				v2 = map.local(it->first[1]);

				assert(v1 != INVALID_INDEX);
				assert(v2 != INVALID_INDEX);

				const Real e_len = (mesh.point(v1) - mesh.point(v2)).squared_norm();
				len.emplace_back(e_len, it->first);
			}

			std::sort(len.begin(), len.end());
			std::reverse(len.begin(), len.end());

			for(Integer i = 0; i < ne; ++i) {
				const Edge &e = len[i].second;

				Edge local_edge(
					map.local(e[0]),
					map.local(e[1])
					);

				edge_rank_.update(local_edge, i);
			}

			assert(is_valid(mesh));
		}

		const Map &map;
		EdgeNodeMap edge_rank_;
		bool online_update;
	};

}

#endif //MARS_RANKED_EDGE_HPP
