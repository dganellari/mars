#ifndef MARS_EDGE_SPLIT_POOL_HPP
#define MARS_EDGE_SPLIT_POOL_HPP

#include "edge.hpp"
#include "edge_split.hpp"

namespace mars {
	class EdgeSplitPool {
	public:
		explicit EdgeSplitPool(
			const Integer partition_id,
			const Integer n_parts)
		: partition_id(partition_id), n_parts(n_parts)
		{}

		//@return true if the split has been solved globally in terms
		//of ownership and node_id
		bool add_split(const EdgeSplit &edge_split)
		{
			assert(!edge_split.partitions.empty());

			auto &es = get_split(edge_split.edge);
			if(es.midpoint != INVALID_INDEX) {
				//does not need to comunicate
				//is global already
				return true;
			}

			//put the split in the global pool
			const EdgeSplit &gs = add_global_split(partition_id, edge_split);

			if(gs.only_on_partition(partition_id)) {
				return false;
			}
			
			if(!es.is_valid()) {
				to_communicate.push_back(gs);
			} 

			return false;
		}

		void pack(
			std::vector<std::vector<EdgeSplit>> &splits,
			const bool clear_buffers)
		{
			splits.resize(n_parts);

			if(clear_buffers) {
				for(Integer p = 0; p < n_parts; ++p) {
					splits[p].clear();
				}
			}

			for(const auto &s : to_communicate) {
				for(auto p : s.partitions) {
					if(p == partition_id) continue;

					splits[p].push_back(s);
				}
			}

			to_communicate.clear();
		}

		inline void unpack(
			const Integer originator,
			const std::vector<EdgeSplit> &splits)
		{
			for(const auto &s : splits) {
				add_global_split(originator, s);
			}
		}

		const EdgeSplit &add_global_split(
			const Integer originator,
			const EdgeSplit &edge_split)
		{
			assert(!edge_split.partitions.empty());

			auto it = global_splits.find(edge_split.edge);
			if(it == global_splits.end()) {
				auto &es = global_splits[edge_split.edge] = edge_split;
				es.owner = originator;
				return es;
			}

			EdgeSplit &g_split = it->second;
			if(g_split.midpoint != INVALID_INDEX) {
				//ownership and global id already determined
				return g_split;
			}

			//update ids
			if(edge_split.midpoint != INVALID_INDEX) {
				g_split = edge_split;
				return g_split;
			}

			//determine owner
			if(g_split.owner == INVALID_INDEX) {
				g_split.owner = originator;
			} else {
				g_split.owner = std::min(originator, g_split.owner);
			}

			return g_split;
		}

		inline const EdgeSplit &get_split(const Edge &e) const
		{
			static const EdgeSplit null_;

			auto it = global_splits.find(e);
			if(it == global_splits.end()) {
				return null_;
			}

			return it->second;
		}

		inline EdgeSplit &get_split(const Edge &e)
		{
			static EdgeSplit null_;

			auto it = global_splits.find(e);
			if(it == global_splits.end()) {
				return null_;
			}

			return it->second;
		}

		inline Integer owner(const Edge &e) const
		{
			auto it = global_splits.find(e);
			if(it == global_splits.end()) {
				return INVALID_INDEX;
			}

			return it->second.owner;
		}

		inline void resolve_midpoint_id(const Edge &e, const Integer m_id)
		{
			auto &s = get_split(e);
			assert(s.is_valid());
			assert(s.midpoint == INVALID_INDEX);
			s.midpoint = m_id;

			if(!s.only_on_partition(partition_id)) {
				to_communicate.push_back(s);
			}
		}

		inline Integer midpoint(const Edge &e) const
		{
			auto it = global_splits.find(e);
			if(it == global_splits.end()) {
				return INVALID_INDEX;
			}

			return it->second.midpoint;
		}

		template<Integer Dim, Integer ManifoldDim>
		void set_midpoint_ids_to(
			const EdgeNodeMap &enm,
			MeshPartition<Dim, ManifoldDim> &part)
		{
			for(auto gs_it = global_splits.begin(); gs_it != global_splits.end(); ) {
				const auto &gs = gs_it->second;
				Edge e = part.local_edge(gs.edge);
				
				if(!e.is_valid()) {
					++gs_it;
					continue;
				}

				Integer local_mp = enm.get(e);

				if(local_mp == INVALID_INDEX) { ++gs_it; continue; }

				assert(local_mp != INVALID_INDEX);

				part.assign_node_owner(local_mp, gs.owner);

				if(gs.midpoint == INVALID_INDEX) { ++gs_it; continue; }
				
				part.assign_node_local_to_global(local_mp, gs.midpoint);

				global_splits.erase(gs_it++);
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void collect_splits_to_apply(
			const MeshPartition<Dim, ManifoldDim> &part,
			std::vector<EdgeSplit> &splits)
		{
			splits.clear();
			splits.reserve(global_splits.size());

			for(auto gs_it = global_splits.begin(); gs_it != global_splits.end();) {
				const auto &gs = gs_it->second;
				splits.push_back(gs);
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void collect_splits_to_local_edges(
			const MeshPartition<Dim, ManifoldDim> &part,
			std::vector<Edge> &splits)
		{
			splits.clear();
			splits.reserve(global_splits.size());

			for(auto gs_it = global_splits.begin(); gs_it != global_splits.end(); ++gs_it) {
				const auto &gs = gs_it->second;
				auto le = part.local_edge(gs.edge);
				splits.push_back(le);
			}
		}

		inline bool empty() const
		{
			return global_splits.empty() && to_communicate.empty();
		}

		inline std::map<Edge, EdgeSplit>::iterator begin()
		{
			return global_splits.begin();
		}

		inline std::map<Edge, EdgeSplit>::iterator end()
		{
			return global_splits.end();
		}

		inline std::map<Edge, EdgeSplit>::const_iterator begin() const
		{
			return global_splits.begin();
		}

		inline std::map<Edge, EdgeSplit>::const_iterator end() const
		{
			return global_splits.end();
		}

		void describe(std::ostream &os) const
		{
			os << ";;;;;;;;;;;;;;;;;;;;;;;;;;;\n";
			os << "global_splits(" << partition_id << "/" << n_parts << ") empty = " << empty() << "\n"; 
			for(const auto &gs : global_splits) {
				gs.second.describe(os);
			}

			os << "to_communicate\n";
			for(const auto &gs : to_communicate) {
				gs.describe(os);
			}

			os << ";;;;;;;;;;;;;;;;;;;;;;;;;;;\n";
		}


		std::map<Edge, EdgeSplit> global_splits;
		std::vector<EdgeSplit> to_communicate;

		//not really necessary outside debugging
		// std::vector<EdgeSplit> resolved;
		Integer partition_id;
		Integer n_parts;
	};
}

#endif //MARS_EDGE_SPLIT_POOL_HPP
