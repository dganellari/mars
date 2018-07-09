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
			assert(edge_to_split_map.size() == splits_.size());

			auto e_temp = edge_split.edge;
			e_temp.fix_ordering();

			auto it = edge_to_split_map.find(e_temp);
			
			if(it == edge_to_split_map.end()) {
				Integer split_id = splits_.size();
				edge_to_split_map[e_temp] = split_id;
				splits_.push_back(edge_split);
				auto &es = splits_.back();
				es.owner = originator;
				assert(edge_to_split_map.size() == splits_.size());

				return es;
			}

			EdgeSplit &g_split = get_split(it->second);
			if(g_split.midpoint != INVALID_INDEX) {
				//ownership and global id already determined
				return g_split;
			}

			//update ids
			if(edge_split.midpoint != INVALID_INDEX) {
				assert(edge_split.owner != INVALID_INDEX);
				g_split = edge_split;
				return g_split;
			}

			//determine owner
			if(g_split.owner == INVALID_INDEX) {
				g_split.owner = originator;
			} else {
				g_split.owner = std::min(originator, g_split.owner);
			}

			assert(edge_to_split_map.size() == splits_.size());
			return g_split;
		}

		inline const EdgeSplit &get_split(const Edge &e) const
		{
			static const EdgeSplit null_;

			Edge e_temp = e;
			e_temp.fix_ordering();

			auto it = edge_to_split_map.find(e_temp);
			if(it == edge_to_split_map.end()) {
				return null_;
			}

			return get_split(it->second);
		}

		inline const EdgeSplit &get_split(const Integer split_id) const
		{
			assert(split_id >= 0);
			assert(split_id < splits_.size());
			assert(edge_to_split_map.size() == splits_.size());

			return splits_[split_id];
		}

		inline EdgeSplit &get_split(const Edge &e)
		{
			static EdgeSplit null_;


			Edge e_temp = e;
			e_temp.fix_ordering();

			auto it = edge_to_split_map.find(e_temp);
			if(it == edge_to_split_map.end()) {
				return null_;
			}

			return get_split(it->second);
		}

		inline EdgeSplit &get_split(const Integer split_id)
		{
			assert(split_id >= 0);
			assert(split_id < splits_.size());
			assert(edge_to_split_map.size() == splits_.size());
			
			return splits_[split_id];
		}


		inline Integer owner(const Edge &e) const
		{
			auto it = edge_to_split_map.find(e);
			if(it == edge_to_split_map.end()) {
				return INVALID_INDEX;
			}

			return get_split(it->second).owner;
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
			auto it = edge_to_split_map.find(e);
			if(it == edge_to_split_map.end()) {
				return INVALID_INDEX;
			}

			return get_split(it->second).midpoint;
		}

		bool check_consistent() const
		{
			assert(edge_to_split_map.size() == splits_.size());

			if(splits_.empty()) return true;

			std::vector<bool> visited(splits_.size(), false);

			for(auto etoid : edge_to_split_map) {
				auto &s = get_split(etoid.second);
				assert(s.edge == etoid.first);

				if(s.edge != etoid.first) return false;
				visited[etoid.second] = true;
			}


			for(Integer i = 0; i < visited.size(); ++i) {
				if(!visited[i]) {
					std::cerr << "split(" << i << ") not visited: ";
					get_split(i).describe(std::cerr);
					std::cerr << std::endl;
					return false;
				}
			}

			return true;
		}

		template<Integer Dim, Integer ManifoldDim>
		void set_midpoint_ids_to(
			const EdgeNodeMap &enm,
			MeshPartition<Dim, ManifoldDim> &part)
		{
			assert(check_consistent());
			if(splits_.empty()) return;

			Integer remove_index = 0;
			for(auto gs_it = splits_.begin(); gs_it != splits_.end(); ) {
				const auto &gs = *gs_it;
				auto map_it = edge_to_split_map.find(gs.edge);
				map_it->second -= remove_index;

				Edge e = part.local_edge(gs.edge);
				
				if(!e.is_valid()) {
					std::cout << "invalid local edge: ";
					gs.describe(std::cout);
					std::cout << "\n";

					++gs_it;
					continue;
				}

				Integer local_mp = enm.get(e);

				if(local_mp == INVALID_INDEX) { ++gs_it; continue; }

				assert(local_mp != INVALID_INDEX);

				part.assign_node_owner(local_mp, gs.owner);

				if(gs.midpoint == INVALID_INDEX) { ++gs_it; continue; }
				
				part.assign_node_local_to_global(local_mp, gs.midpoint);


				
				
				if(map_it != edge_to_split_map.end()) {
					edge_to_split_map.erase(map_it);
				} else {
					assert(false);
				}

				splits_.erase(gs_it);
				remove_index++;
			}

			assert(edge_to_split_map.size() == splits_.size());
		}

		template<Integer Dim, Integer ManifoldDim>
		void collect_splits_to_apply(
			const MeshPartition<Dim, ManifoldDim> &part,
			std::vector<EdgeSplit> &local_splits)
		{
			assert(edge_to_split_map.size() == splits_.size());

			local_splits.clear();
			local_splits.reserve(splits_.size());

			for(auto gs_it = splits_.begin(); gs_it != splits_.end();) {
				const auto &gs = *gs_it;
				local_splits.push_back(gs);
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void collect_splits_to_local_edges(
			const MeshPartition<Dim, ManifoldDim> &part,
			std::vector<Edge> &local_splits)
		{
			local_splits.clear();
			local_splits.reserve(splits_.size());

			for(auto gs_it = splits_.begin(); gs_it != splits_.end(); ++gs_it) {
				const auto &gs = *gs_it;
				auto le = part.local_edge(gs.edge);
				local_splits.push_back(le);
			}
		}

		inline bool empty() const
		{
			return splits_.empty() && to_communicate.empty();
		}

		inline std::vector<EdgeSplit>::iterator begin()
		{
			return splits_.begin();
		}

		inline std::vector<EdgeSplit>::iterator end()
		{
			return splits_.end();
		}

		inline std::vector<EdgeSplit>::const_iterator begin() const
		{
			return splits_.begin();
		}

		inline std::vector<EdgeSplit>::const_iterator end() const
		{
			return splits_.end();
		}

		void describe(std::ostream &os) const
		{
			os << ";;;;;;;;;;;;;;;;;;;;;;;;;;;\n";
			os << "splits(" << partition_id << "/" << n_parts << ") empty = " << empty() << "\n"; 
			for(const auto &s : splits_) {
				s.describe(os);
				auto it = edge_to_split_map.find(s.edge);
				assert(it != edge_to_split_map.end());
				os << "local_id = " << it->second << "\n";
			}

			os << "to_communicate\n";
			for(const auto &s : to_communicate) {
				s.describe(os);
			}

			os << ";;;;;;;;;;;;;;;;;;;;;;;;;;;\n";
		}

	private:
		std::vector<EdgeSplit>  splits_;
		std::map<Edge, Integer> edge_to_split_map;
		std::vector<EdgeSplit>  to_communicate;
		Integer partition_id;
		Integer n_parts;

		// std::vector<Edge, std::set<Integer> > edge_interface;
	};
}

#endif //MARS_EDGE_SPLIT_POOL_HPP
