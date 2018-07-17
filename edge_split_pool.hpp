#ifndef MARS_EDGE_SPLIT_POOL_HPP
#define MARS_EDGE_SPLIT_POOL_HPP

#include "edge.hpp"
#include "edge_split.hpp"
#include "dof_map.hpp"

namespace mars {
	class EdgeSplitPool {
	public:
		template<typename T>
		using ptr = std::shared_ptr<T>;
	

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

		inline void resolve_midpoint_id(
			const Edge &e,
			const Integer local_m_id,
			const Integer m_id,
			Map &node_map)
		{
			auto &s = get_split(e);
			assert(s.is_valid());
			assert(s.midpoint == INVALID_INDEX);
			s.midpoint = m_id;

			if(!s.only_on_partition(partition_id)) {
				to_communicate.push_back(s);
			}

			update_edge_interface(e, m_id);

			std::vector<Integer> parts;
			edge_interface(e, parts);
			assert(!parts.empty());

			node_map.set_partitions(local_m_id, std::begin(parts), std::end(parts));
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

		template<class Mesh>
		void set_midpoint_ids_to(
			const EdgeNodeMap &enm,
			MeshPartition<Mesh> &part)
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
				update_edge_interface(gs.edge, gs.midpoint);

				std::vector<Integer> parts;
				edge_interface(gs.edge, parts);
				assert(!parts.empty());
				part.node_map().set_partitions(local_mp, std::begin(parts), std::end(parts));

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

		template<class Mesh>
		void collect_splits_to_apply(
			const MeshPartition<Mesh> &part,
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

		template<class Mesh>
		void collect_splits_to_local_edges(
			const MeshPartition<Mesh> &part,
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

		void describe_edge_interface(std::ostream &os) const
		{
			os << ";;;;;;;;;;;;;;;;;;;;;;;;;;;\n";

			os << "edge_interface(" << partition_id << ")\n";

			for(const auto &ei : edge_interface_) {
				ei.first.describe(os);
				os << " ->";
				
				for(auto p : ei.second) {
					os << " " << p;
				}

				os << "\n";
			}

			os << ";;;;;;;;;;;;;;;;;;;;;;;;;;;\n";
		}

		template<class Mesh>
		void build_edge_interface(
			const std::vector< ptr< Bisection<Mesh> > > &b,
			const std::vector< ptr< MeshPartition<Mesh> > > &parts)
		{
			// edge_interface_.clear();

			auto &p    = *parts[partition_id];
			const auto &mesh = p.get_mesh();

			std::vector<Integer> sharing;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				// if(!mesh.is_active(i)) continue;
				const auto &e = mesh.elem(i);


				for(Integer i = 0; i < n_edges(e); ++i) {
					Edge edge;
					e.edge(i, edge.nodes[0], edge.nodes[1]);
					edge.fix_ordering();

					assert(edge.is_valid());

					Edge global_edge(
						p.node_map().global(edge.nodes[0]),
						p.node_map().global(edge.nodes[1])
					);

					auto &inter = edge_interface_[global_edge];
					inter.clear();

					if(!inter.empty()) continue;

					p.node_map().intersect_partitions(
						std::begin(edge.nodes),
						std::end(edge.nodes),
						sharing);

					if(sharing.size() == 1) {
						assert(sharing[0] == partition_id);
						inter.push_back(partition_id);
						continue;
					}

					//FIXME requires communication
					for(auto s : sharing) {
						if(parts[s]->local_edge_exists(global_edge)) {
							auto le = parts[s]->local_edge(global_edge);
							if(!b[s]->edge_element_map().elements(le).empty()) {
								inter.push_back(s);
							}
						}
					}


					assert(!inter.empty());
					assert(p.node_map().is_unique(inter));
				}
			}
		}

		template<class Mesh>
		void update_midpoint_parts(
			const EdgeNodeMap &enm,
			MeshPartition<Mesh> &p
			)
		{
			for(const auto &ei : edge_interface_) {
				const auto edge = p.local_edge(ei.first);
				Integer local_midpoint = enm.get(edge);

				if(local_midpoint != INVALID_INDEX && !p.node_map().has_partitions(local_midpoint)) {
					p.node_map().set_partitions(local_midpoint, std::begin(ei.second), std::end(ei.second));
				}
			}
		}

		void edge_interface(const Edge &global_edge, std::vector<Integer> &partitions) const
		{
			auto it = edge_interface_.find(global_edge);
			if(it == edge_interface_.end()) {
				partitions.clear();
				return;
			}

			partitions = it->second;
		}

		void update_edge_interface(
			const Edge &global_edge,
			const Integer global_midpoint)
		{
			Edge e1(global_edge.nodes[0], global_midpoint);
			Edge e2(global_edge.nodes[1], global_midpoint);

			auto it = edge_interface_.find(global_edge);
			const auto parent_interface = it->second;

			assert(it != edge_interface_.end());
			assert(!it->second.empty());
			// assert(edge_interface_[e1].empty());
			// assert(edge_interface_[e2].empty());
			
			auto &ei1 = edge_interface_[e1];
			if(ei1.empty()) {
				ei1 = parent_interface;
			}

			auto &ei2 = edge_interface_[e2];
			if(ei2.empty()) {
				ei2 = parent_interface;
			}
		}

	private:
		std::vector<EdgeSplit>  splits_;
		std::map<Edge, Integer> edge_to_split_map;
		std::vector<EdgeSplit>  to_communicate;
		Integer partition_id;
		Integer n_parts;

		std::map<Edge, std::vector<Integer> > edge_interface_;
	};
}

#endif //MARS_EDGE_SPLIT_POOL_HPP
