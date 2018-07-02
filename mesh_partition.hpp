#ifndef MARS_MESH_PARTITION_HPP
#define MARS_MESH_PARTITION_HPP

#include "edge_split.hpp"

namespace mars {

	template<Integer Dim, Integer ManifoldDim>
	class MeshPartition {
	public:
		void set_partition_id(const Integer partition_id)
		{
			this->partition_id_ = partition_id;
		}

		Integer partition_id() const
		{
			return partition_id_;
		}

		void mark_partition_boundary(
			const Integer local_element_id,
			const Integer side_num,
			const Integer adj_partition)
		{
			auto &e = mesh.elem(local_element_id);
			e.side_tags[side_num] = partition_id_to_tag(adj_partition);
		}

		Integer add_and_index_elem(const Simplex<Dim, ManifoldDim> &element)
		{
			const Integer local_element_id = mesh.add_elem(element);
			assert(local_element_id == global_elem_id_.size());
			global_elem_id_.push_back(element.id);
			global_to_local_elem_[element.id] = local_element_id;
			return local_element_id;
		}

		Simplex<Dim, ManifoldDim> &global_elem(const Integer element_id)
		{
			auto it = global_to_local_elem_.find(element_id);
			assert(it != global_to_local_elem_.end());
			return mesh.elem(it->second);
		}

		const Simplex<Dim, ManifoldDim> &global_elem(const Integer element_id) const
		{
			auto it = global_to_local_elem_.find(element_id);
			assert(it != global_to_local_elem_.end());
			return mesh.elem(it->second);
		}

		void describe(std::ostream &os) const
		{
			os << "===============================\n";
			os << "partition: " << partition_id() << std::endl;
			os << "elements:\n";
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;
				os << "[" << i << ", " << global_elem_id_[i] << "]";

				for(auto n : mesh.elem(i).nodes) {
					os << " " << n;
				}

				os << " ->";

				for(auto n : mesh.elem(i).nodes) {
					os << " " << global_node_id_[n];
				}

				if(is_elem_interface(i)) {
					os << " | interfaces with: [";

					for(auto st : mesh.elem(i).side_tags) {
						if(is_partition_tag(st)) {
							os << " " << tag_to_partition_id(st);	
						} else {
							os << " -";
						}
					}

					os << " ]";
				}

				os << "\n";
			}

			os << "nodes:\n";
			for(Integer i = 0; i < mesh.n_nodes(); ++i) {
				os << "[" << i << ", " << global_node_id_[i] << "] " << node_partition_id_[i] << "\n";
			}
		}

		bool is_elem_interface(const Integer local_element_id) const
		{
			const auto &e = mesh.elem(local_element_id);
			for(auto st : e.side_tags) {
				if(is_partition_tag(st)) return true;
			}

			return false;
		}

		void add_and_index_nodes(
			const Mesh<Dim, ManifoldDim> &parent_mesh, 
			const std::vector<Integer> &node_partitioning,
			std::vector<Integer> &visited)
		{
			if(visited.empty()) {
				visited.resize(node_partitioning.size(), INVALID_INDEX);
			}

			const auto n_nodes = node_partitioning.size();

			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				assert(mesh.is_active(i));
				auto &e = mesh.elem(i);
				const auto &parent_e = parent_mesh.elem(global_elem_id_[i]);

				for(Integer k = 0; k < parent_e.nodes.size(); ++k) {
					const Integer n = parent_e.nodes[k];
					Integer local_n = INVALID_INDEX;

					if(visited[n] == INVALID_INDEX) {
						local_n = mesh.add_point(parent_mesh.point(n));
						global_node_id_.push_back(n);
						node_partition_id_.push_back(node_partitioning[n]);
						global_to_local_node_[n] = local_n;
						
						visited[n] = local_n;
					} else {
						local_n = visited[n];
					}

					assert(e.nodes[k] == parent_e.nodes[k]);
					e.nodes[k] = local_n;

					// std::cout << local_n << " -> " << n << " == " <<  visited[n] << " -> " << global_node_id_[local_n] << std::endl;
				}
			}

			for(auto i : global_node_id_) {
				visited[i] = INVALID_INDEX;
			}

			mesh.update_dual_graph();
		}

		const Mesh<Dim, ManifoldDim> &get_mesh() const
		{
			return mesh;
		}

		Mesh<Dim, ManifoldDim> &get_mesh()
		{
			return mesh;
		}

		bool edge_interfaces(
			const EdgeElementMap &eem,
			const Edge &e,
			std::vector<Integer> &partitions)
		{
			if(global_node_id_[e.nodes[0]] == INVALID_INDEX || 
			   global_node_id_[e.nodes[1]] == INVALID_INDEX)
			{
				// assert(false);
				return false;
			}

			Simplex<Dim, ManifoldDim-1> side;

			const auto &neighs = eem.elements(e);
			for(auto n : neighs) {
				if(n == INVALID_INDEX) continue;
				if(!is_elem_interface(n)) continue;

				const auto &n_e = mesh.elem(n);

				for(Integer s = 0; s < n_sides(n_e); ++s) {
					if(!is_partition_tag(n_e.side_tags[s])) continue;

					n_e.side(s, side);

					Integer n_matched = 0;
					for(auto node : side.nodes) {
						n_matched += node == e.nodes[0];
						n_matched += node == e.nodes[1];
					}

					assert(n_matched <= 2);
					
					if(n_matched == 2) {
						partitions.push_back(tag_to_partition_id(n_e.side_tags[s]));
					}
				}
			}
			
			return !partitions.empty();
		}

		void append_separate_interface_edges(
			const EdgeElementMap &eem,
			const EdgeNodeMap &enm,
			const std::vector<Edge> &local_edges,
			std::vector< std::vector<EdgeSplit> > &global_edges)
		{
			node_partition_id_.resize(mesh.n_nodes(), INVALID_INDEX);
			global_node_id_.resize(mesh.n_nodes(), INVALID_INDEX);


			// for(auto &ge : global_edges) 
			// {
			// 	ge.clear();
			// }

			for(auto &e : local_edges) {
				std::vector<Integer> partitions;

				if(edge_interfaces(eem, e, partitions)) {
					assert(e.nodes[0] < global_node_id_.size());
					assert(e.nodes[1] < global_node_id_.size());

					Edge global_edge(global_node_id_[e.nodes[0]], global_node_id_[e.nodes[1]]);
					// global_edge.fix_ordering

					//find midpoint owner and global_id
					Integer local_midpoint = enm.get(e);
					assert(local_midpoint != INVALID_INDEX);

					assert(local_midpoint < mesh.n_nodes());

					Integer midpoint = global_node_id_[local_midpoint];
					Integer owner = node_partition_id_[local_midpoint];

					if(owner == INVALID_INDEX) {
						//try to claim ownership
						owner = partition_id();
					} 

					EdgeSplit global_edge_split(global_edge, midpoint, owner);

					global_edge_split.partitions.insert(partition_id());
					global_edge_split.partitions.insert(partitions.begin(), partitions.end());

					for(auto p : global_edge_split.partitions) {
						global_edges[p].push_back(global_edge_split);
					}

					// std::cout << "[" << partition_id() << "] ";
					// global_edge_split.describe(std::cout);

					// global_edges[partition_id()].push_back(global_edge_split);
				} else {
					Integer local_midpoint = enm.get(e);
					node_partition_id_[local_midpoint] = partition_id();
				}
			}
		}

		// Integer update_ownership_of_midpoints(
		// 	const EdgeNodeMap &enm,
		// 	const std::vector<Edge> &refined_edges
		// 	)
		// {
		// 	Integer n_new_local_nodes = 0;

		// 	for(auto e : refined_edges) {
		// 		Integer local_midpoint_id = enm.get(e);

		// 		if(node_partition_id_.size() <= local_midpoint_id) {
		// 			node_partition_id_.resize(local_midpoint_id + 1, INVALID_INDEX);
		// 		}

		// 		if(node_partition_id_[local_midpoint_id] == INVALID_INDEX) {
		// 			auto p0 = node_partition_id_[e.nodes[0]];
		// 			auto p1 = node_partition_id_[e.nodes[1]];

		// 			assert(p0 != INVALID_INDEX);
		// 			assert(p1 != INVALID_INDEX);

		// 			if(p0 < p1) {
		// 				node_partition_id_[local_midpoint_id] = p0;
		// 			} else {
		// 				node_partition_id_[local_midpoint_id] = p1;
		// 			}

		// 			if(node_partition_id_[local_midpoint_id] == partition_id_) {
		// 				++n_new_local_nodes;
		// 			}
		// 		}
		// 	}

		// 	return n_new_local_nodes;
		// }

		Edge local_edge(const Edge &global_edge) const
		{
			auto it_0 = global_to_local_node_.find(global_edge.nodes[0]);
			auto it_1 = global_to_local_node_.find(global_edge.nodes[1]);
			
			// if(it_0 == global_to_local_node_.end() ||
			//    it_1 == global_to_local_node_.end()) {
			// 	return Edge();
			// }

			assert(it_0 != global_to_local_node_.end());
			assert(it_1 != global_to_local_node_.end());

			return Edge(it_0->second, it_1->second);
		}

		void update_midpoint_ids(
			const EdgeNodeMap &enm,
			const std::vector<EdgeSplit> &global_edges)
		{
			node_partition_id_.resize(mesh.n_nodes(), INVALID_INDEX);

			for(auto es : global_edges) {
				auto local_e = local_edge(es.edge);
				if(!local_e.is_valid()) {
					std::cerr << "[Error] met invalid edge for ";
					es.edge.describe(std::cerr);
					std::cout << std::endl;
					continue;
				}

				Integer local_mp = enm.get(local_e);
				if(es.midpoint != INVALID_INDEX) {
					assign_node_local_to_global(local_mp, es.midpoint);
				} 

				assert(es.owner != INVALID_INDEX);
				// assert(local_mp != INVALID_INDEX);
				
				if(local_mp != INVALID_INDEX) {
					if(node_partition_id_.size() <= local_mp) {
						node_partition_id_.resize(local_mp, INVALID_INDEX);
					}

					node_partition_id_[local_mp] = es.owner;
				}
			}
		}

		void assign_node_local_to_global(
			const Integer local_id,
			const Integer global_id)
		{
			global_node_id_[local_id]        = global_id;
			global_to_local_node_[global_id] = local_id;
		}

		void localize_edges(
			const std::vector<Edge> &global_edges,
			std::vector<Edge> &local_edges)
		{
			local_edges.clear();
			local_edges.reserve(global_edges.size());

			for(auto e : global_edges) {
				auto local_e = local_edge(e);
				local_edges.push_back(local_e);
			}
		}

		inline Integer max_gobal_node_id() const
		{
			return *std::max_element(global_node_id_.begin(), global_node_id_.end());
		}

		inline Integer max_gobal_elem_id() const
		{
			return *std::max_element(global_elem_id_.begin(), global_elem_id_.end());
		}

		Integer n_owned_nodes() const
		{
			Integer ret = 0;
			for(auto n : node_partition_id_) {
				ret += (n == partition_id());
			}

			return ret;
		}

		void assign_global_node_ids(
			const Integer begin,
			const Integer end)
		{
			global_node_id_.resize(mesh.n_nodes(), INVALID_INDEX);

			if(begin == end) return;

			Integer current_id = begin;
			
			Integer count = 0;
			for(Integer i = 0; i < mesh.n_nodes(); ++i) {
				if(node_partition_id_[i] != partition_id_ || 
				   global_node_id_[i] != INVALID_INDEX) continue;

				    // global_node_id_[i] = current_id++;

					assign_node_local_to_global(i, current_id++);
					count++;
			}

			assert(count == end - begin);
		}


		void assign_global_elem_ids(
			const Integer begin,
			const Integer end)
		{
			global_elem_id_.resize(mesh.n_elements(), INVALID_INDEX);

			Integer current_id = begin;
			
			Integer count = 0;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(global_elem_id_[i] == INVALID_INDEX) {
					global_elem_id_[i] = current_id++;
					count++;
				}
			}

			assert(count == end - begin);
		}

		const std::vector<Integer> &global_elem_id() const
		{
			return global_elem_id_;
		}
		const std::vector<Integer> &global_node_id() const
		{
			return global_node_id_;
		}

		const Integer global_node_id(const Integer local_node_id) const
		{
			return global_node_id_[local_node_id];
		}

	private:
		static bool is_partition_tag(const Integer tag)
		{
			return tag < INVALID_INDEX;
		}

		static Integer partition_id_to_tag(const Integer adj_partition)
		{
			return -(adj_partition + 2);
		}

		static Integer tag_to_partition_id(const Integer tag)
		{
			return -(tag + 2);
		}

		Mesh<Dim, ManifoldDim> mesh;
		Integer partition_id_;
		std::vector<Integer> global_elem_id_;
		std::vector<Integer> global_node_id_;
		std::vector<Integer> node_partition_id_;

		//global to local
		std::map<Integer, Integer> global_to_local_elem_;
		std::map<Integer, Integer> global_to_local_node_;
	};
}

#endif //MARS_MESH_PARTITION_HPP
