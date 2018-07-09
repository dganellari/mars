#ifndef MARS_MESH_PARTITION_HPP
#define MARS_MESH_PARTITION_HPP

#include "edge_split.hpp"
#include "edge_split_pool.hpp"
#include "dof_map.hpp"

namespace mars {

	template<Integer Dim, Integer ManifoldDim>
	class MeshPartition {
	public:
		MeshPartition(
			const Integer partition_id,
			const Integer n_partitions)
		: node_map_(partition_id, n_partitions),
		  elem_map_(partition_id, n_partitions)
		{

		}

		Integer partition_id() const
		{
			return node_map_.partition_id();
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
			// assert(local_element_id == global_elem_id_.size());
			// global_elem_id_.push_back(element.id);
			// global_to_local_elem_[element.id] = local_element_id;
			elem_map_.add(local_element_id, element.id, partition_id());
			return local_element_id;
		}

		Simplex<Dim, ManifoldDim> &global_elem(const Integer element_id)
		{
			auto index = elem_map_.local(element_id);
			assert(index != INVALID_INDEX);
			return mesh.elem(index);
		}

		const Simplex<Dim, ManifoldDim> &global_elem(const Integer element_id) const
		{
			auto index = elem_map_.local(element_id);
			assert(index != INVALID_INDEX);
			return mesh.elem(index);
		}

		void describe_element(const Integer element_id, std::ostream &os) const
		{
			os << "[" << element_id << ", " << elem_map_.global(element_id) << "]";

			for(auto n : mesh.elem(element_id).nodes) {
				os << " " << n;
			}

			os << " ->";

			for(auto n : mesh.elem(element_id).nodes) {
				os << " " << node_map_.global(n);
			}

			if(is_elem_interface(element_id)) {
				os << " | interfaces with: [";

				for(auto st : mesh.elem(element_id).side_tags) {
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

		void describe(std::ostream &os) const
		{
			os << "===============================\n";
			os << "partition: " << partition_id() << std::endl;
			os << "elements:\n";
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;
				describe_element(i, os);
			}

			os << "nodes:\n";
			for(Integer i = 0; i < mesh.n_nodes(); ++i) {
				os << "[" << i << ", " << node_map_.global(i) << "] " << node_map_.owner(i) << "\n";
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
				const auto &parent_e = parent_mesh.elem(elem_map_.global(i));

				for(Integer k = 0; k < parent_e.nodes.size(); ++k) {
					const Integer n = parent_e.nodes[k];
					Integer local_n = INVALID_INDEX;

					if(visited[n] == INVALID_INDEX) {
						local_n = mesh.add_point(parent_mesh.point(n));
						// global_node_id_.push_back(n);
						// node_partition_id_.push_back(node_partitioning[n]);
						// global_to_local_node_[n] = local_n;

						node_map_.add(
							local_n,
							n,
							node_partitioning[n]
						);
						
						visited[n] = local_n;
					} else {
						local_n = visited[n];
					}

					assert(e.nodes[k] == parent_e.nodes[k]);
					e.nodes[k] = local_n;
				}
			}

			for(Integer i = 0; i < mesh.n_nodes(); ++i) {
				visited[node_map_.global(i)] = INVALID_INDEX;
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
			partitions.clear();
			
			if(node_map_.global(e.nodes[0]) == INVALID_INDEX || 
			   node_map_.global(e.nodes[1]) == INVALID_INDEX)
			{
				assert(false);
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


		void update_edge_split_pool(
			const EdgeElementMap &eem,
			const EdgeNodeMap &enm,
			const std::vector<Edge> &local_edges,
			EdgeSplitPool &pool)
		{
			for(auto &e : local_edges) {
				std::vector<Integer> partitions;

				if(edge_interfaces(eem, e, partitions)) {
					Edge global_edge(
						node_map_.global(e.nodes[0]),
						node_map_.global(e.nodes[1])
					);

					//find midpoint owner and global_id
					Integer local_midpoint = enm.get(e);

					assert(local_midpoint != INVALID_INDEX);
					assert(local_midpoint < mesh.n_nodes());

					Integer midpoint = node_map_.insert_global(local_midpoint);
					Integer owner    = node_map_.insert_owner(local_midpoint);

					EdgeSplit global_edge_split(global_edge, midpoint, owner);
					global_edge_split.partitions.insert(partitions.begin(), partitions.end());

					pool.add_split(global_edge_split);
				} else {
					Integer local_midpoint = enm.get(e);
					node_map_.set_owner(local_midpoint, partition_id());
				}
			}
		}

		void read_from_edge_pool(
			const EdgeNodeMap &enm,
			EdgeSplitPool &pool)
		{
			pool.set_midpoint_ids_to(enm, *this);
		}

		void write_to_edge_pool(
			const EdgeNodeMap &enm,
			EdgeSplitPool &pool) const
		{
			for(auto it = pool.begin(); it != pool.end(); ++it) {
				auto &e_split = *it;
				auto le = local_edge(e_split.edge);
				const Integer midpoint = enm.get(le);

				if(midpoint == INVALID_INDEX) continue;
				if(node_map_.global(midpoint) == INVALID_INDEX) continue;
				pool.resolve_midpoint_id(e_split.edge, node_map_.global(midpoint));
			}
		}

		Integer local_midpoint(
			const EdgeNodeMap &enm,
			const Edge &global_edge) const
		{
			Edge e = local_edge(global_edge);
			if(!e.is_valid()) return INVALID_INDEX;

			return enm.get(e);
		}

		Integer global_midpoint(
			const EdgeNodeMap &enm,
			const Edge &global_edge) const
		{
			auto local = local_midpoint(enm, global_edge);
			if(local == INVALID_INDEX) return local;

			return node_map_.global(local);
		}

		Edge local_edge(const Edge &global_edge) const
		{
			auto n1 = node_map_.local(global_edge.nodes[0]);
			auto n2 = node_map_.local(global_edge.nodes[1]);

			if(n1 == INVALID_INDEX || n2 == INVALID_INDEX) {
				//retruning invalid edge
				std::cerr << "[Error] invalid edge ";
				global_edge.describe(std::cerr);
				std::cerr << "\n"; 
				// return Edge();
			}
			
			assert(n1 != INVALID_INDEX);
			assert(n2 != INVALID_INDEX);
			return Edge(n1, n2);
		}

		void update_midpoint_ids(
			const EdgeNodeMap &enm,
			const std::vector<EdgeSplit> &global_edges)
		{
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
				
				if(local_mp != INVALID_INDEX) {
					node_map_.set_owner(local_mp, es.owner);
				}
			}
		}

		void assign_node_local_to_global(
			const Integer local_id,
			const Integer global_id)
		{
			node_map_.set_global(local_id, global_id);
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

		inline Integer max_global_node_id() const
		{
			return node_map_.max_id();
		}

		inline Integer max_global_elem_id() const
		{
			return elem_map_.max_id();
		}

		Integer n_owned_nodes() const
		{
			return node_map_.n_owned();
		}

		void assign_node_owner(
			const Integer local_node_id,
			const Integer owner)
		{
			node_map_.set_owner(local_node_id, owner);
		}

		void assign_global_node_ids(
			const Integer begin,
			const Integer end)
		{
			node_map_.assign_global_ids(begin, end);
		}


		void assign_global_elem_ids(
			const Integer begin,
			const Integer end)
		{
			elem_map_.assign_global_ids(begin, end);
		}

		const Integer global_node_id(const Integer local_node_id) const
		{
			return node_map_.global(local_node_id);
		}

		inline const Map &node_map() const
		{
			return node_map_;
		}
		
		inline const Map &elem_map() const
		{
			return elem_map_;
		}

		inline Map &node_map()
		{
			return node_map_;
		}
		
		inline Map &elem_map()
		{
			return elem_map_;
		}

		inline bool is_valid()
		{
			bool ret = true;
			for(Integer i = 0; i < mesh.n_nodes(); ++i) {
				auto g = node_map_.global(i);
				auto o = node_map_.owner(i);

				if(g == INVALID_INDEX) {
					std::cerr << "global(" << i << ") = INVALID_INDEX" << std::endl;
					ret = false;
				}

				if(o == INVALID_INDEX) {
					std::cerr << "owner(" << i << ")  = INVALID_INDEX" << std::endl;
					ret = false;
				}
			}

			for(Integer i = 1; i < mesh.n_nodes(); ++i) {
				auto g      = node_map_.global(i);
				auto g_prev = node_map_.global(i-1);
				if(g == g_prev) {
					std::cerr << "global(" << (i-1) << ")  == global(" << (i) << ") == " << g << std::endl;
					ret = false;
				}

			}

			return ret;
		}

		void update_maps()
		{
			elem_map().resize(get_mesh().n_elements(), partition_id());
			node_map().resize(get_mesh().n_nodes(), INVALID_INDEX);
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
		Map node_map_;
		Map elem_map_;

	};
}

#endif //MARS_MESH_PARTITION_HPP
