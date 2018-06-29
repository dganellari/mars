#ifndef MARS_MESH_PARTITION_HPP
#define MARS_MESH_PARTITION_HPP

namespace mars {

	template<Integer Dim, Integer ManifoldDim>
	class MeshPartition {
	public:
		void set_partition_id(const Integer partition_id)
		{
			this->partition_id_ = partition_id;
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
			return local_element_id;
		}

		void describe(std::ostream &os) const
		{
			os << "===============================\n";
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
				os << "[i, " << global_node_id_[i] << "] " << node_partition_id_[i] << "\n";
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
						
						visited[n] = local_n;
					} else {
						local_n = visited[n];
					}

					assert(e.nodes[k] == parent_e.nodes[k]);
					e.nodes[k] = local_n;

					std::cout << local_n << " -> " << n << " == " <<  visited[n] << " -> " << global_node_id_[local_n] << std::endl;
				}
			}

			for(auto i : global_node_id_) {
				visited[i] = INVALID_INDEX;
			}
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

	private:
		Mesh<Dim, ManifoldDim> mesh;
		Integer partition_id_;
		std::vector<Integer> global_elem_id_;
		std::vector<Integer> global_node_id_;
		std::vector<Integer> node_partition_id_;
	};
}

#endif //MARS_MESH_PARTITION_HPP
