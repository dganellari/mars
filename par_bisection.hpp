#ifndef MARS_PAR_BISECTION_HPP
#define MARS_PAR_BISECTION_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer Dim, Integer ManifoldDim>
	class Bisection;

	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer Dim, Integer ManifoldDim>
	class ParBisection {
	public:
		template<typename T>
		using ptr = std::shared_ptr<T>;
		using B = Bisection<Dim, ManifoldDim>;

		ParBisection(std::vector<MeshPartition<Dim, ManifoldDim>> &parts)
		: parts(parts)
		{
			for(auto &m : parts) {
				bisection.push_back(std::make_shared<B>(m.get_mesh()));
			}
		}

		void update_offsets()
		{
			node_offsets.resize(parts.size() + 1, 0);
			elem_offsets.resize(parts.size() + 1, 0);
			
			max_node_id = 0;
			max_elem_id = 0;

			node_offsets[0] = 0;
 			elem_offsets[0] = 0;

			for(Integer k = 0; k < parts.size(); ++k) {
				auto b_ptr = bisection[k];
				max_node_id = std::max(max_node_id, parts[k].max_gobal_node_id());
				max_elem_id = std::max(max_elem_id, parts[k].max_gobal_elem_id());

				auto n_nodes    = parts[k].n_owned_nodes();
				auto n_elements = b_ptr->get_mesh().n_elements();

				delta_node_offsets[k+1] = n_nodes    - node_offsets[k+1];
				delta_elem_offsets[k+1] = n_elements - elem_offsets[k+1];

				node_offsets[k+1] = n_nodes;
				elem_offsets[k+1] = n_elements;
			}

			delta_node_offsets[0] = max_node_id + 1;
			std::partial_sum(delta_node_offsets.begin(), delta_node_offsets.end(), delta_node_offsets.begin());

			delta_elem_offsets[0] = max_elem_id + 1;
			std::partial_sum(delta_elem_offsets.begin(), delta_elem_offsets.end(), delta_elem_offsets.begin());
		}

		void refine(std::vector<std::vector<mars::Integer>> &elements)
		{
			std::cout << "------------------------------\n";

			std::vector< std::vector<EdgeSplit> > global_refined_edges(parts.size());

			update_offsets();

			//local parallel refinement
			for(Integer k = 0; k < parts.size(); ++k) {
				if(elements[k].empty()) continue;

				auto b_ptr = bisection[k];
				b_ptr->refine(elements[k]);
			}

			//post-local refinement synchronization
			bool complete = false;
			Integer synchronization_loops = 0;
			
			while(!complete) {
				for(Integer k = 0; k < parts.size(); ++k) {
					auto b_ptr = bisection[k];
					parts[k].update_ownership_of_midpoints(
						b_ptr->edge_node_map(),
						b_ptr->bisected_edges()
					);
				}

				update_offsets();

				for(Integer k = 0; k < parts.size(); ++k) {
					parts[k].assign_global_node_ids(
						delta_node_offsets[parts[k].partition_id()],
						delta_node_offsets[parts[k].partition_id() + 1]
					);

					parts[k].assign_global_elem_ids(
						delta_elem_offsets[parts[k].partition_id()],
						delta_elem_offsets[parts[k].partition_id() + 1]
					);

					parts[k].append_separate_interface_edges(
						bisection[k]->edge_element_map(),
						bisection[k]->bisected_edges(),
						global_refined_edges
					);

					auto b_ptr = bisection[k];
					b_ptr->clear_bisected_edges();
				}

				complete = true;
				for(Integer k = 0; k < parts.size(); ++k) {
					if(!global_refined_edges[k].empty()) {
						complete = false;
						auto b_ptr = bisection[k];

						for(auto ges: global_refined_edges[k]) {	
							std::vector<Edge> global_edge = {ges.edge};
							std::vector<Edge> local_edges;

							parts[k].localize_edges(global_edge, local_edges);
							b_ptr->if_exist_refine_edges(local_edges);

							parts[k].update_midpoint_ids(
								b_ptr->edge_node_map(),
								{ges});

							complete = false;
						}
					} 
				}

				write_mesh_partitions("parts_" + std::to_string(synchronization_loops) + ".eps", parts, PLOT_UNIFORM);
				++synchronization_loops;
			}

			std::cout << "synchronization_loops: " << synchronization_loops << std::endl;


			write_mesh_partitions("parts_final.eps", parts, PLOT_UNIFORM);

			Integer p = 0;
			for(auto &m : parts) {
				std::cout << "p=[" <<  m.partition_id() << "]-------------------------------\n";
				m.describe(std::cout);
				// write_mesh("mesh_2_p" + std::to_string(p++) + ".eps", m.get_mesh(), 10., PLOT_ID);
				std::cout << "-------------------------------\n";
			}
		}

		std::vector<MeshPartition<Dim, ManifoldDim>> &parts;
		std::vector< ptr<B> > bisection;

		Integer max_node_id = 0;
		Integer max_elem_id = 0;
		std::vector<Integer> node_offsets;
		std::vector<Integer> elem_offsets;

		std::vector<Integer> delta_node_offsets;
		std::vector<Integer> delta_elem_offsets;
	};
}

#endif //MARS_PAR_BISECTION_HPP
