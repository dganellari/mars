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
		using P = MeshPartition<Dim, ManifoldDim>;

		ParBisection(std::vector<ptr<P>> &parts)
		: parts(parts)
		{
			for(auto &m : parts) {
				bisection.push_back(std::make_shared<B>(m->get_mesh()));
				
				edge_split_pool_.push_back(
					std::make_shared<EdgeSplitPool>(
						m->partition_id(),
						parts.size())
				);
			}
		}

		void update_offsets()
		{
			node_offsets.resize(parts.size() + 1, 0);
			elem_offsets.resize(parts.size() + 1, 0);
			delta_node_offsets.resize(parts.size() + 1, 0);
			delta_elem_offsets.resize(parts.size() + 1, 0);
			
			max_node_id = 0;
			max_elem_id = 0;

			node_offsets[0] = 0;
 			elem_offsets[0] = 0;

			for(Integer k = 0; k < parts.size(); ++k) {
				auto b_ptr = bisection[k];
				max_node_id = std::max(max_node_id, parts[k]->max_global_node_id());
				max_elem_id = std::max(max_elem_id, parts[k]->max_global_elem_id());

				auto n_nodes    = parts[k]->n_owned_nodes();
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

			verify_offsets(delta_node_offsets);
			verify_offsets(delta_elem_offsets);
		}

		static bool verify_offsets(const std::vector<Integer> &o)
		{
			for(std::size_t i = 1; i < o.size(); ++i) {
				if(o[i-1] > o[i]) {
					assert(o[i-1] <= o[i]);
					return false;
				}
			}

			return true;
		}

		//1) (parallel)
		void refine_elements(std::vector<std::vector<mars::Integer>> &elements)
		{
			if(verbose) {
				std::cout << "------------------------------\n";
			}

			if(elements.empty()) {
				std::cerr << "refinement for elements is empty" << std::endl;
				return;
			}

			//local parallel refinement
			for(Integer k = 0; k < parts.size(); ++k) {
				if(elements[k].empty()) continue;

				auto b_ptr = bisection[k];
				b_ptr->refine(elements[k]);
			}
		}

		void exchange_edge_pools()
		{
			//communicate edges
			//fake commmunication
			{
				for(Integer k = 0; k < parts.size(); ++k) {
					std::vector<std::vector<EdgeSplit>> splits;
					edge_split_pool_[k]->pack(splits, false);

					for(Integer j = 0; j < parts.size(); ++j) {
						if(j == k) continue;

						std::cout << "edge_exchange(" << k << ", " << j << ") " << splits[j].size() << std::endl;
						edge_split_pool_[j]->unpack(k, splits[j]);
					}
				}
			}

		}

		//2) (synchronized)
		void update_edge_pools()
		{
			//parallel for
			for(Integer k = 0; k < parts.size(); ++k) {
				parts[k]->update_edge_split_pool(
					bisection[k]->edge_element_map(),
					bisection[k]->edge_node_map(),
					bisection[k]->bisected_edges(),
					*edge_split_pool_[k]
				);

				bisection[k]->clear_bisected_edges();
			}

			exchange_edge_pools();
		}

		//3) (synchronized)
		void update_global_ids()
		{
			//parallel for
			for(Integer k = 0; k < parts.size(); ++k) {
				parts[k]->read_from_edge_pool(
					bisection[k]->edge_node_map(),
					*edge_split_pool_[k]
				);

				parts[k]->elem_map().resize(parts[k]->get_mesh().n_elements(), k);
			}


			// std::cout << "before " << std::endl;
			// for(Integer k = 0; k < parts.size(); ++k) {
			// 	parts[k]->elem_map().describe(std::cout);
			// }

			// update global ids

			//communicate offsets
			//fake commmunication
			{
				std::vector<Integer> offsets;

				for(Integer k = 0; k < parts.size(); ++k) {
					parts[k]->node_map().pack_for_global(offsets, false);
				}

				for(Integer k = 0; k < parts.size(); ++k) {
					parts[k]->node_map().unpack_for_global(offsets);
				}
			}

			//fake commmunication
			{
				std::vector<Integer> offsets;
				
				for(Integer k = 0; k < parts.size(); ++k) {
					parts[k]->elem_map().pack_for_global(offsets, false);
				}

				for(Integer k = 0; k < parts.size(); ++k) {
					parts[k]->elem_map().unpack_for_global(offsets);
				}
			}

			//parallel for
			for(Integer k = 0; k < parts.size(); ++k) {
				parts[k]->write_to_edge_pool(
					bisection[k]->edge_node_map(),
					*edge_split_pool_[k]
				);
			}

			exchange_edge_pools();

			//parallel for
			for(Integer k = 0; k < parts.size(); ++k) {
				parts[k]->read_from_edge_pool(
					bisection[k]->edge_node_map(),
					*edge_split_pool_[k]
				);
			}

			// std::cout << "after " << std::endl;
			for(Integer k = 0; k < parts.size(); ++k) {
				bool ok = parts[k]->is_valid();
				
				if(!ok) {
					// parts[k]->node_map().describe(std::cout);
					// edge_split_pool_[k]->describe(std::cout);

					for(Integer j = 0; j < parts.size(); ++j) {
						parts[j]->node_map().describe(std::cout);
						edge_split_pool_[j]->describe(std::cout);
					}
				}

				assert( ok );
				// parts[k]->elem_map().describe(std::cout);
			}
		}	

		//4) (parallel)
		bool conform_interfaces()
		{
			//refine edges and go to 2)
			bool has_more = false;
			for(Integer k = 0; k < parts.size(); ++k) {
				// std::vector<EdgeSplit> splits;
			 // 	edge_split_pool_[k]->collect_splits_to_apply(
				// 		*parts[k],
				// 		splits
				// );

				std::vector<Edge> splits;
			 	edge_split_pool_[k]->collect_splits_to_local_edges(
						*parts[k],
						splits
				);

			 	//refine edges
				bisection[k]->if_exist_refine_edges(splits);
				
				has_more = has_more || !edge_split_pool_[k]->empty();

				if(!edge_split_pool_[k]->empty()) {
					edge_split_pool_[k]->describe(std::cout);
				}
			}

			return !has_more;
		}

		void refine(std::vector<std::vector<mars::Integer>> &elements)
		{
			//1)
			refine_elements(elements);

			bool complete = false;

			Integer max_loops = 10;
			Integer loops = 0;
			while(!complete) {
				//2)
				update_edge_pools();
				
				//3)
				update_global_ids();

				//4)
				complete = conform_interfaces();

				std::cout << "loop " << loops++ << " done. complete = " << complete << std::endl;

				// for(auto &es : edge_split_pool_) {
				// 	es->describe(std::cout);
				// }

				// for(auto &p : parts) {
				// 	p->describe(std::cout);
				// }

				write_mesh_partitions(
					"par2" + std::to_string(loops) + ".eps",
					parts,
					PLOT_UNIFORM);

				if(loops >= max_loops) {
					assert(false);
				}
			}
		}


		void refine_old(std::vector<std::vector<mars::Integer>> &elements)
		{
			if(verbose) {
				std::cout << "------------------------------\n";
			}

			if(elements.empty()) {
				std::cerr << "refinement for elements is empty" << std::endl;
				return;
			}

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

			std::vector<std::vector<EdgeSplit>> midpoint_id_data(parts.size());
			
			while(!complete) {
				for(Integer k = 0; k < parts.size(); ++k) {
					auto b_ptr = bisection[k];
				}

				update_offsets();

				//update midpoint ids
				for(Integer k = 0; k < parts.size(); ++k) {
					parts[k]->assign_global_node_ids(
						delta_node_offsets[parts[k]->partition_id()],
						delta_node_offsets[parts[k]->partition_id() + 1]
					);

					for(auto &mpd : midpoint_id_data) {
						for(auto &es : mpd) {
							if(es.owner == k) {
								Integer midpoint = bisection[k]->edge_node_map().get(parts[k]->local_edge(es.edge));
								assert(midpoint != INVALID_INDEX);

								Integer global_midpoint = parts[k]->global_node_id(midpoint);
								assert(global_midpoint != INVALID_INDEX);
								es.midpoint = global_midpoint;
							}
						}
					}
				}

				for(auto &gre : global_refined_edges)
				{
					gre.clear();
				}

				for(Integer k = 0; k < parts.size(); ++k) {
					for(auto &ges : midpoint_id_data[k]) {
						// std::cout << "[" << k << "] "; ges.describe(std::cout);

						parts[k]->update_midpoint_ids(
								bisection[k]->edge_node_map(),
								{ges}
						);
					}

					parts[k]->assign_global_elem_ids(
						delta_elem_offsets[parts[k]->partition_id()],
						delta_elem_offsets[parts[k]->partition_id() + 1]
					);

					parts[k]->append_separate_interface_edges(
						bisection[k]->edge_element_map(),
						bisection[k]->edge_node_map(),
						bisection[k]->bisected_edges(),
						global_refined_edges
					);

					auto b_ptr = bisection[k];
					b_ptr->clear_bisected_edges();
				}

				complete = true;

				midpoint_id_data.clear();
				midpoint_id_data.resize(parts.size());

				for(Integer k = 0; k < parts.size(); ++k) {

					if(!global_refined_edges[k].empty()) {

						complete = false;
						auto b_ptr = bisection[k];

						std::map<Edge, EdgeSplit> edge_2_owner;
						
						for(auto ges: global_refined_edges[k]) {
							if(verbose) {
								std::cout << "[" << k << "] "; ges.describe(std::cout);
							}

							auto it = edge_2_owner.find(ges.edge);
							
							if(it == edge_2_owner.end()) {
								edge_2_owner[ges.edge] = ges;
							} else {
								if(ges.midpoint != INVALID_INDEX) {
									assert(it->second.midpoint == INVALID_INDEX ||
										   it->second.midpoint == ges.midpoint);

									it->second.midpoint = ges.midpoint;
									it->second.owner    = ges.owner;
								} else if(it->second.midpoint == INVALID_INDEX)  {
									it->second.owner = std::min(k, it->second.owner);
									assert(it->second.owner >= 0);
								}
							}
						}

						for(auto ges: global_refined_edges[k]) {
							auto it = edge_2_owner.find(ges.edge);

							if(it->second.owner == k) {
								for(auto p : it->second.partitions) {
									if(p != k) {
										midpoint_id_data[p].push_back(it->second);
									}
								}
							}
						}

						for(auto p_ges: edge_2_owner) {	
							auto ges = p_ges.second;
						
							std::vector<Edge> global_edge = {ges.edge};
							std::vector<Edge> local_edges;

							parts[k]->localize_edges(global_edge, local_edges);
							b_ptr->if_exist_refine_edges(local_edges);

							parts[k]->update_midpoint_ids(
								b_ptr->edge_node_map(),
								{ges});

							complete = false;
						}
					} 
				}
				
				if(verbose) {
					write_mesh_partitions("parts_"  + std::to_string(n_refinements) + "_" + std::to_string(synchronization_loops) + ".eps", parts, PLOT_UNIFORM);
				}

				++synchronization_loops;
			}

			if(verbose) {
				std::cout << "synchronization_loops: " << synchronization_loops << std::endl;
				write_mesh_partitions("parts_" + std::to_string(n_refinements) + "_final.eps", parts, PLOT_UNIFORM);


				Integer p = 0;
				for(auto &m : parts) {
					std::cout << "p=[" <<  m->partition_id() << "]-------------------------------\n";
					m->describe(std::cout);
					// write_mesh("mesh_2_p" + std::to_string(p++) + ".eps", m.get_mesh(), 10., PLOT_ID);
					std::cout << "-------------------------------\n";
				}
			}
			++n_refinements;
		}

		std::vector< ptr<P> > &parts;
		std::vector< ptr<B> > bisection;
		std::vector< ptr<EdgeSplitPool> > edge_split_pool_;

		Integer max_node_id = 0;
		Integer max_elem_id = 0;
		std::vector<Integer> node_offsets;
		std::vector<Integer> elem_offsets;

		std::vector<Integer> delta_node_offsets;
		std::vector<Integer> delta_elem_offsets;
		bool verbose = false;
		Integer n_refinements = 0;
	};
}

#endif //MARS_PAR_BISECTION_HPP
