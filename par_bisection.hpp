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


		void set_edge_select(const std::shared_ptr<EdgeSelect<Dim, ManifoldDim>> &edge_select)
		{
			for(auto &b : bisection) {
				b->set_edge_select(edge_select);
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

						// std::cout << "edge_exchange(" << k << ", " << j << ") " << splits[j].size() << std::endl;
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

				parts[k]->update_maps();
			}

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

				// if(!edge_split_pool_[k]->empty()) {
				// 	edge_split_pool_[k]->describe(std::cout);
				// }
			}

			return !has_more;
		}

		void refine(std::vector<std::vector<mars::Integer>> &elements)
		{
			//1)
			refine_elements(elements);

			bool complete = false;

			Integer max_loops = 5;
			Integer loops = 0;
			while(!complete) {
				//2)
				update_edge_pools();
				
				//3)
				update_global_ids();

				//4)
				complete = conform_interfaces();

				std::cout << "loop " << loops++ << " done. complete = " << complete << std::endl;

				if(loops >= max_loops && !complete) {
					for(Integer k = 0; k < parts.size(); ++k) {
						if(!edge_split_pool_[k]->empty()) {
							edge_split_pool_[k]->describe(std::cout);
							parts[k]->node_map().describe(std::cout);
						}
					}
					
					write_mesh_partitions(
						"par_" + std::to_string(ManifoldDim) + "_" + std::to_string(loops) + ".eps",
						parts,
						PLOT_UNIFORM);

					assert(false);
				}
			}

			if(!complete) print_all();
		}

		void  __attribute__ ((used))  print_all() const
		{
			for(Integer k = 0; k < parts.size(); ++k) {
				edge_split_pool_[k]->describe(std::cout);
				parts[k]->node_map().describe(std::cout);
			}
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
