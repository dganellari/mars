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

			build_edge_interface();

			//local parallel refinement
			for(Integer k = 0; k < parts.size(); ++k) {
				if(elements[k].empty()) continue;

				auto b_ptr = bisection[k];
				b_ptr->refine(elements[k]);
			}
		}

		void build_edge_interface()
		{
			for(auto &m : parts) {
				if(bisection[m->partition_id()]->edge_element_map().empty()) {
					bisection[m->partition_id()]->edge_element_map().update(m->get_mesh());
				}
			}


			for(auto &m : parts) {
				edge_split_pool_[m->partition_id()]->build_edge_interface(
					bisection,
					parts
				);

				// m->describe(std::cout);
				// edge_split_pool_[m->partition_id()]->describe_edge_interface(std::cout);
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

				edge_split_pool_[k]->update_midpoint_parts(bisection[k]->edge_node_map(), *parts[k]);
			}

			build_edge_interface();

			// std::cout << "after " << std::endl;
			for(Integer k = 0; k < parts.size(); ++k) {
				bool ok = parts[k]->is_valid();
				
				if(!ok) {
					// parts[k]->node_map().describe(std::cout);
					// edge_split_pool_[k]->describe(std::cout);

					// for(Integer j = 0; j < parts.size(); ++j) {
					// 	parts[j]->node_map().describe(std::cout);
					// 	edge_split_pool_[j]->describe(std::cout);
					// }

					print_analysis(k);
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

		void print_analysis(const Integer partition_id)
		{
			const auto &part   = *parts[partition_id];
			const auto &e_pool = *edge_split_pool_[partition_id];
			auto &b 	       = *bisection[partition_id];
			auto &mesh         = part.get_mesh();
			const auto &enm    = b.edge_node_map();

			std::cout << "print_analysis(" << partition_id << ")" << std::endl;
			
			std::vector<EdgeSplit> inconsistent_edge_splits;

			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;
				const auto &e = mesh.elem(i);

				std::vector<Integer> invalid_nodes;
				for(auto n : e.nodes) {
					auto g = part.node_map().global(n);
					auto o = part.node_map().owner(n);

					if(g == INVALID_INDEX || o == INVALID_INDEX) {
						invalid_nodes.push_back(n);
					}
				}

				if(!invalid_nodes.empty()) {
					auto &parent = mesh.elem(e.parent_id);

					for(Integer k = 0; k < n_edges(parent); ++k) {
						Edge edge;
						parent.edge(k, edge.nodes[0], edge.nodes[1]);
						edge.fix_ordering();

						Integer midpoint = enm.get(edge);

						if(midpoint == INVALID_INDEX) continue;

						for(auto in : invalid_nodes) {
							if(in != midpoint) continue;
							Edge global_edge(
								part.node_map().global(edge.nodes[0]),
								part.node_map().global(edge.nodes[1])
							);

							std::cout << midpoint << " split(" 
									  << global_edge.nodes[0] << ","
									  << global_edge.nodes[1] << ")\n";

							auto &es = e_pool.get_split(global_edge);

							if(es.is_valid()) {
								std::cout << "found split ";
								es.describe(std::cout);

								inconsistent_edge_splits.push_back(es);
							}
						}
					}
				}
			}

			for(auto e_it = e_pool.begin(); e_it != e_pool.end(); ++e_it) {
				const auto &e_split = *e_it;
				const auto local_edge = part.local_edge(e_split.edge);

				if(!local_edge.is_valid()) {
					std::cout << "INVALID EDGE for ";
					e_split.describe(std::cout);
					std::cout << std::endl;
				} else {
					Integer midpoint = enm.get(local_edge);
					if(midpoint != INVALID_INDEX) continue;

					std::cout << "unapplied edge splitting for ";
					e_split.describe(std::cout);
					std::cout << std::endl;
					std::cout << "local_edge = ";
					local_edge.describe(std::cout);
					std::cout << std::endl;

					const auto &incident = b.edge_element_map().elements(local_edge);
					if(incident.empty()) {
						std::cout << "no incidence for edge\n";

						b.edge_element_map().update(mesh);

						const auto &incident = b.edge_element_map().elements(local_edge);
						if(!incident.empty()) {
							std::cout << "you should update the edge_element_map" << std::endl;
						} else {
							std::vector<Integer> elems;
							mesh.find_elements_by_nodes(
								local_edge.nodes.begin(),
								local_edge.nodes.end(),
								elems,
								false,
								true);

							for(auto e : elems) {
								mesh.describe_element(e, std::cout, false);
							}

							if(elems.empty()) {
								std::cout << "no elements found\n";
							}
						}
					}

					for(auto a : incident) {
						if(mesh.is_active(a)) {
							std::cout << "splitting OK for "    << a  << std::endl;
						}  else {
							std::cout << "splitting not OK for " << a << std::endl;
						}
					}
				}

			}

			edge_split_pool_[partition_id]->describe(std::cout);
			parts[partition_id]->node_map().describe(std::cout);

			for(auto es : inconsistent_edge_splits) {
				auto owner = es.owner;
				if(owner == INVALID_INDEX) continue;

				auto b = bisection[owner];
				auto p = parts[owner];

				Integer lmp = p->local_midpoint(b->edge_node_map(),  es.edge);
				Integer mp  = p->global_midpoint(b->edge_node_map(), es.edge);
				es.describe(std::cout);
				std::cout << "found midpoint in part " << owner <<  " with lid " << lmp << " and gid " << mp << std::endl;

				const auto le = p->local_edge(es.edge);
				const auto &incident = b->edge_element_map().elements(le);

				for(auto i : incident) {
					p->describe_element(i, std::cout);
				}

				std::vector<Integer> edge_interface;
				edge_split_pool_[partition_id]->edge_interface(
					es.edge,
					edge_interface);

				for(auto i : edge_interface) {
					std::cout << i << " ";
				}

				std::cout << " == ";

				edge_split_pool_[partition_id]->edge_interface(
					 es.edge,
					 edge_interface);

				for(auto i : edge_interface) {
					std::cout << i << " ";
				}

				std::cout << std::endl;

			}
		}


		void uniform_refine(const Integer n_levels) {
			
			for(Integer l = 0; l < n_levels; ++l) {
				if(verbose) {
					std::cout << "------------------------------\n";
				}

				build_edge_interface();

				//local parallel refinement
				for(Integer k = 0; k < parts.size(); ++k) {
					auto b_ptr = bisection[k];
					b_ptr->uniform_refine(1);
				}

				green_refinement();
			}
		}

		void green_refinement()
		{
			bool complete = false;

			Integer max_loops = 20;
			Integer loops = 0;
			while(!complete) {
				//2)
				update_edge_pools();
				
				//3)
				update_global_ids();

				//4)
				complete = conform_interfaces();

				std::cout << "loop " << loops++ << " done. complete = " << complete << std::endl;

				// for(Integer k = 0; k < parts.size(); ++k) {
				// 	if(!edge_split_pool_[k]->empty()) {
				// 		// edge_split_pool_[k]->describe(std::cout);
				// 		// parts[k]->node_map().describe(std::cout);
				// 		print_analysis(parts[k]->partition_id());
				// 	}
				// }

				if(loops >= max_loops && !complete) {
					for(Integer k = 0; k < parts.size(); ++k) {
						if(!edge_split_pool_[k]->empty()) {
							print_analysis(parts[k]->partition_id());
						}
					}
					
					// write_mesh_partitions(
					// 	"par_" + std::to_string(ManifoldDim) + "_" + std::to_string(loops) + ".eps",
					// 	parts,
					// 	PLOT_UNIFORM);

					assert(false);
				}
			}

			if(!complete) print_all();
		}

		void refine(std::vector<std::vector<mars::Integer>> &elements)
		{
			//1)
			refine_elements(elements);
			green_refinement();
			
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
		bool verbose = false;
	};
}

#endif //MARS_PAR_BISECTION_HPP
