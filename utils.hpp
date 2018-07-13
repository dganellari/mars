#ifndef MARS_UTILS_HPP
#define MARS_UTILS_HPP


#include <cstdlib>
#include <iostream>
#include <cassert>
#include "simplex.hpp"
#include "lagrange_element.hpp"
#include "mesh.hpp"
#include "bisection.hpp"
#include "vtk_writer.hpp"
#include "quality.hpp"

namespace mars {

	template<Integer Dim, Integer ManifoldDim>
	class MeshPartition;



	template<Integer Dim, Integer ManifoldDim>
	void mark_hypersphere_for_refinement(
		const Mesh<Dim, ManifoldDim> &mesh,
		const Vector<Real, Dim> &center,
		const Real &radius,
		std::vector<Integer> &elements)
	{
		elements.clear();

		std::vector<Vector<Real, Dim>> points;
		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i)) continue;
			
			mesh.points(i, points);

			bool inside = false;
			bool outside = false;


			for(const auto &p : points) {
				auto dir = p - center;
				auto d = dir.norm();
				if(d < radius) {
					inside = true;
				} else if(d > radius) {
					outside = true;
				} else if(std::abs(d) < 1e-16) {
					inside = true;
					outside = true;
					break;
				}
			}

			if(inside && outside) {
				elements.push_back(i);
			}
		}
	}

	template<Integer Dim, Integer ManifoldDim>
	void mark_boundary(Mesh<Dim, ManifoldDim> &mesh)
	{
		using namespace mars;

		Simplex<Dim, ManifoldDim-1> side;
		mesh.update_dual_graph();
		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i) || !mesh.is_boundary(i)) continue;
			auto &e = mesh.elem(i);

			std::fill(e.side_tags.begin(), e.side_tags.end(), INVALID_INDEX);

			auto &adj = mesh.dual_graph().adj(i);

			for(Integer k = 0; k < n_sides(e); ++k) {
				if(adj[k] == INVALID_INDEX) {
					e.side(k, side);
					auto n = normal(side, mesh.points());

					Integer tag = 0;
					for(Integer d = 0; d < Dim; ++d) {
						if(std::abs(n(d) - 1.) < 1e-8) {
							tag = (d+1);
						} else if(std::abs(n(d) + 1.) < 1e-8) {
							tag = Dim + (d+1);
						}
					}

					e.side_tags[k] = tag;
				}
			}
		}
	}

	template<Integer Dim, Integer ManifoldDim>
	void export_mesh(
		const Mesh<Dim, ManifoldDim> &mesh,
		std::ostream &os)
	{

		os << "MFEM mesh v1.0\n\n";
		os << "dimension\n" << ManifoldDim << "\n";
		os << "elements\n" << mesh.n_active_elements() << "\n";

		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i)) continue;

			os << "-1 " << i << " ";

			const auto &e = mesh.elem(i);
			for(auto n : e.nodes) {
				os << " " << n;
			}

			os << "\n";
		}

		os << "\n";

		os << "vertices\n" << mesh.n_nodes() << "\n";
		os << Dim << "\n";

		for(Integer n = 0; n < mesh.n_nodes(); ++n) {
			
			for(Integer d = 0; d < Dim; ++d) {
				os << mesh.point(n)(d);

				if(d < Dim -1) { os << " ";  }
			}

			os << "\n";
		}
	}



	template<Integer Dim, Integer ManifoldDim>
	void export_elems_with_bad_tags(
		const Mesh<Dim, ManifoldDim> &mesh,
		std::ostream &os,
		const bool export_parent = false,
		const bool export_neighs = false)
	{
		std::set<Integer> elems;
		std::set<Integer> nodes;

		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i) || !mesh.is_boundary(i)) continue;
			const auto &e = mesh.elem(i);
			const auto &adj = mesh.dual_graph().adj(i);

			for(Integer k = 0; k < n_sides(e); ++k) {
				if(adj[k] == INVALID_INDEX) {
					if(e.side_tags[k] == INVALID_INDEX) {

						Integer e_id = i;
						if(export_parent) {
							e_id = mesh.elem(i).parent_id;
						}
						
						elems.insert(e_id);

						for(auto n : mesh.elem(e_id).nodes) {
							nodes.insert(n);
						}


						if(export_neighs) {
							const auto &e_adj = mesh.dual_graph().adj(e_id);

							for(Integer k = 0; k < n_sides(mesh.elem(e_id)); ++k) {
								if(e_adj[k] != INVALID_INDEX) {
									elems.insert(e_adj[k]);

									for(auto n : mesh.elem(e_adj[k]).nodes) {
										nodes.insert(n);
									}
								}
							}
						}

						break;
					}
				}
			}
		}
		std::map<Integer, Integer> global_to_local;

		Integer local_id = 0;
		for(auto s : nodes) {
			global_to_local[s] = local_id++;
		}

		os << "MFEM mesh v1.0\n\n";
		os << "dimension\n" << ManifoldDim << "\n";
		os << "elements\n" << elems.size() << "\n";

		for(auto e_id : elems) {
			os << "-1 " << e_id << " ";

			const auto &e = mesh.elem(e_id);
			for(auto n : e.nodes) {
				auto it = global_to_local.find(n);
				assert(it != global_to_local.end());

				os << " " << it->second;
			}

			os << "\n";
		}

		os << "\n";

		os << "vertices\n" << nodes.size() << "\n";
		os << Dim << "\n";

		for(auto n : nodes) {
			
			for(Integer d = 0; d < Dim; ++d) {
				os << mesh.point(n)(d);

				if(d < Dim -1) { os << " ";  }
			}

			os << "\n";
		}
	}

	template<Integer Dim, Integer ManifoldDim>
	void print_boundary_tags(
		const Simplex<Dim, ManifoldDim> &e
	) 
	{
		for(auto t : e.side_tags) {
			if(t == INVALID_INDEX) {
				std::cout << "- ";
			} else {
				std::cout << t << " ";
			}
		}

		std::cout << "\n";
	}

	template<Integer Dim, Integer ManifoldDim>
	void print_boundary_points(
		const Mesh<Dim, ManifoldDim> &mesh,
		std::ostream &os = std::cout
		const bool not_on_unit_cube = false)
	{
		std::vector<bool> is_boundary(mesh.n_nodes(), false);

		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			Simplex<Dim, ManifoldDim-1> side;
			auto &e   = mesh.elem(i);
			auto &adj = mesh.dual_graph().adj(i);

			for(Integer k = 0; k < n_sides(e); ++k) {
				if(adj[k] == INVALID_INDEX) {
					e.side(k, side);

					if(not_on_unit_cube) {

					} else {
						for(auto n : side.nodes) {
							is_boundary[n] = true;
						}
					}
				}
			}
		}

		for(Integer i = 0; i < mesh.n_nodes(); ++i) {
			if(is_boundary[i]) {
				os << "[" << i << "] " << mesh.point(i);
			}
		}
	}


	template<Integer Dim, Integer ManifoldDim>
	void print_boundary_info(
		const Mesh<Dim, ManifoldDim> &mesh, 
		const Integer i,
		const bool only_bad_tags, const bool sort_nodes = true)
	{
		Simplex<Dim, ManifoldDim-1> side;

		auto &e = mesh.elem(i);
		auto &adj = mesh.dual_graph().adj(i);


		if(!only_bad_tags) {
			std::cout << "[" << i << "]\n"; 
		}

		for(Integer k = 0; k < n_sides(e); ++k) {
			if(adj[k] == INVALID_INDEX) {
				if(e.side_tags[k] == INVALID_INDEX) {
					if(only_bad_tags) {
						std::cout << "[" << i << "]\n"; 
					}
					
					std::cerr << "+++++++++ bad boundary tag ++++++++++++++\n";
				}

				if(only_bad_tags && e.side_tags[k] != INVALID_INDEX) continue;

				std::cout << "\ttag(" << k << ") = " << e.side_tags[k] << " ( ";
				e.side(k, side);

				if(sort_nodes) {
					auto sorted_nodes = side.nodes;
					std::sort(sorted_nodes.begin(), sorted_nodes.end());
					for(auto n : sorted_nodes) {
						std::cout << n << " ";
					}

				} else {

					for(auto n : side.nodes) {
						std::cout << n << " ";
					}
				}

				std::cout << ")\n";

				if(e.side_tags[k] == INVALID_INDEX) {
					// if(e.parent_id != INVALID_INDEX) {
					// 	std::cout << "parent:\n";
					// 	mesh.describe_element(e.parent_id, std::cout);
					// 	std::cout << "children:\n";
					// 	for(auto c : mesh.elem(e.parent_id).children) {
					// 		mesh.describe_element(c, std::cout);
					// 	}
					// }

					std::cerr << "+++++++++++++++++++++++++++++++++++++++\n";
				}
			}
		}
	}

	template<Integer Dim, Integer ManifoldDim>
	void print_boundary_info(const Mesh<Dim, ManifoldDim> &mesh, const bool only_bad_tags, const bool sort_nodes = true)
	{
		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i) || !mesh.is_boundary(i)) continue;
			print_boundary_info(mesh, i, only_bad_tags, sort_nodes);
		}
	}

	template<Integer Dim, Integer ManifoldDim>
	void parition_mesh(
		Mesh<Dim, ManifoldDim> &mesh,
		const Integer n_partitions,
		const std::vector<Integer> &partitioning,
		std::vector<std::shared_ptr<MeshPartition<Dim, ManifoldDim>>> &parts)
	{
		assert(partitioning.size() == mesh.n_elements());

		mesh.update_dual_graph();

		parts.clear();
		parts.resize(n_partitions);

		for(Integer i = 0; i < n_partitions; ++i) {
			parts[i] = std::make_shared< MeshPartition<Dim, ManifoldDim> >(i, n_partitions);
		}

		std::vector<Integer> node_partitioning(mesh.n_nodes(), INVALID_INDEX);

		std::vector<std::set<Integer>> shared_nodes(mesh.n_nodes());

		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i)) continue;
			const Integer partition_id = partitioning[i];
			const auto &e = mesh.elem(i);

			auto &p = parts[partition_id];

			const Integer local_element_id = p->add_and_index_elem(e);

			const auto &adj = mesh.dual_graph().adj(i);
			for(Integer f = 0; f < adj.size(); ++f) {
				if(adj[f] == INVALID_INDEX) continue;

				const Integer adj_partition_id = partitioning[adj[f]];
				if(adj_partition_id != partition_id) {

					p->mark_partition_boundary(
						local_element_id, f, adj_partition_id
						);
				}
			}

			for(Integer k = 0; k < n_nodes(e); ++k) {
				if(node_partitioning[e.nodes[k]] == INVALID_INDEX) {
					node_partitioning[e.nodes[k]] = partition_id;
				}

				shared_nodes[e.nodes[k]].insert(partition_id);
			}
		}

		std::vector<Integer> visited;
		for(auto &p : parts) {
			p->add_and_index_nodes(mesh, node_partitioning, visited);

			auto &nm = p->node_map();
			for(Integer i = 0; i < p->get_mesh().n_nodes(); ++i) {
				const auto global_n = nm.global(i);
				
				for(auto pid : shared_nodes[global_n]) {
					nm.add_partition(i, pid);
				}
			}
		}

		for(const auto &p : parts) {
			// p->describe(std::cout);
			std::cout << p->partition_id() << " n_active_elements: " << p->get_mesh().n_active_elements() << std::endl;
		}
	}

	template<Integer Dim>
	bool write_mesh(
		const std::string &path,
		const Mesh<Dim, 4> &mesh
	)
	{			
		std::ofstream os(path + ".MFEM");
		if(!os.good()) return false;
		export_mesh(mesh, os);
		os.close();
		return false;
	}

	template<Integer Dim>
	bool write_mesh(
		const std::string &path,
		const Mesh<Dim, 5> &mesh
	)
	{	
		std::ofstream os(path + ".MFEM");
		if(!os.good()) return false;
		export_mesh(mesh, os);
		os.close();
		return false;
	}

	template<Integer Dim>
	bool write_mesh(
		const std::string &path,
		const Mesh<Dim, 6> &mesh
	)
	{	
		std::ofstream os(path + ".MFEM");
		if(!os.good()) return false;
		export_mesh(mesh, os);
		os.close();
		return false;
	}

}


#endif //MARS_UTILS_HPP