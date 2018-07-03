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
	void print_boundary_info(const Mesh<Dim, ManifoldDim> &mesh)
	{
		Simplex<Dim, ManifoldDim-1> side;
		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i) || !mesh.is_boundary(i)) continue;
			auto &e = mesh.elem(i);
			auto &adj = mesh.dual_graph().adj(i);

			std::cout << "[" << i << "]\n"; 
			for(Integer k = 0; k < n_sides(e); ++k) {
				if(adj[k] == INVALID_INDEX) {
					if(e.side_tags[k] == INVALID_INDEX) {
						std::cerr << "+++++++++ bad boundary tag ++++++++++++++\n";
					}

					std::cout << "\ttag(" << k << ") = " << e.side_tags[k] << " ( ";
					e.side(k, side);

					for(auto n : side.nodes) {
						std::cout << n << " ";
					}

					if(e.side_tags[k] == INVALID_INDEX) {
						std::cerr << "+++++++++++++++++++++++++++++++++++++++";
					}

					std::cout << ")\n";
				}
			}

			std::cout << std::endl;
		}
	}

	template<Integer Dim, Integer ManifoldDim>
	void parition_mesh(
		Mesh<Dim, ManifoldDim> &mesh,
		const Integer n_partitions,
		const std::vector<Integer> &partitioning,
		std::vector<MeshPartition<Dim, ManifoldDim>> &meshes)
	{
		assert(partitioning.size() == mesh.n_elements());

		mesh.update_dual_graph();

		meshes.clear();
		meshes.resize(n_partitions);

		std::vector<Integer> node_partitioning(mesh.n_nodes(), INVALID_INDEX);

		for(Integer i = 0; i < mesh.n_elements(); ++i) {
			if(!mesh.is_active(i)) continue;
			const Integer partition_id = partitioning[i];
			const auto &e = mesh.elem(i);

			auto &p = meshes[partition_id];
			p.set_partition_id(partition_id);

			const Integer local_element_id = p.add_and_index_elem(e);

			const auto &adj = mesh.dual_graph().adj(i);
			for(Integer f = 0; f < adj.size(); ++f) {
				if(adj[f] == INVALID_INDEX) continue;

				const Integer adj_partition_id = partitioning[adj[f]];
				if(adj_partition_id != partition_id) {

					p.mark_partition_boundary(
						local_element_id, f, adj_partition_id
						);
				}
			}

			for(Integer k = 0; k < n_nodes(e); ++k) {
				if(node_partitioning[e.nodes[k]] == INVALID_INDEX) {
					node_partitioning[e.nodes[k]] = partition_id;
				}
			}
		}

		std::vector<Integer> visited;
		for(auto &p : meshes) {
			p.add_and_index_nodes(mesh, node_partitioning, visited);
		}

		// for(const auto &p : meshes) {
		// 	p.describe(std::cout);
		// }
	}

}


#endif //MARS_UTILS_HPP