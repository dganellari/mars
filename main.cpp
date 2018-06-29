#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "simplex.hpp"
#include "lagrange_element.hpp"
#include "mesh.hpp"
#include "bisection.hpp"
#include "vtk_writer.hpp"
#include "quality.hpp"
#include "utils.hpp"
#include "mesh_partition.hpp"

void test_bisection_2D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<2, 2> mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/square_2_def.MFEM", mesh);

	Quality<2, 2> q(mesh);
	q.compute();

	mark_boundary(mesh);
	// print_boundary_info(mesh);
	write_mesh("mesh_2_bisect_0.eps", mesh, 10., PLOT_ID);

	Bisection<2, 2> b(mesh);
	b.set_limit_to_1_level_diff(false);
	b.uniform_refine(1);

	Integer n_levels = 16;
	for(Integer i = 0; i < n_levels; ++i) {
		std::vector<mars::Integer> elements;
		
		mark_hypersphere_for_refinement(
			mesh,
			{0.5, 0.5},
			0.25,
			elements
			);

		std::cout << "n_marked(" << i << "/" << n_levels << ") : " << elements.size() << std::endl;

		b.refine(elements);
		q.compute();

		write_mesh("mesh_2_bisect_" + std::to_string(i+1) + ".eps", mesh, 10., PLOT_NUMERIC_TAG);
	}

	mesh.update_dual_graph();
	// print_boundary_info(mesh); 

	q.report.normalize_data_points();
	q.save_report("quality2.svg");
}

void test_bisection_3D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<3, 3> mesh;
	read_mesh("../data/cube_6.MFEM", mesh, true);

	Quality<3, 3> q(mesh);
	q.compute();

	
	mark_boundary(mesh);
	// print_boundary_info(mesh);

	std::cout << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

	Bisection<3, 3> b(mesh);
	b.uniform_refine(1);
	b.set_limit_to_1_level_diff(false);

	Integer n_levels = 10;
	for(Integer i = 0; i < n_levels; ++i) {
		std::vector<mars::Integer> elements;
		
		mark_hypersphere_for_refinement(
			mesh,
			{0.5, 0.5, 0.5},
			0.25,
			elements
			);

		std::cout << "n_marked(" << i << "/" << n_levels << ") : " << elements.size() << std::endl;

		b.refine(elements);
		q.compute();
	}

	VTKMeshWriter<Mesh<3, 3>> w;
	w.write("mesh_bisect_refined.vtu", mesh);

	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

	mesh.update_dual_graph();
	// print_boundary_info(mesh);

	q.report.normalize_data_points();
	q.save_report("quality3.svg");

	std::cout << "======================================\n";
}

void test_bisection_4D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<4, 4> mesh;
	read_mesh("../data/cube4d_24.MFEM", mesh);
	
	Quality<4, 4> q(mesh);
	q.compute();


	mark_boundary(mesh);
	// print_boundary_info(mesh);

	std::cout << "volume: " << mesh.volume() << std::endl;
	Bisection<4, 4> b(mesh);
	b.uniform_refine(1);
	b.set_limit_to_1_level_diff(false);

	Integer n_levels = 8;
	for(Integer i = 0; i < n_levels; ++i) {
		std::vector<mars::Integer> elements;
		
		mark_hypersphere_for_refinement(
			mesh,
			{0.5, 0.5, 0.5, 0.5},
			0.25,
			elements
			);

		std::cout << "n_marked(" << i << "/" << n_levels << ") : " << elements.size() << std::endl;
		
		b.refine(elements);
		q.compute();
	}
	
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	mesh.update_dual_graph();
	// print_boundary_info(mesh);

	q.report.normalize_data_points();
	q.save_report("quality4.svg");
}

namespace mars {
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

		for(const auto &p : meshes) {
			p.describe(std::cout);
		}
	}

	template<typename T>
	using ptr = std::shared_ptr<T>;


	template<Integer Dim, Integer ManifoldDim>
	void test_bisection_on(std::vector<MeshPartition<Dim, ManifoldDim>> &meshes)
	{
		using B = Bisection<Dim, ManifoldDim>;

		std::vector< ptr<B> > bisection;
		for(auto &m : meshes) {
			// m.describe(std::cout);
			bisection.push_back(std::make_shared<B>(m.get_mesh()));
		}

		std::cout << "------------------------------\n";

		Integer n_levels = 1;
		for(Integer i = 0; i < n_levels; ++i) {
			bool complete = false;

			Integer synchronization_loops = 0;
			while(!complete) {

				//add midpoint global-id if any
				std::vector< std::vector<Edge> > global_refined_edges(meshes.size());
				std::vector<Integer> nodes_offsets(meshes.size() + 1, 0);
				std::vector<Integer> elem_offsets(meshes.size() + 1, 0);

				//parallel step
				Integer max_node_id = 0;
				Integer max_elem_id = 0;
				for(Integer k = 0; k < meshes.size(); ++k) {
					auto b_ptr = bisection[k];

					max_node_id = std::max(max_node_id, meshes[k].max_gobal_node_id());
					max_elem_id = std::max(max_elem_id, meshes[k].max_gobal_elem_id());

					Integer prev_n_elem = b_ptr->get_mesh().n_elements();

					std::vector<mars::Integer> elements;
					if(meshes[k].partition_id() != 1) {
						mark_hypersphere_for_refinement(
							b_ptr->get_mesh(),
							{0.5, 0.5},
							0.25,
							elements
							);

						b_ptr->refine(elements);
						std::cout << "n_marked(" << i << "/" << n_levels << ") : " << elements.size() << std::endl;
					}

					nodes_offsets[k+1] = meshes[k].update_ownership_of_midpoints(
						b_ptr->edge_node_map(),
						b_ptr->bisected_edges()
					);

					elem_offsets[k+1] = b_ptr->get_mesh().n_elements() - prev_n_elem;

					meshes[k].append_separate_interface_edges(
						b_ptr->edge_element_map(),
						b_ptr->bisected_edges(),
						global_refined_edges);
				}
				
				//sync step
				nodes_offsets[0] = max_node_id + 1;
				std::partial_sum(nodes_offsets.begin(), nodes_offsets.end(), nodes_offsets.begin());

				elem_offsets[0] = max_elem_id + 1;
				std::partial_sum(elem_offsets.begin(), elem_offsets.end(), elem_offsets.begin());

				for(Integer k = 0; k < meshes.size(); ++k) {
					meshes[k].assign_global_node_ids(
						nodes_offsets[meshes[k].partition_id()],
						nodes_offsets[meshes[k].partition_id() + 1]
					);

					meshes[k].assign_global_elem_ids(
						elem_offsets[meshes[k].partition_id()],
						elem_offsets[meshes[k].partition_id() + 1]
					);

				}

				complete = true;
				for(Integer k = 0; k < meshes.size(); ++k) {
					if(!global_refined_edges[k].empty()) {
						complete = false;

						//TODO
						//refine edges
						auto b_ptr = bisection[k];

						std::vector<Edge> local_edges;
						meshes[k].localize_edges(global_refined_edges[k], local_edges);
						b_ptr->if_exist_refine_edges(local_edges);

						write_mesh(
							"mesh_2_inter_" + std::to_string(synchronization_loops) + "_" + std::to_string(k) + ".eps",
							meshes[k].get_mesh(), 10., PLOT_ID);
					} 
				}

				++synchronization_loops;
			}

			std::cout << "synchronization_loops: " << synchronization_loops << std::endl;
		}

		Integer p = 0;
		for(auto &m : meshes) {
			// m.get_mesh().describe(std::cout);
			m.describe(std::cout);
			write_mesh("mesh_2_p" + std::to_string(p++) + ".eps", m.get_mesh(), 10., PLOT_ID);
		}
	}

}


void test_partition()
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<2, 2> mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/square_2_def.MFEM", mesh);

	Bisection<2, 2> b(mesh);
	b.uniform_refine(1);

	std::vector<MeshPartition<2, 2>> partitions;
	parition_mesh(mesh, 2, {0, 1, 0, 1, 0, 1}, partitions);

	write_mesh("mesh_2_p.eps", mesh, 10., PLOT_ID);

	test_bisection_on(partitions);

}

int main(const int argc, const char *argv[])
{
	using namespace mars;
	// test_bisection_2D();
	// test_bisection_3D();
	// test_bisection_4D();

	test_partition();
	return EXIT_SUCCESS;
}


