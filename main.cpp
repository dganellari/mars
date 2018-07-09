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
#include "par_bisection.hpp"
#include "benchmark.hpp"

void test_bisection_2D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<2, 2> mesh;
	read_mesh("../data/square_2.MFEM", mesh);
	// read_mesh("../data/square_2_def.MFEM", mesh);

	Quality<2, 2> q(mesh);
	q.compute();

	mark_boundary(mesh);
	// print_boundary_info(mesh);
	

	Bisection<2, 2> b(mesh);
	auto edge_select = std::make_shared<NewestVertexEdgeSelect<2, 2>>();
	// auto edge_select = std::make_shared<LongestEdgeSelect<2, 2>>();
	// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<2, 2>>();
	edge_select->set_recursive(true);
	// edge_select->set_recursive(false);
	b.uniform_refine(1);
	b.set_edge_select(edge_select);
	// b.uniform_refine(3);

	write_mesh("mesh_2_bisect_0.eps", mesh, 10., PLOT_ID);

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
	// auto edge_select = std::make_shared<NewestVertexEdgeSelect<3, 3>>();
	// auto edge_select = std::make_shared<LongestEdgeSelect<3, 3>>();
	auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<3, 3>>();
	
	b.uniform_refine(7);

	edge_select->set_recursive(true);
	b.set_edge_select(edge_select);

	Integer n_levels = 12;
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
	
	// template<typename T>
	// using ptr = std::shared_ptr<T>;


	template<Integer Dim, Integer ManifoldDim>
	void test_bisection(
		const Integer n_levels,
		std::vector<std::shared_ptr<MeshPartition<Dim, ManifoldDim>>> &parts)
	{
		ParBisection<Dim, ManifoldDim> b(parts);
		// auto edge_select = std::make_shared<NewestVertexEdgeSelect<Dim, ManifoldDim>>();
		// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<Dim, ManifoldDim>>(true, true);
		auto edge_select = std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>(true, true);
		edge_select->set_recursive(true);
		// edge_select->set_recursive(false);
		b.set_edge_select(edge_select);
		std::vector<std::vector<mars::Integer>> elements(parts.size());

		// for(Integer k = 0; k < parts.size(); ++k) {
		// 	parts[k]->node_map().describe(std::cout);
		// }
		
		for(Integer i = 0; i < n_levels; ++i) {
			std::cout << "xxxxxxxxxxxxxxxxxxxxxx\n";
			
			for(Integer k = 0; k < parts.size(); ++k) {
				// if(k % 2 == 1) {
					Vector<Real, Dim> center;
					center.set(0.5);
					mark_hypersphere_for_refinement(
						parts[k]->get_mesh(),
						center,
						0.25,
						elements[k]
						);
				// }
			}

			b.verbose = i == n_levels-1;
			// b.verbose = true;
			b.refine(elements);

			std::cout << "xxxxxxxxxxxxxxxxxxxxxx\n";
		}
	}

}


void test_partition_2D()
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<2, 2> mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/square_2_def.MFEM", mesh);

	Bisection<2, 2> b(mesh);
	b.uniform_refine(2);

	std::vector<Integer> partitioning(mesh.n_elements());

	Integer n_parts = 2;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		partitioning[i] = i % n_parts;
	}

	std::vector<std::shared_ptr<MeshPartition<2, 2>>> partitions;
	parition_mesh(mesh, n_parts, partitioning, partitions);

	write_mesh("mesh_2_p.eps", mesh, 10., PLOT_ID);

	

	test_bisection(5, partitions);

	write_mesh_partitions(
		"par2.eps",
		partitions,
		PLOT_UNIFORM);

	// export_parts("part2", partitions);	

}

void test_partition_3D()
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<3, 3> mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/cube_6.MFEM", mesh, true);

	Bisection<3, 3> b(mesh);
	b.uniform_refine(1);

	std::vector<Integer> partitioning(mesh.n_elements());

	Integer n_parts = 6;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		partitioning[i] = i % n_parts;
	}

	std::vector<std::shared_ptr<MeshPartition<3, 3>>> partitions;
	parition_mesh(mesh, n_parts, partitioning, partitions);

	write_mesh_partitions(
		"before_par3_",
		partitions,
		PLOT_UNIFORM);

	// write_mesh("mesh_3", mesh, 10., PLOT_ID);

	test_bisection(1, partitions);

	write_mesh_partitions(
		"after_par3_",
		partitions,
		PLOT_UNIFORM);
}

void run_benchmarks()
{	
	using namespace mars;

	Benchmark<2, 2> b2;
	Mesh<2, 2> m2;
	read_mesh("../data/square_2_def.MFEM", m2);

	b2.run(14, m2, "b2");

	Benchmark<3, 3> b3;
	Mesh<3, 3> m3;
	read_mesh("../data/cube_6.MFEM", m3);

	b3.run(15, m3, "b3");

	Benchmark<4, 4> b4;
	Mesh<4, 4> m4;
	read_mesh("../data/cube4d_24.MFEM", m4);

	b4.run(6, m4, "b4");
}

int main(const int argc, const char *argv[])
{
	using namespace mars;
	// test_bisection_2D();
	// test_bisection_3D();
	// test_bisection_4D();

	// run_benchmarks();

	test_partition_2D();
	// test_partition_3D();
	return EXIT_SUCCESS;
}


