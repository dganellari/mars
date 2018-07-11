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


	Bisection<2, 2> b(mesh);
	b.uniform_refine(1);
	auto edge_select = std::make_shared<NewestVertexEdgeSelect<2, 2>>();
	// auto edge_select = std::make_shared<LongestEdgeSelect<2, 2>>();
	// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<2, 2>>();
	b.set_edge_select(edge_select);

	write_mesh("mesh_2_bisect_0.eps", mesh, 10., PLOT_ID);

	Integer n_levels = 14;
	for(Integer i = 0; i < n_levels; ++i) {
		std::vector<mars::Integer> elements;
		
		mark_hypersphere_for_refinement(
			mesh,
			{0.5, 0.5},
			0.25,
			elements
			);

		std::cout << "n_marked(" << (i+1) << "/" << n_levels << ") : " << elements.size() << std::endl;

		b.refine(elements);
		q.compute();

		write_mesh("mesh_2_bisect_" + std::to_string(i+1) + ".eps", mesh, 10., PLOT_NUMERIC_TAG);

		mesh.update_dual_graph();
		print_boundary_info(mesh, true);
	}

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
	b.uniform_refine(8);
	// auto edge_select = std::make_shared<NewestVertexEdgeSelect<3, 3>>();
	auto edge_select = std::make_shared<LongestEdgeSelect<3, 3>>(true);
	// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<3, 3>>();
	b.set_edge_select(edge_select);

	Integer n_levels = 8;
	for(Integer i = 0; i < n_levels; ++i) {
		std::vector<mars::Integer> elements;
		
		mark_hypersphere_for_refinement(
			mesh,
			{0.5, 0.5, 0.5},
			0.25,
			elements
			);

		std::cout << "n_marked(" << (i+1) << "/" << n_levels << ") : " << elements.size() << std::endl;

		b.refine(elements);
		q.compute();

		mesh.update_dual_graph();
		print_boundary_info(mesh, true);
	}

	VTKMeshWriter<Mesh<3, 3>> w;
	w.write("mesh_bisect_refined.vtu", mesh);

	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;


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

		std::cout << "n_marked(" << (i+1) << "/" << n_levels << ") : " << elements.size() << std::endl;
		
		b.refine(elements);
		q.compute();

		mesh.update_dual_graph();
		print_boundary_info(mesh, true);
	}
	
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

	q.report.normalize_data_points();
	q.save_report("quality4.svg");
}

namespace mars {
	
	// template<typename T>
	// using ptr = std::shared_ptr<T>;


	template<Integer Dim, Integer ManifoldDim>
	void test_bisection(
		const Integer n_levels,
		std::vector<std::shared_ptr<MeshPartition<Dim, ManifoldDim>>> &parts,
		const bool uniform_refine = false)
	{
		ParBisection<Dim, ManifoldDim> b(parts);
		auto edge_select = std::make_shared<NewestVertexEdgeSelect<Dim, ManifoldDim>>();
		// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<Dim, ManifoldDim>>(true, false);
		// auto edge_select = std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>(true, true);
		edge_select->set_recursive(true);
		// edge_select->set_recursive(false);
		b.set_edge_select(edge_select);
		std::vector<std::vector<mars::Integer>> elements(parts.size());

		// for(Integer k = 0; k < parts.size(); ++k) {
		// 	parts[k]->node_map().describe(std::cout);
		// }
		
		for(Integer i = 0; i < n_levels; ++i) {
			std::cout << "xxxxxxxxxxxxxxxxxxxxxx\n";
			std::cout << "level " << (i+1) << "/" << n_levels << std::endl;
			b.verbose = i == n_levels-1;
			// b.verbose = true;
			if(uniform_refine) {
				b.uniform_refine(1);
			} else {
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
				
				b.refine(elements);
			}

			for(auto p : parts) {
				p->get_mesh().update_dual_graph();
				print_boundary_info(p->get_mesh(), true);
			}

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
	mark_boundary(mesh);

	Bisection<2, 2> b(mesh);
	b.uniform_refine(3);

	std::vector<Integer> partitioning(mesh.n_elements());

	Integer n_parts = 7;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		partitioning[i] = i % n_parts;
	}

	std::vector<std::shared_ptr<MeshPartition<2, 2>>> parts;
	parition_mesh(mesh, n_parts, partitioning, parts);

	write_mesh_partitions(
		"par2_in.eps",
		parts,
		PLOT_UNIFORM);

	test_bisection(12, parts);

	write_mesh_partitions(
		"par2.eps",
		parts,
		PLOT_UNIFORM);

	for(const auto &p : parts) {
		std::cout << p->partition_id() << " n_active_elements: " << p->get_mesh().n_active_elements() << std::endl;
		p->get_mesh().update_dual_graph();
		print_boundary_info(p->get_mesh(), true);
	}	
}

void test_partition_3D()
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<3, 3> mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/cube_6.MFEM", mesh, true);
	mark_boundary(mesh);

	Bisection<3, 3> b(mesh);
	b.uniform_refine(2);

	std::vector<Integer> partitioning(mesh.n_elements(), 0);

	Integer n_parts = mesh.n_active_elements();
	Integer element_index = 0;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		if(mesh.is_active(i)) {
			partitioning[i] = (element_index++) % n_parts;
		}
	}

	std::vector<std::shared_ptr<MeshPartition<3, 3>>> parts;
	parition_mesh(mesh, n_parts, partitioning, parts);

	write_mesh_partitions(
		"before_par3_",
		parts,
		PLOT_UNIFORM);

	test_bisection(14, parts, false);

	write_mesh_partitions(
		"after_par3_",
		parts,
		PLOT_UNIFORM);

	for(const auto &p : parts) {
		std::cout << "---------------------\n";
		std::cout << p->partition_id() << " n_active_elements: " << p->get_mesh().n_active_elements() << std::endl;
		p->get_mesh().update_dual_graph();
		print_boundary_info(p->get_mesh(), true);
	}
}


void test_partition_4D()
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<4, 4> mesh;
	read_mesh("../data/cube4d_24.MFEM", mesh);
	mark_boundary(mesh);

	Bisection<4, 4> b(mesh);
	b.uniform_refine(1);

	std::vector<Integer> partitioning(mesh.n_elements());

	Integer n_parts = mesh.n_active_elements();
	Integer element_index = 0;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		if(mesh.is_active(i)) {
			partitioning[i] = (element_index++) % n_parts;
		}
	}

	std::vector<std::shared_ptr<MeshPartition<4, 4>>> parts;
	parition_mesh(mesh, n_parts, partitioning, parts);

	test_bisection(8, parts, false);

	for(const auto &p : parts) {
		std::cout << p->partition_id() << " n_active_elements: " << p->get_mesh().n_active_elements() << std::endl;
		p->get_mesh().update_dual_graph();
		print_boundary_info(p->get_mesh(), true);
	}
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
	// test_partition_2D();
	// test_partition_3D();
	test_partition_4D();
	return EXIT_SUCCESS;
}


