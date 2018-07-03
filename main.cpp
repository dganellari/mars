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
	
	// template<typename T>
	// using ptr = std::shared_ptr<T>;


	template<Integer Dim, Integer ManifoldDim>
	void test_bisection(
		const Integer n_levels,
		std::vector<std::shared_ptr<MeshPartition<Dim, ManifoldDim>>> &parts)
	{
		ParBisection<Dim, ManifoldDim> b(parts);
		std::vector<std::vector<mars::Integer>> elements(parts.size());

		;
		for(Integer i = 0; i < n_levels; ++i) {
			std::cout << "xxxxxxxxxxxxxxxxxxxxxx\n";
			
			for(Integer k = 0; k < parts.size(); ++k) {
				if(k % 2 == 1) {
					Vector<Real, Dim> center;
					center.set(0.5);
					mark_hypersphere_for_refinement(
						parts[k]->get_mesh(),
						center,
						0.25,
						elements[k]
						);
				}
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

	Integer n_parts = 4;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		partitioning[i] = i % n_parts;
	}

	std::vector<std::shared_ptr<MeshPartition<2, 2>>> partitions;
	parition_mesh(mesh, n_parts, partitioning, partitions);

	write_mesh("mesh_2_p.eps", mesh, 10., PLOT_ID);

	test_bisection(10, partitions);

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

	Integer n_parts = 2;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		partitioning[i] = i % n_parts;
	}

	std::vector<std::shared_ptr<MeshPartition<3, 3>>> partitions;
	parition_mesh(mesh, n_parts, partitioning, partitions);

	write_mesh_partitions(
		"par3_",
		partitions,
		PLOT_NUMERIC_TAG);

	// write_mesh("mesh_3", mesh, 10., PLOT_ID);

	test_bisection(2, partitions);

}

int main(const int argc, const char *argv[])
{
	using namespace mars;
	// test_bisection_2D();
	// test_bisection_3D();
	// test_bisection_4D();

	test_partition_2D();
	// test_partition_3D();
	return EXIT_SUCCESS;
}


