#include <cstdlib>
#include <iostream>
#include <cassert>
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

int main(const int argc, const char *argv[])
{
	using namespace mars;
	test_bisection_2D();
	test_bisection_3D();
	test_bisection_4D();
	return EXIT_SUCCESS;
}


