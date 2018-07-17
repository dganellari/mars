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
#include "test.hpp"
#include "ranked_edge.hpp"
#include "oldest_edge.hpp"
#include "longest_edge.hpp"


void test_bisection_2D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<2, 2> mesh;
	read_mesh("../data/square_2.MFEM", mesh);
	// read_mesh("../data/square_2_def.MFEM", mesh);

	Bisection<2, 2> b(mesh);
	b.uniform_refine(3);
	b.clear();

	mesh.clean_up();
	mesh.reorder_nodes();

	Quality<2, 2> q(mesh);
	q.compute();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	auto edge_select = std::make_shared<ImplicitOrderEdgeSelect<2, 2>>();
	b.set_edge_select(edge_select);
	b.uniform_refine(1);

	// write_mesh("mesh_2_bisect_0.eps", mesh, 10., PLOT_ID);

	Integer n_levels = 1;
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

		// write_mesh("mesh_2_bisect_" + std::to_string(i+1) + ".eps", mesh, 10., PLOT_NUMERIC_TAG);

		mesh.update_dual_graph();
		print_boundary_info(mesh, true);
	}

	// q.report.normalize_data_points();
	q.save_report("quality2.svg");


	b.tracking_begin();
	b.uniform_refine(2);
	b.tracking_end();
	
	write_mesh("mesh2_before.eps", mesh, 10., PLOT_NUMERIC_TAG);
	
	//use tracking information
	b.undo();
	write_mesh("mesh2_after.eps", mesh, 10., PLOT_NUMERIC_TAG);
}

void test_bisection_3D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<3, 3> mesh;
	read_mesh("../data/cube_6.MFEM", mesh, true);
	mesh.renumber_nodes();
	mesh.reorder_nodes();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	Bisection<3, 3> b(mesh);
	b.uniform_refine(2); b.clear();
	print_boundary_points(mesh, std::cout, true);
	
	mesh.clean_up();
	mesh.update_dual_graph();

	Quality<3, 3> q(mesh);
	q.compute();

	std::cout << "n_boundary_sides: " << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

	// auto edge_select = std::make_shared<NewestVertexEdgeSelect<3, 3>>();
	auto edge_select = std::make_shared<LongestEdgeSelect<3, 3>>(true);

	// auto edge_select = std::make_shared<ImplicitOrderEdgeSelect<3, 3>>();
	// auto edge_select = std::make_shared<UniqueLongestEdgeSelect<3, 3>>(true);	
	// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<3, 3>>();
	b.set_edge_select(edge_select);
	b.uniform_refine(1);

	// mesh.describe(std::cout);

	Integer n_levels = 10;
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


	// q.report.normalize_data_points();
	q.save_report("quality3.svg");

	std::cout << "======================================\n";
}

void test_bisection_4D()
{	
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<4, 4> mesh;
	read_mesh("../data/cube4d_24.MFEM", mesh);
	mesh.renumber_nodes();
	
	Quality<4, 4> q(mesh);
	q.compute();
	mark_boundary(mesh);

	std::cout << "volume: " << mesh.volume() << std::endl;
	Bisection<4, 4> b(mesh);
	b.set_edge_select(std::make_shared<UniqueLongestEdgeSelect<4, 4>>());
	b.uniform_refine(1);
	print_boundary_info(mesh, true);

	Integer n_levels = 5;

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

	// q.report.normalize_data_points();
	q.save_report("quality4.svg");
}

namespace mars {


	template<class GlobalEdgeSelect, Integer Dim, Integer ManifoldDim>
	bool test_incomplete(
		Mesh<Dim, ManifoldDim> &mesh,
		Map &map,
		const std::shared_ptr<GlobalEdgeSelect> &edge_select,
		const bool use_uniform_refinement = false)
	{
		Integer each_element = 11;
		bool bypass_incomplete = true;

		if(bypass_incomplete) {
			map.resize(mesh.n_nodes(), 0);
			map.identity();
		} else {
			Integer element_index = 0;
			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;

				if(((element_index++ % each_element) == 0) && 
					!edge_select->can_refine(mesh, i)) {
					//repair global index
					for(auto n : mesh.elem(i).nodes) {
						map.set_global(n, n);
					}
				}
			}
		}

		Bisection<Dim, ManifoldDim> b(mesh);
		b.set_edge_select(edge_select);
		b.tracking_begin();

		edge_select->update(mesh);

		if(use_uniform_refinement) {
			b.uniform_refine(1);
		} else {
			Vector<Real, Dim> center;
			center.set(0.5);
			std::vector<Integer> elements;
			mark_hypersphere_for_refinement(
				mesh,
				center,
				0.25,
				elements
				);

			b.refine(elements);
		}

		write_mesh("m2_inc.eps", mesh);

		///////////////////////////////////////////////////

		Integer max_iter = 20;
		for(Integer i = 0; i < max_iter; ++i) {
			std::cout << "iter: " << (i + 1) << "/" << max_iter << std::endl;
			map.resize(mesh.n_nodes(), 0);
			map.identity();

			b.set_fail_if_not_refine(i == max_iter -1);

			edge_select->update(mesh);
			// edge_select->describe(std::cout);

			std::cout << "n_nodes " << mesh.n_nodes() << std::endl;

			if(b.refine_incomplete()) {
				break;
			}
		}

		write_mesh("m2.eps", mesh);

		print_boundary_info(mesh, true, true);
		
		if(ManifoldDim <= 4) {
			print_boundary_points(mesh, std::cout, true);
		}

		b.tracking_end();
		if(!mesh.is_conforming()) {
			b.undo(); 
			std::cerr << "[Warning] encountered non-conforming mesh undoing refinement" << std::endl;
			return false;
		}

		return true;
	}

	template<class GlobalEdgeSelect, Integer Dim, Integer ManifoldDim>
	bool test_incomplete(
		Mesh<Dim, ManifoldDim> &mesh,
		const bool use_uniform_refinement = false)
	{
		Map map(0, 1);
		map.resize(mesh.n_nodes(), 0);
		auto edge_select = std::make_shared<GlobalEdgeSelect>(map);
		return test_incomplete(mesh, map, edge_select, use_uniform_refinement);
	}

	template<Integer Dim, Integer ManifoldDim>
	bool test_incomplete_with_edge_rank(
		Mesh<Dim, ManifoldDim> &mesh,
		const Integer n_tests,
		const bool use_uniform_refinement = false,
		const bool online_update = true,
		const std::shared_ptr<NodeRank> &node_rank = std::make_shared<NodeRank>())
	{
		Map map(0, 1);
		map.resize(mesh.n_nodes(), 0);
		map.identity();

		Quality<Dim, ManifoldDim> q(mesh);
		q.compute();

		// auto edge_select = std::make_shared<RankedEdgeSelect<Dim, ManifoldDim>>(map, online_update);
		auto edge_select = std::make_shared<OldestEdgeSelect<Dim, ManifoldDim>>(map, node_rank);
		
		for(Integer i = 0; i < n_tests; ++i) {
			std::cout << "test_incomplete : " << (i+1) << "/" << n_tests << std::endl; 
			const bool ok = test_incomplete(mesh, map, edge_select, use_uniform_refinement);
			if(!ok) {
				assert(false);
				return false;
			}
			
			mesh.clean_up();
			q.compute();
			std::cout << "n_active_elements: " << mesh.n_active_elements() << "/" << mesh.n_elements() << std::endl;
		}

		// edge_select->describe(std::cout);

		q.save_csv("edge_rank", std::to_string(ManifoldDim) + "D_er.csv", true);
		q.save_report(std::to_string(ManifoldDim) + "D_er.svg");
		return true;
	}


	template<Integer Dim, Integer ManifoldDim>
	void test_incomplete_ND(
		Mesh<Dim, ManifoldDim> &mesh,
		const Integer n_tests,
		const bool use_edge_rank)
	{
		using NVES = mars::GlobalNewestVertexEdgeSelect<Dim, ManifoldDim>;
		using LEES = mars::GloballyUniqueLongestEdgeSelect<Dim, ManifoldDim>;

		std::cout << "======================================\n";
		mesh.renumber_nodes();
		mark_boundary(mesh);

		Integer n_serial_ref = 3;
		Bisection<Dim, ManifoldDim> b(mesh);
		b.uniform_refine(n_serial_ref);

		Map map(0, 1);
		map.resize(mesh.n_nodes(), 0);

		auto edge_select = std::make_shared<LEES>(map);
		auto node_rank   = std::make_shared<NodeRank>();
		edge_select->set_node_rank(node_rank);

		if(use_edge_rank) {
			test_incomplete_with_edge_rank(mesh, n_tests, false, true);
		} else {
			Quality<Dim, ManifoldDim> q(mesh);
			q.compute();

			for(Integer i = 0; i < n_tests; ++i) {
				std::cout << "test_incomplete : " << (i+1) << "/" << n_tests << std::endl; 
				
				if(!test_incomplete<LEES>(mesh, map, edge_select, false)) {
					std::cout << "using edge_rank" << std::endl;
					if(!test_incomplete_with_edge_rank(mesh, 1, false, true, node_rank)) {
						assert(false);
						std::cout << "edge_rank failed" << std::endl;
					}
				}

				mesh.clean_up();
				std::cout << "n_active_elements: " << mesh.n_active_elements() << "/" << mesh.n_elements() << std::endl;
				q.compute();
			}

			q.save_csv("ref", std::to_string(ManifoldDim) + "D.csv", true);
			q.save_report(std::to_string(ManifoldDim) + "D.svg");
		}
	}


	template<Integer Dim, Integer ManifoldDim>
	void test_bisection(
		const Integer n_levels,
		std::vector<std::shared_ptr<MeshPartition<Dim, ManifoldDim>>> &parts,
		const bool uniform_refine = false)
	{
		ParBisection<Dim, ManifoldDim> b(parts);
		// auto edge_select = std::make_shared<NewestVertexEdgeSelect<Dim, ManifoldDim>>();
		// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<Dim, ManifoldDim>>(true, false);
		// auto edge_select = std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>(true, true);
		auto edge_select = std::make_shared<UniqueLongestEdgeSelect<Dim, ManifoldDim>>();
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
	b.uniform_refine(1);
	b.set_edge_select(std::make_shared<UniqueLongestEdgeSelect<3, 3>>());

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

	test_bisection(5, parts, false);

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

	test_bisection(3, parts, false);

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

void test_incomplete_2D()
{
	using namespace mars;

	Mesh<2, 2> mesh(true);
	read_mesh("../data/square_2.MFEM", mesh);
	// read_mesh("../data/square_2_def.MFEM", mesh);
	test_incomplete_ND(mesh, 8, true);
}

void test_incomplete_3D()
{
	using namespace mars;

	Mesh<3, 3> mesh(true);
	read_mesh("../data/cube_6.MFEM", mesh);
	// test_incomplete_ND(mesh, 10, false);
	test_incomplete_ND(mesh, 10, true);
}

void test_incomplete_4D()
{
	using namespace mars;
	
	Mesh<4, 4> mesh(true);
	read_mesh("../data/cube4d_24.MFEM", mesh);
	test_incomplete_ND(mesh, 8, false);
}

void test_incomplete_5D()
{
	using namespace mars;
	using NVES = mars::GlobalNewestVertexEdgeSelect<5, 5>;
	using LEES = mars::GloballyUniqueLongestEdgeSelect<5, 5>;

	std::cout << "======================================\n";
	Mesh<5, 5> mesh(true);
	read_mesh("../data/hexateron_1.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<5, 5> q(mesh);
	q.compute();
	mark_boundary(mesh);

	Bisection<5, 5> b(mesh);
	b.uniform_refine(2);

	Integer n_tests = 5;
	for(Integer i = 0; i < n_tests; ++i) {
		std::cout << "-----------------\n";
		std::cout << "test_incomplete : " << (i+1) << "/" << n_tests << std::endl; 
		test_incomplete<NVES>(mesh, true);
		q.compute();
		std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
		std::cout << "-----------------\n";
	}

	q.save_csv("glob_long_edge", "5D.csv", true);
	q.save_report("5D.svg");

	std::ofstream m_os("mesh5.MFEM");
	export_mesh(mesh, m_os);
	m_os.close();
}


void test_incomplete_6D()
{
	using namespace mars;
	using NVES = mars::GlobalNewestVertexEdgeSelect<6, 6>;
	using LEES = mars::GloballyUniqueLongestEdgeSelect<6, 6>;

	std::cout << "======================================\n";
	Mesh<6, 6> mesh(true);
	read_mesh("../data/uniform_polypeton_1.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<6, 6> q(mesh);
	q.compute();
	mark_boundary(mesh);

	Bisection<6, 6> b(mesh);
	b.uniform_refine(2);

	Integer n_tests = 5;
	for(Integer i = 0; i < n_tests; ++i) {
		std::cout << "-----------------\n";
		std::cout << "test_incomplete : " << (i+1) << "/" << n_tests << std::endl; 
		test_incomplete<NVES>(mesh, true);
		q.compute();
		std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
		std::cout << "-----------------\n";
	}

	q.save_csv("glob_long_edge", "6D.csv", true);
	q.save_report("6D.svg");

	std::ofstream m_os("mesh6.MFEM");
	export_mesh(mesh, m_os);
	m_os.close();
	// mesh.describe(std::cout);
}

void test_incomplete_bad_4D()
{
	using namespace mars;
	using NVES = mars::GlobalNewestVertexEdgeSelect<4, 4>;
	using LEES = mars::GloballyUniqueLongestEdgeSelect<4, 4>;

	std::cout << "======================================\n";
	Mesh<4, 4> mesh(true);
	// read_mesh("../data/bad_mesh_p_wn.MFEM", mesh);
	read_mesh("../data/big_mesh_2.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<4, 4> q(mesh);
	q.compute();
	mark_boundary(mesh);

	// mesh.describe(std::cout);

	Bisection<4, 4> b(mesh);
	// b.uniform_refine(2);

	Integer n_tests = 4;
	for(Integer i = 0; i < n_tests; ++i) {
		std::cout << "-----------------\n";
		std::cout << "test_incomplete : " << (i+1) << "/" << n_tests << std::endl; 
		test_incomplete<LEES>(mesh, true);
		q.compute();
		std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
		std::cout << "-----------------\n";
	}

	q.save_csv("glob_long_edge_big", "4D_b.csv", true);
	q.save_report("glob_long_edge_big.svg");

	std::ofstream m_os("big_mesh_3.MFEM");
	export_mesh(mesh, m_os);
	m_os.close();

	std::ofstream os("bad_mesh.MFEM");
	export_elems_with_bad_tags(mesh, os);
	os.close();

	std::ofstream os_p("bad_mesh_p.MFEM");
	export_elems_with_bad_tags(mesh, os_p, true);
	os_p.close();
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
	// test_partition_4D();
	// test_incomplete_2D();
	// test_incomplete_3D();
	test_incomplete_4D();
	// test_incomplete_5D();
	// test_incomplete_6D();
	// test_incomplete_bad_4D();
	// run_tests();
	return EXIT_SUCCESS;
}
