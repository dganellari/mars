#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "mars_simplex.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_mesh.hpp"
#include "mars_bisection.hpp"
#include "mars_vtk_writer.hpp"
#include "mars_quality.hpp"
#include "mars_utils.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_benchmark.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_test.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_longest_edge.hpp"
#include "generation/mars_memory.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include <err.h>


#include "generation/mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "generation/mars_test_kokkos.hpp"
#endif //WITH_KOKKOS

#include "mars_entity.hpp"
#include "mars_connectivity.hpp"
#include "mars_examples.hpp"
#include <Kokkos_Crs.hpp>
// #include <Kokkos_MultiVector.hpp>
// #include <Kokkos_DefaultNode.hpp>
// typedef Kokkos::DefaultNode::DefaultNodeType                 KokkosNode;
// typedef KokkosExamples::DummySparseKernel<Node>         KokkosSparseOps;
// typedef Kokkos::CrsMatrix< float,int,KokkosNode,KokkosSparseOps>     FloatMat;

#ifdef WITH_PAR_MOONOLITH
#include "mars_moonolith_test.hpp"
#include <mpi.h>
#endif //WITH_PAR_MOONOLITH


#ifdef WITH_MPI
#include "mars_par_bisection.hpp"
#include "mars_par_mesh.hpp"
#endif //WITH_MPI
#include <chrono>

using namespace std::chrono;

mars::Mesh1 test_mars_mesh_generation_1D(const int x) {

	using namespace mars;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Mesh1 mesh;
	generate_line(mesh, x);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast < seconds > (t2 - t1).count();

	std::cout << "Generation took: "<< duration<<" seconds."<<std::endl;

	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;


	if (x < 100) {
		VTKMeshWriter<Mesh1> w;
		w.write("build_line" + std::to_string(x) + ".vtu", mesh);
	}

	return mesh;
}

mars::Mesh2 test_mars_mesh_generation_2D(const int x,
		const int y) {

	using namespace mars;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Mesh2 mesh;
	generate_square(mesh, x, y);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast < seconds > (t2 - t1).count();

	std::cout << "Generation took: "<< duration<<" seconds."<<std::endl;

	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;


	if (x <= 1000) {

		VTKMeshWriter<Mesh2> w;
		w.write("build_square" + std::to_string(x) + std::to_string(y) + ".vtu", mesh);
	}

	return mesh;
}

mars::Mesh3 test_mars_mesh_generation_3D(const int x,
		const int y, const int z) {

	using namespace mars;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Mesh3 mesh;
	generate_cube<3, 3>(mesh, x, y, z);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast < seconds > (t2 - t1).count();

	std::cout << "Generation took: "<< duration<<" seconds."<<std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;

	if (z < 102) {
		std::cout<<"build_cube" + std::to_string(x) + std::to_string(y) + ".vtu"<<std::endl;
		VTKMeshWriter<Mesh3> w;
		w.write("build_cube" + std::to_string(x) + std::to_string(y) + ".vtu", mesh);
	}
	return mesh;
}

void test_mars_mesh_generation_unit_cube() {

	using namespace mars;

	Mesh3 mesh = generate_unit_cube();

	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;

	VTKMeshWriter < Mesh3 > w;
	w.write("build_unit_cube_tetrakis.vtu", mesh);
}

void test_read_write_3D(const std::string filename)
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh3 mesh;
	read_mesh(filename, mesh, true);

	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;
	VTKMeshWriter<Mesh3> w;
	w.write("read_write.vtu", mesh);
}

void test_write_3D(const mars::Mesh3 mesh)
{
	using namespace mars;
	std::cout << "======================================\n";
	/*Mesh3 mesh;
	read_mesh(filename, mesh, true);*/

	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;
	VTKMeshWriter<Mesh3> w;
	w.write("write.vtu", mesh);
}

void test_uniform_bisection_3D(const int level, const std::string filename)
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh3 mesh;
	read_mesh(filename, mesh, true);
	mesh.renumber_nodes();
	mesh.reorder_nodes();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	Bisection<Mesh3> b(mesh);
	b.uniform_refine(level); b.clear();
	print_boundary_points(mesh, std::cout, true);

	mesh.clean_up();
	mesh.update_dual_graph();

	Quality<Mesh3> q(mesh);
	q.compute();

	std::cout << "n_boundary_sides: " << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes()<< std::endl;
	VTKMeshWriter<Mesh3> w;
						w.write("cube_bisect_"+ std::to_string(mesh.Dim) +".vtu", mesh);
}

void test_uniform_bisection_3D(const int level, mars::Mesh3 mesh)
{
	using namespace mars;
	std::cout << "======================================\n";
	/*Mesh3 mesh;
	read_mesh(filename, mesh, true);*/

	mesh.renumber_nodes();
	mesh.reorder_nodes();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	Bisection<Mesh3> b(mesh);
	b.uniform_refine(level); b.clear();

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast < seconds > (t2 - t1).count();

	std::cout << "Refinment took: "<< duration<<" seconds."<<std::endl;
	print_boundary_points(mesh, std::cout, true);

	mesh.clean_up();
	mesh.update_dual_graph();

	Quality<Mesh3> q(mesh);
	q.compute();

	std::cout << "n_boundary_sides: " << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes()<< std::endl;
	VTKMeshWriter<Mesh3> w;
						w.write("cube_bisect_"+ std::to_string(mesh.Dim) +".vtu", mesh);
}

void test_uniform_bisection_2D(const int level, mars::Mesh2 mesh)
{
	using namespace mars;
	std::cout << "======================================\n";
	/*Mesh2 mesh;
	read_mesh(filename, mesh, true);*/
	mesh.renumber_nodes();
	mesh.reorder_nodes();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	Bisection<Mesh2> b(mesh);
	b.uniform_refine(level); b.clear();
	print_boundary_points(mesh, std::cout, true);

	mesh.clean_up();
	mesh.update_dual_graph();

	Quality<Mesh2> q(mesh);
	q.compute();

	std::cout << "n_boundary_sides: " << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes()<< std::endl;
	VTKMeshWriter<Mesh2> w;
						w.write("cube_bisect_"+ std::to_string(mesh.Dim) +".vtu", mesh);
}

void test_uniform_bisection_2D(const int level, const std::string filename)
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh<3,2> mesh;
	read_mesh(filename, mesh, true);
	mesh.renumber_nodes();
	mesh.reorder_nodes();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	Bisection<Mesh<3,2>> b(mesh);
	b.uniform_refine(level); b.clear();
	print_boundary_points(mesh, std::cout, true);

	mesh.clean_up();
	mesh.update_dual_graph();

	Quality<Mesh<3,2>> q(mesh);
	q.compute();

	std::cout << "n_boundary_sides: " << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << mesh.n_nodes()<< std::endl;
	VTKMeshWriter<Mesh<3,2>> w;
						w.write("cube_bisect_"+ std::to_string(mesh.Dim) +".vtu", mesh);
}

void test_bisection_2D()
{
	using namespace mars;
	std::cout << "======================================\n";
	Mesh2 mesh;
	read_mesh("../data/square_2.MFEM", mesh);
	// read_mesh("../data/square_2_def.MFEM", mesh);

	Bisection<Mesh2> b(mesh);
	b.uniform_refine(3);
	b.clear();

	mesh.clean_up();
	mesh.reorder_nodes();

	Quality<Mesh2> q(mesh);
	q.compute();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	auto edge_select = std::make_shared<ImplicitOrderEdgeSelect<Mesh2>>();
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
	Mesh3 mesh;
	read_mesh("../data/write/cube.MFEM", mesh, true);
	mesh.renumber_nodes();
	mesh.reorder_nodes();

	mesh.update_dual_graph();
	mark_boundary(mesh);

	Bisection<Mesh3> b(mesh);
	b.uniform_refine(2); b.clear();
	print_boundary_points(mesh, std::cout, true);

	mesh.clean_up();
	mesh.update_dual_graph();

	Quality<Mesh3> q(mesh);
	q.compute();

	std::cout << "n_boundary_sides: " << mesh.n_boundary_sides() << std::endl;
	std::cout << "volume: " << mesh.volume() << std::endl;
	std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

	VTKMeshWriter<Mesh3> w;
	w.write("cube_bisect_"+ std::to_string(mesh.Dim) +".vtu", mesh);
	// auto edge_select = std::make_shared<NewestVertexEdgeSelect<3, 3>>();
	auto edge_select = std::make_shared<LongestEdgeSelect<Mesh3>>(true);

	// auto edge_select = std::make_shared<ImplicitOrderEdgeSelect<3, 3>>();
	// auto edge_select = std::make_shared<UniqueLongestEdgeSelect<3, 3>>(true);
	// auto edge_select = std::make_shared<NewestVertexAndLongestEdgeSelect<3, 3>>();
	b.set_edge_select(edge_select);
	b.uniform_refine(1);

	// mesh.describe(std::cout);

	Integer n_levels = 1;
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
	Mesh4 mesh;
	read_mesh("../data/cube4d_24.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<Mesh4> q(mesh);
	q.compute();
	mark_boundary(mesh);

	std::cout << "volume: " << mesh.volume() << std::endl;
	Bisection<Mesh4> b(mesh);
	b.set_edge_select(std::make_shared<UniqueLongestEdgeSelect<Mesh4>>());
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
		using MeshD = mars::Mesh<Dim, ManifoldDim>;

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

		Bisection<MeshD> b(mesh);
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
		using MeshD = mars::Mesh<Dim, ManifoldDim>;

		Map map(0, 1);
		map.resize(mesh.n_nodes(), 0);
		map.identity();

		Quality<MeshD> q(mesh);
		q.compute();

		// auto edge_select = std::make_shared<RankedEdgeSelect<MeshD>>(map, online_update);
		auto edge_select = std::make_shared<OldestEdgeSelect<MeshD>>(map, node_rank);

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
		using MeshD = mars::Mesh<Dim, ManifoldDim>;
		using NVES = mars::GlobalNewestVertexEdgeSelect<MeshD>;
		using LEES = mars::GloballyUniqueLongestEdgeSelect<MeshD>;

		std::cout << "======================================\n";
		mesh.renumber_nodes();
		mark_boundary(mesh);

		Integer n_serial_ref = 3;
		Bisection<MeshD> b(mesh);
		b.uniform_refine(n_serial_ref);

		Map map(0, 1);
		map.resize(mesh.n_nodes(), 0);

		auto edge_select = std::make_shared<LEES>(map);
		auto node_rank   = std::make_shared<NodeRank>();
		edge_select->set_node_rank(node_rank);

		if(use_edge_rank) {
			test_incomplete_with_edge_rank(mesh, n_tests, false, true);
		} else {
			Quality<MeshD> q(mesh);
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


	template<class Mesh>
	void test_bisection(
		const Integer n_levels,
		std::vector<std::shared_ptr<MeshPartition<Mesh>>> &parts,
		const bool uniform_refine = false)
	{
		static const Integer Dim = Mesh::Dim;
		using Point = typename Mesh::Point;

		PartitionedBisection<Mesh> b(parts);
		auto edge_select = std::make_shared<UniqueLongestEdgeSelect<Mesh>>();
		edge_select->set_recursive(true);
		b.set_edge_select(edge_select);
		std::vector<std::vector<mars::Integer>> elements(parts.size());

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
	Mesh2 mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/square_2_def.MFEM", mesh);
	mark_boundary(mesh);

	Bisection<Mesh2> b(mesh);
	b.uniform_refine(3);

	std::vector<Integer> partitioning(mesh.n_elements());

	Integer n_parts = 7;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		partitioning[i] = i % n_parts;
	}

	std::vector<std::shared_ptr<MeshPartition<Mesh2>>> parts;
	partition_mesh(mesh, n_parts, partitioning, parts);

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
	using Mesh = mars::Mesh<3, 3>;

	std::cout << "======================================\n";
	Mesh mesh;
	// read_mesh("../data/square_2.MFEM", mesh);
	read_mesh("../data/cube_6.MFEM", mesh, true);
	mark_boundary(mesh);

	Bisection<Mesh> b(mesh);
	b.uniform_refine(1);
	b.set_edge_select(std::make_shared<UniqueLongestEdgeSelect<Mesh>>());

	std::vector<Integer> partitioning(mesh.n_elements(), 0);

	Integer n_parts = mesh.n_active_elements();
	Integer element_index = 0;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		if(mesh.is_active(i)) {
			partitioning[i] = (element_index++) % n_parts;
		}
	}

	std::vector<std::shared_ptr<MeshPartition<Mesh>>> parts;
	partition_mesh(mesh, n_parts, partitioning, parts);

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
	using Mesh = mars::Mesh<4, 4>;

	std::cout << "======================================\n";
	Mesh mesh;
	read_mesh("../data/cube4d_24.MFEM", mesh);
	mark_boundary(mesh);

	Bisection<Mesh> b(mesh);
	b.uniform_refine(1);

	std::vector<Integer> partitioning(mesh.n_elements());

	Integer n_parts = mesh.n_active_elements();
	Integer element_index = 0;
	for(Integer i = 0; i < mesh.n_elements(); ++i) {
		if(mesh.is_active(i)) {
			partitioning[i] = (element_index++) % n_parts;
		}
	}

	std::vector<std::shared_ptr<MeshPartition<Mesh>>> parts;
	partition_mesh(mesh, n_parts, partitioning, parts);

	test_bisection(3, parts, false);

	for(const auto &p : parts) {
		std::cout << p->partition_id() << " n_active_elements: " << p->get_mesh().n_active_elements() << std::endl;
		p->get_mesh().update_dual_graph();
		print_boundary_info(p->get_mesh(), true);
	}
}

void run_benchmarks(int level)
{
	using namespace mars;
/*
	Benchmark<Mesh2> b;
	Mesh2 m;
	read_mesh("../data/square_2_def.MFEM", m);

	b.run(level, m, "b");

	LeppBenchmark<Mesh2> lb;
	lb.run(level, m, "lb");*/

	Benchmark<Mesh3> b3;
	Mesh3 m3;
	read_mesh("../data/cube_6.MFEM", m3);

	b3.run(level, m3, "b3");

	LeppBenchmark<Mesh3> lb3;
	lb3.run(level, m3, "lb3");

	/*ParallelMesh2 pMesh;
	generate_square(pMesh, 97, 10);

	Mesh2 sMesh;
	convert_parallel_mesh_to_serial(sMesh, pMesh);

	std::cout << "n_active_elements: " << sMesh.n_active_elements()
			<< std::endl;
	std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

	VTKMeshWriter<Mesh2> w;
	w.write("build_square_parallel.vtu", sMesh);

	Benchmark<Mesh2> b1;
	//b1.run(level, sMesh, "b1");

	Mesh2 Mesh;
	generate_square(Mesh, 97, 10);

	std::cout << "n_active_elements: " << Mesh.n_active_elements() << std::endl;
	std::cout << "n_nodes: " << Mesh.n_nodes() << std::endl;

	VTKMeshWriter<Mesh2> w2;
	w2.write("build_square_serial.vtu", Mesh);

	LeppBenchmark<Mesh2> b2;
	b2.run(level, Mesh, "b2");*/

	/*ParallelMesh3 pMesh3;
	generate_cube(pMesh3, 4, 2, 3);

	Mesh3 sMesh3;
	convert_parallel_mesh_to_serial(sMesh3, pMesh3);

	std::cout << "n_active_elements: " << sMesh3.n_active_elements()
			<< std::endl;
	std::cout << "n_nodes: " << sMesh3.n_nodes() << std::endl;

	VTKMeshWriter<Mesh3> w3;
	w3.write(
			"build_cube_parallel" + std::to_string(1) + std::to_string(1)
					+ ".vtu", sMesh3);

	Benchmark<Mesh3> b3;
	b3.run(level, sMesh3, "b3");

	LeppBenchmark<Mesh3> lb3;
		lb3.run(level, sMesh3, "lb3");*/

/*	Benchmark<Mesh4> b4;
	Mesh4 m4;
	read_mesh("../data/cube4d_24.MFEM", m4);

	b4.run(level, m4, "b4");

	LeppBenchmark<Mesh4> lb4;
	lb4.run(level, m4, "lb4");*/
}

void test_incomplete_2D()
{
	using namespace mars;

	Mesh2 mesh(true);
	read_mesh("../data/square_2.MFEM", mesh);
	// read_mesh("../data/square_2_def.MFEM", mesh);
	test_incomplete_ND(mesh, 8, true);
}

void test_incomplete_3D()
{
	using namespace mars;

	Mesh3 mesh(true);
	read_mesh("../data/cube_6.MFEM", mesh);
	// test_incomplete_ND(mesh, 10, false);
	test_incomplete_ND(mesh, 10, true);
}

void test_incomplete_4D()
{
	using namespace mars;

	Mesh4 mesh(true);
	read_mesh("../data/cube4d_24.MFEM", mesh);
	// test_incomplete_ND(mesh, 8, false);
	test_incomplete_ND(mesh, 8, true);
}

void test_incomplete_5D()
{
	using namespace mars;
	using Mesh = mars::Mesh5;
	using NVES = mars::GlobalNewestVertexEdgeSelect<Mesh>;
	using LEES = mars::GloballyUniqueLongestEdgeSelect<Mesh>;


	std::cout << "======================================\n";
	Mesh mesh(true);
	read_mesh("../data/hexateron_1.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<Mesh> q(mesh);
	q.compute();
	mark_boundary(mesh);

	Bisection<Mesh> b(mesh);
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
	using Mesh = mars::Mesh6;
	using NVES = mars::GlobalNewestVertexEdgeSelect<Mesh6>;
	using LEES = mars::GloballyUniqueLongestEdgeSelect<Mesh6>;


	std::cout << "======================================\n";
	Mesh mesh(true);
	read_mesh("../data/uniform_polypeton_1.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<Mesh> q(mesh);
	q.compute();
	mark_boundary(mesh);

	Bisection<Mesh> b(mesh);
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
	using Mesh = mars::Mesh4;
	using NVES = mars::GlobalNewestVertexEdgeSelect<Mesh>;
	using LEES = mars::GloballyUniqueLongestEdgeSelect<Mesh>;

	std::cout << "======================================\n";
	Mesh mesh(true);
	// read_mesh("../data/bad_mesh_p_wn.MFEM", mesh);
	read_mesh("../data/big_mesh_2.MFEM", mesh);
	mesh.renumber_nodes();

	Quality<Mesh> q(mesh);
	q.compute();
	mark_boundary(mesh);

	// mesh.describe(std::cout);

	Bisection<Mesh> b(mesh);
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

void read_file()
{
	using namespace mars;
	std::cout << "======================================\n";
	mars::Mesh<3, 3> mesh;
	// mars::Mesh<2, 2> mesh;
	read("../data/test_3D.msh", mesh, true);
	// read("../data/square_2.MFEM", mesh);
	// read_msh("../data/test_3D.msh", mesh, true);
	mesh.describe(std::cout);
	std::cout << mesh.n_elements() << std::endl; 
	std::cout << mesh.n_nodes() << std::endl; 
}


void write_file()
{
	using namespace mars;
	std::cout << "======================================\n";
	mars::Mesh<3, 3> mesh;
	// mars::Mesh<2, 2> mesh;
	read("../data/test_3D.msh", mesh, true);
	// read("../data/square_2.MFEM", mesh, true);
	mesh.describe(std::cout);
	write("../data/test_write_tetra3.msh", mesh);
	// write("../data/test_write_square2.MFEM", mesh); 

}

void par_mesh_test()
{
#ifdef WITH_MPI
	using namespace mars;

	using ParMesh2 = mars::ParMesh<2, 2>;

	Mesh2 serial_mesh(true);
	read_mesh("../data/square_2.MFEM", serial_mesh);
	mark_boundary(serial_mesh);

	std::vector<Integer> partitioning = {0, 1};

	Communicator world; //world.set_verbose(true);
	ParMesh2 mesh(world);
	mesh.init(serial_mesh, partitioning);

	assert(mesh.is_conforming());

	// ParBisection<ParMesh2> b(mesh);
	// b.uniform_refine(1);

#endif //WITH_MPI
}


using std::cout;
using std::endl;
int main(int argc, char *argv[])
{
	using namespace mars;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
#else
#ifdef WITH_PAR_MOONOLITH
	MPI_Init(&argc, &argv);
#endif
#endif 

	// test_bisection_2D();
	 //test_bisection_3D(atoi(argv[1]));
	// test_bisection_4D();



// functionspaces_example2D();
// assembly_example();
// boundary_example();

// simplex_reference_normals_3D();

// normals_example2D();
// normals_example3D();
//normals_example4D();
//connectivity_example4D();

	//run_benchmarks();
	//run_benchmarks(atoi(argv[1]));


	// test_partition_2D();
	// test_partition_3D();
	// test_partition_4D();
	// test_incomplete_2D();
	// test_incomplete_3D();
	// test_incomplete_4D();
	// test_incomplete_5D();
	// test_incomplete_6D();
	// test_incomplete_bad_4D();
	// read_file();
	// write_file(); 

	// int level = 1;
	// std::string filename = "../data/write/tetrakis.MFEM";
	// if (argc > 1) {
	// 	char *end_ptr = argv[1];
	// 	level = strtol(argv[1], &end_ptr, 10);
	// 	if (*end_ptr != '\0' || end_ptr == argv[1])
	// 		warnx("'%s' could not be (completely) converted to long", argv[1]);
	// } else
	// 	std::cout
	// 			<< "No level of refinement was specified. Setting the default to 1!"
	// 			<< std::endl;



	/*if (argc > 2) {
		filename = argv[2];
	} else
		std::cout
				<< "No file name was specified. Setting the default to tetrakis!"
				<< std::endl;*/

	//run_tests(level,filename);

	//test_uniform_bisection_2D(level,filename);
	//test_read_write_3D(filename);


	// if(level <100){

		// test_mars_mesh_generation_2D(level,level);

		//test_mars_mesh_generation_1D(level);
	// }

	/*test_mars_mesh_generation_3D(100,100,100);
	test_mars_mesh_generation_3D(78,100,80);*/

	//test_mars_mesh_generation_3D(150,150,150);
	//test_mars_mesh_generation_3D(200,200,200);

	//equivalent 3D generation using refinement and libmesh like mesh generation technique.

	/*test_uniform_bisection_3D(3, test_mars_mesh_generation_3D(1,1,1));*/
	//parallel with kokkos.
#ifdef WITH_KOKKOS
	Kokkos::initialize(argc,argv);
	{

Integer M=4;
Integer N=4;
Kokkos::View<Real **[3][4]> aaa("aaa",M,N);
assembly_example();
		// run_benchmarks(level);

		//test_mars_mesh_generation_kokkos_2D(2,4);

		//test_mars_mesh_generation_kokkos_2D(level + 4,level);
		//test_mars_mesh_generation_kokkos_3D(level,level,level);
		//test_mars_mesh_generation_kokkos_1D(level);
	}

	Kokkos::finalize();

#endif


#ifdef WITH_PAR_MOONOLITH
	run_mars_moonolith_test();
#endif //WITH_PAR_MOONOLITH

	// if (level < 100) {
		//test_mars_mesh_generation_3D(level , level, level );

		//test_mars_mesh_generation_2D(level,level);
		//test_mars_mesh_generation_1D(level);
	// }


#ifdef WITH_MPI
	// par_mesh_test();
	return MPI_Finalize();
#else
	#ifdef WITH_PAR_MOONOLITH
		return MPI_Finalize();
	#else 
		return 0;
	#endif
#endif //WITH_MPI
}
