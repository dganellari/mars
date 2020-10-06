#include "mars_mesh_kokkos.hpp"

void test_mars_mesh_generation_kokkos_1D(const int level) {

	using namespace mars;

	Kokkos::Timer timer;

	ParallelMesh1 pMesh;
	generate_cube(pMesh, level, 0, 0);

	Kokkos::fence();

	double time = timer.seconds();

	std::cout << "Generation 1D Kokkos took: " << time << " seconds." << std::endl;

	if (level < 100) {

		Mesh1 sMesh;
		convert_parallel_mesh_to_serial(sMesh, pMesh);

		std::cout << "n_active_elements: " << sMesh.n_active_elements()
				<< std::endl;
		std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

		VTKMeshWriter<Mesh1> w;
		w.write("build_line_parallel" + std::to_string(level) + ".vtu", sMesh);
	}

}

void test_mars_mesh_generation_kokkos_2D(const int x, const int y) {

	using namespace mars;

	/*Kokkos::initialize();
	 {*/
	Kokkos::Timer timer;

	ParallelMesh2 pMesh;
	generate_cube(pMesh, x, y, 0);

	//Kokkos::fence();

	double time = timer.seconds();

	std::cout << "Generation 2D kokkos took: " << time << " seconds." << std::endl;

	if (x < 100) {

		Mesh2 sMesh;
		convert_parallel_mesh_to_serial(sMesh, pMesh);

		std::cout << "n_active_elements: " << sMesh.n_active_elements()
				<< std::endl;
		std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

		VTKMeshWriter<Mesh2> w;
		w.write(
				"build_square_parallel" + std::to_string(x) + std::to_string(y)
						+ ".vtu", sMesh);

		sMesh.repair(true);
	}
	/*}

	 Kokkos::finalize();*/
}


void test_mars_nonsimplex_mesh_generation_kokkos_2D(const int x, const int y) {

	using namespace mars;


	Kokkos::Timer timer;

	ParallelQuad4Mesh pMesh;
	generate_cube(pMesh, x, y, 0);

	double time = timer.seconds();

	std::cout << "Generation 2D kokkos took: " << time << " seconds." << std::endl;

	std::cout << "n_active_elements pmesh: " << pMesh.n_active_elements(pMesh.get_view_elements().extent(0))
			<< std::endl;
	if (x < 100) {

		Quad4_Mesh sMesh;
		convert_parallel_mesh_to_serial(sMesh, pMesh);

		std::cout << "n_active_elements: " << sMesh.n_active_elements()
				<< std::endl;
		std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

		VTKMeshWriter<Quad4_Mesh> w;
		w.write(
				"build_quad4_parallel" + std::to_string(x) + std::to_string(y)
						+ ".vtu", sMesh);
	}
	
}

void test_mars_nonsimplex_mesh_generation_kokkos_3D(const int x, const int y, const int z) {

	using namespace mars;

	Kokkos::Timer timer;

	ParallelHex8Mesh pMesh;
	generate_cube(pMesh, x, y, z);

	double time = timer.seconds();

	std::cout << "Generation 3D kokkos took: " << time << " seconds." << std::endl;
	std::cout << "n_active_elements pmesh: " << pMesh.n_active_elements(pMesh.get_view_elements().extent(0))
				<< std::endl;
				
	if (x < 100) {

		Hex8_Mesh sMesh;
		convert_parallel_mesh_to_serial(sMesh, pMesh);

		std::cout << "n_active_elements: " << sMesh.n_active_elements()
				<< std::endl;
		std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

		VTKMeshWriter<Hex8_Mesh> w;
		w.write(
				"build_hex8_parallel" + std::to_string(x) + std::to_string(y)
						+ ".vtu", sMesh);
	}
	
}
void test_mars_mesh_generation_kokkos_3D(const int x, const int y, const int z) {

	using namespace mars;

	Kokkos::Timer timer;

	ParallelMesh3 pMesh;
	generate_cube(pMesh, x, y, z);

	//Kokkos::fence();

	double time = timer.seconds();

	std::cout << "Generation 3D kokkos took: " << time << " seconds." << std::endl;

	if (x < 100) {

		Mesh3 sMesh;
		convert_parallel_mesh_to_serial(sMesh, pMesh);

		std::cout << "n_active_elements: " << sMesh.n_active_elements()
				<< std::endl;
		std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

		std::cout<<"build_cube_parallel" + std::to_string(x) + std::to_string(y) + ".vtu"<<std::endl;
		VTKMeshWriter<Mesh3> w;
		w.write(
				"build_cube_parallel" + std::to_string(x) + std::to_string(y)
						+ ".vtu", sMesh);

		sMesh.repair(true);

	}
}
