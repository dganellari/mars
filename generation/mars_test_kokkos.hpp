
#include "generation/mars_mesh_kokkos.hpp"
#include "generation/mars_mesh_generation_kokkos.hpp"
#include "generation/mars_utils_kokkos.hpp"

void test_mars_mesh_generation_kokkos_1D(const int level){

	using namespace mars;


	Kokkos::initialize();
		{
			Kokkos::Timer timer;

			generation::kokkos::Parallel_Mesh<1, 1> pMesh;
			generation::kokkos::generate_cube(pMesh, level, 0, 0);

			//fence();

			double time = timer.seconds();

			std::cout<< "Generation took: "<<time<<" seconds."<<std::endl;

			if(level<100){

				Mesh<1,1> sMesh;
				generation::convertParallelMeshToSerial(sMesh,pMesh);

				std::cout << "n_active_elements: " << sMesh.n_active_elements() << std::endl;
				std::cout << "n_nodes: " << sMesh.n_nodes() << std::endl;

				VTKMeshWriter<Mesh1> w;
				w.write("build_line_parallel" + to_string(level) + ".vtu", sMesh);
			}
		}

		Kokkos::finalize();
}



