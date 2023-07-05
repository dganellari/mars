#include <iostream>

#include "mars_base.hpp"
#include "mars_err.hpp"

#include "mars.hpp"
#include "mars_env.hpp"

#ifdef MARS_ENABLE_AMR_BACKEND
#include "mars_kokkos_generate_cube.hpp"
#endif

// void test_mars_mesh_generation(const int x) {
//     using namespace mars;

//     high_resolution_clock::time_point t1 = high_resolution_clock::now();

//     Mesh1 mesh;
//     generate_line(mesh, x);

//     high_resolution_clock::time_point t2 = high_resolution_clock::now();
//     auto duration = duration_cast<seconds>(t2 - t1).count();

//     std::cout << "Generation took: " << duration << " seconds." << std::endl;

//     std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
//     std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;

//     if (x < 100) {
//         VTKMeshWriter<Mesh1> w;
//         w.write("build_line" + std::to_string(x) + ".vtu", mesh);
//     }

//     return mesh;
// }

#ifdef MARS_ENABLE_AMR_BACKEND
void test_mars_mesh_generation(const int x) {
#ifdef MARS_ENABLE_KOKKOS
    Kokkos::Timer timer;

    mars::ParallelMesh3 pMesh;
    mars::generate_cube(pMesh, x, x, x);

    double time = timer.seconds();

    std::cout << "Generation 3D kokkos took: " << time << " seconds." << std::endl;
#endif  // MARS_ENABLE_KOKKOS
}

int main(int argc, char *argv[]) {
    mars::Env env(argc, argv);
    test_mars_mesh_generation(100);
    return env.exit_code();
}
#endif
