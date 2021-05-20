#include <adios2.h>
#include <iostream>
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_mesh_writer.hpp"

// pass to interpolate.
// class ST3Analitcal {
// public:
//     MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_exact(p); }
// };

// class ST3RHS {
// public:
//     MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_spacetime(p); }
// };

// template <class Mesh>
// class ST3BC {
// public:
//     /* BC --> zero dirichlet + natural neumann on upper bound */
//     static const int Dim = Mesh::Dim;

//     MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
//         if (is_boundary(p)) {
//             val = ex3_st_exact(p);
//         }
//     }

//     MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
//         bool ret = false;
//         for (int d = 0; d < Dim; ++d) {
//             if (p[d] <= 1e-14) {
//                 ret = true;
//                 break;
//             }

//             if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
//                 ret = true;
//                 break;
//             }
//         }

//         return ret;
//     }
// };

/*
 * Write a structured image.
 */

void write_image() {
    Settings settings("solutionImage.bp");
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");

    ImageWriter main_image(settings, io_main);

    // Create example data vector
    main_image.new_data(30, 40, 2);
    // Open the writer, (which now is in write mode), with the settings found.
    main_image.open(settings.output);
    // Write with the following steps
    main_image.write(1);
    // Close writer
    main_image.close();
}

void write_mesh() {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");
    mars::ParallelQuad4Mesh parMesh;

    MeshWriter<mars::ParallelQuad4Mesh> writer(parMesh, io_main);
    writer.open("example.bp");
    writer.generate_data_cube();
    // writer.interpolate()
    writer.write();
    writer.close();
}

int main(int argc, char *argv[]) {
// write_image();
#ifdef WITH_KOKKOS
    Kokkos::initialize();
    write_mesh();
    Kokkos::finalize();
#endif
}