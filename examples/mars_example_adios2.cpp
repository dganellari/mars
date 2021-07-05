#include <adios2.h>
#include <iostream>
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_spacetime_ex.hpp"

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

void write_mesh(const std::string fileName) {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");
    // Switch to figure out which mesh is given
    mars::ParallelMesh2 parMesh;

    MeshWriter<mars::ParallelMesh2> writer(parMesh, io_main);
    writer.open(fileName);
    writer.generate_data_cube(2);
    // writer.interpolate()
    writer.write();
    writer.close();
}

void read_mesh(const std::string inputFile) {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationInput");
    adios2::Engine reader;
    reader = io_main.Open(inputFile, adios2::Mode::Read);
    std::vector<double> data;

    // adios2::Variable<double> &data = io_main.DefineVariable<double>("U", {}, {}, {n_nodes});

    reader.BeginStep();

    reader.Get("U", data.data(), adios2::Mode::Deferred);
    reader.EndStep();
    reader.Close();
}
void read_image(const std::string inputFile) {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationInput");
    adios2::Engine reader;
    reader = io_main.Open(inputFile, adios2::Mode::Read);
    std::vector<double> data;

    reader.BeginStep();

    reader.EndStep();
}
int main(int argc, char *argv[]) {
    // write_image();
    // std::cout << type;
#ifdef WITH_KOKKOS
    // mars::ParallelMesh2 parMesh;
    // mars::ParallelMesh3 parMesh3;
    // mars::ParallelQuad4Mesh parQuadMesh;
    Kokkos::initialize();
    // write_mesh("example.bp");
    Kokkos::finalize();
#endif
    read_image(argv[1]);
}