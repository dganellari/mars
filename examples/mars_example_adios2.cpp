#include <adios2.h>
#include <iostream>
#include "cxxopts.hpp"
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_spacetime_ex.hpp"

/**
 * Run the writing operation of an imageusing adios2.
 *
 * @param fileName name of .bp file we want to create.
 **/
void write_image(const std::string fileName) {
    Settings settings(fileName);
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

/**
 * Run the writing operatino of a mesh using adios2.
 *
 * @param fileName name of .bp file we want to create.
 **/
void write_mesh(const std::string fileName) {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");
    // Switch to figure out which mesh is given
    mars::ParallelMesh2 parMesh;

    MeshWriter<mars::ParallelMesh2> writer(parMesh, io_main);
    writer.open(fileName);
    writer.generate_data_cube(2);
    writer.write();
    writer.close();
}

/**
 * Given a .bp file of a mesh we read the adios variables and
 * attributes defined for the mesh. Ultimately print out the
 * values to terminal.
 *
 * @param name of input file
 **/

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

/**
 * Given a .bp file of an image we read the adios variables and
 * attributes defined for the image. Ultimately print out the
 * values to terminal.
 *
 * @param name of input file
 **/
void read_image(const std::string inputFile) {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationInput");
    adios2::Engine reader;
    reader = io_main.Open(inputFile, adios2::Mode::Read);
    std::vector<double> data;

    reader.BeginStep();
    reader.Get("U", data.data(), adios2::Mode::Deferred);
    reader.EndStep();
    reader.Close();
}

/**
 * Given the arguments of the command line execution
 * choose which method to run. Read/Write, Mesh/Image.
 *
 **/
void run(cxxopts::ParseResult &args) {
    bool values[4] = {
        args["image"].as<bool>(), args["mesh"].as<bool>(), args["write"].as<bool>(), args["read"].as<bool>()};
    std::string fileName = args["file"].as<std::string>();

    if (values[0] && values[2]) {
        fileName = fileName + "-image.bp";
        write_image(fileName);
    }

    if (values[1] && values[2]) {
        fileName = fileName + "-mesh.bp";
        Kokkos::initialize();
        write_mesh(fileName);
        Kokkos::finalize();
    }
}

int main(int argc, char *argv[]) {
    using namespace mars;
    using namespace cxxopts;

    Options options("./adios_example", "Run M.A.R.S with Adios2");
    options.add_options()("i,image", "Write/Read Image")                               // bool param
        ("m,mesh", "Write/Read Mesh")("r,read", "Read Mode")("w,write", "Write Mode")  // bool param
        ("f,file", "File Name", value<std::string>());

    auto result = options.parse(argc, argv);

    run(result);
}