#include <adios2.h>
#include <iostream>
#include "mars_config.hpp"
#ifdef WITH_CXXOPTS
#include "cxxopts.hpp"
#include "mars_base.hpp"
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_io.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_spacetime_ex.hpp"

// void testReadWriteImage(const std::vector read, const std::vector write) {}

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
    main_image.new_data(3000, 4000, 4);
    // Open the writer, (which now is in write mode), with the settings found.
    main_image.open(settings.output);
    // Write with the following steps
    main_image.write(1);
    // Close writer
    main_image.close();
}

// /**
//  * Run the writing operatino of a mesh using adios2.
//  *
//  * @param fileName name of .bp file we want to create.
//  **/
// void write_mesh(const std::string fileName) {
//     adios2::ADIOS adios(adios2::DebugON);
//     adios2::IO io_main = adios.DeclareIO("SimulationOutput");
//     // Switch to figure out which mesh is given
//     mars::ParallelMesh2 parMesh;

//     MeshWriter<mars::ParallelMesh2> writer(parMesh, io_main);
//     writer.open(fileName);
//     writer.generate_data_cube(2);
//     writer.write();
//     writer.close();
// }

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

    reader.BeginStep();

    adios2::Variable<double> uVar = io_main.InquireVariable<double>("U");
    adios2::Attribute<std::string> uAttr = io_main.InquireAttribute<std::string>("format/mars_mesh");
    adios2::Variable<uint32_t> numOfElements = io_main.InquireVariable<uint32_t>("NumOfElements");
    if (uVar) {
        std::vector<double> example;
        std::cout << "Got it\n";
        reader.Get(uVar, example, adios2::Mode::Deferred);
        std::cout << example.size() << std::endl;
    }

    // TODO: reader not finding the Get methods for some reason.
    // if (numOfElements) {
    //     uint32_t numOfLocalElements;
    //     reader.Get(numOfLocalElements, numOfElements, adios2::Mode::Deferred);
    //     // std::cout << numOfLocalElements;
    // }
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
    adios2::Variable<double> uVar = io_main.InquireVariable<double>("U");
    if (uVar)  // it exists
    {
        size_t n = 1;
        std::cout << "Got it ";
        for (auto i : uVar.Shape()) {
            n *= i;
        }
        data.resize(n);
        reader.Get(uVar, data.data());
    }
    reader.EndStep();
    reader.Close();

    for (auto i = data.begin(); i != data.end(); ++i) {
        std::cout << *i << '\n';
    }
}

/**
 * Given the arguments of the command line execution
 * choose which method to run. Read/Write, Mesh/Image.
 *
 **/
void run(cxxopts::ParseResult &args) {
    const std::string help =
        "Usage is as follows: \n"
        "        -i : Create an image.\n"
        "        -m : Create a mesh.\n"
        "        -w : Write the object.\n"
        "        -r : Read the .bp file.\n"
        "        -f : Name of the file ";

    bool values[4] = {
        args["image"].as<bool>(), args["mesh"].as<bool>(), args["write"].as<bool>(), args["read"].as<bool>()};
    std::string fileName = args["file"].as<std::string>();
    if (args["help"].as<bool>()) {
        std::cout << help;
    }

    if (values[0] && values[2]) {
        fileName = fileName + "-image.bp";
        write_image(fileName);
    }

    if (values[1] && values[2]) {
        fileName = fileName + "-mesh.bp";
        mars::Mesh_IO<mars::ParallelMesh3> io;
        io.write(fileName);
    }

    if (values[0] && values[3]) {
        read_image(fileName);
    }

    if (values[1] && values[3]) {
        read_image(fileName);
    }
}

int main(int argc, char *argv[]) {
    using namespace mars;
    using namespace cxxopts;

    Env env(argc, argv);

    Options options("./adios_example", "Run M.A.R.S with Adios2");
    options.add_options()("i,image", "Write/Read Image")  // bool param
        ("m,mesh", "Write/Read Mesh")("r,read", "Read Mode")("w,write", "Write Mode")("h,help",
                                                                                      "Print usage")  // bool param
        ("f,file", "File Name", value<std::string>()->default_value("solution"));

    auto result = options.parse(argc, argv);

    run(result);
}
#endif