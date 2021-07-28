#ifndef MARS_MESH_IO_HPP
#define MARS_MESH_IO_HPP
#include "adios2.h"
#include "mars_mesh_writer.hpp"

/**
 * Class for handling Reading and writing.
 **/
namespace mars {

    template <class Mesh>
    class Mesh_IO {
    public:
        Mesh_IO() : writer_(mesh, io_){};
        ~Mesh_IO(){};

        void write(const std::string fileName) {
            writer_.open(fileName);
            writer_.generate_data_cube(3);
            writer_.write();
            writer_.close();
        };

        void read(const std::string inputFile) {
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
                // }
            }
        }

    private:
        Mesh mesh;
        adios2::ADIOS adios_;
        adios2::IO io_ = adios_.DeclareIO("SimulationOutput");
        MeshWriter<Mesh> writer_;
    };
}  // namespace mars

#endif