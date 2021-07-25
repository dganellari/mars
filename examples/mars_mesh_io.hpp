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

    private:
        Mesh mesh;
        adios2::ADIOS adios_;
        adios2::IO io_ = adios_.DeclareIO("SimulationOutput");
        MeshWriter<Mesh> writer_;
    };
}  // namespace mars

#endif