#ifndef MARS_MESH_IO_HPP
#define MARS_MESH_IO_HPP
#include "adios2.h"
#include "mars_mesh_writer.hpp"

/**
 * Class for handling Reading and writing.
 **/
class Mesh_IO {
public:
    Mesh_IO();
    ~Mesh_IO();

private:
    adios2::IO io_;
    adios2::Engine engine_;
    MeshWriter<class Mesh> writer;
};

#endif