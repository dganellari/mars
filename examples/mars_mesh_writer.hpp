#ifndef MARS_MESH_WRITER_HPP
#define MARS_MESH_WRITER_HPP

#include <iostream>
#include <vector>
#include "adios2.h"
// #include "mars_mesh_kokkos.hpp"

class MeshWriter {
public:
    MeshWriter(adios2::IO io, const std::string engineType);
    std::string VTKSchema();
    void open(const std::string& fname);
    void write(int step);
    void close();

protected:
    adios2::IO io;
    adios2::Engine engine;
};
#endif