#ifndef MARS_MESH_WRITER_HPP
#define MARS_MESH_WRITER_HPP

#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <iostream>
#include <vector>
#include "adios2.h"

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
#include "mars_env.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_longest_edge.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_generation.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_prelepp_benchmark.hpp"
#include "mars_quality.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_simplex.hpp"
#include "mars_test.hpp"
#include "mars_utils.hpp"
#include "mars_vtk_writer.hpp"

class MeshWriter {
public:
    static const int Dim;
    static const int ManifoldDim;

    mars::Mesh<3, 3> mesh;

    MeshWriter(adios2::IO io);
    std::string VTKSchema();
    void open(const std::string& fname);
    void write(int step);
    void close();

protected:
    adios2::IO io;
    adios2::Engine engine;
};
#endif