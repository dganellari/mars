#include <cmath>
#include <complex>
#include <iostream>
#include "adios2.h"
#include "cxxopts.hpp"
#include "mars_globals.hpp"
#include "mars_mandelbrot.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_generation.hpp"
#include "mars_test_mesh.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_mesh_kokkos.hpp"
#endif

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

void read(std::string filename) {
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationInput");
    adios2::Engine reader;
    reader = io_main.Open(filename, adios2::Mode::Read);
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
    reader.EndStep();
}

int main(int argc, char const *argv[]) {
    std::string inputFileName;
    if (argc == 2) {
        inputFileName = argv[1];
    }
}