#include <cmath>
#include <complex>
#include <iostream>
#include <map>
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

    const std::map<std::string, adios2::Params> variables = io_main.AvailableVariables();

    for (const auto variablePair : variables) {
        std::cout << "Name: " << variablePair.first;

        for (const auto &parameter : variablePair.second) {
            std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
        }
    }

    adios2::Variable<double> uVar = io_main.InquireVariable<double>("U");
    adios2::Attribute<std::string> uAttr = io_main.InquireAttribute<std::string>("format/mars_mesh");
    adios2::Variable<uint32_t> numOfElements = io_main.InquireVariable<uint32_t>("NumOfElements");
    if (uVar) {
        std::vector<double> myU;
        std::cout << "Got it\n";

        reader.Get(uVar, myU, adios2::Mode::Deferred);
        std::cout << myU.size() << std::endl;
    }
    reader.EndStep();
}

int main(int argc, char const *argv[]) {
    std::string inputFileName;
    if (argc == 2) {
        inputFileName = argv[1];
    }

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    read(inputFileName);

    MPI_Finalize();
}