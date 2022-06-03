#include "mars_mesh_reader.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include "adios2.h"

#include "mars_globals.hpp"

#include <mpi.h>

void read(std::string filename, int rank, std::vector<double>& uVector) {
    // INITIAL DECLERATIONS OF ADIOS MODE, ADIOS IO, ADIOS ENGINE
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io = adios.DeclareIO("SimulationInput");
    adios2::Engine reader;

    // SET ENGINE TO OPEN FILE IN READ MODE AND CALL BeginStep() on engine.
    reader = io.Open(filename, adios2::Mode::Read);
    reader.BeginStep();

    // AvaiableVariables on io returns a map of string, adios2::Params, which contatains all the
    // variables declared in the .bp file.
    const std::map<std::string, adios2::Params> variables = io.AvailableVariables();

    // Print out the name and value of those variables.
    if (rank == 0) {
        for (const auto variablePair : variables) {
            std::cout << "Name: " << variablePair.first;

            for (const auto& parameter : variablePair.second) {
                std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
            }
        }
    }

    adios2::Variable<double> uVar = io.InquireVariable<double>("U");
    if (uVar) {
        std::vector<double> myU;
        std::cout << "Got variable U\n";

        reader.Get(uVar, myU, adios2::Mode::Sync);
        for (const auto number : myU) {
            std::cout << number << " ";
        }

        if (myU == uVector) {
            std::cout << "Same vector" << std::endl;
        }
    }
    reader.EndStep();

    reader.Close();
}

int main(int argc, char* argv[]) {
    //     std::string inputFileName;
    //     if (argc == 2) {
    //         inputFileName = argv[1];
    //     }
    //     MPI_Init(&argc, &argv);
    //     int rank, size;
    //     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //     MPI_Comm_size(MPI_COMM_WORLD, &size);
    // #if ADIOS2_USE_MPI
    //     read(inputFileName, rank);
    // #endif

    //     MPI_Finalize();
}