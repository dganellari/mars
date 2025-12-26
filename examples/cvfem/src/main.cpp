#include <chrono>
#include <iomanip>
#include <iostream>

#include "cxxopts.hpp"
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>

// #include <span>

using scalar = double;
// using scalar = float; // TODO: not currently meaningful, stk mesh uses
// hardcoded doubles sometimes
using label = int;

constexpr unsigned BLOCKSIZE = 1;
constexpr unsigned SPATIAL_DIM = 3;

#include <CRSMatrix.h>
#include <CRSNodeGraph.h>
#include <myNodeGraph.h>

#include <meshTestTools.h>

#include <assemblerFactory.h>

#include <fenv.h>

int mpi_rank = 0;
int mpi_size = 0;

int n_reps = 1;

template <size_t BLOCKSIZE>
using Coefficients = linearSolver::coefficients<BLOCKSIZE>;

#include <mainFuncs.h>

int main(int argc, char* argv[])
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    // add default assembler list
    std::set<std::string> assemblerNames;
    // assemblerNames.insert("cvfem"); 
    // assemblerNames.insert("cvfemTeamScratch");
    assemblerNames.insert("cvfemTeamBalance");
    // assemblerNames.insert("edge");
    assemblerNames.insert("cvfemInlineSF");

    std::string baseName;
    int verbose = 0;

    // read runtime options
    parseOptions(argc, argv, baseName, assemblerNames, verbose, n_reps);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    Kokkos::initialize(argc, argv);
    {
        // create mesh pointers
        std::unique_ptr<stk::mesh::BulkData> bulkDataPtr_ = nullptr;
        std::unique_ptr<stk::io::StkMeshIoBroker> stkIoPtr_ = nullptr;
        stk::ParallelMachine pm = MPI_COMM_WORLD;
        stkIoPtr_ = std::make_unique<stk::io::StkMeshIoBroker>(pm);

        // read mesh from exodus file
        read_mesh_from_file(
            baseName + ".exo", pm, bulkDataPtr_, stkIoPtr_, verbose);

        // get primary mesh info references
        stk::mesh::MetaData& metaData = bulkDataPtr_->mesh_meta_data();
        const stk::mesh::BulkData& bulkData = *bulkDataPtr_.get();

        // read nodegraph from file
        accel::myNodeGraph nodegraph(
            baseName + ".bin",
            MPI_COMM_WORLD,
            linearSolver::GraphLayout::ColumnIndexOrder__Local,
            verbose);

        // add fields to registry
        addFieldsToMesh(metaData);

        // populate fields
        fill_phi_gradPhi<scalar>(bulkData, verbose);
        set_field_on_host<scalar>(
            bulkData, "gamma", stk::topology::NODE_RANK, scalar(0.1), verbose);
        set_field_on_host<scalar>(
            bulkData, "beta", stk::topology::NODE_RANK, scalar(1.234), verbose);
        set_field_on_host<scalar>(
            bulkData, "mdot", stk::topology::ELEMENT_RANK, scalar(0), verbose);
        set_field_on_host<scalar>(
            bulkData, "mdotEdge", stk::topology::EDGE_RANK, scalar(0), verbose);

        // populate areaVector fields
        calculateEdgeAreaVectors<scalar, label, BLOCKSIZE, SPATIAL_DIM>(
            bulkData);
        preCalcAreaV<scalar, label, BLOCKSIZE, SPATIAL_DIM>(bulkData);

        // setup coefficients for A, x, b, res
        Coefficients<BLOCKSIZE> coeffs(&nodegraph);

        // evaluate data transfer size for performance eval
        size_t total_bytes = 0;
        getDataSize<scalar, BLOCKSIZE>(total_bytes, coeffs, bulkData, verbose);

        // create assemblers
        std::vector<
            std::unique_ptr<AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>>>
            assemblerPtrs;
        for (const std::string& asmName : assemblerNames)
        {
            assemblerPtrs.push_back(
                createAssemblerInstance<scalar, label, BLOCKSIZE, SPATIAL_DIM>(
                    asmName, bulkData));
        }

        // run assemblers
        for (const auto& asmPtr : assemblerPtrs)
        {
            runAssembler(
                asmPtr, coeffs, metaData, bulkData, total_bytes, verbose);
        }

        // write rows/columns/values to files
        if (verbose > 99)
        {
            coeffs.getAMatrix().writeMatrix("matrixDump");
        }
    }

    Kokkos::finalize();

    MPI_Finalize();
}
