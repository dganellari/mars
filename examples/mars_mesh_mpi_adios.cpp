#include <complex>
#include <iostream>
#include "adios2.h"
#include "cxxopts.hpp"
#include "mars_globals.hpp"
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

std::string VTKSchema() {
    std::string vtkSchema = R"(
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.2" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="NumOfVertices" NumberOfCells="NumOfElements">
      <Points>
        <DataArray Name="vertices" />)";

    vtkSchema += R"(
      </Points>
          <CellData>
        <DataArray Name="attribute" />
    </CellData>
      <Cells>
        <DataArray Name="connectivity" />
        <DataArray Name="types" />
      </Cells>
      <PointData>)";

    // if (point_data_variables.empty()) {
    //     vtkSchema += "\n";
    // } else {
    //     for (const std::string& point_datum : point_data_variables) {
    //         vtkSchema += "        <DataArray Name=\"" + point_datum + "\"/>\n";
    //     }
    // }

    // vtkSchema += "        <DataArray Name=\"TIME\">\n";
    // vtkSchema += "          TIME\n";
    // vtkSchema += "        </DataArray>\n";

    vtkSchema += R"(
       </PointData>
       </Piece>
     </UnstructuredGrid>
   </VTKFile>)";

    return vtkSchema;
}

int main(int argc, char* argv[]) {
    int rank, size;
    int my_grid_rank;

#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm gridComm;

    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("BPFile_SZ");
    adios2::Engine engine;

    Kokkos::initialize();

    // Use generate unit cube, then translate by rank.
    // Will form long mesh, 6,1,1 with 6 processors.
    // Create global id's for nodes.(Calculate offset for each node and element).
    // Export rank.
    // Calculate function for nodes (interpolate mandelbrot function).

    // mars::Mesh<3> mesh = mars::generate_unit_cube();
    {
        mars::ParallelMesh3 mesh;

        mars::generate_cube(mesh, 1, 1, 1);

        auto elementsmesh = mesh.get_view_elements();

        size_t n_nodes = mesh.n_nodes();
        // // auto points = mesh.get_view_points();
        const size_t spaceDim = static_cast<size_t>(mesh.Dim);
        const uint32_t vtktype = 9;
        // // auto elements l= mesh.get_view_elements();
        size_t nelements = mesh.n_elements();
        size_t element_nvertices = 0;

        element_nvertices = elementsmesh.extent(1);

        // Create vector for global index of each node.

        // For each node we want to translate on the x axis by rank.
        for (size_t i = 0; i < n_nodes; ++i) {
            auto point = mesh.point(i)[0];
            mesh.point(i)[0] = point + rank;
        }

        if (rank == 0) {
            std::cout << "SOME INFO ABOUT THE MESH:\n";
            std::cout << "Number of elements: " << nelements << std::endl;
            std::cout << "Number of nodes: " << n_nodes << std::endl;
            std::cout << "Space dim: " << spaceDim << std::endl;
            std::cout << "Type of mesh: " << mesh.type() << std::endl;
            std::cout << "Active elements: " << mesh.n_active_elements() << std::endl;
        }

        engine = io.Open("adios-mesh-mpi.bp", adios2::Mode::Write);
        engine.BeginStep();
        io.DefineVariable<uint32_t>("NumOfElements", {adios2::LocalValueDim});
        io.DefineVariable<uint32_t>("NumOfVertices", {adios2::LocalValueDim});
        io.DefineVariable<double>("vertices", {}, {}, {n_nodes, spaceDim});
        io.DefineVariable<uint64_t>("connectivity", {}, {}, {nelements, element_nvertices + 1});
        io.DefineVariable<int32_t>("attribute", {}, {}, {nelements});
        io.DefineVariable<uint32_t>("types");

        std::string mesh_type = "Mars Unstructured Mesh";
        std::vector<std::string> viz_tools;
        viz_tools.push_back("Paraview: ADIOS2VTXReader");
        viz_tools.push_back("VTK: vtkADIOS2VTXReader.h");

        io.DefineAttribute<std::string>("format/mars_mesh", mesh_type);
        io.DefineAttribute<std::string>("format/viz_tools", viz_tools.data(), viz_tools.size());

        engine.Put("NumOfElements", static_cast<uint32_t>(mesh.n_elements()));

        adios2::Variable<int32_t> varElementAttribute = io.InquireVariable<int32_t>("attribute");
        adios2::Variable<int32_t>::Span spanElementAttribute = engine.Put<int32_t>(varElementAttribute);
        adios2::Variable<uint32_t> varTypes = io.InquireVariable<uint32_t>("types");
        adios2::Variable<uint64_t> varConnectivity = io.InquireVariable<uint64_t>("connectivity");
        adios2::Variable<uint64_t>::Span spanConnectivity = engine.Put<uint64_t>(varConnectivity);

        auto elements = mesh.get_view_elements();
        size_t elementPosition = 0;
        for (int e = 0; e < mesh.n_elements(); ++e) {
            spanElementAttribute[e] = static_cast<int32_t>(e);
            spanConnectivity[elementPosition] = mars::ParallelMesh3::Elem::NNodes;
            for (int v = 0; v < mars::ParallelMesh3::Elem::NNodes; ++v) {
                spanConnectivity[elementPosition + v + 1] = elements(e, v);
            }
            elementPosition += mars::ParallelMesh3::Elem::NNodes + 1;
        }

        engine.Put("NumOfVertices", static_cast<uint32_t>(mesh.n_nodes()));
        adios2::Variable<double> varVertices = io.InquireVariable<double>("vertices");
        // zero-copy access to adios2 buffer to put non-contiguous to contiguous memory
        adios2::Variable<double>::Span spanVertices = engine.Put<double>(varVertices);

        auto points = mesh.get_view_points();

        // For each of the nodes (in this case 16 for cube 3,3,3) we iterate through the size of
        // the space_dim which is in the case of ParallelQuad4Mesh 2.
        for (int v = 0; v < mesh.n_nodes(); ++v) {
            const int space_dim = static_cast<int>(mesh.Dim);
            for (int coord = 0; coord < space_dim; ++coord) {
                spanVertices[v * space_dim + coord] = points(v, coord);
            }
        }
        engine.Put(varTypes, vtktype);
        io.DefineAttribute<std::string>("vtk.xml", VTKSchema(), {}, {});

        engine.EndStep();
        engine.Close();
    }
    Kokkos::finalize();
    MPI_Finalize();
#endif
    return 0;
}