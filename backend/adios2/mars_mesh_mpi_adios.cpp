#include <complex>
#include <iostream>
#include "adios2.h"
#include "mars_config.hpp"
#include "mars_globals.hpp"
#include "mars_interpolate.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_generation.hpp"
#include "mars_mesh_reader.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_mesh_kokkos.hpp"
#endif

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

using VectorReal = mars::ViewVectorType<mars::Real>;

/**
 * Print to stdout some info about the singular mesh each process creates.
 * @param nElements, number of elements
 * @param nNodes, number of nodes.
 * @param space, Space dimention of mesh.
 * @param rank, rank of process.
 **/
void printMeshInfo(const size_t nElements, const size_t nNodes, const size_t space, const int rank) {
    if (rank == 0) {
        std::cout << "SOME INFO ABOUT THE MESH:\n";
        std::cout << "Number of elements: " << nElements << std::endl;
        std::cout << "Number of nodes: " << nNodes << std::endl;
        std::cout << "Space dim: " << space << std::endl;
    }
}

/**
 * .xml attribute for ParaView visualization.
 * Must contain, types, connectivity, vertices.
 * @param std::set<std::string> Set of strings to add to this .xml.
 **/
std::string VTKSchema(std::set<std::string> point_data_variables) {
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

    if (point_data_variables.empty()) {
        vtkSchema += "\n";
    } else {
        for (const std::string& point_datum : point_data_variables) {
            vtkSchema += "        <DataArray Name=\"" + point_datum + "\"/>\n";
        }
    }

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
    std::set<std::string> point_data_variables;

#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm gridComm;

    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("BPFile_SZ");
    adios2::Engine engine;

    // Initialize Kokkos
    Kokkos::initialize();

    // Create a scope for the mesh, without it it will not delete from memory.
    {
        // Mesh creation and generation.
        mars::ParallelMesh3 mesh;
        mars::generate_cube(mesh, 10, 10, 10);

        auto elementsmesh = mesh.get_view_elements();
        size_t n_nodes = mesh.n_nodes();
        const size_t spaceDim = static_cast<size_t>(mesh.Dim);
        const uint32_t vtktype = 10;
        size_t nelements = mesh.n_elements();
        size_t element_nvertices = static_cast<size_t>(elementsmesh.extent(1));

        // For each node translate on the x axis by rank.
        // This create a long polygon.
        for (size_t i = 0; i < n_nodes; ++i) {
            auto point = mesh.point(i)[0];
            mesh.point(i)[0] = point + rank;
        }

        // Print info about the mesh.
        printMeshInfo(nelements, n_nodes, spaceDim, rank);

        // Open engine for writing, define .bp file-name.
        engine = io.Open("adios-mesh-mpi.bp", adios2::Mode::Write);
        engine.BeginStep();

        // Define all adios variables
        io.DefineVariable<uint32_t>("NumOfElements", {adios2::LocalValueDim});
        io.DefineVariable<uint32_t>("NumOfVertices", {adios2::LocalValueDim});
        io.DefineVariable<double>("vertices", {}, {}, {n_nodes, spaceDim});
        io.DefineVariable<uint64_t>("connectivity", {}, {}, {nelements, element_nvertices + 1});
        io.DefineVariable<int32_t>("attribute", {}, {}, {nelements});

        // Vtk type of mesh:
        // 1: POINT, 3:SEGMENT, 5:TRIANGLE, 9:SQUARE, 10:TETRAHEDRON, 12:CUBE, 13:PRISM
        io.DefineVariable<uint32_t>("types");

        // Where the solution to the function is written.
        io.DefineVariable<double>("U", {}, {}, {n_nodes});
        point_data_variables.insert("U");

        // More info written in .bp file.
        std::string mesh_type = "Mars Unstructured Mesh";
        std::vector<std::string> viz_tools;
        viz_tools.push_back("Paraview: ADIOS2VTXReader");
        viz_tools.push_back("VTK: vtkADIOS2VTXReader.h");

        io.DefineAttribute<std::string>("format/mars_mesh", mesh_type);
        io.DefineAttribute<std::string>("format/viz_tools", viz_tools.data(), viz_tools.size());

        engine.Put("NumOfElements", static_cast<uint32_t>(mesh.n_elements()));

        // Getting the previously defined variables and storing in new ones.
        // Could be done with one single line.
        adios2::Variable<int32_t> varElementAttribute = io.InquireVariable<int32_t>("attribute");
        adios2::Variable<int32_t>::Span spanElementAttribute = engine.Put<int32_t>(varElementAttribute);
        adios2::Variable<uint32_t> varTypes = io.InquireVariable<uint32_t>("types");
        adios2::Variable<uint64_t> varConnectivity = io.InquireVariable<uint64_t>("connectivity");
        adios2::Variable<uint64_t>::Span spanConnectivity = engine.Put<uint64_t>(varConnectivity);

        // Iterate through the elements of the mesh, each element has 4 vertices(meshNNodes).
        // Then set up the connectivity between the elements.
        size_t elementPosition = 0;
        for (int e = 0; e < nelements; ++e) {
            const int meshNNodes = mars::ParallelMesh3::Elem::NNodes;
            spanElementAttribute[e] = static_cast<int32_t>(e);
            spanConnectivity[elementPosition] = meshNNodes;
            for (int v = 0; v < meshNNodes; ++v) {
                spanConnectivity[elementPosition + v + 1] = elementsmesh(e, v);
            }
            elementPosition += meshNNodes + 1;
        }

        engine.Put("NumOfVertices", static_cast<uint32_t>(mesh.n_nodes()));
        adios2::Variable<double> varVertices = io.InquireVariable<double>("vertices");
        // zero-copy access to adios2 buffer to put non-contiguous to contiguous memory
        adios2::Variable<double>::Span spanVertices = engine.Put<double>(varVertices);

        auto points = mesh.get_view_points();
        // For each of the nodes (in this case 16 for cube 3,3,3) we iterate through the size of
        // the space_dim which is in the case of ParallelQuad4Mesh 2.
        for (int v = 0; v < n_nodes; ++v) {
            const int space_dim = static_cast<int>(mesh.Dim);
            for (int coord = 0; coord < space_dim; ++coord) {
                spanVertices[v * space_dim + coord] = points(v, coord);
            }
        }

        int maxiter{500};
        int outofbounds{3};
        VectorReal x = VectorReal("x", n_nodes);
        mars::Interpolate<mars::ParallelMesh3> interpMesh(mesh);
        interpMesh.apply(
            x, MARS_LAMBDA(const mars::Real* p)->mars::Real {
                // std::complex<double> c{p[0], p[1]};
                // std::complex<double> z = c;
                // int i = 0;
                // while (i < maxiter) {
                //     double n = std::norm(z);
                //     if (n > outofbounds) {
                //         break;
                //     }

                //     z = z * z + c;
                //     i++;
                // }

                // return i;
                return (p[0] * p[1] * p[2]);
            });

        // Variables for U.
        adios2::Variable<double> varU = io.InquireVariable<double>("U");
        // zero-copy access to adios2 buffer to put non-contiguous to contiguous memory
        adios2::Variable<double>::Span spanU = engine.Put<double>(varU);
        for (int v = 0; v < n_nodes; ++v) {
            spanU[v] = x(v);
        }

        engine.Put(varTypes, vtktype);
        io.DefineAttribute<std::string>("vtk.xml", VTKSchema(point_data_variables), {}, {});

        engine.EndStep();
        engine.Close();
    }
    Kokkos::finalize();

    // read("adios-mesh-mpi.bp", rank, x);

#endif
    return 0;
}