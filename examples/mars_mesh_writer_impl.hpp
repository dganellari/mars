#include <adios2.h>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <iostream>
#include <typeinfo>
#include "mars_mesh_writer.hpp"

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
#include "mars_interpolate.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_prelepp_benchmark.hpp"
#include "mars_quality.hpp"
#include "mars_simplex.hpp"
#include "mars_utils.hpp"
#include "mars_vtk_writer.hpp"

#include "mars_longest_edge.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_test.hpp"

#include "mars_env.hpp"

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_mesh_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#endif  // WITH_KOKKOS

using VectorReal = mars::ViewVectorType<mars::Real>;

// Vtk schema for visualizing in ParaView.
template <class Mesh>
std::string MeshWriter<Mesh>::VTKSchema() {
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

template <class Mesh>
MeshWriter<Mesh>::MeshWriter(Mesh& mesh, adios2::IO io) : mesh_(mesh), io_(io) {}

// Here we open the IO of adios and define what kind of operation we want to do, in this case a write operation.
template <class Mesh>
void MeshWriter<Mesh>::open(const std::string& fname) {
    engine_ = io_.Open(fname, adios2::Mode::Write);
    engine_.BeginStep();
}

template <class Mesh>
void MeshWriter<Mesh>::generate_data_cube(const int space) {
    // Generate 3,3,3 cube on mesh_.

    switch (space) {
        case 2:
            mars::generate_cube(mesh_, 4, 3, 0);
        case 3:
            mars::generate_cube(mesh_, 4, 3, 3);
    }

    // Setup intial values of vertices, nelements and Dim.
    size_t nelements = 0;
    size_t element_nvertices = 0;
    size_t n_nodes = mesh_.n_nodes();
    assert(n_nodes > 0);
    auto points = mesh_.get_view_points();
    const size_t spaceDim = static_cast<size_t>(mesh_.Dim);
    const size_t dimension = static_cast<uint32_t>(mesh_.Dim);

    // Get new values for nelements.
    auto elements = mesh_.get_view_elements();
    nelements = mesh_.n_elements();

    assert(nelements == elements.extent(0));
    element_nvertices = elements.extent(1);

    // Print out some info of the mesh.
    std::cout << "n_active_elements: " << mesh_.n_active_elements() << std::endl;
    std::cout << "n_nodes: " << mesh_.n_nodes() << std::endl;
    std::cout << "elements extent(0): " << elements.extent(0) << std::endl;
    std::cout << "elements extent(1): " << elements.extent(1) << std::endl;
    std::cout << "elements extent(2): " << elements.extent(2) << std::endl;
    std::cout << "Space dim:" << mesh_.Dim << std::endl;

    // TODO: Write the type of Mesh that this was so that we can instantiate one when we read.
    std::string mesh_type = "Mars Unstructured Mesh";
    std::vector<std::string> viz_tools;
    viz_tools.push_back("Paraview: ADIOS2VTXReader");
    viz_tools.push_back("VTK: vtkADIOS2VTXReader.h");

    io_.DefineAttribute<std::string>("format/mars_mesh", mesh_type);
    io_.DefineAttribute<std::string>("format/viz_tools", viz_tools.data(), viz_tools.size());

    // Define the multiple variables that we need, such as vertices, types and elements.
    io_.DefineVariable<uint64_t>("connectivity", {}, {}, {nelements, element_nvertices + 1});
    io_.DefineVariable<uint32_t>("types");

    // LocalValueDim: values local to the MPI process
    // adios2::LocalValueDim is an enum
    io_.DefineVariable<uint32_t>("NumOfElements", {adios2::LocalValueDim});
    io_.DefineVariable<uint32_t>("NumOfVertices", {adios2::LocalValueDim});
    io_.DefineVariable<double>("vertices", {}, {}, {n_nodes, spaceDim});
    io_.DefineVariable<int32_t>("attribute", {}, {}, {nelements});
    io_.DefineVariable<uint32_t>("dimension", {}, {}, {dimension});
    io_.DefineVariable<double>("U", {}, {}, {n_nodes});

    point_data_variables.insert("U");
}

// Here we want to apply a function to the mesh.
template <class Mesh>
void MeshWriter<Mesh>::interpolate() {
    const mars::Integer n_nodes = mesh_.n_nodes();
    VectorReal x = VectorReal("x", n_nodes);
    mars::Interpolate<Mesh> interpMesh(mesh_);
    interpMesh.apply(
        x, MARS_LAMBDA(const mars::Real* p)->mars::Real { return p[0] * p[0] * p[0]; });
}

// Write step, write tthe results to the variables we defined before.
template <class Mesh>
void MeshWriter<Mesh>::write() {
    // Still don't really understand what this does, but on mFEM
    // TODO: FIX THIS.
    const uint32_t vtktype = 9;
    adios2::Variable<uint64_t> varConnectivity = io_.InquireVariable<uint64_t>("connectivity");
    adios2::Variable<uint32_t> varTypes = io_.InquireVariable<uint32_t>("types");
    engine_.Put(varTypes, vtktype);

    //###############################Attribute###########################
    // For MFEM attribute is specifying material properties
    adios2::Variable<int32_t> varElementAttribute = io_.InquireVariable<int32_t>("attribute");
    adios2::Variable<uint64_t>::Span spanConnectivity = engine_.Put<uint64_t>(varConnectivity);

    // // zero-copy access to adios2 buffer to put non-contiguous to contiguous memory
    adios2::Variable<int32_t>::Span spanElementAttribute = engine_.Put<int32_t>(varElementAttribute);

    auto elements = mesh_.get_view_elements();
    size_t elementPosition = 0;
    for (int e = 0; e < mesh_.n_elements(); ++e) {
        spanElementAttribute[e] = static_cast<int32_t>(e);
        spanConnectivity[elementPosition] = Mesh::Elem::NNodes;

        for (int v = 0; v < Mesh::Elem::NNodes; ++v) {
            spanConnectivity[elementPosition + v + 1] = elements(e, v);
        }
        elementPosition += Mesh::Elem::NNodes + 1;
    }

    //###################################################################

    //############################Vertices###############################
    engine_.Put("NumOfElements", static_cast<uint32_t>(mesh_.n_elements()));
    // engine_.Put("vertices", static_cast<double>(mesh_.n_nodes()));
    engine_.Put("NumOfVertices", static_cast<uint32_t>(mesh_.n_nodes()));

    std::vector<mars::Integer> e_nodes;
    // mars::IElem element;
    adios2::Variable<double> varVertices = io_.InquireVariable<double>("vertices");
    // zero-copy access to adios2 buffer to put non-contiguous to contiguous memory
    adios2::Variable<double>::Span spanVertices = engine_.Put<double>(varVertices);

    auto points = mesh_.get_view_points();

    // For each of the nodes (in this case 16 for cube 3,3,3) we iterate through the size of
    // the space_dim which is in the case of ParallelQuad4Mesh 2.
    for (int v = 0; v < mesh_.n_nodes(); ++v) {
        const int space_dim = static_cast<int>(mesh_.Dim);
        for (int coord = 0; coord < space_dim; ++coord) {
            spanVertices[v * space_dim + coord] = points(v, coord);
        }
    }

    //#################################Interpolate#################################

    // Find out how many nodes we have.
    const mars::Integer n_nodes = mesh_.n_nodes();
    // Create a vector of the size of the nodes.
    VectorReal x = VectorReal("x", n_nodes);
    mars::Interpolate<Mesh> interpMesh(mesh_);
    interpMesh.apply(
        x, MARS_LAMBDA(const mars::Real* p)->mars::Real { return p[0] * p[0] * p[0]; });

    adios2::Variable<double> varU = io_.InquireVariable<double>("U");
    // zero-copy access to adios2 buffer to put non-contiguous to contiguous memory
    adios2::Variable<double>::Span spanU = engine_.Put<double>(varU);

    for (int v = 0; v < mesh_.n_nodes(); ++v) {
        spanU[v] = x(v);
    }
}

template <class Mesh>
void MeshWriter<Mesh>::close() {
    // Define attribute vtk.xml which is how we visualize in ParaView.
    io_.DefineAttribute<std::string>("vtk.xml", VTKSchema(), {}, {});
    engine_.EndStep();
    engine_.Close();
}