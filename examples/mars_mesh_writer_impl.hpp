#include <adios2.h>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <iostream>
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
#include "mars_mesh_writer.hpp"
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

template <class Mesh>
void example_func() {}

// Vtk schema for visualizing in ParaView.
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
void MeshWriter<Mesh>::generate_data_cube() {
    // Setup intial values of vertices, nelements and Dim.
    size_t nelements = 0;
    size_t element_nvertices = 0;
    size_t n_nodes = mesh_.n_nodes();
    auto points = mesh_.get_view_points();
    const size_t spaceDim = static_cast<size_t>(mesh_.Dim);

    // Generate 3,3,3 cube on mesh_.
    mars::generate_cube(mesh_, 3, 3, 3);

    // Get new values for nelements.
    auto elements = mesh_.get_view_elements();
    nelements = mesh_.n_elements();

    assert(nelements == elements.extent(0));
    element_nvertices = elements.extent(1);

    // Print out some info on the mesh.
    std::cout << "n_active_elements: " << mesh_.n_active_elements() << std::endl;
    std::cout << "n_nodes: " << mesh_.n_nodes() << std::endl;
    std::cout << "elements extent(0): " << elements.extent(0) << std::endl;
    std::cout << "elements extent(1): " << elements.extent(1) << std::endl;
    std::cout << "Space dim:" << mesh_.Dim << std::endl;

    // Define the multiple variables that we need, such as vertices, types and elements.
    io_.DefineVariable<uint64_t>("connectivity", {}, {}, {nelements, element_nvertices + 1});
    io_.DefineVariable<uint32_t>("types");
    io_.DefineVariable<uint32_t>("NumOfElements", {adios2::LocalValueDim});
    io_.DefineVariable<double>("vertices", {}, {}, {n_nodes, spaceDim});
}

// Here we want to apply a function to the mesh.
template <class Mesh>
void MeshWriter<Mesh>::interpolate(VectorReal& x) {
    x = VectorReal("x", mesh_.n_nodes());
    mars::Interpolate<Mesh> interpMesh(mesh_);
    // interpMesh.apply(x, );
}

// Write step, write tthe results to the variables we defined before.
template <class Mesh>
void MeshWriter<Mesh>::write() {
    adios2::Variable<uint64_t> varConnectivity = io_.InquireVariable<uint64_t>("connectivity");
    adios2::Variable<uint32_t> varTypes = io_.InquireVariable<uint32_t>("types");

    engine_.Put("NumOfElements", static_cast<uint32_t>(mesh_.n_active_elements()));
    engine_.Put("vertices", static_cast<double>(mesh_.n_nodes()));
    engine_.Put<uint64_t>(varConnectivity);
    // engine_.Put("types", varTypes);
}

template <class Mesh>
void MeshWriter<Mesh>::close() {
    // Need to write the ouput in vtk schema:
    // MFEM does like this:
    // Define attribute of vtk schema
    io_.DefineAttribute<std::string>("vtk.xml", VTKSchema(), {}, {});
    engine_.EndStep();
    engine_.Close();
}