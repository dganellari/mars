#include "mars_mesh_writer.hpp"
#include <adios2.h>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <iostream>

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
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
#include "mars_test_kokkos.hpp"
#endif  // WITH_KOKKOS

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

MeshWriter::MeshWriter(adios2::IO io, const std::string engineType) {
    adios2::ADIOS adios(adios2::DebugON);
    io = adios.DeclareIO("SimulationOutput");
    io.SetEngine(engineType);
}

void MeshWriter::open(const std::string& fname) {
    engine = io.Open(fname, adios2::Mode::Write);
    engine.BeginStep();
}
void MeshWriter::write(Mesh& mesh, int step) {
    size_t nelements = 0;
    size_t element_nvertices = 0;
    auto elements = mesh.get_view_elements();
    std::cout << elements.extent(0);
    std::cout << elements.extent(1);
    nelements = mesh.n_elements();
    assert(nelements == elements.extent(0));
    element_nvertices = elements.extent(1);
    int dim = Mesh::Dim;
    size_t n_nodes = mesh.n_nodes();
    auto points = mesh.get_view_points();

    // io.DefineVariable<uint64_t>("connectivity", {}, {}, {nelements, element_nvertices + 1});
}
void MeshWriter::close() {
    // Need to write the ouput in vtk schema:
    // MFEM does like this:
    // SafeDefineAttribute<std::string>(io, "vtk.xml", VTKSchema() );
    // For us would be DefineAttribute...
    engine.EndStep();
    engine.Close();
}