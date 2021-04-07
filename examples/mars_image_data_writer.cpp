#include "adios2.h"
#include "mars_base.hpp"
#include "mars_boundary_conditions.hpp"
#include "mars_copy_operator.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_gradient_recovery.hpp"
#include "mars_identity_operator.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_interpolate.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_serial_mesh_type.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"

void define_bpvtk_attribute(const Settings &s, adios2::IO &io) {
    auto lf_VTKImage = [](const Settings &s, adios2::IO &io) {
        const std::string extent =
            "0 " + std::to_string(s.L) + " " + "0 " + std::to_string(s.L) + " " + "0 " + std::to_string(s.L);

        const std::string imageData = R"(
        <?xml version="1.0"?>
        <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
          <ImageData WholeExtent=")" + extent +
                                      R"(" Origin="0 0 0" Spacing="1 1 1">
            <Piece Extent=")" + extent +
                                      R"(">
              <CellData Scalars="U">
                  <DataArray Name="U" />
                  <DataArray Name="V" />
                  <DataArray Name="TIME">
                    step
                  </DataArray>
              </CellData>
            </Piece>
          </ImageData>
        </VTKFile>)";

        io.DefineAttribute<std::string>("vtk.xml", imageData);
    };

    if (s.mesh_type == "image") {
        lf_VTKImage(s, io);
    } else if (s.mesh_type == "structured") {
        throw std::invalid_argument(
            "ERROR: mesh_type=structured not yet "
            "   supported in settings.json, use mesh_type=image instead\n");
    }
    // TODO extend to other formats e.g. structured
}

ImageWriter::ImageWriter(const Settings &settings, adios2::IO io) {
    io.DefineAttribute<double>("Du", settings.Du);

    if (!settings.mesh_type.empty()) {
        define_bpvtk_attribute(settings, io);
    }

    var_u = io.DefineVariable<double>("U", {settings.L, settings.L, settings.L});

    // if (settings.adios_memory_selection) {
    //     var_u.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    //     var_v.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    // }

    var_step = io.DefineVariable<int>("step");
}

ImageWriter::open(const std::string &fname) { writer = io.Open(fname, adios2::Mode::Write); }

void ImageWriter::write(int step, double &data[]) {
    writer.BeginStep();
    writer.Put<double>(var_u, data.data());
    writer.EndStep();
}
void ImageWriter::close() { writer.Close(); }