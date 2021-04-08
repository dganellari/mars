#include "mars_image_data_writer.hpp"
#include <iostream>
#include "adios2.h"
#include "mars_image_data_writer_settings.hpp"

double simple_func(const double &x, const double &y, const double &z) { return x * y * z; }

void ImageWriter::new_data(const unsigned long &Nx, const unsigned long &Ny, const unsigned long &Nz) {
    unsigned long size = Nx * Ny * Nz;
    double H[3];
    H[0] = 0.5;
    H[1] = 0.25;
    H[2] = 1;
    data.resize(size);
    double x;
    double y;
    double z;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                x = i * H[0];
                y = j * H[1];
                z = k * H[2];
                data[i * Ny * Nz + j * Nz + k] = simple_func(x, y, z);
            }
        }
    }
    var_data = io.DefineVariable<double>("U", {Nx, Ny, Nz}, {0UL, 0UL, 0UL}, {Nx, Ny, Nz});
}

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

ImageWriter::ImageWriter(const Settings &settings, adios2::IO io) : io(io) {
    // if (!settings.mesh_type.empty()) {
    //     define_bpvtk_attribute(settings, io);
    // }

    // var_u = io.DefineVariable<double>("U", {settings.L, settings.L, settings.L});

    // // if (settings.adios_memory_selection) {
    // //     var_u.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    // //     var_v.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    // // }

    // var_step = io.DefineVariable<int>("step");
}

void ImageWriter::open(const std::string &fname) { writer = io.Open(fname, adios2::Mode::Write); }

void ImageWriter::write(int step) {
    writer.BeginStep();
    writer.Put<double>(var_data, data.data());
    writer.EndStep();
}
void ImageWriter::close() { writer.Close(); }