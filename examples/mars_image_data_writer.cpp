#include "mars_image_data_writer.hpp"
#include <cmath>
#include <iostream>
#include "adios2.h"
#include "mars_base.hpp"
#include "mars_globals.hpp"
#include "mars_image_data_writer_settings.hpp"

double simple_func(const double &x, const double &y, const double &z) { return std::sqrt(x * x + y * y + z * z); }

// Example data vector creation.
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
    std::cout << "Var_Data:" << var_data.Type() << std::endl;
    std::cout << Nx << std::endl << Ny << std::endl << Nz << std::endl;
    std::cout << "ImageWriter::\n";
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

        io.DefineAttribute<std::string>("adios2.xml", imageData);
    };

    if (s.mesh_type == "image") {
        lf_VTKImage(s, io);
    }
}

ImageWriter::ImageWriter(const Settings &settings, adios2::IO io) : io(io) { define_bpvtk_attribute(settings, io); }

void ImageWriter::open(const std::string &fname) { writer = io.Open(fname, adios2::Mode::Write); }

/*
 *Writing step, Begin part, put, End part
 */
void ImageWriter::write(int step) {
    writer.BeginStep();
    writer.Put<double>(var_data, data.data());
    writer.EndStep();
}
void ImageWriter::close() { writer.Close(); }