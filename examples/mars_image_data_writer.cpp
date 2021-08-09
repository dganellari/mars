#include "mars_image_data_writer.hpp"
#include <cmath>
#include <iostream>
#include "adios2.h"
#include "mars_globals.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mpi.h"

/**
 * Simple function to perform on the vector values.
 *
 * @param x,y,z values of 3-dimentional vector.
 **/
double simple_func(const double &x, const double &y, const double &z) { return std::sqrt(x * x + y * y + z * z); }

/**
 * Given Nx,Ny,Nz, which are the number of cells we want to create.
 * Apply a function simple_func(x,y,z) to each of the cells and then store these values
 * in a adios2 Variable called "U".
 *
 * @param Nx, Ny, Nz size of the image.
 **/
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

    // TODO: Need to parallelize this somehow.
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                x = i * H[0];
                y = j * H[1];
                // Create global index for k.
                z = k * H[2];
                data[i * Ny * Nz + j * Nz + k] = simple_func(x, y, z);
            }
        }
    }
    // Indicizzazione local to global.
    // NZlocal, NZGlobal.
    // Remember remainder of div.
    // data = NX * Ny * NZlocal
    var_data = io.DefineVariable<double>("U", {Nx, Ny, Nz}, {0UL, 0UL, 0UL}, {Nx, Ny, Nz});

    // std::cout<< var_data.Sizeof()
    //     std::cout << "Var_Data:" << var_data.Type() << std::endl;
    //     std::cout << Nx << std::endl << Ny << std::endl << Nz << std::endl;
    //     std::cout << "ImageWriter::\n";
}

/**
 * Define the adios2.xml attribute which contains this xml string
 * which can be read by ParaView for visualizing the data files.
 *
 * @param s Settings object which contains details on how to create the xml.
 * @param io IO for Defining attributes,variables and other functionalities.
 **/
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

/**
 * Constructor of class ImageWriter. For now it calls the define_bpvtk_attribute() method.
 *
 * @param settings
 * @param io IO of adios2
 **/
ImageWriter::ImageWriter(const Settings &settings, adios2::IO io) : io(io) { define_bpvtk_attribute(settings, io); }

/**
 * Tells the engine(writer) to open fname for writing.
 * @param fname, name of file.
 **/
void ImageWriter::open(const std::string &fname) { writer = io.Open(fname, adios2::Mode::Write); }

/**
 * Writing step: Begin, Put data from the data vector into var_data adios variable.
 * @param step, which is useless right now.
 **/
void ImageWriter::write(int step) {
    writer.BeginStep();
    writer.Put<double>(var_data, data.data());
    writer.EndStep();
}
/**
 * Close the engine.
 **/
void ImageWriter::close() { writer.Close(); }