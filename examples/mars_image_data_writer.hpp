#ifndef MARS_IMAGE_DATA_WRITER_HPP
#define MARS_IMAGE_DATA_WRITER_HPP

#include "adios2.h"
#include "array"
#include "mars_image_data_writer_settings.hpp"
#include "vector"

class ImageWriter {
public:
    ImageWriter(const Settings &settings, adios2::IO io);
    void new_data(const unsigned long &Nx, const unsigned long &Ny, const unsigned long &Nz);
    void open(const std::string &fname);
    void write(int step);
    void close();

protected:
    Settings settings;
    adios2::IO io;
    adios2::Engine writer;
    adios2::Variable<double> var_data;
    std::vector<double> data;
};
#endif