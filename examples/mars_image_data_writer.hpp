#ifndef MARS_IMAGE_DATA_WRITER_HPP
#define MARS_IMAGE_DATA_WRITER_HPP

#include "adios2.h"
#include "mars_image_data_writer_settings.hpp"

class ImageWriter {
public:
    Settings settings;
    adios2::IO io;
    adios2::Engine writer;
    adios2::Variable<double> var_u;
    adios2::Variable<double> var_v;
    adios2::Variable<int> var_step;

    ImageWriter();
    ImageWriter(const Settings &settings, adios2::IO io);
    void open(const std::string &fname);
    void write(int step, double data[]);
    void close();
};
#endif