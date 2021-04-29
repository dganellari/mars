#include "mars_image_data_writer_settings.hpp"
#include <iostream>

Settings::Settings(const std::string &fname) {
    L = 128;
    step = 2000;
    Du = 0.05;
    output = fname;
    adios_config = "adios2.xml";
    adios_span = false;
    adios_memory_selection = false;
    mesh_type = "image";
}

Settings::Settings() {
    L = 128;
    step = 2000;
    Du = 0.05;
    output = "solution.bp";
    adios_config = "adios2.xml";
    adios_span = false;
    adios_memory_selection = false;
    mesh_type = "image";
}
