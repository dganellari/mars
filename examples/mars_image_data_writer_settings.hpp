#ifndef MARS_IMAGE_DATA_WRITER_SETTINGS
#define MARS_IMAGE_DATA_WRITER_SETTINGS

#include <string>

struct Settings {
    size_t L;
    int step;
    double Du;
    std::string output;
    std::string adios_config;
    bool adios_span;
    bool adios_memory_selection;
    std::string mesh_type;

    Settings();
    Settings(const std::string &fname);
};
#endif