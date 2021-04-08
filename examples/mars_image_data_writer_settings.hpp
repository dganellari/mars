#ifndef MARS_IMAGE_DATA_WRITER_SETTINGS
#define MARS_IMAGE_DATA_WRITER_SETTINGS

#include <string>

struct Settings {
    size_t L;
    int steps;
    int plotgap;
    double F;
    double k;
    double dt;
    double Du;
    double Dv;
    double noise;
    std::string output;
    bool checkpoint;
    int checkpoint_freq;
    std::string checkpoint_output;
    std::string adios_config;
    bool adios_span;
    bool adios_memory_selection;
    std::string mesh_type;

    Settings();
};
#endif