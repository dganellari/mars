#include "mars_image_data_writer_settings.hpp"

Settings::Settings() {
    L = 128;
    steps = 20000;
    plotgap = 200;
    F = 0.04;
    k = 0.06075;
    dt = 0.2;
    Du = 0.05;
    Dv = 0.1;
    noise = 0.0;
    output = "solution.bp";
    checkpoint = false;
    checkpoint_freq = 2000;
    checkpoint_output = "gs_ckpt.bp";
    adios_config = "adios2.xml";
    adios_span = false;
    adios_memory_selection = false;
    mesh_type = "image";
}
