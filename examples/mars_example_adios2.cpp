#include <adios2.h>
#include <iostream>
#include "mars_image_data_writer.hpp"
#include "mars_image_data_writer_settings.hpp"

int main(int argc, char *argv[]) {
    Settings settings;
    adios2::ADIOS adios(adios2::DebugON);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");

    ImageWriter main_image(settings, io_main);

    // Create example data vector
    main_image.new_data(500, 600, 1);
    // Open the writer, (which now is in write mode), with the settings found.
    main_image.open(settings.output);
    // Write with the following steps
    main_image.write(1);
    // Close writer
    main_image.close();
}