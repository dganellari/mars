#ifndef MARS_IMAGE_DATA_WRITER_HPP
#define MARS_IMAGE_DATA_WRITER_HPP

#include "adios2.h"
#include "mars_base.hpp"
#include "mars_boundary_conditions.hpp"
#include "mars_copy_operator.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_gradient_recovery.hpp"
#include "mars_identity_operator.hpp"
#include "mars_image_data_writer_settings.hpp"
#include "mars_interpolate.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_serial_mesh_type.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"

void define_bpvtk_attribute(const Settings &s, adios2::IO &io);

class ImageWriter {
public:
    Settings settings;
    adios2::IO io;
    adios2::Engine writer;
    adios2::Variable<double> var_u;
    adios2::Variable<double> var_v;
    adios2::Variable<int> var_step;

    ImageWriter(const Settings &settings, adios2::IO &io);
    ~ImageWriter();
    void open(const std::string &fname);
    void write(int step);
    void close();
};
#endif