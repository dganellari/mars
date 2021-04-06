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
                  <DataArray Name="V" />
                  <DataArray Name="TIME">
                    step
                  </DataArray>
              </CellData>
            </Piece>
          </ImageData>
        </VTKFile>)";

        io.DefineAttribute<std::string>("vtk.xml", imageData);
    };

    if (s.mesh_type == "image") {
        lf_VTKImage(s, io);
    } else if (s.mesh_type == "structured") {
        throw std::invalid_argument(
            "ERROR: mesh_type=structured not yet "
            "   supported in settings.json, use mesh_type=image instead\n");
    }
    // TODO extend to other formats e.g. structured
}

// TODO add PMesh as GrayScott sim

ImageWriter::ImageWriter(const Settings &settings, PMesh &mesh, adios2::IO &io) {
    // io.DefineAttribute<double>("F", settings.F);
    // io.DefineAttribute<double>("k", settings.k);
    // io.DefineAttribute<double>("dt", settings.dt);
    io.DefineAttribute<double>("Du", settings.Du);
    // io.DefineAttribute<double>("Dv", settings.Dv);
    // io.DefineAttribute<double>("noise", settings.noise);

    if (!settings.mesh_type.empty()) {
        define_bpvtk_attribute(settings, io);
    }

    var_u = io.DefineVariable<double>("U",
                                      {settings.L, settings.L, settings.L},
                                      {sim.offset_z, sim.offset_y, sim.offset_x},
                                      {sim.size_z, sim.size_y, sim.size_x});

    var_v = io.DefineVariable<double>("V",
                                      {settings.L, settings.L, settings.L},
                                      {sim.offset_z, sim.offset_y, sim.offset_x},
                                      {sim.size_z, sim.size_y, sim.size_x});

    // TODO: figure out what this does...
    // size_x, size_z and size_y are dimentions of local array
    // What is sim? PMesh?
    if (settings.adios_memory_selection) {
        var_u.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
        var_v.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    }

    var_step = io.DefineVariable<int>("step");
}

ImageWriter::open(const std::string &fname) { writer = io.Open(fname, adios2::Mode::Write); }

void ImageWriter::write(int step) {  // std::cout << "Writing results to disk..." << std::endl;

    // Integer n_nodes = mesh.n_nodes();

    // VectorReal::HostMirror x_host("x_host", n_nodes);
    // VectorReal::HostMirror rhs_host("rhs_host", n_nodes);
    // Kokkos::deep_copy(x_host, x);
    // Kokkos::deep_copy(rhs_host, rhs);

    // SMesh serial_mesh;
    // convert_parallel_mesh_to_serial(serial_mesh, mesh);

    // VTUMeshWriter<SMesh> w;

    // std::cout << "Input Vector: " << std::endl;

    // if (!w.write("solution.vtu", serial_mesh, x_host)) {
    //     return false;
    // }

    // if (!w.write("rhs.vtu", serial_mesh, rhs_host)) {
    //     return false;
    // }

    // return true;
    int rank = 0;
    adios2::ADIOS adios;
    const int NSTEPS = 5;
    // random size per process, a different size at each step
    unsigned int Nelems;
    const size_t Nglobal = 6;
    // Application variables for output
    // random size per process, 5..10 each
    // v1 has different size on each process (but fixed over time)
    const unsigned int Nx = rand() % 6 + 5;

    std::vector<double> v0(Nglobal);
    std::vector<double> v1(Nx);
    // Local array, size is changing over time on each process
    std::vector<double> v2;

    // Local array, size is changing over time on each process
    // Also, random number of processes will write it at each step
    std::vector<double> &v3 = v2;

    // using SMesh = typename SerialMeshType<PMesh>::Type;
    // std::cout << "Writing results to disk..." << std::endl;

    // Integer n_nodes = mesh.n_nodes();

    // VectorReal::HostMirror x_host("x_host", n_nodes);
    // Kokkos::deep_copy(x_host, x);

    // SMesh serial_mesh;
    // convert_parallel_mesh_to_serial(serial_mesh, mesh);

    // w.write(path, serial_mesh, x_host)

    // Get io settings from the config file or
    // create one with default settings here
    adios2::IO io = adios.DeclareIO("Output");
    io.SetEngine("BP3");
    io.SetParameters({{"verbose", "4"}});

    /*
     * Define local array: type, name, local size
     * Global dimension and starting offset must be an empty vector
     * Here the size of the local array is the same on every process
     */
    adios2::Variable<double> varV0 = io.DefineVariable<double>("v0", {}, {}, {Nglobal});

    /*
     * v1 is similar to v0 but on every process the local size
     * is a different value
     */
    adios2::Variable<double> varV1 = io.DefineVariable<double>("v1", {}, {}, {Nx});

    /*
     * Define local array: type, name
     * Global dimension and starting offset must be an empty vector
     * but local size CANNOT be an empty vector.
     * We can use {adios2::UnknownDim} for this purpose or any number
     * actually since we will modify it before writing
     */
    adios2::Variable<double> varV2 = io.DefineVariable<double>("v2", {}, {}, {adios2::UnknownDim});

    /*
     * v3 is just like v2
     */
    adios2::Variable<double> varV3 = io.DefineVariable<double>("v3", {}, {}, {adios2::UnknownDim});

    adios2::Engine writer = io.Open("solution.bp", adios2::Mode::Write);

    for (int step = 0; step < NSTEPS; step++) {
        writer.BeginStep();
        // Write Step

        // // v0
        // for (size_t i = 0; i < Nglobal; i++) {
        //     v0[i] = rank * 1.0 + step * 0.1;
        // }
        // writer.Put<double>(varV0, v0.data());

        // // v1
        // for (size_t i = 0; i < Nx; i++) {
        //     v1[i] = rank * 1.0 + step * 0.1;
        // }
        // writer.Put<double>(varV1, v1.data());

        // v2

        // random size per process per step, 5..10 each
        Nelems = x.size();
        x.reserve(Nelems);
        for (size_t i = 0; i < Nelems; i++) {
            v2[i] = rank * 1.0 + step * 0.1;
        }

        // Set the size of the array now because we did not know
        // the size at the time of definition
        varV2.SetSelection(adios2::Box<adios2::Dims>({}, {Nelems}));
        writer.Put<double>(varV2, x.data());

        // v3

        // // random chance who writes it
        // unsigned int chance = rand() % 100;
        // /*if (step == 2)
        // {
        //     chance = 0;
        // }*/
        // bool doWrite = (chance > 60);
        // if (doWrite) {
        //     varV3.SetSelection(adios2::Box<adios2::Dims>({}, {Nelems}));
        //     writer.Put<double>(varV3, v3.data());
        // }

        writer.EndStep();
    }
    writer.Close();
    return true;
}
void Writer::close() { writer.Close(); }