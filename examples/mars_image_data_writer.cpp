#include "mars_base.hpp"
#include "mars_boundary_conditions.hpp"
#include "mars_copy_operator.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_gradient_recovery.hpp"
#include "mars_identity_operator.hpp"
#include "mars_interpolate.hpp"
#include "mars_invert.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_serial_mesh_type.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_laplace.hpp"

template <class PMesh>
class ImageWriter {
public:
    using VectorReal = mars::ViewVectorType<Real>;
    ImageWriter(PMesh &mesh) {}

    bool write(VectorReal &x) {  // std::cout << "Writing results to disk..." << std::endl;

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

        using SMesh = typename SerialMeshType<PMesh>::Type;
        std::cout << "Writing results to disk..." << std::endl;

        Integer n_nodes = mesh.n_nodes();

        VectorReal::HostMirror x_host("x_host", n_nodes);
        Kokkos::deep_copy(x_host, x);

        SMesh serial_mesh;
        convert_parallel_mesh_to_serial(serial_mesh, mesh);

        // w.write(path, serial_mesh, x_host)

        // Get io settings from the config file or
        // create one with default settings here
        adios2::IO io = adios.DeclareIO("Output");
        io.SetEngine("BP3");
        io.SetParameters({{"verbose", "4"}});

        std::string extent =
            std::to_string()

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
}