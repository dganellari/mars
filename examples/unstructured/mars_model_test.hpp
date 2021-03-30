#ifndef MARS_MODEL_TEST_HPP
#define MARS_MODEL_TEST_HPP

#include <err.h>

#include <adios2.h>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>
#include <KokkosBlas1_sum.hpp>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>

#include <external/cxxopts/include/cxxopts.hpp>
// #include <Kokkos_sort.hpp>

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
#include "vtu_writer.hpp"

#include "mars_mesh_kokkos.hpp"

using namespace std::chrono;

namespace mars {

    template <class PMesh, class Op, class BC, class RHS, class AnalyticalFun>
    class ModelTest {
    public:
        static const int Dim = PMesh::Dim;
        static const int ManifoldDim = PMesh::ManifoldDim;
        static const int NFuns = ManifoldDim + 1;

        using SMesh = typename SerialMeshType<PMesh>::Type;
        using VectorReal = mars::ViewVectorType<Real>;
        using VectorInt = mars::ViewVectorType<Integer>;
        using VectorBool = mars::ViewVectorType<bool>;

        class Problem {
        public:
            Problem(PMesh &mesh) : mesh(mesh), op(mesh), bc(mesh), rhs("rhs", mesh.n_nodes()), write_output(true) {}

            bool solve(VectorReal &x) {
                x = VectorReal("x", mesh.n_nodes());
                bc.apply(x, bc_fun);

                // auto prec_ptr = op.preconditioner();
                auto prec_ptr = std::make_shared<CopyOperator>();

                Integer num_iter = 0;

                auto start = steady_clock::now();

                bool ok = bcg_stab(op, *prec_ptr, rhs, rhs.extent(0), x, num_iter);

                auto end = steady_clock::now();
                auto diff = end - start;

                std::cout << "bcg_stab/time: " << std::chrono::duration<double>(diff).count() << " s" << std::endl;

                // FILE *file_time = fopen("timing.csv", "a");
                // fprintf(file_time, "%ld, %.4fs\n", mesh.n_nodes(), std::chrono::duration<double>(diff).count());
                // fclose(file_time);

                // if (ok) {
                if (write_output) {
                    return write(x) && ok;
                } else {
                    return ok;
                }

                // } else {
                // std::cerr << "No OK" << std::endl;
                // return false;
                // }
            }

            bool write(VectorReal &x) {
                // std::cout << "Writing results to disk..." << std::endl;

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

                std::string extent = std::to_string()

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

            bool measure_actual_error(VectorReal &x) {
                /////////////////////////////////////////////////////////////////////////////
                // Compute Error
                Integer n_nodes = mesh.n_nodes();

                VectorReal x_exact("X_exact", n_nodes);
                VectorReal diff("Diff", n_nodes);
                VectorReal Ax("Ax", n_nodes);

                Interpolate<PMesh> interp(mesh);
                interp.apply(x_exact, an_fun);

                Kokkos::deep_copy(diff, x_exact);

                KokkosBlas::axpy(-1.0, x, diff);
                Real err = KokkosBlas::nrm2(diff) / KokkosBlas::nrm2(x_exact);

                std::cout << "=====================================" << std::endl;

                std::cout << "n_elements:        " << mesh.n_elements() << std::endl;
                std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
                std::cout << "n_nodes:           " << mesh.n_nodes() << std::endl;

                std::cout << "err : " << err << std::endl;

                ////////////////////////////////////////////////////////////////////////////
                // Compute op (x_exact) - rhs

                op.apply(x_exact, diff);

                Real sum_rhs = KokkosBlas::sum(rhs);

                Real norm_lapl_x = KokkosBlas::sum(diff);
                std::cout << "sum(Op(x)): " << norm_lapl_x << " == " << sum_rhs << std::endl;

                // Diff
                KokkosBlas::axpy(-1.0, rhs, diff);
                err = KokkosBlas::nrm2(diff) / KokkosBlas::nrm2(rhs);
                std::cout << "Op(x_exact) - rhs : " << err << std::endl;

                std::cout << "=====================================" << std::endl;

                // if (write_output) {
                VTUMeshWriter<SMesh> w;
                VectorReal::HostMirror x_exact_host("x_exact_host", mesh.n_nodes());
                Kokkos::deep_copy(x_exact_host, x_exact);
                return w.write("x_exact.vtu", serial_mesh, x_exact_host);
                // }

                return true;
            }

            void init() {
                op.init();
                auto id = std::make_shared<IdentityOperator>();
                id->init(mesh, bc_fun);
                op.set_identity(id);
                op.assemble_rhs(rhs, rhs_fun);
                bc.apply(rhs, bc_fun);

                if (write_output) {
                    convert_parallel_mesh_to_serial(serial_mesh, mesh);
                }
            }

            BC bc_fun;
            RHS rhs_fun;
            AnalyticalFun an_fun;

            PMesh &mesh;
            Op op;
            BoundaryConditions<PMesh> bc;
            VectorReal rhs;

            SMesh serial_mesh;
            bool write_output;
        };

        bool interpolate(PMesh &mesh, FEValues<PMesh> &values, const Integer &new_elem_offset, VectorReal &x) {
            VectorReal refined_x("refined_x", mesh.n_nodes());
            interpolate(mesh, values, new_elem_offset, x, refined_x);
            x = refined_x;
        }

        bool interpolate(PMesh &mesh,
                         FEValues<PMesh> &values,
                         const Integer &new_elem_offset,
                         const VectorReal &x,
                         VectorReal &refined_x) {
            Integer n_elements = mesh.n_elements();
            Integer n_new = n_elements - new_elem_offset;
            if (n_new == 0) return false;

            // auto parents = mesh.parent_map();

            VectorInt parents("parents", n_elements, -1);

            Kokkos::parallel_for(
                new_elem_offset, MARS_LAMBDA(const Integer &i) {
                    auto c = mesh.get_children(i);
                    if (c.is_valid()) {
                        assert(c(0) < n_elements);
                        assert(c(1) < n_elements);
                        assert(c(0) > 0);
                        assert(c(1) > 0);

                        parents(c(0)) = i;
                        parents(c(1)) = i;
                    }
                });

            refined_x = VectorReal("refined_x", mesh.n_nodes());

            Kokkos::parallel_for(
                x.extent(0), MARS_LAMBDA(const Integer &i) { refined_x(i) = x(i); });

            {
                auto elems = mesh.get_view_elements();
                auto points = mesh.get_view_points();
                auto J_inv = values.J_inv();

                Kokkos::parallel_for(
                    n_new, MARS_LAMBDA(const Integer &iter) {
                        const Integer elem_id = new_elem_offset + iter;
                        assert(elem_id < parents.extent(0));

                        Integer parent = parents(elem_id);
                        assert(parent != -1);

                        Real J_inv_e[Dim * Dim];
                        Real tr[Dim], p[Dim], p_ref[Dim];
                        Integer idx[NFuns];
                        Real u_e[NFuns];

                        for (int k = 0; k < Dim * Dim; ++k) {
                            assert(parent < J_inv.extent(0));

                            J_inv_e[k] = J_inv(parent, k);
                            assert(J_inv_e[k] == J_inv_e[k]);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            assert(parent < elems.extent(0));

                            idx[k] = elems(parent, k);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            u_e[k] = x(idx[k]);
                        }

                        Integer n0 = elems(parent, 0);
                        assert(n0 < points.extent(0));

                        for (int d = 0; d < Dim; ++d) {
                            tr[d] = points(n0, d);
                        }

                        for (int k = 0; k < NFuns; ++k) {
                            assert(elem_id < elems.extent(0));

                            Integer c_n_id = elems(elem_id, k);

                            for (int d = 0; d < Dim; ++d) {
                                assert(c_n_id < points.extent(0));

                                p[d] = points(c_n_id, d) - tr[d];
                            }

                            Algebra<Dim>::mv_mult(J_inv_e, p, p_ref);

                            Real val = FESimplex<Dim>::fun(p_ref, u_e);
                            refined_x(c_n_id) = val;
                        }
                    });
            }
        }

        bool refine(PMesh &mesh, FEValues<PMesh> &values, const Real tol, VectorReal &x, const bool write_output) {
            ParallelBisection<PMesh> bisection(&mesh);

            VectorReal error;
            GradientRecovery<PMesh> grad_rec;
            grad_rec.estimate(values, x, error);

            if (KokkosBlas::sum(error) < tol) {
                std::cout << "[Status] Reached desired tolerance" << std::endl;
                return false;
            }

            const Real alpha = 0.3;
            const Real max_error = KokkosBlas::nrminf(error);
            const Real low_bound_err = alpha * max_error;

            Integer old_n_elements = mesh.n_elements();

            VectorBool marked("marked", mesh.n_elements());

            Kokkos::parallel_for(
                mesh.n_elements(), MARS_LAMBDA(const Integer &i) { marked(i) = (error(i) >= low_bound_err) * i; });

            VectorInt marked_list = mark_active(marked);

            Kokkos::parallel_for(
                marked_list.extent(0), MARS_LAMBDA(const Integer &i) {
                    assert(mesh.is_active(marked_list(i)));
                    assert(marked_list(i) < mesh.n_elements());
                });

            bisection.refine(marked_list);

            Integer n_new = mesh.n_elements() - old_n_elements;
            if (n_new == 0) return false;

            // x = VectorReal("x", mesh.n_nodes());

            // FIXME
            interpolate(mesh, values, old_n_elements, x);

            std::cout << "n_elements:        " << mesh.n_elements() << std::endl;
            std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
            std::cout << "n_nodes:           " << mesh.n_nodes() << std::endl;

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (write_output) {
                VectorReal::HostMirror refined_x_host("refined_x_host", mesh.n_nodes());
                Kokkos::deep_copy(refined_x_host, x);

                SMesh serial_mesh;
                convert_parallel_mesh_to_serial(serial_mesh, mesh);

                VTUMeshWriter<SMesh> w;
                w.write("mesh_refined.vtu", serial_mesh, refined_x_host);
            }

            return true;
        }

        bool solve(PMesh &mesh) {
            const Integer n_nodes = mesh.n_nodes();
            VectorReal x("X", n_nodes);

            Integer refs{0};

            // do {
            Problem problem(mesh);
            problem.write_output = write_output;
            problem.init();

            if (problem.solve(x)) {
                problem.measure_actual_error(x);
            } else {
                return false;
            }

            // if (max_refinements == 0 || !use_adaptive_refinement) return true;
            /*


                            if (use_adaptive_refinement) {
                                refined = refine(mesh, problem.op.values(), tol, x, write_output);
                                mesh.clean_up();
                            } else {
                                Integer elem_offset = mesh.n_elements();
                                ParallelBisection<PMesh>(&mesh).uniform_refine(1);
                                interpolate(mesh, problem.op.values(), elem_offset, x);
                                refined = true;
                            }

                            if (reset_x_to_zero) {
                                // reset x to 0
                                Kokkos::parallel_for(
                                    x.extent(0), MARS_LAMBDA(const Integer &i) { x(i) = 0.0; });
                            }

                            if (!refined) {
                                std::cout << "Did not refine" << std::endl;
                                break;
                            } */

            // } while (++refs <= max_refinements);
            return true;
        }

        void run(cxxopts::ParseResult &args) {
            Integer ns[4] = {
                args["nx"].as<Integer>(), args["ny"].as<Integer>(), args["nz"].as<Integer>(), args["nt"].as<Integer>()};

            use_adaptive_refinement = args["adaptive"].as<bool>();
            write_output = args["output"].as<bool>();
            max_refinements = args["refine_level"].as<Integer>();

            PMesh mesh;

            for (int d = Dim; d < 3; ++d) {
                ns[d] = 0;
            }

            if (Dim <= 3) {
                generate_cube(mesh, ns[0], ns[1], ns[2]);
            } else {
                assert(false);
                // SMesh smesh;
                // read_mesh("../data/cube4d_24.MFEM", smesh);
                // convert_serial_mesh_to_parallel(mesh, smesh);

                // ParallelBisection<PMesh>(&mesh).uniform_refine(3);
                /* ParallelBisection<PMesh>(&mesh).uniform_refine(4); */
                /* write_output = false; */
            }

            // ParallelBisection<PMesh> bisection(&mesh);
            // bisection.uniform_refine(5);

            solve(mesh);
        }

    private:
        Real tol{1e-6};
        Integer max_refinements{0};

        bool refined{false};
        bool use_adaptive_refinement{false};
        bool reset_x_to_zero{false};
        bool write_output{true};

    };  // namespace mars
}  // namespace mars

#endif  // MARS_MODEL_TEST_HPP
