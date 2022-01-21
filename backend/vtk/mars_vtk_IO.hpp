#ifndef MARS_VTK_IO_HPP
#define MARS_VTK_IO_HPP

#include <fstream>
#include <string>

#include <vtkSmartPointer.h>
// #include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "mars_base.hpp"
#include "mars_config.hpp"
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_octant.hpp"
#include "mars_globals.hpp"
#include "mars_sfc_code.hpp"
#include "mars_sfc_generation.hpp"

// http://rotorbit.blogspot.com/2017/02/how-to-write-vtk-files-in-parallel.html

namespace mars {

    template <class DofHandler, class TpetraMultiView, typename LocalView>
    void owned_to_local(const DofHandler &dof_handler, const TpetraMultiView &owned_data, LocalView &local_data) {
        using namespace Kokkos;

        assert(dof_handler.get_owned_dof_size() == owned_data.extent(0));
        // const Integer size = dof_handler.get_global_dof_enum().get_elem_size();
        const Integer size = dof_handler.get_owned_dof_size();

        ViewVectorType<Integer> global_to_sfc = dof_handler.get_global_dof_enum().get_view_elements();
        // ViewVectorType<Integer> sfc_to_local = dof_handler.get_local_dof_enum().get_view_sfc_to_local();

        Kokkos::parallel_for(
            "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                const Integer base = dof_handler.compute_base(i);
                const Integer component = dof_handler.compute_component(i);
                const Integer sfc = global_to_sfc(base);

                // const Integer local = sfc_to_local(sfc);
                const Integer local = dof_handler.sfc_to_local(sfc);
                local_data(dof_handler.compute_block_index(local, component)) = owned_data(i, 0);
            });
    }

    namespace vtk {
        template <class FEDofMap>
        class IO {
        private:
            using DofHandler = typename FEDofMap::DHandler;
            static constexpr Integer Type = FEDofMap::ElemType;

            static const int VTU_TRIANGLE = 5;
            static const int VTU_QUAD = 9;
            static const int VTU_HEXAHEDRON = 12;

        public:
            IO(const FEDofMap &fe) : dm(fe.get_dof_handler()), fe(fe) {}

            template <class View>
            bool write(const std::string &path, const View &data) {
                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);

                if (comm_size == 1) {
                    if (write_vtu(path, dm, fe, data)) {
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    return write_pvtu(path, dm, fe, data);
                }
            }

        private:
            const DofHandler &dm;
            const FEDofMap &fe;

            std::string main_path, piece_path;

            void create_paths(const std::string &path) {
                main_path = path;
                piece_path = path;

                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);
                int rank = ::mars::rank(ctx);

                if (comm_size > 1) {
                    std::size_t index = path.find_last_of('.');
                    std::string extension = path.substr(index + 1);
                    std::string prefix = piece_path.substr(0, index) + std::to_string(comm_size);
                    assert(extension == "vtu");
                    piece_path = prefix;
                    piece_path += "_" + std::to_string(rank) + "." + "vtu";

                    main_path = prefix;
                    main_path += ".pvtu";

                    // std::cout << piece_path << "\n";
                }
            }

            void export_pvtu(vtkSmartPointer<vtkUnstructuredGrid> &unstructuredGrid) {
                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);
                int rank = ::mars::rank(ctx);

                if (comm_size > 1 && rank == 0) {
                    // auto pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

                    vtkNew<vtkXMLPUnstructuredGridWriter> pwriter;

                    pwriter->EncodeAppendedDataOff();
                    pwriter->SetFileName(main_path.c_str());
                    pwriter->SetNumberOfPieces(comm_size);
                    // pwriter->SetStartPiece(0);
                    // pwriter->SetEndPiece(comm_size - 1);

                    // #if VTK_MAJOR_VERSION <= 5
                    //                     pwriter->SetInput(unstructuredGrid);
                    // #else
                    pwriter->SetInputData(unstructuredGrid);
                    // #endif

                    // pwriter->Write();
                    pwriter->Update();
                }
            }

            // path = /dir/filename.vtu
            // in parallel it will export /dir/filename_<rank>.vtu and /dir/filename.pvtu

            bool write_vtu(const std::string &path) {
                create_paths(path);

                vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
                    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

                vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

                convert_points(dm, *unstructuredGrid);
                convert_cells(dm, fe, *unstructuredGrid);

                writer->SetFileName(piece_path.c_str());

#if VTK_MAJOR_VERSION <= 5
                writer->SetInput(unstructuredGrid);
#else
                writer->SetInputData(unstructuredGrid);
#endif
                writer->Write();

                export_pvtu(unstructuredGrid);
                return true;
            }

            template <class View>
            bool write_vtu(const std::string &path, const DofHandler &dm, const FEDofMap &fe, const View &data) {
                create_paths(path);

                vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

                convert_points(dm, *unstructuredGrid);
                convert_cells(dm, fe, *unstructuredGrid);
                add_data(*unstructuredGrid, dm, fe, data);

                // auto fun = vtkSmartPointer<vtkDoubleArray>::New();

                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);

                // Integer n;
                // if (comm_size > 1) {
                //     ViewVectorType<Real> local_data("local_data", dm.get_dof_size());
                //     ::mars::owned_to_local(dm, data, local_data);
                //     ::mars::gather_ghost_data(dm, local_data);

                //     typename ViewVectorType<Real>::HostMirror data_host = Kokkos::create_mirror_view(local_data);

                //     n = data_host.extent(0);
                //     fun->SetNumberOfValues(n);

                //     for (Integer i = 0; i < n; ++i) {
                //         fun->SetValue(i, data_host(i));
                //     }

                // } else {
                //     typename View::HostMirror data_host = Kokkos::create_mirror_view(data);

                //     n = data_host.extent(0);
                //     fun->SetNumberOfValues(n);

                //     for (Integer i = 0; i < n; ++i) {
                //         fun->SetValue(i, data_host(i, 0));
                //     }
                // }

                // fun->SetName("u");
                // unstructuredGrid->GetPointData()->SetScalars(fun);

                add_meta_data(*unstructuredGrid);

                {
                    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
                        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

                    writer->SetFileName(piece_path.c_str());
#if VTK_MAJOR_VERSION <= 5
                    writer->SetInput(unstructuredGrid);
#else
                    writer->SetInputData(unstructuredGrid);
#endif
                    writer->Write();
                }

                ctx->distributed->barrier();

                export_pvtu(unstructuredGrid);
                return true;
            }

            template <class View>
            void add_data(vtkUnstructuredGrid &unstructuredGrid,
                          const DofHandler &dm,
                          const FEDofMap &fe,
                          const View &data) {
                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);

                auto fun = vtkSmartPointer<vtkDoubleArray>::New();

                int block_size = dm.get_block();

                ViewVectorType<Real> local_data("local_data", dm.get_dof_size());
                ::mars::owned_to_local(dm, data, local_data);

                if (comm_size > 1) {
                    ::mars::gather_ghost_data(dm, local_data);
                }

                typename ViewVectorType<Real>::HostMirror data_host = Kokkos::create_mirror_view(local_data);

                Integer n = data_host.extent(0);
                fun->SetNumberOfValues(n);

                if (block_size == 1) {
                    fun->SetName("u");

                } else {
                    Integer n_blocks = n / block_size;
                    fun->SetNumberOfComponents(block_size);
                    // std::cout << "n_blocks: " << n_blocks << " block_size: " << block_size << "\n";
                    fun->SetName("vec");
                }

                for (Integer i = 0; i < n; ++i) {
                    fun->SetValue(i, data_host(i));
                }

                if (block_size == 1) {
                    unstructuredGrid.GetPointData()->SetScalars(fun);
                } else {
                    unstructuredGrid.GetPointData()->AddArray(fun);
                }
            }

            template <class View>
            bool write_pvtu(const std::string &path, const DofHandler &dm, const FEDofMap &fe, const View &data) {
                create_paths(path);

                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);
                int rank = ::mars::rank(ctx);

                auto writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
                writer->SetNumberOfPieces(comm_size);
                writer->SetStartPiece(rank);
                writer->SetEndPiece(rank);
                writer->SetWriteSummaryFile(rank == 0);

                vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

                convert_points(dm, *unstructuredGrid);
                convert_cells(dm, fe, *unstructuredGrid);

                add_data(*unstructuredGrid, dm, fe, data);

                add_meta_data(*unstructuredGrid);

                writer->SetFileName(main_path.c_str());
                writer->SetInputData(unstructuredGrid);

                writer->Write();
                return true;
            }

            void add_meta_data(vtkUnstructuredGrid &unstructuredGrid) {
                auto &ctx = dm.get_context();
                int comm_size = ::mars::num_ranks(ctx);
                int rank = ::mars::rank(ctx);

                auto n_elements = dm.get_elem_size();

                auto ranks = vtkSmartPointer<vtkIntArray>::New();
                ranks->SetNumberOfValues(n_elements);

                for (int i = 0; i < n_elements; ++i) {
                    ranks->SetValue(i, rank);
                }

                ranks->SetName("ranks");
                unstructuredGrid.GetCellData()->AddArray(ranks);
                unstructuredGrid.GetCellData()->SetActiveScalars("ranks");
            }

            void convert_points(const DofHandler &dm, vtkUnstructuredGrid &unstructuredGrid) {
                vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

                SFC<Type> dof = dm.get_local_dof_enum();

                auto global_dof_array = vtkSmartPointer<vtkIntArray>::New();
                global_dof_array->SetNumberOfValues(dof.get_elem_size());
                global_dof_array->SetName("global_dof_id");

                // const Integer base = dof_handler.compute_base(i);
                // const Integer component = dof_handler.compute_component(i);
                // const Integer sfc = global_to_sfc(base);

                Integer n_nodes = dof.get_elem_size();
                Integer block_size = dm.get_block();

                // std::cout << "n_nodes: " << n_nodes << "\n";

                Kokkos::parallel_for("for", n_nodes, [&](const int i) {
                    // const Integer base = dm.compute_base(i * block_size);
                    const Integer sfc_elem = dm.local_to_sfc(i * block_size);

                    // const Integer sfc_elem = dm.local_to_sfc(i);
                    const Integer global_dof = dm.local_to_global(i);

                    double point[3];
                    get_vertex_coordinates_from_sfc<Type>(
                        sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

                    points->InsertNextPoint(point[0], point[1], point[2]);

                    global_dof_array->SetValue(i, global_dof);
                });

                unstructuredGrid.SetPoints(points);

                unstructuredGrid.GetPointData()->AddArray(global_dof_array);
            }

            void convert_cells(const DofHandler &dm, const FEDofMap &fe, vtkUnstructuredGrid &unstructuredGrid) {
                SFC<Type> dof = dm.get_local_dof_enum();

                vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();

                int block_size = dm.get_block();

                dm.elem_iterate([&](const Integer elem_index) {
                    vtkSmartPointer<vtkCell> cell;

                    int n_nodes = 0;

                    if (DofHandler::Dim == 2) {
                        cell = vtkSmartPointer<vtkQuad>::New();
                        n_nodes = 4;
                    } else {
                        cell = vtkSmartPointer<vtkHexahedron>::New();
                        n_nodes = 8;
                    }

                    // std::cout << "[" << elem_index << "]: ";
                    for (int i = 0; i < n_nodes; i++) {
                        const Integer local_dof = fe.get_elem_local_dof(elem_index, i * block_size);
                        const Integer base = dm.compute_base(local_dof);

                        // std::cout << base << " ";

                        assert(i < 8);
                        cell->GetPointIds()->SetId(i, base);
                    }

                    // std::cout << "\n";

                    cell_array->InsertNextCell(cell);
                });

                if (DofHandler::Dim == 2) {
                    unstructuredGrid.SetCells(VTK_QUAD, cell_array);
                } else {
                    unstructuredGrid.SetCells(VTK_HEXAHEDRON, cell_array);
                }
            }
        };
    }  // namespace vtk
}  // namespace mars

#endif  // MARS_VTK_IO_HPP