#ifndef MARS_VTK_IO_HPP
#define MARS_VTK_IO_HPP

#include <fstream>
#include <string>

#include <vtkSmartPointer.h>
// #include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "mars_base.hpp"
#include "mars_config.hpp"
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_octant.hpp"
#include "mars_globals.hpp"
#include "mars_sfc_code.hpp"
#include "mars_sfc_generation.hpp"

namespace mars {

    namespace vtk {
        template <class DM, class FEM, Integer Type = FEM::ElemType>
        class IO {
        private:
            static const int VTU_TRIANGLE = 5;
            static const int VTU_QUAD = 9;
            static const int VTU_HEXAHEDRON = 12;

        public:
            IO(const DM &dm, const FEM &fe) : dm(dm), fe(fe) {}

            template <class View>
            bool write(const std::string &path, const View &data) {
                return write_vtu(path, dm, fe, data);
            }

        private:
            const DM &dm;
            const FEM &fe;

            bool write_vtu(const std::string &path) {
                vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
                    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

                vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

                convert_points(dm, *unstructuredGrid);
                convert_cells(dm, fe, *unstructuredGrid);

                writer->SetFileName(path.c_str());
                writer->SetInputData(unstructuredGrid);
                writer->Write();
                return true;
            }

            template <class View>
            bool write_vtu(const std::string &path, const DM &dm, const FEM &fe, const View &data) {
                vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
                    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

                vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

                convert_points(dm, *unstructuredGrid);
                convert_cells(dm, fe, *unstructuredGrid);

                auto fun = vtkSmartPointer<vtkDoubleArray>::New();

                Integer n = data.extent(0);
                fun->SetNumberOfValues(n);

                typename View::HostMirror data_host = Kokkos::create_mirror_view(data);

                for (Integer i = 0; i < n; ++i) {
                    fun->SetValue(i, data_host(i, 0));
                }

                unstructuredGrid->GetPointData()->SetScalars(fun);

                writer->SetFileName(path.c_str());
                writer->SetInputData(unstructuredGrid);
                writer->Write();
                return true;
            }

            void convert_points(const DM &dm, vtkUnstructuredGrid &unstructuredGrid) {
                vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

                SFC<Type> dof = dm.get_local_dof_enum();

                Kokkos::parallel_for("for", dof.get_elem_size(), [&](const int i) {
                    const Integer sfc_elem = dm.local_to_sfc(i);
                    const Integer global_dof = dm.local_to_global(i);

                    double point[3];
                    get_vertex_coordinates_from_sfc<Type>(
                        sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

                    points->InsertNextPoint(point[0], point[1], point[2]);
                });

                unstructuredGrid.SetPoints(points);
            }

            void convert_cells(const DM &dm, const FEM &fe, vtkUnstructuredGrid &unstructuredGrid) {
                SFC<Type> dof = dm.get_local_dof_enum();

                vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();

                dm.elem_iterate([&](const Integer elem_index) {
                    vtkSmartPointer<vtkCell> cell;

                    if (DM::Dim == 2) {
                        cell = vtkSmartPointer<vtkQuad>::New();
                    } else {
                        cell = vtkSmartPointer<vtkHexahedron>::New();
                    }

                    for (int i = 0; i < FEM::elem_nodes; i++) {
                        const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                        Dof d = dm.local_to_global_dof(local_dof);

                        const Integer g_id = d.get_gid();
                        cell->GetPointIds()->SetId(i, g_id);
                    }

                    cell_array->InsertNextCell(cell);
                });

                if (DM::Dim == 2) {
                    unstructuredGrid.SetCells(VTK_QUAD, cell_array);
                } else {
                    unstructuredGrid.SetCells(VTK_HEXAHEDRON, cell_array);
                }
            }
        };
    }  // namespace vtk
}  // namespace mars

#endif  // MARS_VTK_IO_HPP