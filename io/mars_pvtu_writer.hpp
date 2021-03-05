#ifndef MARS_PVTU_WRITER
#define MARS_PVTU_WRITER

#include <fstream>
#include <string>

#include <vtkSmartPointer.h>
// #include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
// #include <vtkXMLUnstructuredGridReader.h>
// #include <vtkDataSetMapper.h>
// #include <vtkActor.h>
// #include <vtkRenderer.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
// #include <vtkVertexGlyphFilter.h>

#include "mars_base.hpp"
#include "mars_config.hpp"
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_octant.hpp"
#include "mars_globals.hpp"
#include "mars_sfc_code.hpp"
#include "mars_sfc_generation.hpp"

namespace mars {

    template <class DM, Integer Type, class FEM>
    class PVTUMeshWriter {
    private:
        static const int VTU_TRIANGLE = 5;
        static const int VTU_QUAD = 9;

    public:
        bool write_vtu(const std::string &path, const DM &dm, const FEM &fe) {
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

            vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

            convert_points(dm, *unstructuredGrid);
            convert_cells(dm, fe, *unstructuredGrid);

            writer->SetFileName(path.c_str());
            writer->SetInputData(unstructuredGrid);
            writer->Write();
            return true;
        }

        bool write_vtu(const std::string &path, const DM &dm, const FEM &fe, const ViewVectorType<Real> &data) {
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

            vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

            convert_points(dm, *unstructuredGrid);
            convert_cells(dm, fe, *unstructuredGrid);

            auto fun = vtkSmartPointer<vtkDoubleArray>::New();

            // SFC<Type> dof = dm.get_local_dof_enum();

            Integer n = data.extent(0);
            fun->SetNumberOfValues(n);

            Kokkos::parallel_for("write_vtu", n, [&](const int i) { fun->SetValue(i, data(i)); });

            // Kokkos::parallel_for("for", dof.get_elem_size(), [&](const int i) {
            //     const Integer sfc_elem = dm.local_to_sfc(i);
            //     const Integer global_dof = dm.local_to_global(i);

            //     double point[3];
            //     get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(), dof.get_YDim(),
            //     dof.get_ZDim());

            //     fun->SetValue(global_dof, point[0] * point[1]);
            // });

            unstructuredGrid->GetPointData()->SetScalars(fun);
            // unstructuredGrid->GetPointData()->SetFi

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
                get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

                points->InsertNextPoint(point[0], point[1], point[2]);
            });

            unstructuredGrid.SetPoints(points);
        }

        void convert_cells(const DM &dm, const FEM &fe, vtkUnstructuredGrid &unstructuredGrid) {
            SFC<Type> dof = dm.get_local_dof_enum();

            vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();

            dm.elem_iterate([&](const Integer elem_index) {
                vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

                for (int i = 0; i < FEM::elem_nodes; i++) {
                    const Integer local_dof = fe.get_elem_local_dof(elem_index, i);
                    Dof d = dm.local_to_global_dof(local_dof);

                    const Integer g_id = d.get_gid();
                    quad->GetPointIds()->SetId(i, g_id);
                }

                cell_array->InsertNextCell(quad);
            });

            unstructuredGrid.SetCells(VTK_QUAD, cell_array);
        }
    };
}  // namespace mars

#endif MARS_PVTU_WRITER
