#ifndef MARS_PVTU_WRITER
#define MARS_PVTU_WRITER

namespace mars{


	template<class DM,  Integer Type>
	class PVTUMeshWriter {

	private: 
		static const int VTU_TRIANGLE = 5;
		static const int VTU_QUAD = 9;

	public:

	    bool write_vtu(const std::string &path, const DM &dm) {
	        std::ofstream os;
	        os.open(path.c_str());
	        if (!os.good()) {
	            os.close();
	            return false;
	        }

	        // writeHeader(mesh, os);
	        
			// writeCellData(mesh, sol, 1, os);
			// writePointData(mesh, mesh.points(), sol, 1, os, coord);
		
	        writePoints(dm, os);
	        writeCells(dm, os);

	        writeFooter(dm, os);
	        os.close();
	        // clear();
	        return true;
    	}


    	void writePoints(const DM &dm, std::ostream &os) {
    		os << "vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();" << "\n";

    		SFC<Type> dof = dm.get_local_dof_enum();

	    	Kokkos::parallel_for(
		      "for", dof.get_elem_size(), [&](const int i) {
		        const Integer sfc_elem = dm.local_to_sfc(i);
		        const Integer global_dof = dm.local_to_global(i);

		        double point[3];
		        get_vertex_coordinates_from_sfc<Type>(sfc_elem, point, dof.get_XDim(),
		                                              dof.get_YDim(), dof.get_ZDim());

		        // printf("dof: %li - gdof: %li --- (%lf, %lf) - rank: %i\n", i,
		        //        global_dof, point[0], point[1], rank);

		        // std::cout << "points->InsertNextPoint(" << point[0] << "," << point[1] << std::endl;
		        
		        os << "points->InsertNextPoint(" << point[0] << "," << point[1] << "," << point[2] << ")" <<";\n";

		    });

    	}


    	void writeCells(const DM &dm, std::ostream &os) {

    		SFC<Type> dof = dm.get_local_dof_enum();
    		dm.elem_iterate([&](const Integer elem_index) {
    			os <<"vtkSmartPointer<vtkQuad> quad"<<std::to_string(elem_index)<<"= vtkSmartPointer<vtkQuad>::New();" << "\n";

    			for (int i = 0; i < DMQ2::elem_nodes; i++) {
    				const Integer local_dof = dm.get_elem_local_dof(elem_index, i);
    				Dof d = dm.local_to_global_dof(local_dof);

    				int c = elem_index+i;

    				os << "quad"+std::to_string(elem_index)+"->GetPointIds()->SetId("<< c <<","<< d.get_gid() <<");"<< "\n";

    			}

    		});

    		os << "vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();\n";

    		dm.elem_iterate([&](const Integer elem_index) {

    			os<< "cellArray->InsertNextCell(quad" << std::to_string(elem_index) << ");\n";
    		});


    	}


    	void writeFooter(const DM &dm, std::ostream &os) {

    		os << "vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();\n";

    		os << "unstructuredGrid->SetPoints(points);\n";

    		os << "unstructuredGrid->SetCells(VTK_QUAD, cellArray);\n";





    	}





	};
}


#endif MARS_PVTU_WRITER

