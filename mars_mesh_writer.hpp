#ifndef MARS_MESH_WRITER_HPP
#define MARS_MESH_WRITER_HPP

#include "mars_simplex.hpp"
#include "mars_mesh.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <algorithm> 

namespace mars {
    template<Integer Dim, Integer ManifoldDim>
    class MeshWriter {
    public: 
        virtual ~MeshWriter() {}
        virtual bool write(const std::string &path, const Mesh<Dim, ManifoldDim> &mesh) = 0; 
    };


    template<Integer Dim, Integer ManifoldDim>
	class MFEMWriter final : public MeshWriter<Dim, ManifoldDim> {
	public:
        bool write(const std::string &path, const Mesh<Dim, ManifoldDim> &mesh) override
        {
            std::ofstream os;
    		os.open(path.c_str());
    		if (!os.good()) {
    			os.close();
    			return false;
    		}

			write_header(mesh, os);
			write_elements(mesh, os);
			write_vertices(mesh, os);

			os.close();
		//	clear();
			return true;

        }

        //////////////// Format ////////////////////////

        // template<Integer Dim, Integer ManifoldDim>
	    void write_header(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) 
        {
            Integer dim = mesh.ManifoldDim;
            os << "MFEM mesh v1.0\n\ndimension\n";
            os << dim<<"\n\n";
	    }  

        ////////////////// Elements ////////////////////////

        // template<Integer Dim, Integer ManifoldDim>
	    void write_elements(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) 
        {
            os << "elements\n";
            os << mesh.n_elements() << "\n";

            for (Integer k = 0; k < mesh.n_elements(); ++k) {
                if (!mesh.is_active(k))
                    continue;

                if(mesh.tags().size() == mesh.n_elements())
                    os<<mesh.tags()[k]<< " " << mesh.elem(k).type()<<" ";
                else
                    os<<INVALID_INDEX<< " " << INVALID_INDEX<<" ";

                for (Integer i = 0; i < ManifoldDim + 1; ++i) {
                    const Integer v = mesh.elem(k).nodes[i];
                    os << v;
                    if (i < ManifoldDim) {
                        os << " ";
                    }
                }
                os << "\n";
            }
            os << "\n";
	    }

        ////////////////// Nodes ////////////////////////

        // template<Integer Dim, Integer ManifoldDim>
	    void write_vertices(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) 
        {
            Integer dim = mesh.Dim;
            os << "vertices\n";

            const std::vector<Vector<Real, Dim>> points = mesh.points();

            os << points.size() << "\n";
            os << dim << "\n";

            for (Integer i = 0; i < points.size(); ++i) {
                for (Integer d = 0; d < Dim; ++d) {
                    os << points[i](d);
                    if (d < Dim - 1) {
                        os << " ";
                    }
                }
                os << "\n";
            }

	    }


    }; 

    template<Integer Dim, Integer ManifoldDim>
	class MSHWriter final : public MeshWriter<Dim, ManifoldDim> {
	public: 
        bool write(const std::string &path, const Mesh<Dim, ManifoldDim> &mesh) override
        {
            std::cout << mesh.n_nodes() << "\n";
            std::ofstream os;
    		os.open(path.c_str());
    		if (!os.good()) {
    			os.close();
    			return false;
    		}

            write_format(mesh, os);
			write_nodes(mesh, os);
			write_elements_msh(mesh, os);

			os.close();
		//	clear();
			return true;
        }

        //////////////// Format ////////////////////////

        // template<Integer Dim, Integer ManifoldDim>
	    void write_format(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) 
        {
            os << "$MeshFormat\n";
            os << "4.1 0 8\n";
            os << "$EndMeshFormat\n";

	    }

        ////////////////// Nodes ////////////////////////

        // template<Integer Dim, Integer ManifoldDim>
	    void write_nodes(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os) 
        {
            Integer dim = mesh.ManifoldDim;
            
            const std::vector<Vector<Real, Dim>> points = mesh.points();
            
            os << "$Nodes\n";
            //first line 
            os << "1" << " " << points.size() << " " << "1" << " " <<  points.size() << "\n";
            // second line 
            os << dim << " " << "1" << " " << "0" << " " <<   points.size() << "\n"; 

            // nodes id 
            for(Integer i = 0; i <  points.size(); ++i){
                
                os << i+1 << "\n";  
            }

            for (Integer i = 0; i <  points.size(); ++i) {
                for (Integer d = 0; d < Dim; ++d) {
                    os << points[i](d);
                    if (d < Dim - 1) {
                        os << " ";
                    }
                }
                // os << " " << 0  << "\n";

                for(Integer k = Dim; k < 3; ++k) {
                                os << " " << 0;
                            }


                os << "\n";
            }
            os << "$EndNodes\n";

	    }

    ////////////////// Elements ////////////////////////

        // template<Integer Dim, Integer ManifoldDim>
        void write_elements_msh(const Mesh<Dim, ManifoldDim> &mesh, std::ostream &os)
        {
            os << "$Elements\n";
            os << "1"  <<  " " << mesh.n_elements() << " " <<  "1" << " " << mesh.n_elements() << "\n";

            int t = -1; 
            if (ManifoldDim == 1) {
                t = 1;
            } else if (ManifoldDim == 2){
                t = 2; 

            }else if (ManifoldDim == 3){
                t = 4; 
            }

            os << ManifoldDim << " " << "1" << " " << t << " " << mesh.n_elements() << "\n";

            for (Integer k = 0; k < mesh.n_elements(); ++k) {
                if (!mesh.is_active(k))
                    continue;

                os << k+1 << " "; 

                for (Integer i = 0; i < ManifoldDim + 1; ++i) {
                    const Integer v = mesh.elem(k).nodes[i];
                    os << v+1;
                    if (i < ManifoldDim) {
                        os << " ";
                    }
                }
                



                os << "\n";
            }
            os << "$EndElements\n";
        }
    }; 



    template<Integer Dim, Integer ManifoldDim>
    bool write(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false) 
        {
            auto ext = path.substr(path.find_last_of(".") + 1);

            if (ext == "MFEM"){
                MFEMWriter<Dim, ManifoldDim> w2; 
                w2.write(path, mesh); 
            } else if (ext == "msh" ){
                MSHWriter<Dim, ManifoldDim> w1; 
                w1.write(path,mesh);  
            } else {
                std::cout << "ERROR: Only MFEM and msh formats are supported" << std::endl; 
            }
            

        return true;
        }


}

#endif //MARS_MESH_WRITER_HPP
