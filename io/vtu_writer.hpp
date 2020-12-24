#ifndef VTU_WRITER_H
#define VTU_WRITER_H

#include <KokkosSparse_spmv.hpp>
#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "mars_base.hpp"
#include "mars_mesh.hpp"

#include <cassert>

namespace mars {

    template <class Mesh>
    class VTUMeshWriter {
    private:
        using KokkosVector = mars::ViewVectorType<Real>::HostMirror;
        // typedef typename Mesh::Matrix Matrix;
        static const Integer Dim = Mesh::Dim;
        static const Integer ManifoldDim = Mesh::ManifoldDim;
        static const Integer NNodes = Mesh::Elem::NNodes;

        static const int VTU_TETRA = 10;
        static const int VTU_TRIANGLE = 5;
        static const int VTU_QUAD = 9;
        static const int VTU_HEXAHEDRON = 12;
        static const int VTU_POLYGON = 7;
        static const int VTU_LINE = 3;

    public:
        inline static int VTKTagVolume(const Integer nVertices) {
            switch (nVertices) {
                case 4:
                    return VTU_TETRA;
                case 8:
                    return VTU_HEXAHEDRON;
                default:
                    assert(
                        false);  // element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
                    return -1;
            }
        }

        inline static int VTKTag(const Integer nVertices) {
            switch (nVertices) {
                case 2:
                    return VTU_LINE;
                default:
                    std::cerr << "[Error] " << nVertices << " not supported" << std::endl;
                    assert(
                        false);  // element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
                    return -1;
            }
        }

        inline static int VTKTagPlanar(const Integer nVertices) {
            if (nVertices > 4) {
                return VTU_POLYGON;
            }

            switch (nVertices) {
                case 3:
                    return VTU_TRIANGLE;
                case 4:
                    return VTU_QUAD;
                default:
                    std::cerr << "[Error] " << nVertices << " not supported" << std::endl;
                    assert(
                        false);  // element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
                    return -1;
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////

        bool write(const std::string &path,
                   const Mesh &mesh,
                   const KokkosVector &sol,
                   bool cell = false,
                   const bool &coord = true) {
            std::ofstream os;
            os.open(path.c_str());
            if (!os.good()) {
                os.close();
                return false;
            }

            writeHeader(mesh, os);

            if (cell) {
                writeCellData(mesh, sol, 1, os);
            } else {
                writePointData(mesh, mesh.points(), sol, 1, os, coord);
            }

            writePoints(mesh, os);
            writeCells(mesh, os);

            writeFooter(mesh, os);
            os.close();
            // clear();
            return true;
        }

        void writeHeader(const Mesh &mesh, std::ostream &os) {
            Integer nCells = mesh.n_active_elements();
            os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
            os << "<UnstructuredGrid>\n";
            os << "<Piece NumberOfPoints=\"" << mesh.n_nodes() << "\" NumberOfCells=\"" << nCells << "\">\n";
        }

        void writePointData(const Mesh &mesh,
                            const std::vector<typename Mesh::Point> &points,
                            const KokkosVector &sol,
                            const Integer &nc,
                            std::ostream &os,
                            const bool coord) {
            Integer nNodes = mesh.n_nodes();
            os << "<PointData Scalars=\"scalars\">\n";

            if (nc == 1) {
                os << "<DataArray type=\"Float32\"  Name=\"solution\" Format=\"ascii\">\n";
                for (int i = 0; i < nNodes; ++i) {
                    if (std::abs(sol(i)) < 1e-30) {
                        os << 0 << "\n";
                    } else {
                        os << sol(i) << "\n";
                    }
                    // os << sol(i) << "\n";
                }
                os << "</DataArray>\n";

                if (coord == true) {
                    os << "<DataArray type=\"Float32\"  Name=\"X_coord\" Format=\"ascii\">\n";
                    for (Integer i = 0; i < points.size(); ++i) {
                        os << points[i](0) << "\n";

                        // os << sol(i) << "\n";
                    }
                    os << "</DataArray>\n";

                    os << "<DataArray type=\"Float32\"  Name=\"Y_coord\" Format=\"ascii\">\n";
                    for (Integer i = 0; i < points.size(); ++i) {
                        os << points[i](1) << "\n";

                        // os << sol(i) << "\n";
                    }
                    // os  << "</DataArray>\n";

                    // os << "<DataArray type=\"Float32\"  Name=\"Z_coord\" Format=\"ascii\">\n";
                    // for (Integer i = 0; i < points.size(); ++i) {
                    // 	os << points[i](2) << "\n";

                    // 	// os << sol(i) << "\n";
                    // }
                }

            } else {
                os << "<DataArray type=\"Float32\"  Name=\"solution\" NumberOfComponents=\"" << nc << "\">\n";
                for (int i = 0; i < nNodes; ++i) {
                    for (int j = 0; j < nc; ++j) {
                        os << sol(i) << " ";
                    }
                    os << "\n";
                }
            }

            os << "</DataArray>\n";
            os << "</PointData>\n";
        }

        void writeCellData(const Mesh &mesh, const KokkosVector &sol, const Integer &nc, std::ostream &os) {
            Integer nCells = mesh.n_elements();
            os << "<CellData Scalars=\"scalars\">\n";

            if (nc == 1) {
                os << "<DataArray type=\"Float32\"  Name=\"solution_cell\" Format=\"ascii\">\n";
                for (int i = 0; i < nCells; ++i) {
                    if (std::abs(sol(i)) < 1e-30) {
                        os << 0 << "\n";
                    } else {
                        os << sol(i) << "\n";
                    }
                }
            } else {
                os << "<DataArray type=\"Float32\" NumberOfComponents=\"" << nc << "\" Format=\"ascii\">\n";
                for (int i = 0; i < nCells; ++i) {
                    for (int j = 0; j < nc; ++j) {
                        os << sol(i) << " ";
                    }
                    os << "\n";
                }
            }
            os << "</DataArray>\n";
            os << "</CellData>\n";
        }

        void writePoints(const std::vector<typename Mesh::Point> &points, std::ostream &os) {
            os << "<Points>\n";
            os << "<DataArray type=\"Float32\" NumberOfComponents=\"" << ((Dim < 3) ? 3 : Dim)
               << "\" format=\"ascii\">\n";
            for (Integer i = 0; i < points.size(); ++i) {
                for (Integer d = 0; d < Dim; ++d) {
                    os << points[i](d);
                    if (d < Dim - 1) {
                        os << " ";
                    } else if (Dim == 2) {  // padding for paraview vtu format visualisation.
                        os << " ";
                        os << 0;

                    } else if (Dim == 1) {
                        os << " ";
                        os << 0;
                        os << " ";
                        os << 0;
                    }
                }
                os << "\n";
            }

            os << "</DataArray>\n";
            os << "</Points>\n";
        }

        void writePoints(const Mesh &mesh, std::ostream &os) { writePoints(mesh.points(), os); }

        void writeCells(const Mesh &mesh, std::ostream &os) {
            const auto n_active_elements = mesh.n_active_elements();

            os << "<Cells>\n";
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            os << "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n";

            for (Integer k = 0; k < mesh.n_elements(); ++k) {
                if (!mesh.is_active(k)) continue;

                for (Integer i = 0; i < NNodes; ++i) {
                    const Integer v = mesh.elem(k).nodes[i];
                    os << v;
                    if (i < NNodes - 1) {
                        os << " ";
                    }
                }
                os << "\n";
            }

            os << "</DataArray>\n";
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            int minTag, maxTag;
            if (ManifoldDim == 1) {
                minTag = VTKTag(NNodes);
                maxTag = VTKTag(NNodes);
            } else if (ManifoldDim == 2) {
                minTag = VTKTagPlanar(NNodes);
                maxTag = VTKTagPlanar(NNodes);
            } else {
                minTag = VTKTagVolume(NNodes);
                maxTag = VTKTagVolume(NNodes);
            }

            os << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\" RangeMin=\"" << minTag << "\" RangeMax=\""
               << maxTag << "\">\n";
            for (Integer i = 0; i < n_active_elements; ++i) {
                if (ManifoldDim == 3) {
                    os << VTKTagVolume(NNodes) << "\n";
                } else if (ManifoldDim == 2) {
                    os << VTKTagPlanar(NNodes) << "\n";
                } else
                    os << VTKTag(NNodes) << "\n";
            }

            os << "</DataArray>\n";

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            os << "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"" << (NNodes)
               << "\" RangeMax=\"" << (n_active_elements * (NNodes)) << "\">\n";

            for (Integer i = 0, offset = (NNodes); i < n_active_elements; ++i, offset += (NNodes)) {
                os << offset << "\n";
            }

            os << "</DataArray>\n";
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            os << "</Cells>\n";
        }

        void writeFooter(const Mesh &mesh, std::ostream &os) {
            os << "</Piece>\n";
            os << "</UnstructuredGrid>\n";
            os << "</VTKFile>\n";
        }

        virtual ~VTUMeshWriter() {
            std::ofstream os;
            os.close();
        }
    };
}  // namespace mars

#endif  // VTU_WRITER_H
