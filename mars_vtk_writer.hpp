#ifndef GEOM_VTK_WRITER_H
#define GEOM_VTK_WRITER_H

#include "mars_base.hpp"
#include "mars_mesh.hpp"

#include <cassert>

namespace mars {
 //    template<class Matrix>
	// class VTKDataNode {
	// private:
	// 	std::string _name;
 //        std::string _algebricType; ///Scalars/Vectors/Tensors
 //        std::string _numericType; ///Float32/
 //        Matrix _data;
 //        Integer _nComponents;

 //    public:

 //    	VTKDataNode() { }

 //    	VTKDataNode(const std::string &name, const std::string &numericType, const Matrix &data = Matrix(),
 //    		const Integer nComponents = 1)
 //    	: _name(name), _numericType(numericType), _data(data), _nComponents(nComponents) { }


 //    	inline Matrix &data() { return _data; }

 //    	void initialize(const std::string &name, const std::string &numericType, const Matrix &data,
 //    		const Integer nComponents = 1) {
 //    		_name = name;
 //    		_numericType = numericType;
 //    		_data = data;
 //    		_nComponents = nComponents;
 //    	}

 //    	void write(std::ostream &os) const {
 //    		os << "<DataArray type=\"" << _numericType << "\" Name=\"" << _name << "\" NumberOfComponents=\"" <<
 //    		_nComponents << "\" format=\"ascii\">\n";
 //    		if (!_data.isVector()) {
 //    			std::cerr << "Warning: writing matrix in vtk file (check ordering of values)" << std::endl;
 //    		}
 //    		os << _data;
 //    		os << "</DataArray>\n";
 //    	}

 //    	inline bool empty() const { return _data.isZero(); }
 //    };


    template<class Mesh>
    class VTKMeshWriter {
    private:
    	// typedef typename Mesh::Matrix Matrix;
    	static const Integer Dim = Mesh::Dim;
    	static const Integer ManifoldDim = Mesh::ManifoldDim;

    	static const int VTK_TETRA = 10;
    	static const int VTK_TRIANGLE = 5;
    	static const int VTK_QUAD = 9;
    	static const int VTK_HEXAHEDRON = 12;
    	static const int VTK_POLYGON = 7;

    	// std::vector<VTKDataNode<Matrix> > _pointData;
    	// std::vector<VTKDataNode<Matrix> > _cellData;
    	std::string _currentScalarPointData, _currentVectorPointData;

    	// typedef typename std::vector<VTKDataNode<Matrix> >::const_iterator DataIter;


    	// void writePolyData(const Matrix &polyIn, const Integer stride, std::ostream &os) {

    	// 	const std::string type = "Lines";
    	// 	Matrix poly;
    	// 	if(poly.rows() < 3) {
    	// 		poly.resize(3, polyIn.columns());
    	// 		poly.allSet(0);
    	// 		poly.rowRange(0, polyIn.rows()) = polyIn;
    	// 	} else {
    	// 		poly = polyIn;
    	// 	}

    	// 	const Integer nPoints = poly.columns();
    	// 	const Integer nPolys = nPoints / stride;
    	// 	assert(nPoints % stride == 0);

    	// 	os << "<PolyData>\n";
    	// 	os << "<Piece NumberOfPoints=\"" << nPoints << "\" NumberOf" << type << "=\"" << nPolys << "\">\n";
    	// 	writePoints(poly, os);


    	// 	os << "<"<< type << ">\n";
    	// 	os << "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n";

    	// 	for (Integer i = 0; i < nPoints; ++i) {
    	// 		os << i << "\n";
    	// 	}

    	// 	os << "</DataArray>\n";
     //        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    	// 	os << "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n";
    	// 	for (Integer i = 1; i <= nPolys; ++i) {
    	// 		os << i * stride << "\n";
    	// 	}

    	// 	os << "</DataArray>\n";
     //        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    	// 	os << "</" << type << ">\n";
    	// 	os << "</Piece>\n";
    	// 	os << "</PolyData>\n";
    	// }

    	// void writePointData(const Mesh &mesh, std::ostream &os) {
    	// 	if (_currentScalarPointData.empty() && _currentVectorPointData.empty())
    	// 		return;

    	// 	os << "<PointData ";
    	// 	if (!_currentScalarPointData.empty())
    	// 		os << "Scalars=\"" << _currentScalarPointData << "\" ";
    	// 	if (!_currentVectorPointData.empty())
    	// 		os << "Vectors=\"" << _currentVectorPointData << "\" ";
    	// 	os << ">\n";

    	// 	for (DataIter it = _pointData.begin(); it != _pointData.end(); ++it) {
    	// 		it->write(os);
    	// 	}

    	// 	VTKDataNode<express::Matrix<Integer> > ids("id", "UInt64"), gIds("globalId", "UInt64");
    	// 	auxWriteIds(mesh.verticesBegin(), mesh.verticesEnd(), ids.data(), gIds.data());
    	// 	ids.write(os);
    	// 	if (!gIds.empty())
    	// 		gIds.write(os);
    	// 	os << "</PointData>\n";
    	// }

     //    template<class Iterator, class Matrix>
    	// void auxWriteIds(const Iterator &begin, const Iterator &end, Matrix &ids, Matrix &gIds) {
    	// 	Integer nIds = std::distance(begin, end);
    	// 	ids.resize(nIds, 1);
    	// 	gIds.resize(nIds, 1);

    	// 	for (Iterator it = begin; it != end; ++it) {
    	// 		ids[it->handle()] = it->handle();
    	// 		gIds[it->handle()] = it->gId();
    	// 	}

    	// 	if (gIds[0] == -1) {
    	// 		gIds.data().nullify();
    	// 	}
    	// }


    public:
    	void clear() {
    		// _pointData.clear();
    		// _cellData.clear();
    	}

     //    template<class MatrixT>
    	// void addScalarPointData(const std::string &name, const MatrixT &data) {
    	// 	_pointData.push_back(VTKDataNode<Matrix>());
    	// 	_pointData.back().initialize(name, "Float32", data);
    	// 	_currentScalarPointData = name;
    	// }

     //    template<class MatrixT>
    	// void addVectorPointData(const std::string &name, const Integer nComponents, const MatrixT &data) {
    	// 	_pointData.push_back(VTKDataNode<Matrix>());

     //        //FIXME?
    	// 	if (nComponents == 2) {
    	// 		express::BlockMatrix<typename MatrixT::Real> data3((data.rows() * 3) / 2, data.columns(), 3,
    	// 			data.columns());
    	// 		data3.allSet(0);
    	// 		for (Integer i = 0; i < data3.nBlockRows(); ++i) {
    	// 			data3.setBlockAt(i, 0, data.rowRange(i * nComponents, (i + 1) * nComponents));
    	// 		}
    	// 		_pointData.back().initialize(name, "Float32", data3, 3);
    	// 	} else
    	// 	_pointData.back().initialize(name, "Float32", data, nComponents);


    	// 	_currentVectorPointData = name;

    	// }

    	// bool writePolyLine(const std::string &path, const Matrix &poly) {
    	// 	return writePolyLine(path, poly, poly.columns());
    	// }

    	// bool writeMorphedMesh(const std::string &path, const Mesh &mesh, IMorphing<Matrix> &morphing)
    	// {
    	// 	using namespace express;
    	// 	Integer nPointsXEdge = 20;

    	// 	Mesh copied = mesh;
    	// 	copied.createEdges();

    	// 	Matrixd pts(copied.nDims(), copied.nEdges() * nPointsXEdge);
    	// 	CUTK_DEBUG(pts.allSet(-666));

    	// 	Integer index = 0;
    	// 	for(typename Mesh::ConstEdgeIter eIt = copied.edgesBegin(); eIt != copied.edgesEnd(); ++eIt) {

    	// 		const Matrixd p1 = copied.point(eIt->v1());
    	// 		const Matrixd p2 = copied.point(eIt->v2());

    	// 		for(Integer i = 0; i < nPointsXEdge; ++i) {
    	// 			const double t = double(i)/(nPointsXEdge-1);
    	// 			pts.col(index++) =  (1-t)*p1 + t*p2;
    	// 		}
    	// 	}

    	// 	assert(index == pts.columns());

    	// 	Matrix warped;
    	// 	morphing.morph(pts, warped);

    	// 	return writePolyLine(path, warped, nPointsXEdge);
    	// }

    	// bool writePolyLine(const std::string &path, const Matrix &poly, const Integer stride) {

    	// 	std::ofstream os;
    	// 	os.open(path.c_str());
    	// 	if (!os.good()) {
    	// 		os.close();
    	// 		std::cerr << "Unable to open file at " << path << std::endl;
    	// 		return false;
    	// 	}

    	// 	os << "<?xml version=\"1.0\" ?>\n";
    	// 	os << "<VTKFile type=\"PolyData\" version=\"0.1\">\n";
    	// 	writePolyData(poly, stride, os);
    	// 	os << "</VTKFile>\n";
    	// 	os.close();
    	// 	return true;
    	// }


    	bool write(const std::string &path, const Mesh &mesh) {
    		
    		std::ofstream os;
    		os.open(path.c_str());
    		if (!os.good()) {
    			os.close();
    			return false;
    		}

    		writeHeader(mesh, os);
    		writePoints(mesh, os);
    		// writePointData(mesh, os);
    		writeCells(mesh, os);

    		writeFooter(mesh, os);
    		os.close();
    		clear();
    		return true;
    	}

    	void writeHeader(const Mesh &mesh, std::ostream &os) {
    		Integer nCells = mesh.n_active_elements();
    		os << "<VTKFile type=\"UnstructuredGrid\">\n";
    		os << "<UnstructuredGrid>\n";
    		os << "<Piece NumberOfPoints=\"" << mesh.n_nodes() << "\" NumberOfCells=\"" << nCells << "\">\n";

    	}

    	void writeFooter(const Mesh &mesh, std::ostream &os) {
    		os << "</Piece>\n";
    		os << "</UnstructuredGrid>\n";
    		os << "</VTKFile>\n";
    	}


    	inline static int VTKTagVolume(const Integer nVertices) {
    		switch (nVertices) {
    			case 4:
    			return VTK_TETRA;
    			case 8:
    			return VTK_HEXAHEDRON;
    			default:
    			assert(
                            false);//element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
    			return -1;
    		}
    	}

    	inline static int VTKTagPlanar(const Integer nVertices) {

    		if (nVertices > 4) {
    			return VTK_POLYGON;
    		}

    		switch (nVertices) {
    			case 3:
    			return VTK_TRIANGLE;
    			case 4:
    			return VTK_QUAD;
    			default:
    			std::cerr << "[Error] " << nVertices << " not supported" << std::endl;
    			assert(
                            false);//element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
    			return -1;
    		}
    	}

    	void writeCells(const Mesh &mesh, std::ostream &os) {
    		const auto n_active_elements = mesh.n_active_elements();

    		os << "<Cells>\n";
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    		os << "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n";

    		for (Integer k = 0; k < mesh.n_elements(); ++k) {
    			if(!mesh.is_active(k)) continue;

    			for (Integer i = 0; i < ManifoldDim+1; ++i) {
    				const Integer v = mesh.elem(k).nodes[i];
    				os << v;
    				if (i < ManifoldDim) {
    					os << " ";
    				}
    			}
    			os << "\n";
    		}

    		os << "</DataArray>\n";
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    		int minTag, maxTag;
    		if (ManifoldDim == 2) {
    			minTag = VTKTagPlanar(ManifoldDim+1);
    			maxTag = VTKTagPlanar(ManifoldDim+1);
    		} else {
    			minTag = VTKTagVolume(ManifoldDim+1);
    			maxTag = VTKTagVolume(ManifoldDim+1);
    		}

    		

    		os << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"" << minTag <<
    		"\" RangeMax=\"" << maxTag << "\">\n";
    		for (Integer i = 0; i < n_active_elements; ++i) {
    			if (ManifoldDim == 3) {
    				os << VTKTagVolume(ManifoldDim+1) << "\n";
    			} else {
    				os << VTKTagPlanar(ManifoldDim+1) << "\n";

    			}
    		}

    		os << "</DataArray>\n";

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    		
    		os << "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"" << (ManifoldDim+1) <<
    		"\" RangeMax=\"" << (n_active_elements * (ManifoldDim+1)) << "\">\n";

    		for (Integer i = 0, offset = (ManifoldDim+1); i < n_active_elements; ++i, offset += (ManifoldDim+1)) {
    			os << offset << "\n";
    		}

    		os << "</DataArray>\n";
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    		os << "</Cells>\n";

    	}

    	void writePoints(const std::vector<Vector<Real, Dim>> &points, std::ostream &os) {

    		os << "<Points>\n";
    		os << "<DataArray type=\"Float32\" NumberOfComponents=\"" << ((Dim <3) ? 3 : Dim) << "\" format=\"ascii\">\n";
    		for (Integer i = 0; i < points.size(); ++i) {
    			for (Integer d = 0; d < Dim; ++d) {
    				os << points[i](d);
    				if (d < Dim - 1) {
    					os << " ";
    				}else if(Dim == 2){ //padding for paraview vtu format visualisation.
    					os << " ";
    					os << 0;

    				}else if(Dim == 1){
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

    	void writePoints(const Mesh &mesh, std::ostream &os) {
    		writePoints(mesh.points(), os);
    	}

    	virtual ~VTKMeshWriter() { }
    };
}

#endif //GEOM_VTK_WRITER_H
