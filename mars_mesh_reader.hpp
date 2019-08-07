#ifndef MARS_MESH_READER_HPP
#define MARS_MESH_READER_HPP

#include "mars_simplex.hpp"
#include "mars_mesh.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

namespace mars {

    template <Integer Dim, Integer ManifoldDim>
    class MeshReader {
    public:
        virtual ~MeshReader() {}
        virtual bool read(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false) = 0;
    };

    template <Integer Dim, Integer ManifoldDim>
    class MFEMReader final : public MeshReader<Dim, ManifoldDim> {
    public:
        bool read(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false) override
        {
            std::ifstream is(path);
            if (!is.good())
            {
                return false;
            }

            int dim = -1;
            int n_elements = -1;
            int n_nodes = -1;
            int n_coords = -1;

            std::string line;
            while (is.good())
            {
                std::getline(is, line);

                if (line == "dimension")
                {
                    std::getline(is, line);
                    dim = atoi(line.c_str());
                    assert(dim == ManifoldDim);
                }
                else if (line == "elements")
                {
                    std::getline(is, line);
                    n_elements = atoi(line.c_str());

                    for (Integer i = 0; i < n_elements; ++i)
                    {
                        assert(is.good());
                        std::getline(is, line);
                        std::stringstream ss(line);
                        int attr, type;

                        std::array<Integer, ManifoldDim + 1> nodes;
                        ss >> attr >> type;

                        for (Integer k = 0; k < ManifoldDim + 1; ++k)
                        {
                            ss >> nodes[k];
                        }

                        mesh.add_elem(nodes);
                    }
                }
                else if (line == "vertices")
                {
                    std::getline(is, line);
                    n_nodes = atoi(line.c_str());
                    std::getline(is, line);
                    n_coords = atoi(line.c_str());
                    assert(n_coords == Dim);

                    Vector<Real, Dim> p;
                    p.zero();
                    for (Integer i = 0; i < n_nodes; ++i)
                    {
                        assert(is.good());

                        for (Integer k = 0; k < n_coords; ++k)
                        {
                            is >> p(k);
                        }

                        mesh.add_point(p);
                    }
                }
            }

            is.close();

            mesh.repair(verbose);
            return true;
        }
    };

    template <Integer Dim, Integer ManifoldDim>
    class MSHReader final : public MeshReader<Dim, ManifoldDim> {
    public:
        bool read(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false) override
        {
            std::ifstream is(path);
            if (!is.good())
            {
                return false;
            }

            int dim = -1;
            int n_elements = -1;
            int n_nodes = -1;
            int n_coords = -1;

            std::string line;

            while (is.good())
            {
                std::getline(is, line);

                if (line == "$MeshFormat")
                {

                    char version[50], ascii[50], size[50];
                    std::getline(is, line);
                    sscanf((line.c_str()), "%s %s %s", version, ascii, size);
                    if (atof(version) != 4.1)
                    {
                        std::cout << "Only 4.1 version supported \n";
                    }
                    if (atoi(ascii) != 0)
                    {
                        std::cout << "Binary not supported \n";
                    }
                }
                else if (line == "$Nodes")
                {

                    //first line
                    char entity[5], tot_nodes[50], min_tag[50], max_tag[50];
                    std::getline(is, line);
                    sscanf((line.c_str()), "%s %s %s %s", entity, tot_nodes, min_tag, max_tag);
                    n_nodes = atoi(tot_nodes);

                    //second line
                    char dimension[5], something[50], parametric[50], tot_nodes2[50];
                    std::getline(is, line);
                    sscanf((line.c_str()), "%s %s %s %s", dimension, something, parametric, tot_nodes2);
                    n_coords = atoi(dimension);
                    dim = atoi(dimension);
                    assert(dim == ManifoldDim);
                    // assert(n_coords == Dim);

                    //nodes id
                    std::vector<int> id_nodes(n_nodes);
                    // std::cout << n_nodes <<std::endl;
                    for (Integer i = 0; i < n_nodes; i++)
                    {
                        // std::cout << i;
                        // std::cout << Dim <<std::endl;
                        // assert(is.good());
                        // is >> id_nodes(i);

                        std::getline(is, line);
                        id_nodes[i] = atoi(line.c_str());

                        // std::cout << id_nodes << std::endl; // printing the id
                    }

                    //nodes coord
                    Vector<Real, Dim> p;
                    p.zero();
                    Real dump;

                    for (Integer i = 0; i < n_nodes; ++i)
                    {
                        assert(is.good());

                        for (Integer k = 0; k < n_coords; ++k)
                        {
                            is >> p(k);
                        }

                        for (Integer k = n_coords; k < 3; ++k)
                        {
                            is >> dump;
                        }

                        mesh.add_point(p);
                    }
                }
                else if (line == "$Elements")
                {
                    //first line
                    char entity[5], tot_elements[50], min_el_tag[50], max_el_tag[50];
                    std::getline(is, line);
                    sscanf((line.c_str()), "%s %s %s %s", entity, tot_elements, min_el_tag, max_el_tag);
                    n_elements = atoi(tot_elements);

                    //second line
                    char dim[5], something[50], type[50], tot_elements2[50];
                    std::getline(is, line);
                    sscanf((line.c_str()), "%s %s %s %s", dim, something, type, tot_elements2);
                    assert(n_coords == ManifoldDim);

                    //elements
                    for (Integer i = 0; i < n_elements; ++i)
                    {
                        assert(is.good());
                        std::getline(is, line);
                        std::stringstream ss(line);
                        int id_element;

                        std::array<Integer, ManifoldDim + 1> nodes;
                        ss >> id_element;

                        for (Integer k = 0; k < ManifoldDim + 1; ++k)
                        {
                            ss >> nodes[k];
                            nodes[k] -= 1;
                        }

                        mesh.add_elem(nodes);
                    }
                }
            }

            is.close();

            mesh.repair(verbose);
            return true;
        }
    };

    template <Integer Dim, Integer ManifoldDim>
    bool read(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false)
    {
        auto ext = path.substr(path.find_last_of(".") + 1);

        if (ext == "MFEM")
        {
            MFEMReader<Dim, ManifoldDim> r2;
            r2.read(path, mesh, verbose);
        }
        else if (ext == "msh")
        {
            MSHReader<Dim, ManifoldDim> r1;
            r1.read(path, mesh, verbose);
        }
        else
        {
            std::cout << "ERROR: Only MFEM and msh formats are supported" << std::endl;
        }

        return true;
    }

} // namespace mars

#endif //MARS_MESH_READER_HPP
