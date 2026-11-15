#include "public_tet_mesh.hpp"
#include <netcdf.h>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
void check(int status)
{
    if (status != NC_NOERR) throw std::runtime_error(nc_strerror(status));
}

void write_fixture(const char* input_path, const char* path)
{
    const auto input = amr_test::read_public_tet(input_path);
    struct Mesh
    {
        std::vector<std::array<double, 3>> points;
        std::vector<std::array<int, 4>> cells;
    } mesh;
    for (size_t i = 0; i < input.mesh.xyz[0].size(); ++i)
        mesh.points.push_back({input.mesh.xyz[0][i], input.mesh.xyz[1][i], input.mesh.xyz[2][i]});
    for (size_t e = 0; e < input.mesh.conn[0].size(); ++e)
        mesh.cells.push_back({int(input.mesh.conn[0][e]), int(input.mesh.conn[1][e]), int(input.mesh.conn[2][e]),
                              int(input.mesh.conn[3][e])});
    constexpr int opposite[4] = {2, 0, 1, 3};
    std::map<std::array<uint64_t, 3>, std::pair<int, int>> exterior;
    std::map<std::array<uint64_t, 3>, int> incidence;
    for (size_t e = 0; e < mesh.cells.size(); ++e)
        for (int f = 0; f < 4; ++f)
        {
            std::array<uint64_t, 3> key;
            for (int v = 0, j = 0; v < 4; ++v)
                if (v != opposite[f]) key[j++] = mesh.cells[e][v];
            std::sort(key.begin(), key.end());
            ++incidence[key];
            exterior[key] = {int(e) + 1, f + 1};
        }
    const int num_sets = int(input.boundary.size());
    std::vector<std::vector<int>> elements(num_sets), sides(num_sets);
    std::vector<int> set_ids;
    int tag = 0;
    for (const auto& [id, faces] : input.boundary)
    {
        set_ids.push_back(id);
        for (auto face : faces)
        {
            std::sort(face.begin(), face.end());
            if (incidence.at(face) != 1) throw std::runtime_error("Label is not an exterior face");
            elements[tag].push_back(exterior.at(face).first);
            sides[tag].push_back(exterior.at(face).second);
        }
        ++tag;
    }
    int file = -1;
    check(nc_create(path, NC_CLOBBER, &file));
    try
    {
        auto dimension = [&](const std::string& name, size_t size)
        {
            int id;
            check(nc_def_dim(file, name.c_str(), size, &id));
            return id;
        };
        auto variable = [&](const std::string& name, nc_type type, std::initializer_list<int> dims)
        {
            int id;
            check(nc_def_var(file, name.c_str(), type, int(dims.size()), dims.begin(), &id));
            return id;
        };
        const std::string title = "Public MFEM Escher Tet4 mesh for labelled-boundary Poisson";
        check(nc_put_att_text(file, NC_GLOBAL, "title", title.size(), title.data()));
        const float version = 7.22f;
        const int word_size = 8, large_model = 1;
        check(nc_put_att_float(file, NC_GLOBAL, "version", NC_FLOAT, 1, &version));
        check(nc_put_att_float(file, NC_GLOBAL, "api_version", NC_FLOAT, 1, &version));
        check(nc_put_att_int(file, NC_GLOBAL, "floating_point_word_size", NC_INT, 1, &word_size));
        check(nc_put_att_int(file, NC_GLOBAL, "file_size", NC_INT, 1, &large_model));
        dimension("num_dim", 3);
        dimension("num_elem", mesh.cells.size());
        int nodes = dimension("num_nodes", mesh.points.size());
        int block = dimension("num_el_blk", 1), cells = dimension("num_el_in_blk1", mesh.cells.size());
        int corners = dimension("num_nod_per_el1", 4), sets = dimension("num_side_sets", num_sets);
        int width  = dimension("len_name", 33);
        int xyz[3] = {variable("coordx", NC_DOUBLE, {nodes}), variable("coordy", NC_DOUBLE, {nodes}),
                      variable("coordz", NC_DOUBLE, {nodes})};
        int conn   = variable("connect1", NC_INT, {cells, corners});
        check(nc_put_att_text(file, conn, "elem_type", 6, "TETRA4"));
        int names    = variable("ss_names", NC_CHAR, {sets, width});
        int block_id = variable("eb_prop1", NC_INT, {block}), ids = variable("ss_prop1", NC_INT, {sets});
        check(nc_put_att_text(file, block_id, "name", 2, "ID"));
        check(nc_put_att_text(file, ids, "name", 2, "ID"));
        int block_status = variable("eb_status", NC_INT, {block}), set_status = variable("ss_status", NC_INT, {sets});
        std::vector<int> element_ids(num_sets), side_ids(num_sets);
        for (int tag = 0; tag < num_sets; ++tag)
        {
            auto suffix      = std::to_string(tag + 1);
            int count        = dimension("num_side_ss" + suffix, elements[tag].size());
            element_ids[tag] = variable("elem_ss" + suffix, NC_INT, {count});
            side_ids[tag]    = variable("side_ss" + suffix, NC_INT, {count});
        }
        check(nc_enddef(file));
        std::vector<double> coordinate(mesh.points.size());
        for (int d = 0; d < 3; ++d)
        {
            for (size_t i = 0; i < mesh.points.size(); ++i)
                coordinate[i] = mesh.points[i][d];
            check(nc_put_var_double(file, xyz[d], coordinate.data()));
        }
        std::vector<int> connectivity;
        for (const auto& cell : mesh.cells)
            for (int node : cell)
                connectivity.push_back(node + 1);
        check(nc_put_var_int(file, conn, connectivity.data()));
        std::vector<std::array<char, 33>> labels(num_sets);
        for (int tag = 0; tag < num_sets; ++tag)
        {
            auto name = "boundary" + std::to_string(tag);
            std::copy(name.begin(), name.end(), labels[tag].begin());
            check(nc_put_var_int(file, element_ids[tag], elements[tag].data()));
            check(nc_put_var_int(file, side_ids[tag], sides[tag].data()));
        }
        std::vector<int> enabled(num_sets, 1);
        const int block_value = 1;
        check(nc_put_var_text(file, names, &labels[0][0]));
        check(nc_put_var_int(file, block_id, &block_value));
        check(nc_put_var_int(file, ids, set_ids.data()));
        check(nc_put_var_int(file, block_status, enabled.data()));
        check(nc_put_var_int(file, set_status, enabled.data()));
        check(nc_close(file));
        file = -1;
    }
    catch (...)
    {
        if (file >= 0) nc_close(file);
        throw;
    }
}
} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc != 3) throw std::invalid_argument("Supply input MFEM and output Exodus paths");
        write_fixture(argv[1], argv[2]);
        std::cout << "Wrote public Escher Tet4 Exodus fixture\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
