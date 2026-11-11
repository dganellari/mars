#include "tet_exodus_fixture.hpp"
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

void write_fixture(const char* path)
{
    const auto mesh = amr_test::exodus_fixture();
    // Exodus Tet4 sides 1..4 are opposite local vertices 3, 1, 2, 4.
    constexpr int opposite[4] = {2, 0, 1, 3};
    std::array<std::vector<int>, 7> elements, sides;
    for (size_t e = 0; e < mesh.cells.size(); ++e)
        for (int face = 0; face < 4; ++face)
            for (int tag = 0; tag < 6; ++tag)
            {
                bool boundary = true;
                for (int c = 0; c < 4; ++c)
                    if (c != opposite[face]) boundary &= mesh.points[mesh.cells[e][c]][tag / 2] == tag % 2;
                if (!boundary) continue;
                elements[tag].push_back(int(e) + 1);
                sides[tag].push_back(face + 1);
            }
    elements[6] = elements[0];
    sides[6]    = sides[0];

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
        const std::string title = "Public synthetic tetrahedral AMR cube with shuffled file numbering";
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
        int corners = dimension("num_nod_per_el1", 4), sets = dimension("num_side_sets", 7);
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
        std::array<int, 7> element_ids, side_ids;
        for (int tag = 0; tag < 7; ++tag)
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
        char labels[7][33]{};
        for (int tag = 0; tag < 7; ++tag)
        {
            auto name = "boundary" + std::to_string(tag);
            std::copy(name.begin(), name.end(), labels[tag]);
            check(nc_put_var_int(file, element_ids[tag], elements[tag].data()));
            check(nc_put_var_int(file, side_ids[tag], sides[tag].data()));
        }
        const int set_ids[7] = {101, 307, 509, 701, 907, 1103, 1301};
        const int enabled[7] = {1, 1, 1, 1, 1, 1, 1}, block_value = 41;
        check(nc_put_var_text(file, names, &labels[0][0]));
        check(nc_put_var_int(file, block_id, &block_value));
        check(nc_put_var_int(file, ids, set_ids));
        check(nc_put_var_int(file, block_status, enabled));
        check(nc_put_var_int(file, set_status, enabled));
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
        if (argc != 2) throw std::invalid_argument("Supply an output Exodus path");
        write_fixture(argv[1]);
        std::cout << "Wrote public Tet4 Exodus fixture\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
