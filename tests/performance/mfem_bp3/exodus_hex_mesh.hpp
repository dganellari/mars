#pragma once
#include "mfem.hpp"
#include <netcdf.h>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cctype>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace bp3_mesh
{
inline void nc_check(int status)
{
    if (status != NC_NOERR) throw std::runtime_error(nc_strerror(status));
}
struct File
{
    int id;
    explicit File(const char* name) { nc_check(nc_open(name, NC_NOWRITE, &id)); }
    ~File() { nc_close(id); }
};
inline size_t dimension(int file, const std::string& name)
{
    int id;
    size_t size;
    nc_check(nc_inq_dimid(file, name.c_str(), &id));
    nc_check(nc_inq_dimlen(file, id, &size));
    return size;
}
inline int variable(int file, const std::string& name, const std::vector<size_t>& shape)
{
    int id, nd;
    nc_check(nc_inq_varid(file, name.c_str(), &id));
    nc_check(nc_inq_varndims(file, id, &nd));
    if (nd != int(shape.size())) throw std::runtime_error("Wrong variable rank: " + name);
    std::vector<int> dims(nd);
    nc_check(nc_inq_vardimid(file, id, dims.data()));
    for (int i = 0; i < nd; ++i)
    {
        size_t size;
        nc_check(nc_inq_dimlen(file, dims[i], &size));
        if (size != shape[i]) throw std::runtime_error("Wrong variable size: " + name);
    }
    return id;
}
// Only offline preprocessing holds the full mesh. Preserve file IDs and geometry;
// homogeneous Dirichlet faces are reconstructed from topology, not side-set names.
inline std::unique_ptr<mfem::Mesh> read(const char* name)
{
    File file(name);
    const int f     = file.id;
    const size_t nv = dimension(f, "num_nodes"), ne = dimension(f, "num_elem"), nb = dimension(f, "num_el_blk");
    if (dimension(f, "num_dim") != 3 || !nv || !ne || nv > size_t(INT_MAX) || ne > size_t(INT_MAX))
        throw std::runtime_error("Expected 3D mesh within MFEM integer limits");
    auto mesh = std::make_unique<mfem::Mesh>(3, int(nv), int(ne), 0, 3);
    int coord[3], packed;
    const bool split = nc_inq_varid(f, "coordx", &coord[0]) == NC_NOERR;
    if (split)
    {
        const char* names[] = {"coordx", "coordy", "coordz"};
        for (int d = 0; d < 3; ++d)
            coord[d] = variable(f, names[d], {nv});
    }
    else
        packed = variable(f, "coord", {3, nv});
    constexpr size_t chunk = 8192;
    std::vector<double> xyz[3];
    for (auto& a : xyz)
        a.resize(chunk);
    for (size_t first = 0; first < nv; first += chunk)
    {
        size_t count = std::min(chunk, nv - first);
        for (int d = 0; d < 3; ++d)
        {
            if (split)
                nc_check(nc_get_vara_double(f, coord[d], &first, &count, xyz[d].data()));
            else
            {
                size_t start[] = {size_t(d), first}, sizes[] = {1, count};
                nc_check(nc_get_vara_double(f, packed, start, sizes, xyz[d].data()));
            }
        }
        for (size_t i = 0; i < count; ++i)
        {
            double x[] = {xyz[0][i], xyz[1][i], xyz[2][i]};
            for (double v : x)
                if (!std::isfinite(v)) throw std::runtime_error("Nonfinite coordinate");
            mesh->AddVertex(x);
        }
    }
    std::vector<int> active(nb, 1);
    int status;
    if (nc_inq_varid(f, "eb_status", &status) == NC_NOERR)
    {
        variable(f, "eb_status", {nb});
        nc_check(nc_get_var_int(f, status, active.data()));
    }
    size_t total = 0;
    std::vector<long long> connectivity(chunk * 8);
    for (size_t b = 0; b < nb; ++b)
    {
        if (!active[b]) continue;
        const auto suffix = std::to_string(b + 1);
        size_t cells      = dimension(f, "num_el_in_blk" + suffix);
        if (dimension(f, "num_nod_per_el" + suffix) != 8 || cells > ne - total)
            throw std::runtime_error("Expected Hex8 blocks");
        int id = variable(f, "connect" + suffix, {cells, 8});
        nc_type type;
        nc_check(nc_inq_vartype(f, id, &type));
        if (type != NC_INT && type != NC_INT64) throw std::runtime_error("Connectivity must use integer node IDs");
        size_t length;
        nc_check(nc_inq_attlen(f, id, "elem_type", &length));
        if (!length || length > 32) throw std::runtime_error("Invalid element type");
        std::string kind(length, '\0');
        nc_check(nc_get_att_text(f, id, "elem_type", &kind[0]));
        while (!kind.empty() && (kind.back() == '\0' || kind.back() == ' '))
            kind.pop_back();
        std::transform(kind.begin(), kind.end(), kind.begin(), [](unsigned char c) { return std::toupper(c); });
        if (kind != "HEX" && kind != "HEX8" && kind != "HEXAHEDRON")
            throw std::runtime_error("Only linear Hex8 supported");
        for (size_t first = 0; first < cells; first += chunk)
        {
            size_t count = std::min(chunk, cells - first), start[] = {first, 0}, sizes[] = {count, 8};
            nc_check(nc_get_vara_longlong(f, id, start, sizes, connectivity.data()));
            for (size_t i = 0; i < count; ++i)
            {
                int ids[8];
                for (int j = 0; j < 8; ++j)
                {
                    auto v = connectivity[8 * i + j];
                    if (v < 1 || v > static_cast<long long>(nv)) throw std::runtime_error("Invalid node ID");
                    ids[j] = int(v - 1);
                    for (int k = 0; k < j; ++k)
                        if (ids[k] == ids[j]) throw std::runtime_error("Repeated Hex8 corner");
                }
                mesh->AddHex(ids);
            }
        }
        total += cells;
    }
    if (total != ne) throw std::runtime_error("Element count mismatch");
    mesh->FinalizeHexMesh(1, 0, false);
    if (mesh->CheckElementOrientation(false)) throw std::runtime_error("Inverted input element");
    return mesh;
}
} // namespace bp3_mesh
