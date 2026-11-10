#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <map>
#include <array>
#include <unordered_set>
#include <utility>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <initializer_list>
#include <cmath>
#include "mars_mesh_io_utils.hpp"
#include "mars_read_exodus_raw.hpp"

#ifdef MARS_HAVE_NETCDF
#include <netcdf.h>
#endif

namespace mars
{

// Defined outside the netCDF guard: code that only names these members (domain.hpp,
// the drivers) must compile in a build without netCDF too.
struct ExodusSideSets
{
    std::map<std::string, std::vector<uint64_t>> nodesByName;   // side-set name -> sorted unique GLOBAL 0-based node IDs
    // side-set name -> per-node coordinate (x,y,z), aligned 1:1 with
    // nodesByName[name]. The driver resolves an Exodus boundary node to its
    // RUNTIME node id by recomputing the SFC key from this coordinate -- because
    // after mesh ingest MARS/cstone index nodes by their SFC-sorted local id,
    // NOT by the Exodus node id. Resolving through the Exodus id (the old g2l)
    // was wrong on >1 rank (the two index spaces differ in order and size).
    std::map<std::string, std::vector<std::array<double, 3>>> nodeCoordsByName;
    // side-set name -> flat list of face triangle coordinates, 3 (x,y,z) per
    // face (size = 3*num_faces), in the same order as trianglesByName. Lets the
    // inlet-normal calc use coordinates directly instead of an index lookup.
    std::map<std::string, std::vector<std::array<double, 3>>> triangleCoordsByName;
    // side-set name -> flat list of face triangles, 3 GLOBAL 0-based node IDs
    // per face (size = 3*num_faces). Used to compute the inlet face normal so
    // an inlet velocity can be applied along the surface normal, not a global
    // axis. Tet faces are triangles; only the tet reader fills this.
    std::map<std::string, std::vector<uint64_t>> trianglesByName;
};

// Batches carry file identity and block provenance until the consumer resolves them.
// Face indices preserve Exodus winding and order, including repeated faces.
struct ExodusSideSetNode
{
    uint64_t id;
    int block;
    double x, y, z;
};

struct ExodusSideSetBatch
{
    std::vector<ExodusSideSetNode> nodes;
    std::vector<unsigned> face_nodes;
    std::vector<uint64_t> face_ids; // Side-set entry indices retain order and repeated faces after routing.
};

#ifdef MARS_HAVE_NETCDF

// Helper macro for netCDF error checking
#define NC_CHECK(call) \
    do { \
        int status = (call); \
        if (status != NC_NOERR) { \
            throw std::runtime_error(std::string("NetCDF error: ") + nc_strerror(status)); \
        } \
    } while(0)

namespace exodus_detail
{
struct File
{
    int id = -1;
    explicit File(const std::string& path) { NC_CHECK(nc_open(path.c_str(), NC_NOWRITE, &id)); }
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File() { if (id >= 0) nc_close(id); }
    void close() { NC_CHECK(nc_close(id)); id = -1; }
};

inline size_t dimension(int file, const std::string& name)
{
    int id;
    size_t size;
    NC_CHECK(nc_inq_dimid(file, name.c_str(), &id));
    NC_CHECK(nc_inq_dimlen(file, id, &size));
    return size;
}

inline int variable(int file, const std::string& name, std::initializer_list<size_t> shape)
{
    int id, ndims;
    NC_CHECK(nc_inq_varid(file, name.c_str(), &id));
    NC_CHECK(nc_inq_varndims(file, id, &ndims));
    if (size_t(ndims) != shape.size())
        throw std::runtime_error("Invalid Exodus variable rank: " + name);
    std::vector<int> dims(ndims);
    NC_CHECK(nc_inq_vardimid(file, id, dims.data()));
    size_t axis = 0;
    for (size_t expected : shape)
    {
        size_t actual;
        NC_CHECK(nc_inq_dimlen(file, dims[axis++], &actual));
        if (actual != expected) throw std::runtime_error("Invalid Exodus variable shape: " + name);
    }
    return id;
}

struct Blocks
{
    std::vector<size_t> sizes;
    size_t total = 0;
};

inline Blocks read_blocks(int file)
{
    Blocks blocks;
    const size_t block_count = dimension(file, "num_el_blk");
    if (block_count > size_t(std::numeric_limits<int>::max()))
        throw std::runtime_error("Exodus block count exceeds local block index capacity");
    blocks.sizes.resize(block_count);
    std::vector<int> active(blocks.sizes.size(), 1);
    int status_id;
    if (nc_inq_varid(file, "eb_status", &status_id) == NC_NOERR)
    {
        status_id = variable(file, "eb_status", {blocks.sizes.size()});
        NC_CHECK(nc_get_var_int(file, status_id, active.data()));
    }
    for (size_t b = 0; b < blocks.sizes.size(); ++b)
    {
        // Null Exodus blocks have no connectivity or element dimensions.
        if (!active[b]) continue;
        const size_t count = dimension(file, "num_el_in_blk" + std::to_string(b + 1));
        if (count > std::numeric_limits<size_t>::max() - blocks.total)
            throw std::runtime_error("Exodus element count overflows size_t");
        blocks.sizes[b] = count;
        blocks.total += count;
    }
    int id;
    if (nc_inq_dimid(file, "num_elem", &id) == NC_NOERR && dimension(file, "num_elem") != blocks.total)
        throw std::runtime_error("Exodus element count disagrees with its blocks");
    return blocks;
}

inline std::pair<size_t, size_t> element_range(size_t total, int rank, int ranks)
{
    if (ranks <= 0 || rank < 0 || rank >= ranks)
        throw std::runtime_error("Invalid Exodus rank partition");
    // Preserve the binary reader's order, including its final-rank remainder.
    const size_t per_rank = total / size_t(ranks);
    return {size_t(rank) * per_rank, rank == ranks - 1 ? total : size_t(rank + 1) * per_rank};
}
} // namespace exodus_detail

// Independent per-rank reads work with serial netCDF as well as parallel builds.
// Only local connectivity and referenced coordinates reach ElementDomain's existing SFC redistribution.
template<int N, typename RealType = double, typename KeyType = unsigned>
inline auto readExodusMeshWithElementPartitioning(const std::string& meshFile, int rank, int numRanks)
{
    static_assert(N == 4 || N == 8, "Exodus reader supports Tet4 and Hex8");
    static_assert(std::is_integral_v<KeyType>);
    using namespace exodus_detail;
    File file(meshFile);
    const size_t total_nodes = dimension(file.id, "num_nodes");
    const auto blocks = read_blocks(file.id);
    const auto [begin, end] = element_range(blocks.total, rank, numRanks);
    const size_t element_count = end - begin;
    if (element_count > std::numeric_limits<size_t>::max() / N)
        throw std::runtime_error("Exodus local connectivity size overflows size_t");

    std::vector<long long> conn(element_count * N);
    size_t block_begin = 0;
    for (size_t b = 0; b < blocks.sizes.size(); ++b)
    {
        const size_t block_end = block_begin + blocks.sizes[b];
        if (blocks.sizes[b])
        {
            // Check every block on every rank, including blocks outside this rank's slice.
            const int id = variable(file.id, "connect" + std::to_string(b + 1), {blocks.sizes[b], N});
            nc_type type;
            NC_CHECK(nc_inq_vartype(file.id, id, &type));
            if (type != NC_INT && type != NC_INT64)
                throw std::runtime_error("Exodus connectivity must contain integer node IDs");
            const size_t lo = std::max(begin, block_begin), hi = std::min(end, block_end);
            if (lo < hi)
            {
                const size_t start[] = {lo - block_begin, 0}, count[] = {hi - lo, N};
                NC_CHECK(nc_get_vara_longlong(file.id, id, start, count, conn.data() + (lo - begin) * N));
            }
        }
        block_begin = block_end;
    }

    std::vector<KeyType> node_ids(conn.size());
    for (size_t k = 0; k < conn.size(); ++k)
    {
        if (conn[k] < 1 || uint64_t(conn[k]) > total_nodes ||
            uint64_t(conn[k] - 1) > uint64_t(std::numeric_limits<KeyType>::max()))
            throw std::runtime_error("Exodus connectivity contains an invalid or unrepresentable node ID");
        node_ids[k] = static_cast<KeyType>(conn[k] - 1);
    }
    auto local_to_global = dedupNodes(node_ids);
    const size_t node_count = local_to_global.size();
    auto connectivity = createNVectors<N, KeyType>(element_count);
    auto local = [&](long long id) {
        return static_cast<KeyType>(std::lower_bound(local_to_global.begin(), local_to_global.end(),
                                                    static_cast<KeyType>(id - 1)) - local_to_global.begin());
    };
    [&]<size_t... I>(std::index_sequence<I...>) {
        for (size_t e = 0; e < element_count; ++e)
            ((std::get<I>(connectivity)[e] = local(conn[e * N + I])), ...);
    }(std::make_index_sequence<N>{});
    std::vector<long long>().swap(conn);
    std::vector<KeyType>().swap(node_ids);

    int coord_ids[3], packed_id;
    const bool split_coords = nc_inq_varid(file.id, "coordx", &coord_ids[0]) == NC_NOERR;
    if (split_coords)
    {
        const char* names[] = {"coordx", "coordy", "coordz"};
        for (int a = 0; a < 3; ++a) coord_ids[a] = variable(file.id, names[a], {total_nodes});
    }
    else packed_id = variable(file.id, "coord", {3, total_nodes});
    std::vector<RealType> x(node_count), y(node_count), z(node_count);
    std::vector<RealType>* coords[] = {&x, &y, &z};
    std::vector<double> buffer;
    // Bound both scratch (8 MiB) and requested values (at most 4x the local nodes).
    // A gap-only rule can read nearly the whole coordinate file on every rank.
    constexpr size_t max_run = size_t(1) << 20;
    buffer.reserve(std::min(max_run / 4, node_count) * 4);
    for (size_t i = 0; i < node_count;)
    {
        const size_t first = size_t(local_to_global[i]);
        size_t j = i + 1;
        while (j < node_count)
        {
            const size_t span = size_t(local_to_global[j]) - first + 1;
            if (span > max_run || (span + 3) / 4 > j - i + 1) break;
            ++j;
        }
        const size_t count = size_t(local_to_global[j - 1]) - first + 1;
        buffer.resize(count);
        for (int a = 0; a < 3; ++a)
        {
            if (split_coords)
                NC_CHECK(nc_get_vara_double(file.id, coord_ids[a], &first, &count, buffer.data()));
            else
            {
                const size_t start[] = {size_t(a), first}, extent[] = {1, count};
                NC_CHECK(nc_get_vara_double(file.id, packed_id, start, extent, buffer.data()));
            }
            for (size_t k = i; k < j; ++k)
                (*coords[a])[k] = static_cast<RealType>(buffer[size_t(local_to_global[k]) - first]);
        }
        i = j;
    }
    file.close();
    if (rank == 0 && std::getenv("MARS_VERBOSE_MESH"))
        std::cout << "Rank 0: " << element_count << " elements, " << node_count << " nodes (from Exodus)\n";
    return std::make_tuple(node_count, element_count, std::move(x), std::move(y), std::move(z),
                           std::move(connectivity), std::move(local_to_global), std::vector<uint8_t>(node_count, 0));
}

// Match the main reader's element order without constructing a global per-element array.
inline std::pair<std::vector<int>, size_t>
readExodusElementBlocks(const std::string& meshFile, int rank, int numRanks)
{
    exodus_detail::File file(meshFile);
    const auto blocks = exodus_detail::read_blocks(file.id);
    const auto [begin, end] = exodus_detail::element_range(blocks.total, rank, numRanks);
    std::vector<int> ids(end - begin);
    size_t block_begin = 0;
    for (size_t b = 0; b < blocks.sizes.size(); ++b)
    {
        const size_t block_end = block_begin + blocks.sizes[b];
        const size_t lo = std::max(begin, block_begin), hi = std::min(end, block_end);
        if (lo < hi) std::fill(ids.begin() + (lo - begin), ids.begin() + (hi - begin), static_cast<int>(b));
        block_begin = block_end;
    }
    file.close();
    return {std::move(ids), std::max(size_t(1), blocks.sizes.size())};
}

// MARS_MESH_DIAG=1 -- SHAREABLE mesh diagnostics, for confidential meshes.
//
// MARS_VERBOSE_MESH prints side-set NAMES and exact face/node counts, which identify the mesh and
// cannot be pasted into a bug report for an NDA geometry. This mode prints only:
//   - topology facts (block count, nodes per element)
//   - side sets by INDEX, never by name
//   - booleans (did it resolve to anything at all)
//   - a hash of the resolved node-id set
// The hash is what makes it useful: identical hashes before and after a change prove the side-set
// resolution is byte-identical, and a 64-bit digest of node IDS reveals nothing about the geometry.
// No coordinates are ever touched here.
inline bool meshDiagEnabled()
{
    static const bool on = std::getenv("MARS_MESH_DIAG") != nullptr;
    return on;
}

// FNV-1a over the raw id bytes. Only needs to be stable within a build, not portable.
inline uint64_t meshDiagHash(const std::vector<uint64_t>& ids)
{
    uint64_t h = 1469598103934665603ull;
    for (uint64_t v : ids)
        for (int b = 0; b < 8; ++b)
        {
            h ^= static_cast<uint64_t>((v >> (b * 8)) & 0xff);
            h *= 1099511628211ull;
        }
    return h;
}

namespace exodus_detail
{
constexpr size_t side_set_batch_faces = 4096;

// Read only requested coordinates. Duplicated IDs on different blocks stay aligned.
inline void read_side_set_coordinates(int file, size_t total_nodes, std::vector<ExodusSideSetNode>& nodes)
{
    int coord_ids[3], packed_id = -1;
    const bool split = nc_inq_varid(file, "coordx", &coord_ids[0]) == NC_NOERR;
    if (split)
    {
        const char* names[] = {"coordx", "coordy", "coordz"};
        for (int a = 0; a < 3; ++a)
            coord_ids[a] = variable(file, names[a], {total_nodes});
    }
    else
        packed_id = variable(file, "coord", {3, total_nodes});
    std::vector<double> buffer;
    for (size_t i = 0; i < nodes.size();)
    {
        const size_t first = nodes[i].id;
        size_t j           = i + 1;
        while (j < nodes.size())
        {
            const size_t span = nodes[j].id - first + 1;
            if (span > (size_t(1) << 20) || (span + 3) / 4 > j - i + 1) break;
            ++j;
        }
        const size_t count = nodes[j - 1].id - first + 1;
        buffer.resize(count);
        for (int a = 0; a < 3; ++a)
        {
            if (split)
                NC_CHECK(nc_get_vara_double(file, coord_ids[a], &first, &count, buffer.data()));
            else
            {
                const size_t start[] = {size_t(a), first}, extent[] = {1, count};
                NC_CHECK(nc_get_vara_double(file, packed_id, start, extent, buffer.data()));
            }
            for (size_t k = i; k < j; ++k)
            {
                const double value = buffer[nodes[k].id - first];
                if (!std::isfinite(value)) throw std::runtime_error("Nonfinite Exodus side-set coordinate");
                if (a == 0)
                    nodes[k].x = value;
                else if (a == 1)
                    nodes[k].y = value;
                else
                    nodes[k].z = value;
            }
        }
        i = j;
    }
}

inline int integer_variable(int file, const std::string& name, std::initializer_list<size_t> shape)
{
    const int id = variable(file, name, shape);
    nc_type type;
    NC_CHECK(nc_inq_vartype(file, id, &type));
    if (type != NC_INT && type != NC_INT64) throw std::runtime_error("Exodus indices must be integer: " + name);
    return id;
}
} // namespace exodus_detail

// No MPI calls. Optional slices divide list entries once across readers.
// Every reader emits the same number of batches, including empty ones, so a
// collective consumer can route each batch without mismatched MPI calls.
// An empty first batch announces each name, including inactive/locally empty sets.
template<int N, class Consumer>
void visitExodusSideSets(const std::string& mesh_file, Consumer&& consume, int rank = 0, int ranks = 1)
{
    static_assert(N == 4 || N == 8);
    if (ranks < 1 || rank < 0 || rank >= ranks) throw std::invalid_argument("Invalid side-set reader partition");
    using namespace exodus_detail;
    constexpr int corners = N == 4 ? 3 : 4, sides = N == 4 ? 4 : 6;
    constexpr int tet_faces[4][4] = {{0, 1, 3, 0}, {1, 2, 3, 0}, {0, 3, 2, 0}, {0, 2, 1, 0}};
    constexpr int hex_faces[6][4] = {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6},
                                     {3, 0, 4, 7}, {0, 3, 2, 1}, {4, 5, 6, 7}};
    const auto* faces             = N == 4 ? tet_faces : hex_faces;
    File file(mesh_file);
    int dim;
    if (nc_inq_dimid(file.id, "num_side_sets", &dim) == NC_EBADDIM) return;
    const size_t sets = dimension(file.id, "num_side_sets");
    if (!sets) return;
    const size_t total_nodes = dimension(file.id, "num_nodes");
    const auto blocks        = read_blocks(file.id);
    std::vector<size_t> ends;
    std::vector<int> connectivity(blocks.sizes.size(), -1);
    size_t offset = 0;
    for (size_t b = 0; b < blocks.sizes.size(); ++b)
    {
        offset += blocks.sizes[b];
        ends.push_back(offset);
        if (blocks.sizes[b])
            connectivity[b] = integer_variable(file.id, "connect" + std::to_string(b + 1), {blocks.sizes[b], N});
    }
    const size_t name_length = dimension(file.id, "len_name");
    if (!name_length) throw std::runtime_error("Invalid Exodus side-set name length");
    const int names_id       = variable(file.id, "ss_names", {sets, name_length});
    const int ids_id         = integer_variable(file.id, "ss_prop1", {sets});
    int status_id            = -1;
    if (nc_inq_varid(file.id, "ss_status", &status_id) == NC_NOERR)
        status_id = integer_variable(file.id, "ss_status", {sets});
    std::vector<char> raw_name(name_length);
    std::unordered_set<std::string> names;
    for (size_t set = 0; set < sets; ++set)
    {
        const size_t name_start[] = {set, 0}, name_count[] = {1, name_length}, one = 1;
        NC_CHECK(nc_get_vara_text(file.id, names_id, name_start, name_count, raw_name.data()));
        size_t length = strnlen(raw_name.data(), raw_name.size());
        while (length && std::strchr(" \t\r\n", raw_name[length - 1]))
            --length;
        const std::string name(raw_name.data(), length);
        if (!names.insert(name).second) throw std::runtime_error("Duplicate Exodus side-set name");
        long long set_id, active = 1;
        NC_CHECK(nc_get_vara_longlong(file.id, ids_id, &set, &one, &set_id));
        if (status_id >= 0) NC_CHECK(nc_get_vara_longlong(file.id, status_id, &set, &one, &active));
        consume(name, ExodusSideSetBatch{});
        if (!active) continue;
        std::string suffix = std::to_string(set + 1);
        if (nc_inq_dimid(file.id, ("num_side_ss" + suffix).c_str(), &dim) != NC_NOERR) suffix = std::to_string(set_id);
        const size_t face_count = dimension(file.id, "num_side_ss" + suffix);
        const int elem_id       = integer_variable(file.id, "elem_ss" + suffix, {face_count});
        const int side_id       = integer_variable(file.id, "side_ss" + suffix, {face_count});
        const size_t quotient = face_count / size_t(ranks), remainder = face_count % size_t(ranks);
        const size_t begin = quotient * size_t(rank) + std::min(size_t(rank), remainder);
        const size_t local_count = quotient + (size_t(rank) < remainder);
        const size_t largest_count = quotient + (remainder != 0);
        for (size_t offset = 0; offset < largest_count; offset += side_set_batch_faces)
        {
            const size_t first = begin + std::min(offset, local_count);
            const size_t count = std::min(side_set_batch_faces, local_count - std::min(offset, local_count));
            if (!count)
            {
                consume(name, ExodusSideSetBatch{});
                continue;
            }
            std::vector<long long> elements(count), side(count);
            NC_CHECK(nc_get_vara_longlong(file.id, elem_id, &first, &count, elements.data()));
            NC_CHECK(nc_get_vara_longlong(file.id, side_id, &first, &count, side.data()));
            std::vector<size_t> element_ids;
            for (size_t f = 0; f < count; ++f)
            {
                if (elements[f] < 1 || uint64_t(elements[f]) > blocks.total || side[f] < 1 || side[f] > sides)
                    throw std::runtime_error("Invalid Exodus side-set element or side index");
                element_ids.push_back(size_t(elements[f] - 1));
            }
            std::sort(element_ids.begin(), element_ids.end());
            element_ids.erase(std::unique(element_ids.begin(), element_ids.end()), element_ids.end());
            std::vector<long long> conn(element_ids.size() * N);
            std::vector<int> element_blocks(element_ids.size());
            for (size_t i = 0; i < element_ids.size();)
            {
                const size_t b = std::upper_bound(ends.begin(), ends.end(), element_ids[i]) - ends.begin();
                size_t j       = i + 1;
                while (j < element_ids.size() && element_ids[j] == element_ids[j - 1] + 1 && element_ids[j] < ends[b])
                    ++j;
                const size_t start[] = {element_ids[i] - (b ? ends[b - 1] : 0), 0}, extent[] = {j - i, N};
                NC_CHECK(nc_get_vara_longlong(file.id, connectivity[b], start, extent, conn.data() + i * N));
                std::fill(element_blocks.begin() + i, element_blocks.begin() + j, int(b));
                i = j;
            }
            using Identity = std::pair<uint64_t, int>;
            std::vector<Identity> face_nodes;
            face_nodes.reserve(count * corners);
            for (size_t f = 0; f < count; ++f)
            {
                const size_t e = std::lower_bound(element_ids.begin(), element_ids.end(), size_t(elements[f] - 1)) -
                                 element_ids.begin();
                for (int c = 0; c < corners; ++c)
                {
                    const long long node = conn[e * N + faces[side[f] - 1][c]];
                    if (node < 1 || uint64_t(node) > total_nodes)
                        throw std::runtime_error("Invalid Exodus side-set node index");
                    face_nodes.emplace_back(uint64_t(node - 1), element_blocks[e]);
                }
            }
            auto nodes = face_nodes;
            std::sort(nodes.begin(), nodes.end());
            nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            ExodusSideSetBatch batch;
            for (size_t f = 0; f < count; ++f) batch.face_ids.push_back(first + f);
            batch.nodes.reserve(nodes.size());
            for (auto [id, block] : nodes)
                batch.nodes.push_back({id, block, 0, 0, 0});
            read_side_set_coordinates(file.id, total_nodes, batch.nodes);
            batch.face_nodes.reserve(face_nodes.size());
            for (auto node : face_nodes)
                batch.face_nodes.push_back(
                    unsigned(std::lower_bound(nodes.begin(), nodes.end(), node) - nodes.begin()));
            consume(name, batch);
        }
    }
    file.close();
}

// Explicit global-return APIs keep their contract for callers that need complete
// boundaries. Only their requested boundary output grows with the global mesh.
template<int N>
ExodusSideSets readExodusSideSets(const std::string& mesh_file, int rank)
{
    ExodusSideSets out;
    std::map<std::string, std::map<uint64_t, std::array<double, 3>>> nodes;
    visitExodusSideSets<N>(mesh_file,
                           [&](const std::string& name, const ExodusSideSetBatch& batch)
                           {
                               auto& set = nodes[name];
                               if constexpr (N == 4)
                               {
                                   auto& triangles = out.trianglesByName[name];
                                   auto& coords    = out.triangleCoordsByName[name];
                                   for (unsigned i : batch.face_nodes)
                                   {
                                       const auto& node = batch.nodes[i];
                                       triangles.push_back(node.id);
                                       coords.push_back({node.x, node.y, node.z});
                                   }
                               }
                               for (const auto& node : batch.nodes)
                                   set.emplace(node.id, std::array<double, 3>{node.x, node.y, node.z});
                           });
    if (rank == 0 && meshDiagEnabled())
        std::cout << "[mesh-diag] side_sets=" << nodes.size() << "\n";
    size_t set_index = 0;
    for (const auto& [name, set] : nodes)
    {
        auto& ids    = out.nodesByName[name];
        auto& coords = out.nodeCoordsByName[name];
        for (const auto& [id, coord] : set)
        {
            ids.push_back(id);
            coords.push_back(coord);
        }
        if (rank == 0 && std::getenv("MARS_VERBOSE_MESH"))
            std::cout << "Exodus side-set [" << name << "]: " << ids.size() << " unique nodes\n";
        if (rank == 0 && meshDiagEnabled())
            std::cout << "[mesh-diag]   sideset[" << set_index << "]: nodes>0=" << (ids.empty() ? "no" : "yes")
                      << " nodes_hash=0x" << std::hex << meshDiagHash(ids) << std::dec << "\n";
        ++set_index;
    }
    return out;
}

inline ExodusSideSets readExodusSideSetsHex8(const std::string& mesh_file, int rank)
{
    return readExodusSideSets<8>(mesh_file, rank);
}

inline ExodusSideSets readExodusSideSetsTet4(const std::string& mesh_file, int rank)
{
    return readExodusSideSets<4>(mesh_file, rank);
}

#undef NC_CHECK

#else // !MARS_HAVE_NETCDF

template<int N, typename RealType = double, typename KeyType = unsigned>
inline auto readExodusMeshWithElementPartitioning(const std::string&, int, int)
{
    static_assert(N == 4 || N == 8, "Exodus reader supports Tet4 and Hex8");
    throw std::runtime_error("Exodus mesh reading requires netCDF. Please rebuild with netCDF support.");
    // Keep the tuple type available to builds that cannot open Exodus files.
    return std::make_tuple(size_t(0), size_t(0), std::vector<RealType>(), std::vector<RealType>(),
                           std::vector<RealType>(), createNVectors<N, KeyType>(0), std::vector<KeyType>(),
                           std::vector<uint8_t>());
}

template<int N, class Consumer>
void visitExodusSideSets(const std::string&, Consumer&&, int = 0, int = 1)
{
    throw std::runtime_error("Exodus side-set reading requires netCDF.");
}

// Stub: without netCDF a .exo is unreadable (the main reader above throws first), so return the
// single-block default. Signature matches the real implementation.
inline std::pair<std::vector<int>, size_t> readExodusElementBlocks(const std::string&, int, int)
{
    return {std::vector<int>(), size_t(1)};
}

inline ExodusSideSets readExodusSideSetsHex8(const std::string&, int)
{
    throw std::runtime_error("readExodusSideSetsHex8 requires netCDF.");
}

inline ExodusSideSets readExodusSideSetsTet4(const std::string&, int)
{
    throw std::runtime_error("readExodusSideSetsTet4 requires netCDF.");
}

#endif // MARS_HAVE_NETCDF

} // namespace mars
