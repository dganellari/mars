#include <mpi.h>
#include <netcdf.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#endif

void require(bool ok, const std::string& message)
{
    if (!ok) throw std::runtime_error(message);
}
void nc(int result) { require(result == NC_NOERR, nc_strerror(result)); }
namespace trace
{
size_t largest = 0, conn_values = 0, coord_values = 0, element_rows = 0, side_rows = 0;
int opened = 0, closed = 0;
int open(const char* path, int mode, int* file)
{
    int rc = nc_open(path, mode, file);
    if (!rc) ++opened;
    return rc;
}
int close(int file)
{
    ++closed;
    return nc_close(file);
}
int full_int(int file, int id, int* values)
{
    char name[NC_MAX_NAME + 1];
    nc(nc_inq_varname(file, id, name));
    require(std::string(name) == "eb_status", "unbounded integer read");
    return nc_get_var_int(file, id, values);
}
int full_i64(int, int, long long*) { throw std::runtime_error("unbounded index read"); }
int full_double(int, int, double*) { throw std::runtime_error("unbounded coordinate read"); }
int partial_i64(int file, int id, const size_t* start, const size_t* count, long long* values)
{
    int dims;
    char name[NC_MAX_NAME + 1];
    nc(nc_inq_varndims(file, id, &dims));
    nc(nc_inq_varname(file, id, name));
    size_t n = count[0] * (dims == 2 ? count[1] : 1);
    largest  = std::max(largest, n);
    if (std::string(name).starts_with("connect")) conn_values += n;
    if (std::string(name).starts_with("elem_ss")) element_rows += n;
    if (std::string(name).starts_with("side_ss")) side_rows += n;
    return nc_get_vara_longlong(file, id, start, count, values);
}
int partial_double(int file, int id, const size_t* start, const size_t* count, double* values)
{
    int dims;
    nc(nc_inq_varndims(file, id, &dims));
    size_t n = count[0] * (dims == 2 ? count[1] : 1);
    largest  = std::max(largest, n);
    coord_values += n;
    return nc_get_vara_double(file, id, start, count, values);
}
void reset()
{
    largest = conn_values = coord_values = element_rows = side_rows = 0;
    opened = closed = 0;
}
} // namespace trace
#define nc_open trace::open
#define nc_close trace::close
#define nc_get_var_int trace::full_int
#define nc_get_var_longlong trace::full_i64
#define nc_get_var_double trace::full_double
#define nc_get_vara_longlong trace::partial_i64
#define nc_get_vara_double trace::partial_double
#include "mars_sideset_lookup.hpp"
#ifdef MARS_MESH_IO_ROUTING
#include "mars_sideset_distribution.hpp"
#endif
#undef nc_open
#undef nc_close
#undef nc_get_var_int
#undef nc_get_var_longlong
#undef nc_get_var_double
#undef nc_get_vara_longlong
#undef nc_get_vara_double

struct Fixture
{
    int arity   = 4;
    bool packed = false, legacy_ids = false;
    std::string fault;
    static constexpr size_t faces = 8201, total_nodes = 200001;
    size_t node(size_t e, int c) const { return ((e * 5 + c * 3) % 37) * 1009; }
    std::array<double, 3> coord(uint64_t id) const { return {double(id), double(id) / 8, -double(id) / 4}; }
    size_t element(size_t f) const { return (f * 3 + 1) % 7; }
    int side(size_t f) const { return int(f % (arity == 4 ? 4 : 6)); }
    std::vector<uint64_t> face(size_t f) const
    {
        // Independent Exodus tables, expressed in one-based local numbering.
        const int tet[4][3] = {{1, 2, 4}, {2, 3, 4}, {1, 4, 3}, {1, 3, 2}};
        const int hex[6][4] = {{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 4, 8, 7}, {4, 1, 5, 8}, {1, 4, 3, 2}, {5, 6, 7, 8}};
        std::vector<uint64_t> out;
        for (int c = 0; c < (arity == 4 ? 3 : 4); ++c)
            out.push_back(node(element(f), (arity == 4 ? tet[side(f)][c] : hex[side(f)][c]) - 1));
        return out;
    }
};

void write_fixture(const std::string& path, const Fixture& f)
{
    int file, nn, ne, nb, dim3, status, coord[3], conn[2], sets, len, names, ids, ss_status, face_dim, elem, side;
    nc(nc_create(path.c_str(), NC_CLOBBER | NC_64BIT_DATA, &file));
    nc(nc_def_dim(file, "num_nodes", f.total_nodes, &nn));
    nc(nc_def_dim(file, "num_elem", 7, &ne));
    nc(nc_def_dim(file, "num_el_blk", 3, &nb));
    nc(nc_def_dim(file, "num_dim", 3, &dim3));
    nc(nc_def_var(file, "eb_status", NC_INT, 1, &nb, &status));
    for (int b = 0; b < 2; ++b)
    {
        const std::string suffix = b ? "3" : "1";
        int n, c;
        nc(nc_def_dim(file, ("num_el_in_blk" + suffix).c_str(), b ? 4 : 3, &n));
        nc(nc_def_dim(file, ("num_nod_per_el" + suffix).c_str(), f.arity, &c));
        int dims[] = {n, c};
        nc(nc_def_var(file, ("connect" + suffix).c_str(), NC_INT64, 2, dims, &conn[b]));
    }
    if (f.packed)
    {
        int dims[] = {dim3, nn};
        nc(nc_def_var(file, "coord", NC_DOUBLE, 2, dims, &coord[0]));
    }
    else
        for (int a = 0; a < 3; ++a)
        {
            const char* names[] = {"coordx", "coordy", "coordz"};
            nc(nc_def_var(file, names[a], NC_DOUBLE, 1, &nn, &coord[a]));
        }
    nc(nc_def_dim(file, "num_side_sets", 2, &sets));
    nc(nc_def_dim(file, "len_name", 16, &len));
    int dims[] = {sets, len};
    nc(nc_def_var(file, "ss_names", NC_CHAR, 2, dims, &names));
    nc(nc_def_var(file, "ss_prop1", NC_INT64, 1, &sets, &ids));
    nc(nc_def_var(file, "ss_status", NC_INT, 1, &sets, &ss_status));
    const std::string suffix = f.legacy_ids ? "101" : "1";
    nc(nc_def_dim(file, ("num_side_ss" + suffix).c_str(), f.faces, &face_dim));
    nc(nc_def_var(file, ("elem_ss" + suffix).c_str(), NC_INT64, 1, &face_dim, &elem));
    nc(nc_def_var(file, ("side_ss" + suffix).c_str(), NC_INT64, 1, &face_dim, &side));
    nc(nc_enddef(file));
    const int active[] = {1, 0, 1}, ss_active[] = {1, 0};
    const long long ss_ids[] = {101, 207};
    nc(nc_put_var_int(file, status, active));
    nc(nc_put_var_int(file, ss_status, ss_active));
    nc(nc_put_var_longlong(file, ids, ss_ids));
    std::array<char, 32> labels{};
    std::copy_n("boundary   ", 11, labels.begin());
    std::copy_n("empty", 5, labels.begin() + 16);
    if (f.fault == "duplicate-name") std::copy_n(labels.begin(), 16, labels.begin() + 16);
    nc(nc_put_var_text(file, names, labels.data()));
    std::vector<double> values(f.total_nodes);
    for (int a = 0; a < 3; ++a)
    {
        for (size_t i = 0; i < values.size(); ++i)
            values[i] = f.coord(i)[a];
        if (f.fault == "nonfinite" && a == 0) values[f.face(0)[0]] = std::numeric_limits<double>::infinity();
        if (f.packed)
        {
            size_t start[] = {size_t(a), 0}, count[] = {1, values.size()};
            nc(nc_put_vara_double(file, coord[0], start, count, values.data()));
        }
        else
            nc(nc_put_var_double(file, coord[a], values.data()));
    }
    for (int b = 0; b < 2; ++b)
    {
        std::vector<long long> v((b ? 4 : 3) * f.arity);
        for (size_t i = 0; i < v.size(); ++i)
            v[i] = 1 + f.node((b ? 3 : 0) + i / f.arity, i % f.arity);
        if (f.fault == "node" && b == 0) v[0] = f.total_nodes + 1;
        nc(nc_put_var_longlong(file, conn[b], v.data()));
    }
    std::vector<long long> e(f.faces), s(f.faces);
    for (size_t i = 0; i < e.size(); ++i)
    {
        e[i] = 1 + f.element(i);
        s[i] = 1 + f.side(i);
    }
    if (f.fault == "element") e[5] = 8;
    if (f.fault == "side") s[5] = f.arity == 4 ? 5 : 7;
    nc(nc_put_var_longlong(file, elem, e.data()));
    nc(nc_put_var_longlong(file, side, s.data()));
    nc(nc_close(file));
}

struct Encode
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    __host__ __device__
#endif
        uint64_t
        operator()(const mars::ExodusSideSetNode& node) const
    {
        return uint64_t(node.x);
    }
};

void check(const std::string& path, const Fixture& f, int rank, int ranks)
{
    using Pair = std::pair<uint64_t, int>;
    std::set<Pair> all;
    std::vector<uint64_t> triangles;
    for (size_t i = 0; i < f.faces; ++i)
        for (uint64_t id : f.face(i))
        {
            all.emplace(id, f.element(i) < 3 ? 0 : 2);
            triangles.push_back(id);
        }
    // Deliberately unrelated to the initial element slice, with shared/ghost-like
    // nodes and one rank with no matches. Both coincident block copies are tested.
    std::vector<uint64_t> keys;
    std::vector<int> blocks;
    for (auto [id, block] : all)
        if (ranks == 1 ||
            (rank != ranks - 1 && ((id / 1009 + size_t(block)) % size_t(ranks - 1) == size_t(rank) || id % 3 == 0)))
        {
            keys.push_back(id);
            blocks.push_back(block);
        }
    const uint64_t* key_data = keys.data();
    const int* block_data    = blocks.data();
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<uint64_t> d_keys(keys);
    thrust::device_vector<int> d_blocks(blocks);
    key_data   = thrust::raw_pointer_cast(d_keys.data());
    block_data = thrust::raw_pointer_cast(d_blocks.data());
#endif
    mars::SideSetBatchLookup<uint64_t, Encode> lookup({key_data, block_data, keys.size(), {}});
    // Conforming blocks share one node identity even when the file tags differ.
    std::vector<uint64_t> conformal_keys = keys;
    conformal_keys.erase(std::unique(conformal_keys.begin(), conformal_keys.end()), conformal_keys.end());
    const uint64_t* conformal_data = conformal_keys.data();
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<uint64_t> d_conformal_keys(conformal_keys);
    conformal_data = thrust::raw_pointer_cast(d_conformal_keys.data());
#endif
    mars::SideSetBatchLookup<uint64_t, Encode> conformal({conformal_data, nullptr, conformal_keys.size(), {}});
    std::map<Pair, mars::SideSetHit<uint64_t>> conformal_hits;
    std::set<Pair> seen, matched;
    std::vector<uint64_t> seen_faces;
    size_t batches = 0, requested_nodes = 0;
    trace::reset();
    auto consume = [&](const std::string& name, const mars::ExodusSideSetBatch& batch)
    {
        require(name == "boundary" || name == "empty", "side-set names changed");
        if (batch.nodes.empty()) return;
        ++batches;
        requested_nodes += batch.nodes.size();
        require(batch.nodes.size() <= 4096 * 4 && batch.face_nodes.size() <= 4096 * 4, "unbounded callback batch");
        for (auto& n : batch.nodes)
        {
            require(std::array<double, 3>{n.x, n.y, n.z} == f.coord(n.id), "wrong batch coordinates");
            seen.emplace(n.id, n.block);
        }
        for (auto i : batch.face_nodes)
            seen_faces.push_back(batch.nodes.at(i).id);
        auto hits = lookup.resolve(batch, true);
        require(hits.size() == batch.nodes.size(), "face lookup lost alignment");
        for (size_t i = 0; i < hits.size(); ++i)
        {
            auto pair    = Pair(batch.nodes[i].id, batch.nodes[i].block);
            int expected = -1;
            for (size_t j = 0; j < keys.size(); ++j)
                if (keys[j] == pair.first && blocks[j] == pair.second)
                {
                    expected = int(j);
                    break;
                }
            require(hits[i].local == expected, "block-aware local/ghost lookup changed");
        }
        for (auto hit : conformal.resolve(batch))
        {
            require(hit.local >= 0 && conformal_keys.at(hit.local) == hit.key, "bad conformal hit");
            conformal_hits.emplace(Pair(hit.id, hit.block), hit);
        }
        for (auto h : lookup.resolve(batch))
        {
            require(h.local >= 0 && keys[h.local] == h.key && blocks[h.local] == h.block, "bad compact hit");
            matched.emplace(h.id, h.block);
        }
    };
    if (f.arity == 4)
        mars::visitExodusSideSets<4>(path, consume);
    else
        mars::visitExodusSideSets<8>(path, consume);
    require(seen == all && seen_faces == triangles && batches == 3, "face winding/order/block provenance changed");
    std::set<Pair> wanted;
    for (size_t i = 0; i < keys.size(); ++i)
        wanted.emplace(keys[i], blocks[i]);
    require(matched == wanted, "post-redistribution boundary labels lost");
    require(conformal_hits.size() == conformal_keys.size(), "conforming block tags duplicated a boundary node");
    for (const auto& [identity, hit] : conformal_hits)
        require(identity.second == 0 && identity.first == conformal_keys.at(hit.local),
                "conforming boundary identity must ignore file block tags");
    require(trace::coord_values <= 12 * requested_nodes && trace::conn_values <= 3 * 7 * size_t(f.arity),
            "read amplification exceeds bounds");
    require(trace::opened == trace::closed && trace::largest <= 4096 * 16, "unbounded scratch or leaked file");
#ifdef MARS_MESH_IO_ROUTING
    // Exercise the real sparse router, including a rank with no matching nodes.
    mars::SideSetPrefix<uint64_t> prefix{10};
    std::vector<uint64_t> registered;
    for (auto key : keys)
        registered.push_back(prefix(key));
    registered.erase(std::unique(registered.begin(), registered.end()), registered.end());
    mars::SideSetDistribution distribution(registered, MPI_COMM_WORLD);
    mars::SideSetBatchPrefixes prefixes(Encode{}, prefix);
    std::vector<int> read_count(f.faces, 0), total_read(f.faces, 0);
    std::set<Pair> routed_nodes;
    std::map<uint64_t, std::vector<uint64_t>> routed_faces, expected_faces;
    for (size_t i = 0; i < f.faces; ++i)
    {
        auto face     = f.face(i);
        bool complete = true;
        for (auto id : face)
            complete &= wanted.contains({id, f.element(i) < 3 ? 0 : 2});
        if (complete) expected_faces.emplace(i, face);
    }
    trace::reset();
    auto distributed = [&](const std::string&, const mars::ExodusSideSetBatch& input)
    {
        const size_t q = f.faces / ranks, rem = f.faces % ranks;
        const size_t first = q * rank + std::min(size_t(rank), rem), count = q + (size_t(rank) < rem);
        for (auto id : input.face_ids)
        {
            require(id >= first && id < first + count, "reader consumed another rank's side-set slice");
            ++read_count.at(id);
        }
        auto meta  = prefixes(input);
        auto batch = f.arity == 4 ? distribution.route<3>(input, meta) : distribution.route<4>(input, meta);
        auto hits  = lookup.resolve(batch, true);
        for (const auto& hit : hits)
            if (hit.local >= 0) routed_nodes.emplace(hit.id, hit.block);
        const int corners = f.arity == 4 ? 3 : 4;
        for (size_t i = 0; i < batch.face_ids.size(); ++i)
        {
            bool complete = true;
            std::vector<uint64_t> face;
            for (int c = 0; c < corners; ++c)
            {
                const auto& hit = hits[batch.face_nodes[i * corners + c]];
                complete &= hit.local >= 0;
                face.push_back(hit.id);
            }
            if (complete) routed_faces.emplace(batch.face_ids[i], std::move(face));
        }
    };
    if (f.arity == 4)
        mars::visitExodusSideSets<4>(path, distributed, rank, ranks);
    else
        mars::visitExodusSideSets<8>(path, distributed, rank, ranks);
    MPI_Allreduce(read_count.data(), total_read.data(), int(f.faces), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    require(std::all_of(total_read.begin(), total_read.end(), [](int n) { return n == 1; }),
            "side-set entries were duplicated or lost between readers");
    uint64_t reads[] = {trace::element_rows, trace::side_rows}, sums[2];
    MPI_Allreduce(reads, sums, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    require(sums[0] == f.faces && sums[1] == f.faces, "netCDF side-set list I/O is replicated");
    require(routed_nodes == wanted && routed_faces == expected_faces,
            "distributed side sets lost block identities, complete faces or winding/order");
    require(trace::opened == trace::closed, "distributed reader leaked a file");
    if (!rank)
        std::cout << "DISTRIBUTED_SIDESET_IO corners=" << f.arity << " entries=" << f.faces
                  << " element_reads=" << sums[0] << " side_reads=" << sums[1] << '\n';
#endif
    auto global = f.arity == 4 ? mars::readExodusSideSetsTet4(path, rank) : mars::readExodusSideSetsHex8(path, rank);
    std::set<uint64_t> expected_ids(triangles.begin(), triangles.end());
    require(global.nodesByName.at("boundary") == std::vector<uint64_t>(expected_ids.begin(), expected_ids.end()),
            "global API nodes changed");
    require(global.nodesByName.contains("empty") && global.nodesByName.at("empty").empty(), "empty set lost");
    if (f.arity == 4)
    {
        require(global.trianglesByName.at("boundary") == triangles, "tet triangle winding changed");
        for (size_t i = 0; i < triangles.size(); ++i)
            require(global.triangleCoordsByName.at("boundary")[i] == f.coord(triangles[i]),
                    "triangle coordinates misaligned");
    }
    require(trace::opened == trace::closed, "global API leaked file");
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        require(argc == 2, "expected fixture directory");
#if defined(__CUDACC__) || defined(__HIPCC__)
        MPI_Comm shared;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared);
        int local, devices = 0;
        MPI_Comm_rank(shared, &local);
#if defined(__HIPCC__)
        require(hipGetDeviceCount(&devices) == hipSuccess && devices > 0, "no HIP devices");
        require(hipSetDevice(local % devices) == hipSuccess, "HIP device selection failed");
#else
        require(cudaGetDeviceCount(&devices) == cudaSuccess && devices > 0, "no CUDA devices");
        require(cudaSetDevice(local % devices) == cudaSuccess, "CUDA device selection failed");
#endif
        MPI_Comm_free(&shared);
#endif
        std::filesystem::create_directories(argv[1]);
        const std::string no_sets = std::string(argv[1]) + "/no-sets.e";
        if (!rank)
        {
            int file;
            nc(nc_create(no_sets.c_str(), NC_CLOBBER, &file));
            nc(nc_close(file));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        trace::reset();
        require(mars::readExodusSideSetsTet4(no_sets, rank).nodesByName.empty(), "absent side sets changed");
        require(trace::opened == trace::closed, "absent side sets leaked a file");
        for (int v = 0; v < 4; ++v)
        {
            Fixture f;
            f.arity          = v % 2 ? 8 : 4;
            f.packed         = v >= 2;
            f.legacy_ids     = v == 3;
            std::string path = std::string(argv[1]) + "/sets-" + std::to_string(v) + ".e";
            if (!rank) write_fixture(path, f);
            MPI_Barrier(MPI_COMM_WORLD);
            check(path, f, rank, ranks);
        }
        for (const char* fault : {"element", "side", "node", "nonfinite", "duplicate-name"})
        {
            Fixture f;
            f.fault          = fault;
            std::string path = std::string(argv[1]) + "/bad-" + fault + ".e";
            if (!rank) write_fixture(path, f);
            MPI_Barrier(MPI_COMM_WORLD);
            trace::reset();
            bool rejected = false;
            try
            {
                mars::readExodusSideSetsTet4(path, rank);
            }
            catch (const std::runtime_error& e)
            {
                rejected         = true;
                std::string text = e.what();
                require(text.find(f.fault == "node"             ? "node index"
                                  : f.fault == "nonfinite"      ? "Nonfinite"
                                  : f.fault == "duplicate-name" ? "Duplicate"
                                                                : "element or side index") != std::string::npos,
                        "wrong rejection reason");
            }
            require(rejected && trace::opened == trace::closed, "bad side set accepted or file leaked");
        }
        if (!rank)
        {
#ifdef MARS_MESH_IO_ROUTING
            std::cout << "PASS: partitioned side-set I/O and sparse routing\n";
#endif
            std::cout << "PASS: bounded side sets, face winding, blocks, local/ghost lookup, empty ranks and refusal "
                         "controls\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "rank " << rank << ": " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
