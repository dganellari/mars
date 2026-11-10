#include <mpi.h>
#include <netcdf.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#endif

namespace io_trace
{
size_t connectivity = 0, coordinates = 0, max_run = 0;
int opened = 0, closed = 0;
void check(int result)
{
    if (result != NC_NOERR) throw std::runtime_error(nc_strerror(result));
}
int open(const char* name, int mode, int* file)
{
    const int result = nc_open(name, mode, file);
    if (result == NC_NOERR) ++opened;
    return result;
}
int close(int file) { ++closed; return nc_close(file); }
int full_int(int file, int var, int* out)
{
    char name[NC_MAX_NAME + 1];
    check(nc_inq_varname(file, var, name));
    if (std::string(name).find("connect") == 0) throw std::runtime_error("full connectivity read");
    return nc_get_var_int(file, var, out);
}
int full_double(int, int, double*) { throw std::runtime_error("full coordinate read"); }
int partial_int64(int file, int var, const size_t* start, const size_t* count, long long* out)
{
    connectivity += count[0] * count[1];
    return nc_get_vara_longlong(file, var, start, count, out);
}
int partial_double(int file, int var, const size_t* start, const size_t* count, double* out)
{
    int ndims;
    check(nc_inq_varndims(file, var, &ndims));
    const size_t n = ndims == 1 ? count[0] : count[0] * count[1];
    coordinates += n;
    max_run = std::max(max_run, n);
    return nc_get_vara_double(file, var, start, count, out);
}
void reset() { connectivity = coordinates = max_run = 0; opened = closed = 0; }
}
#define nc_open io_trace::open
#define nc_close io_trace::close
#define nc_get_var_int io_trace::full_int
#define nc_get_var_double io_trace::full_double
#define nc_get_vara_longlong io_trace::partial_int64
#define nc_get_vara_double io_trace::partial_double
#include "mars_read_exodus_mesh.hpp"
#undef nc_open
#undef nc_close
#undef nc_get_var_int
#undef nc_get_var_double
#undef nc_get_vara_longlong
#undef nc_get_vara_double
#include "mars_read_mesh_binary.hpp"

void require(bool ok, const std::string& message)
{
    if (!ok) throw std::runtime_error(message);
}
void nc(int result) { io_trace::check(result); }

struct Fixture
{
    int arity = 8;
    bool packed = false, int64 = false, omit_total = false, sparse = false;
    std::string fault;
    size_t nodes() const { return sparse ? 700001 : 37; }
    size_t node(size_t element, size_t corner) const
    {
        const size_t base = (element * 5 + corner * 3) % 37;
        return sparse ? base * 17011 : base;
    }
    double coord(size_t node, int axis) const { return (double(node) * 3 + axis) / 8; }
};

void write_fixture(const std::string& path, const Fixture& f)
{
    int file, nn, ne, nb, xyz, corner, b1, b3, status, coord[3], conn[2];
    nc(nc_create(path.c_str(), NC_CLOBBER | (f.int64 ? NC_64BIT_DATA : 0), &file));
    nc(nc_def_dim(file, "num_nodes", f.nodes(), &nn));
    if (!f.omit_total) nc(nc_def_dim(file, "num_elem", f.fault == "total" ? 8 : 7, &ne));
    nc(nc_def_dim(file, "num_el_blk", 3, &nb));
    nc(nc_def_dim(file, "num_dim", 3, &xyz));
    nc(nc_def_dim(file, "num_nod_per_el1", f.fault == "arity" ? f.arity + 1 : f.arity, &corner));
    nc(nc_def_dim(file, "num_el_in_blk1", 3, &b1));
    nc(nc_def_dim(file, "num_el_in_blk3", 4, &b3));
    nc(nc_def_var(file, "eb_status", NC_INT, 1, &nb, &status));
    const char* names[] = {"coordx", "coordy", "coordz"};
    if (f.packed)
    {
        int dims[] = {xyz, nn};
        nc(nc_def_var(file, "coord", NC_FLOAT, 2, dims, &coord[0]));
    }
    else
        for (int a = 0; a < 3; ++a)
        {
            int dims[] = {nn, xyz};
            nc(nc_def_var(file, names[a], NC_DOUBLE, f.fault == "coord-rank" ? 2 : 1, dims, &coord[a]));
        }
    for (int b = 0; b < 2; ++b)
    {
        const std::string name = "connect" + std::to_string(b ? 3 : 1);
        int dims[] = {b ? b3 : b1, corner, xyz};
        nc(nc_def_var(file, name.c_str(), f.int64 ? NC_INT64 : NC_INT,
                     f.fault == "conn-rank" ? 3 : 2, dims, &conn[b]));
    }
    nc(nc_enddef(file));
    const int active[] = {1, 0, 1};
    nc(nc_put_var_int(file, status, active));
    if (f.fault == "total" || f.fault == "arity" || f.fault == "conn-rank")
    {
        nc(nc_close(file));
        return;
    }
    std::vector<double> values(f.nodes());
    for (int a = 0; a < 3 && f.fault != "coord-rank"; ++a)
    {
        for (size_t i = 0; i < values.size(); ++i) values[i] = f.coord(i, a);
        if (f.packed)
        {
            const size_t start[] = {size_t(a), 0}, count[] = {1, f.nodes()};
            nc(nc_put_vara_double(file, coord[0], start, count, values.data()));
        }
        else nc(nc_put_var_double(file, coord[a], values.data()));
    }
    for (int b = 0; b < 2; ++b)
    {
        std::vector<long long> values((b ? 4 : 3) * f.arity);
        for (size_t k = 0; k < values.size(); ++k)
            values[k] = 1 + f.node((b ? 3 : 0) + k / f.arity, k % f.arity);
        // Damage one element only: distributed tests must handle rank-local rejection.
        if (b == 1 && f.fault == "node-zero") values.back() = 0;
        if (b == 1 && f.fault == "node-high") values.back() = f.nodes() + 1;
        nc(nc_put_var_longlong(file, conn[b], values.data()));
    }
    nc(nc_close(file));
}

template<class T> void write_binary(const std::string& path, const std::vector<T>& values)
{
    std::ofstream file(path, std::ios::binary);
    file.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(T));
    require(bool(file), "binary fixture write failed");
}
void write_binary_fixture(const std::string& path, const Fixture& f)
{
    std::filesystem::create_directories(path);
    for (int c = 0; c < f.arity; ++c)
    {
        std::vector<int64_t> values(7);
        for (size_t e = 0; e < 7; ++e) values[e] = f.node(e, c);
        write_binary(path + "/i" + std::to_string(c) + ".int64", values);
    }
    const char* names[] = {"x", "y", "z"};
    for (int a = 0; a < 3; ++a)
    {
        std::vector<double> values(f.nodes());
        for (size_t i = 0; i < values.size(); ++i) values[i] = f.coord(i, a);
        write_binary(path + "/" + names[a] + ".double", values);
    }
    write_binary(path + "/block.uint16", std::vector<uint16_t>{0, 0, 0, 2, 2, 2, 2});
}

template<int N, class Real, class Key>
void check_fixture(const std::string& path, const Fixture& f, int rank, int ranks)
{
    io_trace::reset();
    auto mesh = mars::readExodusMeshWithElementPartitioning<N, Real, Key>(path + ".e", rank, ranks);
    const auto& [nodes, elements, x, y, z, connectivity, globals, boundary] = mesh;
    const size_t begin = size_t(rank) * (7 / ranks), end = rank == ranks - 1 ? 7 : begin + 7 / ranks;
    require(elements == end - begin, "wrong element partition");
    std::vector<Key> expected_nodes;
    for (size_t e = begin; e < end; ++e)
        for (int c = 0; c < N; ++c) expected_nodes.push_back(Key(f.node(e, c)));
    std::sort(expected_nodes.begin(), expected_nodes.end());
    expected_nodes.erase(std::unique(expected_nodes.begin(), expected_nodes.end()), expected_nodes.end());
    require(globals == expected_nodes && nodes == globals.size(), "node ID mapping changed");
    for (size_t i = 0; i < nodes; ++i)
        require(x[i] == Real(f.coord(globals[i], 0)) && y[i] == Real(f.coord(globals[i], 1)) &&
                z[i] == Real(f.coord(globals[i], 2)), "coordinate mapping changed");
    [&]<size_t... I>(std::index_sequence<I...>) {
        for (size_t e = 0; e < elements; ++e)
            (require(globals.at(std::get<I>(connectivity)[e]) == f.node(begin + e, I), "connectivity changed"), ...);
    }(std::make_index_sequence<N>{});
    require(boundary == std::vector<uint8_t>(nodes, 0), "boundary default changed");
    require(io_trace::connectivity == elements * N, "read nonlocal connectivity");
    require(io_trace::coordinates >= 3 * nodes && io_trace::coordinates <= 12 * nodes,
            "unbounded coordinate read amplification");
    require(io_trace::max_run <= (size_t(1) << 20), "coordinate scratch exceeds cap");
    require(io_trace::opened == io_trace::closed, "reader leaked a file");
    const auto [blocks, block_count] = mars::readExodusElementBlocks(path + ".e", rank, ranks);
    require(block_count == 3 && blocks.size() == elements, "block metadata changed");
    for (size_t e = 0; e < elements; ++e) require(blocks[e] == (begin + e < 3 ? 0 : 2), "block order changed");
    const auto binary = mars::readMeshWithElementPartitioning<N, Real, Key>(path + "-binary", rank, ranks);
    require(std::get<0>(binary) == nodes && std::get<1>(binary) == elements && std::get<2>(binary) == x &&
            std::get<3>(binary) == y && std::get<4>(binary) == z && std::get<5>(binary) == connectivity &&
            std::get<6>(binary) == globals, "binary/Exodus ElementDomain inputs differ");
    require(mars::readBinaryElementBlocks(path + "-binary", rank, ranks) == std::make_pair(blocks, block_count),
            "binary/Exodus block metadata differs");
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
#if defined(__HIPCC__) || defined(__CUDACC__)
        MPI_Comm shared;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared);
        int local_rank, devices = 0;
        MPI_Comm_rank(shared, &local_rank);
#if defined(__HIPCC__)
        require(hipGetDeviceCount(&devices) == hipSuccess && devices > 0, "no HIP device");
        require(hipSetDevice(local_rank % devices) == hipSuccess, "HIP device selection failed");
#else
        require(cudaGetDeviceCount(&devices) == cudaSuccess && devices > 0, "no CUDA device");
        require(cudaSetDevice(local_rank % devices) == cudaSuccess, "CUDA device selection failed");
#endif
        MPI_Comm_free(&shared);
#endif
        std::filesystem::create_directories(argv[1]);
        for (int variant = 0; variant < 4; ++variant)
        {
            Fixture f;
            f.arity = variant % 2 ? 4 : 8;
            f.packed = f.int64 = variant == 1;
            f.omit_total = variant == 2;
            f.sparse = variant == 3;
            const std::string path = std::string(argv[1]) + "/case-" + std::to_string(variant);
            if (!rank) { write_fixture(path + ".e", f); write_binary_fixture(path + "-binary", f); }
            MPI_Barrier(MPI_COMM_WORLD);
            if (f.arity == 8) check_fixture<8, double, uint64_t>(path, f, rank, ranks);
            else check_fixture<4, float, unsigned>(path, f, rank, ranks);
        }
        for (const char* fault : {"total", "arity", "conn-rank", "coord-rank", "node-zero", "node-high"})
        {
            Fixture f; f.fault = fault;
            const std::string path = std::string(argv[1]) + "/bad-" + fault + ".e";
            if (!rank) write_fixture(path, f);
            MPI_Barrier(MPI_COMM_WORLD);
            io_trace::reset();
            int rejected = 0, total;
            try { mars::readExodusMeshWithElementPartitioning<8>(path, rank, ranks); }
            catch (const std::runtime_error& error)
            {
                rejected = 1;
                const std::string expected = f.fault == "total" ? "element count disagrees" :
                    f.fault == "arity" ? "variable shape" : f.fault.find("rank") != std::string::npos ?
                    "variable rank" : "invalid or unrepresentable node ID";
                require(std::string(error.what()).find(expected) != std::string::npos,
                        "malformed fixture rejected for the wrong reason");
            }
            require(io_trace::opened == io_trace::closed, "failed reader leaked a file");
            MPI_Allreduce(&rejected, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            require(total > 0, std::string("malformed fixture accepted: ") + fault);
            if (f.fault.find("node-") != 0) require(total == ranks, "schema rejected on only some ranks");
        }
        std::vector<uint64_t> duplicates((size_t(1) << 20) + 7);
        for (size_t i = 0; i < duplicates.size(); ++i) duplicates[i] = (i * 17) % 1009;
        const auto unique = mars::dedupNodes(duplicates);
        require(unique.size() == 1009, "large dedup count changed");
        for (size_t i = 0; i < unique.size(); ++i) require(unique[i] == i, "large dedup IDs changed");
        if (!rank) std::cout << "PASS: Exodus/binary parity, bounded reads, null blocks, empty ranks, rejection and dedup\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "rank " << rank << ": " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
