#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>

bool fail_close = false, fail_write = false;
int checked_close(FILE* f)
{
    const int result = std::fclose(f);
    if (fail_close)
    {
        fail_close = false;
        errno      = ENOSPC;
        return EOF;
    }
    return result;
}
size_t checked_write(const void* data, size_t width, size_t count, FILE* f)
{
    if (fail_write)
    {
        fail_write = false;
        errno      = ENOSPC;
        return 0;
    }
    return std::fwrite(data, width, count, f);
}
#define fclose checked_close
#define fwrite checked_write
#include "mars_sfem_writer.hpp"
#undef fclose
#undef fwrite

void require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}
template<class T>
std::vector<T> read(const std::string& path, size_t count)
{
    std::vector<T> values(count);
    std::ifstream file(path, std::ios::binary);
    file.read(reinterpret_cast<char*>(values.data()), count * sizeof(T));
    require(bool(file) && file.peek() == EOF, "wrong output length");
    return values;
}
template<class Tag>
void run(const std::string& mode, const std::string& dir, int rank, int ranks, bool fp64)
{
    using Writer = mars::fem::SfemWriter<uint64_t, double, Tag>;
    typename Writer::DomainT domain(rank, ranks);
    const int nodes = 2 * Writer::NPE;
    std::vector<double> x(nodes), y(nodes), z(nodes);
    for (int i = 0; i < nodes; ++i)
    {
        x[i] = domain.own[i] ? (nodes - i) + 1.0 / 3 : NAN;
        y[i] = 2 * x[i];
        z[i] = -x[i];
    }
    std::vector<typename Writer::Field> fields{{"p", &x, nullptr, nullptr}, {"u", &x, &y, &z}};
    Writer writer(dir, fp64);
    writer.writeFields(0, 0, domain, fields);
    copied_field_bytes = 0;
    if (mode != "roundtrip" && rank == 0) std::cout << "FAULT_INJECTION_ARMED\n" << std::flush;
    if (mode == "shape" && rank == ranks - 1) x.pop_back();
    if (mode == "null" && rank == ranks - 1) fields[0].x = nullptr;
    if (mode == "vector" && rank == ranks - 1) fields[1].z = nullptr;
    if (mode == "schema" && rank == ranks - 1) fields[0].name = "other";
    if (mode == "alias") fields[0].name = "u.0";
    if (mode == "nul-name") fields[0].name = std::string("p\0q", 3);
    if (mode == "path" && rank == 0)
        std::filesystem::create_directory(dir + "/out/p.000000001." + (fp64 ? "float64" : "float32"));
    MPI_Barrier(MPI_COMM_WORLD);
    if (mode == "copy" && rank == 0) fail_next_copy = true;
    if (mode == "close" && rank == 0) fail_close = true;
    if (mode == "write" && rank == 0) fail_write = true;
    if (mode == "connectivity" || mode == "duplicate" || mode == "orphan")
    {
        Writer invalid(dir + "/invalid", fp64);
        if (rank == 0)
        {
            if (mode == "connectivity") std::get<0>(domain.conn)[0] = nodes;
            if (mode == "duplicate") domain.keys[nodes - 1] = domain.keys[nodes - 3];
            if (mode == "orphan") domain.own[nodes - 1] = 0;
        }
        invalid.writeFields(0, 0, domain, fields);
    }
    else
    {
        const double time = mode == "time" && rank == ranks - 1 ? 2 : 1;
        writer.writeFields(1, time, domain, fields);
    }
    require(mode == "roundtrip", "injected failure was accepted");
    const auto owned = std::count(domain.own.begin(), domain.own.end(), uint8_t(1));
    require(copied_field_bytes == size_t(owned) * sizeof(double) * 4, "field transfer included ghosts");
    Writer resumed(dir, fp64);
    resumed.appendTimes();
    resumed.writeFields(2, 2, domain, fields);
    if (rank != 0) return;
    for (int frame = 0; frame < 3; ++frame)
    {
        for (const auto& name : {"p", "u.0", "u.1", "u.2"})
        {
            const std::string path = dir + "/out/" + name + ".00000000" + std::to_string(frame);
            std::vector<double> values;
            if (fp64)
                values = read<double>(path + ".float64", nodes);
            else
            {
                auto f = read<float>(path + ".float32", nodes);
                values.assign(f.begin(), f.end());
            }
            const double scale = std::string(name) == "u.1" ? 2 : std::string(name) == "u.2" ? -1 : 1;
            for (int i = 0; i < nodes; ++i)
            {
                double expected = scale * ((i + 1) + 1.0 / 3);
                if (!fp64) expected = float(expected);
                require(values[i] == expected, "field owner/numbering/precision mismatch");
            }
        }
    }
    for (int k = 0; k < Writer::NPE; ++k)
    {
        const auto ids = read<int32_t>(dir + "/mesh/i" + std::to_string(k) + ".int32", 2);
        require(ids == std::vector<int32_t>{k, Writer::NPE + k}, "connectivity numbering mismatch");
    }
    const char* xyz        = "xyz";
    const double scales[3] = {1. / 3, 1. / 7, -1. / 11};
    for (int d = 0; d < 3; ++d)
    {
        const auto coords = read<float>(dir + "/mesh/" + xyz[d] + ".float32", nodes);
        for (int i = 0; i < nodes; ++i)
            require(std::abs(coords[i] - float(i * scales[d])) < 1e-6, "mesh mismatch");
    }
    std::ifstream times(dir + "/out/time.txt");
    for (double expected : {0., 1., 2.})
    {
        double value = -1;
        times >> value;
        require(value == expected, "time mismatch");
    }
    std::string extra;
    require(!(times >> extra), "extra time records");
}
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        require(argc == 5, "requires mode directory tet|hex f32|f64");
        if (std::string(argv[1]) == "bounds")
        {
            std::vector<int> offsets;
            using mars::fem::detail::sfemGatherDisplacements;
            require(sfemGatherDisplacements({0, INT_MAX, 0}, offsets) == size_t(INT_MAX), "limit rejected");
            require(offsets == std::vector<int>{0, 0, INT_MAX}, "wrong offsets");
            for (auto counts : {std::vector<int>{INT_MAX, 1}, std::vector<int>{-1}})
            {
                bool rejected = false;
                try
                {
                    sfemGatherDisplacements(counts, offsets);
                }
                catch (const std::runtime_error&)
                {
                    rejected = true;
                }
                require(rejected, "invalid count accepted");
            }
        }
        else if (std::string(argv[3]) == "tet")
            run<mars::TetTag>(argv[1], argv[2], rank, ranks, std::string(argv[4]) == "f64");
        else
            run<mars::HexTag>(argv[1], argv[2], rank, ranks, std::string(argv[4]) == "f64");
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) std::cout << "PASS: SFEM host/MPI check\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "TEST_FAILURE: " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    MPI_Finalize();
}
