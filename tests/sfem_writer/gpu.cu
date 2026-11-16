#include "backend/distributed/unstructured/utils/mars_sfem_writer.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <type_traits>

void require(bool ok, const char* message)
{
    int bad = ok ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (bad) throw std::runtime_error(message);
}
template<class T>
std::vector<T> host(const cstone::DeviceVector<T>& d)
{
    std::vector<T> h(d.size());
    if (!h.empty()) cudaCheckCall(cudaMemcpy(h.data(), d.data(), h.size() * sizeof(T), cudaMemcpyDeviceToHost));
    return h;
}
template<class T>
void upload(const std::vector<T>& h, cstone::DeviceVector<T>& d)
{
    d.resize(h.size());
    if (!h.empty()) cudaCheckCall(cudaMemcpy(d.data(), h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice));
}
template<class T>
std::vector<T> read(const std::string& path, size_t count)
{
    std::vector<T> h(count);
    std::ifstream f(path, std::ios::binary);
    f.read(reinterpret_cast<char*>(h.data()), count * sizeof(T));
    if (!f || f.peek() != EOF) throw std::runtime_error("Wrong raw output length");
    return h;
}
template<class Tag>
void check(const std::string& mode, const std::string& dir, int rank, int ranks)
{
    using Writer = mars::fem::SfemWriter<uint64_t, double, Tag>;
    using Domain = typename Writer::DomainT;
    typename Domain::HostCoordsTuple xyz;
    typename Domain::HostConnectivityTuple conn;
    auto index = [](int x, int y, int z) { return uint64_t(x + 5 * (y + 5 * z)); };
    for (int z = 0; z < 5; ++z)
        for (int y = 0; y < 5; ++y)
            for (int x = 0; x < 5; ++x)
            {
                std::get<0>(xyz).push_back(.1 + .2 * x);
                std::get<1>(xyz).push_back(.2 + .2 * y);
                std::get<2>(xyz).push_back(.3 + .2 * z);
            }
    const int pattern[5][4] = {{1, 2, 4, 7}, {0, 1, 2, 4}, {3, 1, 2, 7}, {5, 1, 4, 7}, {6, 2, 4, 7}};
    int element             = 0;
    auto append             = [&](const auto& nodes)
    {
        if (element++ % ranks != rank) return;
        int a = 0;
        std::apply([&](auto&... v) { (v.push_back(nodes[a++]), ...); }, conn);
    };
    for (int z = 0; z < 4; ++z)
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
            {
                if constexpr (std::is_same_v<Tag, mars::HexTag>)
                {
                    append(std::array<uint64_t, 8>{index(x, y, z), index(x + 1, y, z), index(x + 1, y + 1, z),
                                                   index(x, y + 1, z), index(x, y, z + 1), index(x + 1, y, z + 1),
                                                   index(x + 1, y + 1, z + 1), index(x, y + 1, z + 1)});
                }
                else
                {
                    for (auto& cell : pattern)
                    {
                        std::array<uint64_t, 4> nodes;
                        for (int a = 0; a < 4; ++a)
                        {
                            int c    = cell[a] ^ ((x + y + z) & 1);
                            nodes[a] = index(x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1));
                        }
                        append(nodes);
                    }
                }
            }
    Domain domain(xyz, conn, rank, ranks, 8, false, 1);
    auto x = host(domain.getNodeX()), y = host(domain.getNodeY()), z = host(domain.getNodeZ());
    auto own = host(domain.getNodeOwnershipMap());
    std::vector<double> field(x.size(), NAN);
    for (size_t i = 0; i < x.size(); ++i)
        if (ranks == 1 || own[i] == 1)
            field[i] = double(float(x[i])) + 2 * double(float(y[i])) + 3 * double(float(z[i]));
    cstone::DeviceVector<double> d_field;
    upload(field, d_field);
    if (mode != "roundtrip")
    {
        if (mode == "shape" && rank == ranks - 1) d_field.resize(0);
        if (mode == "path" && rank == 0)
        {
            std::ofstream file(dir);
            file << "blocks output directory";
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) std::cout << "FAULT_INJECTION_ARMED\n" << std::flush;
        Writer writer(dir, true);
        writer.writeFields(0, 0, domain, {{"p", &d_field, nullptr, nullptr}});
        throw std::runtime_error("Injected writer failure was accepted");
    }
    for (bool fp64 : {false, true})
    {
        const auto path = dir + (fp64 ? "/f64" : "/f32");
        Writer writer(path, fp64);
        for (int frame = 0; frame < 2; ++frame)
            writer.writeFields(frame, frame, domain,
                               {{"p", &d_field, nullptr, nullptr}, {"u", &d_field, &d_field, &d_field}});
        bool valid = true;
        if (rank == 0)
        {
            auto hx = read<float>(path + "/mesh/x.float32", 125);
            auto hy = read<float>(path + "/mesh/y.float32", 125);
            auto hz = read<float>(path + "/mesh/z.float32", 125);
            std::set<std::array<float, 3>> points;
            for (size_t i = 0; i < hx.size(); ++i)
                points.insert({hx[i], hy[i], hz[i]});
            valid &= points.size() == 125;
            for (int frame = 0; frame < 2; ++frame)
            {
                for (const auto* name : {"p", "u.0", "u.1", "u.2"})
                {
                    const auto stem = path + "/out/" + name + ".00000000" + std::to_string(frame);
                    std::vector<double> values;
                    if (fp64)
                        values = read<double>(stem + ".float64", 125);
                    else
                    {
                        auto h = read<float>(stem + ".float32", 125);
                        values.assign(h.begin(), h.end());
                    }
                    for (size_t i = 0; i < values.size(); ++i)
                    {
                        double expected = double(hx[i]) + 2 * double(hy[i]) + 3 * double(hz[i]);
                        if (!fp64) expected = float(expected);
                        valid &= values[i] == expected;
                    }
                }
            }
            const size_t elements = Writer::NPE == 4 ? 320 : 64;
            for (int k = 0; k < Writer::NPE; ++k)
                for (auto id : read<int32_t>(path + "/mesh/i" + std::to_string(k) + ".int32", elements))
                    valid &= id >= 0 && id < 125;
        }
        require(valid, "GPU writer changed owners, numbering or field values");
    }
}
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        require(argc == 4, "Expected mode directory tet|hex");
        MPI_Comm local;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
        int local_rank, devices;
        MPI_Comm_rank(local, &local_rank);
        MPI_Comm_free(&local);
        cudaCheckCall(cudaGetDeviceCount(&devices));
        require(devices > 0, "No GPU available");
        cudaCheckCall(cudaSetDevice(local_rank % devices));
        mars::requestFullNodeHaloExchange();
        if (std::string(argv[3]) == "tet")
            check<mars::TetTag>(argv[1], argv[2], rank, ranks);
        else
            check<mars::HexTag>(argv[1], argv[2], rank, ranks);
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: SFEM GPU/MPI check\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "TEST_FAILURE: " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    MPI_Finalize();
}
