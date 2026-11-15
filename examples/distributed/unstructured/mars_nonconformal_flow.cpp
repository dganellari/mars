#include "mars_nonconformal_flow_runtime.hpp"
#include <string>
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    int result = 1;
    try
    {
        std::string path, output, partitions;
        bool plug      = false;
        int iterations = 800;
        double speed = .1, viscosity = .1, density = 1.;
        auto number = [](const std::string& text)
        {
            size_t used = 0;
            double v    = std::stod(text, &used);
            if (used != text.size() || !std::isfinite(v) || v <= 0)
                throw std::runtime_error("expected a positive finite number");
            return v;
        };
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "--mesh" && i + 1 < argc)
                path = argv[++i];
            else if (arg == "--partitioned-mesh" && i + 1 < argc)
                partitions = argv[++i];
            else if (arg == "--output" && i + 1 < argc)
                output = argv[++i];
            else if (arg == "--check-plug")
                plug = true;
            else if (arg == "--iterations" && i + 1 < argc)
            {
                double value = number(argv[++i]);
                if (value != std::floor(value) || value > 100000) throw std::runtime_error("invalid iteration budget");
                iterations = int(value);
            }
            else if (arg == "--inlet-speed" && i + 1 < argc)
                speed = number(argv[++i]);
            else if (arg == "--viscosity" && i + 1 < argc)
                viscosity = number(argv[++i]);
            else if (arg == "--density" && i + 1 < argc)
                density = number(argv[++i]);
            else
                throw std::runtime_error("unknown or incomplete argument");
        }
        if (path.empty() == partitions.empty() || iterations < 1 || iterations > 100000)
            throw std::runtime_error("expected --mesh prepared.mnc or --partitioned-mesh directory");
        auto mesh = partitions.empty() ? read_mesh(path) : read_mesh_partition(partitions);
        for (const auto& x : mesh.x)
            if (x.x < -1e-12 || x.x > 2 + 1e-12 || x.y < -1e-12 || x.y > 1 + 1e-12 || x.z < -1e-12 || x.z > 1 + 1e-12)
                throw std::runtime_error("this example requires the 2 by 1 by 1 channel");
        ncflow::Channel channel(mesh, plug, speed, viscosity, density);
        bool converged = channel.solve(mesh, iterations, 1e-7);
        if (converged) channel.write(mesh, output);
        if (rankId == 0) std::printf("%s: stationary non-conformal channel iteration\n", converged ? "PASS" : "FAIL");
        result = converged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "rank %d: %s\n", rankId, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return result;
}
