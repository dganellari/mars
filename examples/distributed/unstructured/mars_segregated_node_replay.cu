#include "mars_segregated_node.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef MARS_REPLAY_CUDA
#include <cuda_runtime.h>
#define NODE_HD __host__ __device__
#else
#define NODE_HD
#endif

struct Input {
    int stage = 0, blocks = 0, diagonal = 0, consistent = 0;
    std::size_t offset = 0;
    mars::segregated::SteadyMomentumNode node{};
    double factor = 1, volume = 0, lhs[9]{}, rhs[3]{};
};
struct Output { double values[18]{}; bool valid = true; };

NODE_HD void evaluate(const Input& in, const double* rows, Output& out)
{
    using namespace mars::segregated;
    out = Output{};
    if (in.stage == 0) steady_momentum_node(in.node, out.values, out.values+9);
    else if (in.stage == 1) {
        for (int i = 0; i < 9; ++i) out.values[i] = in.lhs[i];
        relax_momentum_diagonal(out.values, in.factor);
    } else if (in.stage == 2) {
        out.valid = momentum_influence(in.volume, rows+in.offset, in.blocks, in.diagonal,
                                       in.consistent != 0, out.values+12, out.values+15);
    } else {
        for (int i = 0; i < 3; ++i) out.values[9+i] = in.rhs[i];
        relax_boundary_rhs(out.values+9, in.factor);
    }
}

template<class T> void read_value(std::istream& stream, T& value)
{
    if (!(stream >> value) || !std::isfinite(double(value)))
        throw std::runtime_error("invalid/truncated node replay data");
}
template<std::size_t N> void read_array(std::istream& stream, double (&values)[N])
{
    for (double& value : values) read_value(stream, value);
}
void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
#ifdef MARS_REPLAY_CUDA
void check_cuda(cudaError_t error)
{
    if (error != cudaSuccess) throw std::runtime_error(cudaGetErrorString(error));
}
__global__ void replay_kernel(const Input* inputs, const double* rows, Output* outputs, int count)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < count) evaluate(inputs[i], rows, outputs[i]);
}
struct DeviceBuffers {
    Input* inputs = nullptr;
    Output* outputs = nullptr;
    double* rows = nullptr;
    ~DeviceBuffers() { cudaFree(inputs); cudaFree(outputs); cudaFree(rows); }
};
#endif

int main(int argc, char** argv)
{
    try {
        require(argc == 2, "usage: mars_segregated_node_replay PUBLIC_NODE_REPLAY_FILE");
        std::ifstream stream(argv[1]);
        std::string magic;
        int count = 0;
        require(bool(stream >> magic >> count) && magic == "MARS_PUBLIC_NODE_REPLAY_V1"
                && count > 0 && count <= 1000000, "invalid node replay header");
        std::vector<Input> inputs(count);
        std::vector<Output> expected(count), actual(count);
        std::vector<unsigned long long> nodes(count);
        std::vector<int> calls(count);
        std::vector<double> rows;
        int stage_count[4]{};
        for (int b = 0; b < count; ++b) {
            auto& in = inputs[b];
            read_value(stream, in.stage); read_value(stream, calls[b]); read_value(stream, nodes[b]);
            require(in.stage >= 0 && in.stage < 4 && calls[b] > 0 && nodes[b] > 0, "invalid node identity");
            ++stage_count[in.stage];
            if (in.stage == 0) {
                auto& n = in.node;
                read_value(stream, n.density); read_value(stream, n.volume);
                read_value(stream, n.pseudo_dt); read_value(stream, n.mass_divergence);
                read_array(stream, n.velocity); read_array(stream, n.pressure_gradient);
                read_array(stream, n.force); read_array(stream, n.source);
                require(n.density > 0 && n.volume > 0 && n.pseudo_dt > 0, "invalid node scales");
            } else if (in.stage == 1) {
                read_value(stream, in.factor); read_array(stream, in.lhs);
            } else if (in.stage == 2) {
                read_value(stream, in.blocks); read_value(stream, in.diagonal);
                read_value(stream, in.consistent); read_value(stream, in.volume);
                require(in.blocks > 0 && in.blocks <= 100000 && in.diagonal >= 0 && in.diagonal < in.blocks
                        && (in.consistent == 0 || in.consistent == 1) && in.volume > 0, "invalid block row");
                in.offset = rows.size();
                require(rows.size()+9*std::size_t(in.blocks) <= 100000000, "node replay too large");
                rows.resize(rows.size()+9*std::size_t(in.blocks));
                for (std::size_t i = in.offset; i < rows.size(); ++i) read_value(stream, rows[i]);
            } else {
                read_value(stream, in.factor); read_array(stream, in.rhs);
            }
            require(in.factor > 0 && in.factor <= 1, "invalid relaxation factor");
            read_array(stream, expected[b].values);
        }
        for (int n : stage_count) require(n > 0, "missing replay stage");
        require(!(stream >> magic), "trailing replay data");
#ifdef MARS_REPLAY_CUDA
        DeviceBuffers device;
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.inputs), inputs.size()*sizeof(Input)));
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.outputs), actual.size()*sizeof(Output)));
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.rows), rows.size()*sizeof(double)));
        check_cuda(cudaMemcpy(device.inputs, inputs.data(), inputs.size()*sizeof(Input), cudaMemcpyHostToDevice));
        check_cuda(cudaMemcpy(device.rows, rows.data(), rows.size()*sizeof(double), cudaMemcpyHostToDevice));
        replay_kernel<<<(count+127)/128,128>>>(device.inputs, device.rows, device.outputs, count);
        check_cuda(cudaGetLastError());
        check_cuda(cudaDeviceSynchronize());
        check_cuda(cudaMemcpy(actual.data(), device.outputs, actual.size()*sizeof(Output), cudaMemcpyDeviceToHost));
        const char* backend = "CUDA";
#else
        for (int b = 0; b < count; ++b) evaluate(inputs[b], rows.data(), actual[b]);
        const char* backend = "host";
#endif
        std::size_t failures = 0, checks = 0;
        double worst[4]{};
        const int offsets[] = {0, 9, 12, 15, 18};
        for (int b = 0; b < count; ++b) {
            require(actual[b].valid, "invalid momentum influence denominator or coefficient");
            for (int field = 0; field < 4; ++field) {
                double scale = 1;
                for (int i = offsets[field]; i < offsets[field+1]; ++i)
                    scale = std::max(scale, std::abs(expected[b].values[i]));
                for (int i = offsets[field]; i < offsets[field+1]; ++i) {
                    ++checks;
                    const double value = actual[b].values[i];
                    const double error = std::abs(value-expected[b].values[i])/scale;
                    worst[inputs[b].stage] = std::max(worst[inputs[b].stage], error);
                    if (!std::isfinite(value) || error > 1e-12) {
                        if (failures++ < 8) std::cerr << "mismatch stage=" << inputs[b].stage
                            << " call=" << calls[b] << " node=" << nodes[b] << " entry=" << i
                            << " expected=" << expected[b].values[i] << " actual=" << value << '\n';
                    }
                }
            }
        }
        for (int s = 0; s < 4; ++s)
            std::cout << "stage=" << s << " nodes=" << stage_count[s] << " worst_scaled=" << worst[s] << '\n';
        require(failures == 0, "node replay mismatch");
        std::cout << "PASS: " << backend << " frozen node replay records=" << count << " scalar_checks=" << checks
                  << "; boundary selection/global assembly/SIMPLE iteration not tested\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
