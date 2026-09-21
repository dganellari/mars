// GPT/Codex, 2026-09-21. Public frozen-input host/CUDA gate; no solver integration.
#include "mars_segregated_tet_interior.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef MARS_REPLAY_CUDA
#include <cuda_runtime.h>
#endif

using mars::segregated::TetInteriorInput;
using mars::segregated::TetInteriorOutput;

template<class T, std::size_t N> void read_array(std::istream& stream, T (&array)[N])
{
    for (auto& value : array)
        if (!(stream >> value) || !std::isfinite(double(value)))
            throw std::runtime_error("invalid/truncated replay data");
}

#ifdef MARS_REPLAY_CUDA
void check_cuda(cudaError_t error)
{
    if (error != cudaSuccess) throw std::runtime_error(cudaGetErrorString(error));
}
__global__ void replay_kernel(const TetInteriorInput* inputs, TetInteriorOutput* outputs, int count)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < count) mars::segregated::tet_interior(inputs[i], outputs[i]);
}
struct DeviceBuffers {
    TetInteriorInput* inputs = nullptr;
    TetInteriorOutput* outputs = nullptr;
    ~DeviceBuffers() { cudaFree(inputs); cudaFree(outputs); }
};
#endif

int main(int argc, char** argv)
{
    try {
        if (argc != 2) throw std::runtime_error("usage: replay PUBLIC_REPLAY_FILE");
        std::ifstream stream(argv[1]);
        std::string magic;
        int count = 0;
        if (!(stream >> magic >> count) || magic != "MARS_PUBLIC_TET_REPLAY_V1"
            || count <= 0 || count > 1000000) throw std::runtime_error("invalid replay header");
        std::vector<TetInteriorInput> inputs(count);
        std::vector<TetInteriorOutput> expected(count), actual(count);
        std::vector<unsigned long long> parents(count);
        std::vector<int> calls(count);
        for (int b = 0; b < count; ++b) {
            auto& x = inputs[b];
            if (!(stream >> x.stage >> calls[b] >> parents[b]) || x.stage < 0 || x.stage > 1
                || calls[b] <= 0 || parents[b] == 0) throw std::runtime_error("invalid block identity");
            read_array(stream, x.edges);
            for (int node : x.edges) if (node < 0 || node >= 4) throw std::runtime_error("bad edge index");
            read_array(stream, x.coordinates); read_array(stream, x.velocity); read_array(stream, x.density);
            read_array(stream, x.velocity_shape); read_array(stream, x.coordinate_shape);
            read_array(stream, x.shape_gradient); read_array(stream, x.area);
            read_array(stream, x.viscosity); read_array(stream, x.velocity_blend);
            read_array(stream, x.velocity_gradient); read_array(stream, x.stored_flux);
            read_array(stream, x.pressure); read_array(stream, x.pressure_gradient);
            read_array(stream, x.influence_lhs); read_array(stream, x.influence_rhs);
            read_array(stream, x.density_blend); read_array(stream, x.density_gradient);
            read_array(stream, expected[b].lhs); read_array(stream, expected[b].rhs);
            read_array(stream, expected[b].flux);
        }
        if (stream >> magic) throw std::runtime_error("trailing replay data");
#ifdef MARS_REPLAY_CUDA
        DeviceBuffers device;
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.inputs), count*sizeof(TetInteriorInput)));
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.outputs), count*sizeof(TetInteriorOutput)));
        check_cuda(cudaMemcpy(device.inputs, inputs.data(), count*sizeof(TetInteriorInput), cudaMemcpyHostToDevice));
        replay_kernel<<<(count+127)/128,128>>>(device.inputs, device.outputs, count);
        check_cuda(cudaGetLastError());
        check_cuda(cudaDeviceSynchronize());
        check_cuda(cudaMemcpy(actual.data(), device.outputs, count*sizeof(TetInteriorOutput), cudaMemcpyDeviceToHost));
        const char* backend = "CUDA";
#else
        for (int b = 0; b < count; ++b) mars::segregated::tet_interior(inputs[b], actual[b]);
        const char* backend = "host";
#endif
        std::size_t checks = 0, failures = 0;
        double worst[2][3] = {};
        int stage_count[2] = {};
        for (int b = 0; b < count; ++b) {
            const int stage = inputs[b].stage, rows = stage ? 12 : 4;
            ++stage_count[stage];
            const double* ref[] = {expected[b].lhs, expected[b].rhs, expected[b].flux};
            const double* got[] = {actual[b].lhs, actual[b].rhs, actual[b].flux};
            const int lengths[] = {rows*rows, rows, 6};
            for (int field = 0; field < 3; ++field) {
                double scale = 1;
                for (int i = 0; i < lengths[field]; ++i) scale = std::max(scale, std::abs(ref[field][i]));
                for (int i = 0; i < lengths[field]; ++i) {
                    ++checks;
                    const double error = std::abs(got[field][i]-ref[field][i])/scale;
                    worst[stage][field] = std::max(worst[stage][field], error);
                    if (!std::isfinite(got[field][i]) || error > 1e-12) {
                        if (failures++ < 8) std::cerr << "mismatch stage=" << stage << " call=" << calls[b]
                            << " parent=" << parents[b] << " field=" << field << " entry=" << i
                            << " expected=" << ref[field][i] << " actual=" << got[field][i] << '\n';
                    }
                }
            }
        }
        for (int stage = 0; stage < 2; ++stage)
            std::cout << (stage ? "momentum" : "pressure") << " blocks=" << stage_count[stage]
                      << " worst_scaled_lhs/rhs/flux=" << worst[stage][0] << '/' << worst[stage][1]
                      << '/' << worst[stage][2] << '\n';
        if (failures) throw std::runtime_error(std::to_string(failures)+" replay mismatches");
        std::cout << "PASS: " << backend << " frozen Tet4 interior replay blocks=" << count
                  << " scalar_checks=" << checks << "; boundary/SIMPLE iteration not tested\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
