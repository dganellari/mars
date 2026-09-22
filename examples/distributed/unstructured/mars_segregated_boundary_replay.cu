#include "mars_segregated_boundary.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef MARS_REPLAY_CUDA
#include <cuda_runtime.h>
#endif

using mars::segregated::BoundaryInput;
using mars::segregated::BoundaryOutput;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
template<class T> void read_value(std::istream& stream, T& value)
{
    require(bool(stream >> value) && std::isfinite(double(value)), "invalid/truncated boundary data");
}
template<class T, std::size_t N> void read_array(std::istream& stream, T (&values)[N])
{
    for (auto& value : values) read_value(stream, value);
}
void validate(const BoundaryInput& x)
{
    const int count = x.stage == 5 ? 3 : 4;
    bool face[4]{}, nearest[4]{};
    for (int i = 0; i < 3; ++i) {
        require(x.face_nodes[i] >= 0 && x.face_nodes[i] < count, "invalid face node");
        require(x.nearest[i] >= 0 && x.nearest[i] < count, "invalid nearest node");
        require(!face[x.face_nodes[i]] && !nearest[x.nearest[i]], "duplicate boundary map");
        face[x.face_nodes[i]] = nearest[x.nearest[i]] = true;
        require(x.reversal[i] == 0 || x.reversal[i] == 1, "invalid reversal flag");
        double magnitude = 0, weight = 0;
        for (int j = 0; j < 3; ++j) {
            magnitude += x.area[3*i+j]*x.area[3*i+j];
            require(x.shape[3*i+j] >= 0, "negative interpolation weight");
            weight += x.shape[3*i+j];
        }
        require(magnitude > 0 && std::isfinite(magnitude), "invalid sample area");
        require(std::abs(weight-1) < 1e-12, "invalid interpolation weights");
        if (x.stage < 3) require(x.density[i] > 0, "invalid density");
        if (x.stage == 3 || x.stage == 4) require(x.viscosity[i] >= 0, "invalid viscosity");
        if (x.stage == 5) {
            require(x.face_nodes[i] == i, "wall block must use native face ordering");
            require(x.wall_coefficient[i] >= 0, "invalid wall coefficient");
        }
    }
    for (int node = 0; node < count; ++node) {
        require(face[node] == nearest[node], "nearest node not on face");
        if (x.stage == 1 || x.stage == 3)
            require(x.bc_multiplier[node] == (face[node] ? 0 : 1), "incorrect boundary column mask");
    }
    if (x.stage != 5)
        for (int opposite : x.opposing)
            require(opposite >= 0 && opposite < 4 && !face[opposite], "invalid opposite node");
}

#ifdef MARS_REPLAY_CUDA
void check_cuda(cudaError_t error)
{
    if (error != cudaSuccess) throw std::runtime_error(cudaGetErrorString(error));
}
__global__ void replay_kernel(const BoundaryInput* inputs, BoundaryOutput* outputs, int count)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < count) mars::segregated::boundary_block(inputs[i], outputs[i]);
}
struct DeviceBuffers {
    BoundaryInput* inputs = nullptr;
    BoundaryOutput* outputs = nullptr;
    ~DeviceBuffers() { cudaFree(inputs); cudaFree(outputs); }
};
#endif

int main(int argc, char** argv)
{
    try {
        require(argc == 2, "usage: mars_segregated_boundary_replay PUBLIC_BOUNDARY_REPLAY_FILE");
        std::ifstream stream(argv[1]);
        std::string magic;
        int count = 0;
        require(bool(stream >> magic >> count) && magic == "MARS_PUBLIC_BOUNDARY_REPLAY_V1"
                && count > 0 && count <= 1000000, "invalid boundary replay header");
        std::vector<BoundaryInput> inputs(count);
        std::vector<BoundaryOutput> expected(count), actual(count);
        std::vector<unsigned long long> faces(count);
        std::vector<int> calls(count);
        int stage_count[6]{}, reversal_count[6]{};
        for (int b = 0; b < count; ++b) {
            auto& x = inputs[b];
            read_value(stream, x.stage); read_value(stream, calls[b]); read_value(stream, faces[b]);
            require(x.stage >= 0 && x.stage < 6 && calls[b] > 0 && faces[b] > 0, "invalid boundary identity");
            ++stage_count[x.stage];
            read_array(stream, x.face_nodes); read_array(stream, x.nearest);
            read_array(stream, x.opposing); read_array(stream, x.reversal);
            read_array(stream, x.area); read_array(stream, x.shape); read_array(stream, x.gradient);
            read_array(stream, x.velocity); read_array(stream, x.boundary_velocity); read_array(stream, x.viscosity);
            read_array(stream, x.density); read_array(stream, x.pressure); read_array(stream, x.pressure_gradient);
            read_array(stream, x.influence_lhs); read_array(stream, x.influence_rhs); read_array(stream, x.bc_multiplier);
            read_array(stream, x.stored_flux); read_array(stream, x.wall_coefficient);
            validate(x);
            for (int flag : x.reversal) reversal_count[x.stage] += flag;
            read_array(stream, expected[b].lhs); read_array(stream, expected[b].rhs);
        }
        for (int n : stage_count) require(n > 0, "missing replay stage");
        require(!(stream >> magic), "trailing replay data");
#ifdef MARS_REPLAY_CUDA
        DeviceBuffers device;
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.inputs), inputs.size()*sizeof(BoundaryInput)));
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&device.outputs), actual.size()*sizeof(BoundaryOutput)));
        check_cuda(cudaMemcpy(device.inputs, inputs.data(), inputs.size()*sizeof(BoundaryInput), cudaMemcpyHostToDevice));
        replay_kernel<<<(count+127)/128,128>>>(device.inputs, device.outputs, count);
        check_cuda(cudaGetLastError());
        check_cuda(cudaDeviceSynchronize());
        check_cuda(cudaMemcpy(actual.data(), device.outputs, actual.size()*sizeof(BoundaryOutput), cudaMemcpyDeviceToHost));
        const char* backend = "CUDA";
#else
        for (int b = 0; b < count; ++b) mars::segregated::boundary_block(inputs[b], actual[b]);
        const char* backend = "host";
#endif
        std::size_t failures = 0, checks = 0;
        double worst[6]{};
        for (int b = 0; b < count; ++b) {
            for (int field = 0; field < 2; ++field) {
                const int length = field == 0 ? 144 : 12;
                const double* value = field == 0 ? actual[b].lhs : actual[b].rhs;
                const double* target = field == 0 ? expected[b].lhs : expected[b].rhs;
                double scale = 1;
                for (int i = 0; i < length; ++i) scale = std::max(scale, std::abs(target[i]));
                for (int i = 0; i < length; ++i) {
                    ++checks;
                    const double error = std::abs(value[i]-target[i])/scale;
                    worst[inputs[b].stage] = std::max(worst[inputs[b].stage], error);
                    if (!std::isfinite(value[i]) || error > 1e-12) {
                        if (failures++ < 8) std::cerr << "mismatch stage=" << inputs[b].stage
                            << " call=" << calls[b] << " face=" << faces[b] << " field=" << field << " entry=" << i
                            << " expected=" << target[i] << " actual=" << value[i] << '\n';
                    }
                }
            }
        }
        for (int s = 0; s < 6; ++s)
            std::cout << "stage=" << s << " faces=" << stage_count[s] << " reversed_samples=" << reversal_count[s]
                      << " worst_scaled=" << worst[s] << '\n';
        require(failures == 0, "boundary replay mismatch");
        std::cout << "PASS: " << backend << " frozen boundary replay records=" << count << " scalar_checks=" << checks
                  << "; trace/reversal updates, global assembly and SIMPLE iteration not tested\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
