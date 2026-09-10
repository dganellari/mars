// GPT/Codex, 2026-09-09. Actual device evaluator with synthetic, public geometry.
#include <cuda_runtime.h>
#include "mars_outlet_boundary_gate.hpp"

namespace {

__global__ void evaluate_cases(const mars::outlet_gate::BoundaryInput* inputs,
                               mars::outlet_gate::BoundaryOutput* outputs, int count)
{
    const int index = blockIdx.x*blockDim.x+threadIdx.x;
    if (index < count) outputs[index] = mars::outlet_gate::evaluate(inputs[index]);
}

bool check_cuda(cudaError_t status, const char* operation)
{
    if (status == cudaSuccess) return true;
    std::fprintf(stderr, "FAIL %s: %s\n", operation, cudaGetErrorString(status));
    return false;
}

} // namespace

int main()
{
    using namespace mars::outlet_gate;
    if (run_host_checks()) return 1;
    std::vector<BoundaryInput> inputs;
    std::vector<ExpectedBoundary> expected;
    boundary_cases(inputs, expected);
    std::vector<BoundaryOutput> outputs(inputs.size());
    BoundaryInput* d_inputs = nullptr;
    BoundaryOutput* d_outputs = nullptr;
    bool ok = check_cuda(cudaMalloc(&d_inputs, inputs.size()*sizeof(BoundaryInput)), "allocate inputs");
    if (ok) ok = check_cuda(cudaMalloc(&d_outputs, outputs.size()*sizeof(BoundaryOutput)), "allocate outputs");
    if (ok) ok = check_cuda(cudaMemcpy(d_inputs, inputs.data(), inputs.size()*sizeof(BoundaryInput),
                                       cudaMemcpyHostToDevice), "publish fixture");
    if (ok) {
        evaluate_cases<<<(inputs.size()+127)/128,128>>>(d_inputs, d_outputs, static_cast<int>(inputs.size()));
        ok = check_cuda(cudaGetLastError(), "launch evaluator") &&
             check_cuda(cudaDeviceSynchronize(), "complete evaluator");
    }
    if (ok) ok = check_cuda(cudaMemcpy(outputs.data(), d_outputs, outputs.size()*sizeof(BoundaryOutput),
                                       cudaMemcpyDeviceToHost), "read gate results");
    if (d_inputs) ok = check_cuda(cudaFree(d_inputs), "free inputs") && ok;
    if (d_outputs) ok = check_cuda(cudaFree(d_outputs), "free outputs") && ok;
    if (!ok) return 1;
    Checks checks;
    check_boundary_outputs(inputs, expected, outputs, checks);
    std::printf("%s: %d CUDA checks of the production outlet evaluator\n",
                checks.failures ? "FAIL" : "PASS", checks.count);
    std::printf("Distributed scatter and one/two/four-rank flow gates remain separate.\n");
    return checks.failures ? 1 : 0;
}
