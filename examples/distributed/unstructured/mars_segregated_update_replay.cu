#include "evaluate.hpp"
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
using update_replay::Input;
using update_replay::Output;
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
template<class T> void read(std::istream& stream, T& value) {
    require(bool(stream >> value) && std::isfinite(double(value)), "invalid/truncated update data");
}
#ifdef MARS_REPLAY_CUDA
void check(cudaError_t error) { if(error!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(error)); }
__global__ void replay(const Input* in, Output* out, int count) {
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<count) update_replay::evaluate(in[i],out[i]);
}
struct Buffers {
    Input* inputs=nullptr; Output* outputs=nullptr;
    ~Buffers() { cudaFree(inputs); cudaFree(outputs); }
};
#endif
int main(int argc,char** argv) {
    try {
        require(argc==2,"usage: mars_segregated_update_replay PUBLIC_UPDATE_FILE");
        std::ifstream stream(argv[1]); std::string magic; int count=0;
        require(bool(stream>>magic>>count) && magic=="MARS_PUBLIC_UPDATE_REPLAY_V1" && count>0
                && count<=1000000,"invalid update header");
        std::vector<Input> inputs(count); std::vector<Output> actual(count),expected(count);
        int counts[10]{}; double worst[10]{};
        for(int i=0;i<count;++i) {
            auto& in=inputs[i]; int iteration,sample; unsigned long long entity;
            read(stream,in.stage); read(stream,iteration); read(stream,entity); read(stream,sample);
            require(in.stage>=0 && in.stage<10 && iteration>0 && entity>0 && sample>=0,"invalid update identity");
            ++counts[in.stage];
            for(int j=0;j<update_replay::input_width[in.stage];++j) read(stream,in.values[j]);
            require(update_replay::valid(in),"invalid update scale, area or flag");
            for(double& value:expected[i].values) read(stream,value);
        }
        for(int n:counts) require(n>0,"missing update stage");
        require(!(stream>>magic),"trailing update data");
#ifdef MARS_REPLAY_CUDA
        Buffers device;
        check(cudaMalloc(reinterpret_cast<void**>(&device.inputs),count*sizeof(Input)));
        check(cudaMalloc(reinterpret_cast<void**>(&device.outputs),count*sizeof(Output)));
        check(cudaMemcpy(device.inputs,inputs.data(),count*sizeof(Input),cudaMemcpyHostToDevice));
        replay<<<(count+127)/128,128>>>(device.inputs,device.outputs,count);
        check(cudaGetLastError()); check(cudaDeviceSynchronize());
        check(cudaMemcpy(actual.data(),device.outputs,count*sizeof(Output),cudaMemcpyDeviceToHost));
        const char* backend="CUDA";
#else
        for(int i=0;i<count;++i) update_replay::evaluate(inputs[i],actual[i]);
        const char* backend="host";
#endif
        std::size_t checks=0,failures=0;
        for(int i=0;i<count;++i) {
            // Scalar/flag outputs stay separately scaled so a large pressure cannot hide a flux defect.
            for(int j=0;j<6;++j) {
                ++checks;
                const double error=std::abs(actual[i].values[j]-expected[i].values[j])
                                   /std::max(1.0,std::abs(expected[i].values[j]));
                worst[inputs[i].stage]=std::max(worst[inputs[i].stage],error);
                if(!std::isfinite(actual[i].values[j]) || error>1e-12) {
                    if(failures++<8) std::cerr<<"mismatch record="<<i<<" stage="<<inputs[i].stage
                        <<" entry="<<j<<" expected="<<expected[i].values[j]<<" actual="<<actual[i].values[j]<<'\n';
                }
            }
        }
        for(int s=0;s<10;++s) std::cout<<"stage="<<s<<" records="<<counts[s]<<" worst_scaled="<<worst[s]<<'\n';
        require(!failures,"update replay mismatch");
        std::cout<<"PASS: "<<backend<<" frozen update replay records="<<count<<" scalar_checks="<<checks
                 <<"; reconstruction, global assembly and full SIMPLE iteration not tested\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
