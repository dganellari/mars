#pragma once
#include "mars_segregated_simple_runtime.hpp"
#ifdef MARS_REPLAY_CUDA
#include <cub/device/device_reduce.cuh>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
namespace mars::segregated::runtime {
struct PackSimpleSums {
    const SimpleSums* partial; double* values;
    __device__ void operator()(int) const {
        SimpleSums a;
        for (int i=0;i<4;++i) a=SimpleSumCombine{}(a,partial[i]);
        values[0]=a.volume; values[1]=a.momentum2; values[2]=a.continuity2; values[3]=a.continuity;
        values[4]=a.velocity_change2; values[5]=a.pressure_change2; values[6]=a.inlet; values[7]=a.outlet;
        values[8]=a.inlet_area; values[9]=a.closed; values[10]=a.changed; values[11]=a.invalid;
        values[12]=a.speed2; values[13]=a.flux_change;
    }
};
struct UnpackSimpleSums {
    const double* values; SimpleSums* result;
    __device__ void operator()(int) const {
        SimpleSums a; a.volume=values[0]; a.momentum2=values[1]; a.continuity2=values[2]; a.continuity=values[3];
        a.velocity_change2=values[4]; a.pressure_change2=values[5]; a.inlet=values[6]; a.outlet=values[7]; a.inlet_area=values[8];
        a.closed=int(values[9]); a.changed=int(values[10]); a.invalid=int(values[11]); a.speed2=values[12]; a.flux_change=values[13];
        *result=a;
    }
};
struct SimpleDeviceReduction {
    thrust::device_vector<unsigned char> scratch;
    Array<SimpleSums> partial{4},result{1};
    Array<double> packed{14};
    int counts[4]={-1,-1,-1,-1};
    size_t sizes[4]={};
    template<class F> void reduce(int count,F f,int slot) {
        auto input=thrust::make_transform_iterator(thrust::make_counting_iterator(0),f);
        if (counts[slot]!=count) {
            assembly_cuda_check(cub::DeviceReduce::Reduce(nullptr,sizes[slot],input,partial.data()+slot,count,SimpleSumCombine{},SimpleSums{}));
            counts[slot]=count;
        }
        size_t bytes=sizes[slot];
        if (scratch.size()<bytes) scratch.resize(bytes);
        assembly_cuda_check(cub::DeviceReduce::Reduce(thrust::raw_pointer_cast(scratch.data()),bytes,input,partial.data()+slot,count,SimpleSumCombine{},SimpleSums{}));
    }
    void finish(MPI_Comm comm) {
        launch(1,PackSimpleSums{partial.data(),packed.data()});
        assembly_cuda_check(cudaStreamSynchronize(nullptr));
        MPI_Request requests[2];
        ensure(MPI_Iallreduce(MPI_IN_PLACE,packed.data(),12,MPI_DOUBLE,MPI_SUM,comm,requests)==MPI_SUCCESS,"SIMPLE diagnostic sum failed");
        ensure(MPI_Iallreduce(MPI_IN_PLACE,packed.data()+12,2,MPI_DOUBLE,MPI_MAX,comm,requests+1)==MPI_SUCCESS,"SIMPLE diagnostic max failed");
        ensure(MPI_Waitall(2,requests,MPI_STATUSES_IGNORE)==MPI_SUCCESS,"SIMPLE diagnostic reductions failed");
        launch(1,UnpackSimpleSums{packed.data(),result.data()});
    }
};
struct SimpleReport { SimpleSums sums; SimpleMetrics metrics; double speed; bool converged; };
struct FinishSimpleReport {
    const SimpleSums* sums; SimpleControls controls; int iterations; double residual,mass,change; SimpleReport* report;
    __device__ void operator()(int) const {
        report->sums=*sums; report->metrics=simple_metrics(*sums,controls); report->speed=sqrt(sums->speed2);
        report->converged=simple_converged(report->metrics,iterations,sums->changed,residual,mass,change);
    }
};
}
#endif
