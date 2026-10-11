#pragma once
#include "mars_segregated_simple_runtime.hpp"
#include <mpi.h>
#include <climits>
#ifdef MARS_REPLAY_CUDA
#include <thrust/binary_search.h>
#include <thrust/sort.h>
#endif
#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_EXPANSION_HD __host__ __device__
#else
#define MARS_EXPANSION_HD
#endif
namespace mars::segregated::runtime {
struct PressureIncidenceEntry {
    SimpleMesh mesh; bool boundary; int *keys,*entries;
    MARS_EXPANSION_HD void operator()(int i) const {
        const int entity=i/(boundary?3:4),local=i%(boundary?3:4);
        const auto f=boundary?mesh.faces[entity]:SimpleFace{entity,0,0};
        keys[i]=mesh.nodes[boundary?tet_face_node(f.ordinal,local):local][f.element]; entries[i]=i;
    }
};
// Persistent complete node stars let each node accumulate a pressure pair without atomics.
struct PressureIncidence {
    Array<int> offsets,entries;
    PressureIncidence(SimpleMesh mesh,bool boundary):offsets(std::size_t(mesh.node_count)+1),entries(count(mesh,boundary)) {
        const int count=int(entries.values.size()); Array<int> keys(count);
        launch(count,PressureIncidenceEntry{mesh,boundary,keys.data(),entries.data()});
#ifdef MARS_REPLAY_CUDA
        thrust::stable_sort_by_key(keys.values.begin(),keys.values.end(),entries.values.begin());
        thrust::lower_bound(keys.values.begin(),keys.values.end(),thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(mesh.node_count+1),offsets.values.begin());
#else
        std::vector<std::pair<int,int>> pairs;
        for (int i=0;i<count;++i) pairs.emplace_back(keys.values[i],entries.values[i]);
        std::sort(pairs.begin(),pairs.end());
        for (int i=0;i<count;++i) { keys.values[i]=pairs[i].first; entries.values[i]=pairs[i].second; }
        for (int i=0;i<=mesh.node_count;++i)
            offsets.values[i]=int(std::lower_bound(keys.values.begin(),keys.values.end(),i)-keys.values.begin());
#endif
    }
    static std::size_t count(SimpleMesh mesh,bool boundary) {
        const auto count=std::size_t(boundary?mesh.face_count:mesh.element_count)*(boundary?3:4);
        ensure(count<=std::size_t(INT_MAX),"pressure incidence exceeds index range"); return count;
    }
};
struct ExpandedGradient {
    SimpleMesh mesh; const int *offsets,*entries;
    const double *high,*low,*volume; double *gradient,*gradient_low; int* error;
    MARS_EXPANSION_HD void operator()(int n) const {
        PressureValue sum[3];
        for (int k=offsets[n];k<offsets[n+1];++k) {
            const int e=entries[k]/4,local=entries[k]%4;
            for (int s=0;s<6;++s) {
                const int left=tet_edge_node(s,0),right=tet_edge_node(s,1);
                if (local!=left && local!=right) continue;
                const auto difference=(PressureValue::load(high,low,mesh.nodes[right][e])
                    -PressureValue::load(high,low,mesh.nodes[left][e]))*.5;
                for (int j=0;j<3;++j) sum[j]=sum[j]+difference*mesh.geometry[e].area[3*s+j];
            }
        }
        // Shifted Tri3 samples are nodal: their incremental boundary numerator is zero.
        if (!(volume[n]>0) || !geometry_finite(volume[n])) { simple_error(error); return; }
        for (int j=0;j<3;++j) (sum[j]/volume[n]).store(gradient,gradient_low,3*n+j);
    }
};
struct ExpandedPressureUpdate {
    double *high,*low; const double *increment,*increment_low; double alpha;
    MARS_EXPANSION_HD void operator()(int i) const {
        (PressureValue::load(high,low,i)+PressureValue::load(increment,increment_low,i)*alpha).store(high,low,i);
    }
};
struct ExpandedVelocityUpdate {
    SimpleState state; const double *gradient,*gradient_low;
    MARS_EXPANSION_HD void operator()(int n) const {
        for (int j=0;j<3;++j) {
            const int i=3*n+j;
            state.velocity[i]=(PressureValue{state.velocity[i],0}
                -PressureValue::load(gradient,gradient_low,i)*state.influence[i]).rounded();
        }
    }
};
struct PressureMoment {
    PressureValue weighted,area;
};
struct PressureMomentAdd {
    MARS_EXPANSION_HD PressureMoment operator()(PressureMoment a,PressureMoment b) const {
        return {a.weighted+b.weighted,a.area+b.area};
    }
};
struct ExpandedMoment {
    SimpleMesh mesh; SimpleState state; const int* owned;
    MARS_EXPANSION_HD PressureMoment operator()(int k) const {
        const int i=owned[k]; const auto f=mesh.faces[i]; PressureMoment result;
        if (f.kind!=1 || (state.reversal && state.reversal[3*i])) return result;
        double a[3]; tet_boundary_area(mesh.geometry[f.element],f.ordinal,a);
        const double area=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        for (int j=0;j<3;++j) {
            const int n=mesh.nodes[tet_face_node(f.ordinal,j)][f.element];
            result.weighted=result.weighted+PressureValue::load(state.pressure,state.pressure_low,n)*area;
            result.area=result.area+PressureValue{area,0};
        }
        return result;
    }
};
class PressureMomentReduction {
    MPI_Datatype type_=MPI_DATATYPE_NULL; MPI_Op op_=MPI_OP_NULL;
    static void add(void* input,void* output,int* count,MPI_Datatype*) {
        const auto* in=static_cast<PressureMoment*>(input); auto* out=static_cast<PressureMoment*>(output);
        for (int i=0;i<*count;++i) out[i]=PressureMomentAdd{}(in[i],out[i]);
    }
public:
    PressureMomentReduction() {
        static_assert(sizeof(PressureMoment)==4*sizeof(double));
        ensure(MPI_Type_contiguous(4,MPI_DOUBLE,&type_)==MPI_SUCCESS,"pressure moment datatype failed");
        ensure(MPI_Type_commit(&type_)==MPI_SUCCESS && MPI_Op_create(add,0,&op_)==MPI_SUCCESS,"pressure moment reduction setup failed");
    }
    ~PressureMomentReduction() { if(op_!=MPI_OP_NULL) MPI_Op_free(&op_); if(type_!=MPI_DATATYPE_NULL) MPI_Type_free(&type_); }
    PressureMomentReduction(const PressureMomentReduction&)=delete;
    PressureValue mean(MPI_Comm comm,int faces,ExpandedMoment f) const {
        PressureMoment local;
#ifdef MARS_REPLAY_CUDA
        local=thrust::transform_reduce(thrust::make_counting_iterator(0),thrust::make_counting_iterator(faces),f,PressureMoment{},PressureMomentAdd{});
#else
        for(int i=0;i<faces;++i) local=PressureMomentAdd{}(local,f(i));
#endif
        // Only four control scalars leave the device, for the paired MPI operation.
        PressureMoment global;
        if (MPI_Allreduce(&local,&global,1,type_,op_,comm)!=MPI_SUCCESS) {
            MPI_Abort(comm,1); throw std::runtime_error("pressure moment reduction failed");
        }
        const double area=global.area.rounded();
        ensure(area>0 && std::isfinite(area) && std::isfinite(global.weighted.high) && std::isfinite(global.weighted.low),
            "all outlet faces closed or nonfinite pressure moment");
        const auto first=global.weighted/area;
        // Account for the low part of area instead of rounding the denominator first.
        const auto remainder=global.weighted-first*global.area.high-first*global.area.low;
        return first+remainder/area;
    }
};
struct ExpandedTrace {
    SimpleMesh mesh; SimpleState state; SimpleControls controls; PressureValue mean;
    MARS_EXPANSION_HD void operator()(int i) const {
        const auto f=mesh.faces[i]; if (f.kind!=1 || (state.reversal && state.reversal[3*i])) return;
        for (int j=0;j<3;++j) {
            const int n=mesh.nodes[tet_face_node(f.ordinal,j)][f.element];
            (PressureValue{controls.pressure_reference,0}+(PressureValue::load(state.pressure,state.pressure_low,n)-mean)*(1-controls.beta))
                .store(state.trace,state.trace_low,3*i+j);
        }
    }
};
struct ExpandedOutletPressure {
    SimpleMesh mesh; SimpleState state; const int *offsets,*entries; const double* area;
    double *high,*low;
    MARS_EXPANSION_HD void operator()(int n) const {
        PressureValue sum;
        for (int k=offsets[n];k<offsets[n+1];++k) {
            const int i=entries[k]/3,j=entries[k]%3; const auto f=mesh.faces[i]; if(f.kind!=1) continue;
            double a[3]; tet_boundary_area(mesh.geometry[f.element],f.ordinal,a);
            const double magnitude=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            sum=sum+PressureValue::load(state.trace,state.trace_low,3*i+j)*magnitude;
        }
        const auto value=area[n]>0?sum/area[n]:PressureValue{};
        value.store(high,low,n);
        if(!geometry_finite(value.high) || !geometry_finite(value.low) || !geometry_finite(area[n]) || area[n]<0) simple_error(state.error);
    }
};
}
#undef MARS_EXPANSION_HD
