#pragma once
#include <mpi.h>
#include <array>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#ifdef MARS_REPLAY_CUDA
#include <cuda_runtime.h>
#endif

namespace mars::segregated::runtime {
// CUDA events measure stream intervals, including idle gaps, not SM busy time.
// Profiling is opt-in; the normal path creates no events and adds no fences.
struct SimpleProfile {
    enum Phase { assembly, diagnostics, momentum, outlet, pressure_assembly, pressure, correction, phase_count };
    struct Total { double wall=0,stream=0; unsigned long long calls=0; };
    std::array<Total,phase_count> totals{};
    std::array<std::array<double,6>,2> linear{};
    bool enabled=false;
    int warmup=10;
    unsigned long long samples=0;
private:
    struct Sample {
        Phase phase=assembly;
        double wall=0;
#ifdef MARS_REPLAY_CUDA
        cudaEvent_t begin=nullptr,end=nullptr;
#endif
    };
    std::array<Sample,32> pending_{};
    int used_=0;
#ifdef MARS_REPLAY_CUDA
    cudaError_t event_error_=cudaSuccess;
    static void check(cudaError_t result) {
        if (result!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
    }
#endif
public:
    SimpleProfile()=default;
    SimpleProfile(const SimpleProfile&)=delete;
    SimpleProfile& operator=(const SimpleProfile&)=delete;
    ~SimpleProfile() {
#ifdef MARS_REPLAY_CUDA
        for (auto& s:pending_) { if (s.begin) cudaEventDestroy(s.begin); if (s.end) cudaEventDestroy(s.end); }
#endif
    }
    void configure(bool on,int excluded_iterations=10) {
        if (used_) throw std::runtime_error("cannot reconfigure active SIMPLE profiling");
        enabled=on; warmup=excluded_iterations;
    }
    struct Scope {
        SimpleProfile* profile;
        int index=-1;
        double start=0;
        Scope(SimpleProfile& p,Phase phase):profile(&p) {
            if (!p.enabled) return;
            if (p.used_==int(p.pending_.size())) throw std::runtime_error("collect SIMPLE timings once per iteration");
            index=p.used_++;
            auto& sample=p.pending_[index]; sample.phase=phase;
#ifdef MARS_REPLAY_CUDA
            if (!sample.begin) check(cudaEventCreate(&sample.begin));
            if (!sample.end) check(cudaEventCreate(&sample.end));
            check(cudaEventRecord(sample.begin));
#endif
            start=MPI_Wtime();
        }
        Scope(const Scope&)=delete;
        Scope& operator=(const Scope&)=delete;
        ~Scope() noexcept {
            if (index<0) return;
            auto& sample=profile->pending_[index]; sample.wall=MPI_Wtime()-start;
#ifdef MARS_REPLAY_CUDA
            const auto result=cudaEventRecord(sample.end);
            if (result!=cudaSuccess) profile->event_error_=result;
#endif
        }
    };
    Scope scope(Phase phase) { return Scope(*this,phase); }
    void collect(int completed) {
        if (!enabled) return;
#ifdef MARS_REPLAY_CUDA
        check(event_error_);
#endif
        const bool keep=completed>warmup;
        if (keep) ++samples;
        for (int i=0;i<used_;++i) {
            const auto& s=pending_[i]; double seconds=0;
#ifdef MARS_REPLAY_CUDA
            // Collect after the fixed-size convergence report has returned from the GPU.
            check(cudaEventSynchronize(s.end));
            float ms=0; check(cudaEventElapsedTime(&ms,s.begin,s.end)); seconds=double(ms)*.001;
#endif
            if (keep) { auto& t=totals[s.phase]; t.wall+=s.wall; t.stream+=seconds; ++t.calls; }
        }
        used_=0;
    }
    template<class Solver> void record_linear(int component,const Solver& solver,int completed) {
        if (!enabled || completed<=warmup) return;
        const auto& t=solver.get_last_timing(); auto& result=linear[component==3?0:1];
        result[0]+=t.prepare_seconds; result[1]+=t.packing_seconds; result[2]+=t.setup_seconds;
        result[3]+=t.solve_seconds; result[4]+=t.finish_seconds; result[5]+=solver.getLastIterations();
    }
    void write(MPI_Comm comm,std::ostream& out) const {
        if (!enabled) return;
        const char* names[]={"assembly","diagnostics","momentum","outlet","pressure_assembly","pressure","correction"};
        double local[2*phase_count],maximum[2*phase_count];
        unsigned long long calls[phase_count],max_calls[phase_count];
        double linear_max[12];
        for (int i=0;i<phase_count;++i) { local[2*i]=totals[i].wall; local[2*i+1]=totals[i].stream; calls[i]=totals[i].calls; }
        if (MPI_Reduce(local,maximum,2*phase_count,MPI_DOUBLE,MPI_MAX,0,comm)!=MPI_SUCCESS ||
            MPI_Reduce(calls,max_calls,phase_count,MPI_UNSIGNED_LONG_LONG,MPI_MAX,0,comm)!=MPI_SUCCESS)
            throw std::runtime_error("SIMPLE profile reduction failed");
        std::array<double,12> linear_local{};
        for (int s=0;s<2;++s) for (int i=0;i<6;++i) linear_local[6*s+i]=linear[s][i];
        if (MPI_Reduce(linear_local.data(),linear_max,12,MPI_DOUBLE,MPI_MAX,0,comm)!=MPI_SUCCESS)
            throw std::runtime_error("SIMPLE linear profile reduction failed");
        int rank; MPI_Comm_rank(comm,&rank);
        if (rank) return;
        out<<std::setprecision(9)<<"[simple-profile] samples="<<samples<<" excluded_iterations="<<warmup
           <<" scope=rank_max_totals stream_intervals_include_idle_gaps=1\n";
        for (int i=0;i<phase_count;++i)
            out<<"[simple-profile] phase="<<names[i]<<" calls="<<max_calls[i]<<" wall_seconds="<<maximum[2*i]
               <<" stream_seconds="<<maximum[2*i+1]<<'\n';
        for (int s=0;s<2;++s)
            out<<"[simple-profile] linear="<<(s?"pressure":"momentum")<<" prepare_seconds="<<linear_max[6*s]
               <<" packing_seconds="<<linear_max[6*s+1]<<" setup_seconds="<<linear_max[6*s+2]
               <<" solve_seconds="<<linear_max[6*s+3]<<" finish_seconds="<<linear_max[6*s+4]
               <<" krylov_iterations="<<linear_max[6*s+5]<<'\n';
    }
};
}
