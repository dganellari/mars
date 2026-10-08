#include "../distributed_simple/channel.hpp"
#include <iostream>
#include <sstream>

using namespace mars::segregated;
using namespace mars::segregated::runtime;

struct HostMatrix {
    std::vector<int> offsets,columns;
    std::vector<double> values;
    void allocate(int rows,int,int nnz) { offsets.resize(rows+1); columns.resize(nnz); values.resize(nnz); }
    int* rowOffsetsPtr() { return offsets.data(); }
    const int* rowOffsetsPtr() const { return offsets.data(); }
    int* colIndicesPtr() { return columns.data(); }
    const int* colIndicesPtr() const { return columns.data(); }
    double* valuesPtr() { return values.data(); }
    const double* valuesPtr() const { return values.data(); }
};

double truth(long long dof) { return .125*(1+dof%17); }

// A known candidate separates the backend verdict from the original CSR check.
template<int C> struct CandidateSolve {
    MPI_Comm comm;
    int mode=0;
    bool configured=false;
    double relative=0,absolute=0;
    std::vector<double> b,x;
    explicit CandidateSolve(MPI_Comm c):comm(c) {}
    void set_tolerances(double r,double a) { configured=true; relative=r; absolute=a; }
    double* rhs(std::size_t n) { b.resize(n); return b.data(); }
    const double* rhs() const { return b.data(); }
    const double* solution() const { return x.data(); }
    std::size_t size() const { return x.size(); }
    template<class System> bool operator()(const System& system) {
        x.resize(system.rows());
        for (std::size_t i=0;i<x.size();++i) x[i]=truth(C*system.first_solver_node()+i);
        int rank=0,ranks=0; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
        if (rank==ranks-1) {
            if (mode==2 || mode==3) x.at(0)+=.25;
            if (mode==4) x.clear();
            if (mode==5) x.at(0)=std::numeric_limits<double>::quiet_NaN();
        }
        return mode==0 || mode==3 || rank!=ranks-1;
    }
};

using Runner=DistributedSimpleRunner<HostMatrix,long long,CandidateSolve>;

template<int C,class System> void cases(Runner& run,System& system,CandidateSolve<C>& solver,bool audit=false) {
    Array<double> blocks(static_cast<std::size_t>(C*C)*run.graph.blocks()),rhs(static_cast<std::size_t>(C)*run.n),increment(static_cast<std::size_t>(C)*run.n);
    blocks.zero(); rhs.zero();
    const auto a=run.graph.template view<C>(blocks.data(),rhs.data());
    for (int row=0;row<run.n;++row) for (int k=a.offsets[row];k<a.offsets[row+1];++k)
        for (int c=0;c<C;++c) {
            const double value=a.columns[k]==row?2+c:-.0625;
            a.values[C*C*k+C*c+c]=value;
            a.rhs[C*row+c]+=value*truth(C*run.solver_node.values[a.columns[k]]+c);
        }
    int rank=0; MPI_Comm_rank(run.comm,&rank);
    for (int mode=0;mode<6;++mode) {
        solver.mode=mode;
        std::ostringstream message;
        auto* previous=std::cerr.rdbuf(message.rdbuf());
        std::string failure;
        try { run.template solve<C>(system,solver,a,increment,false,{}); }
        catch (const std::exception& error) { failure=error.what(); }
        std::cerr.rdbuf(previous);
        simple_collective(run.comm,failure.empty()==(mode==0),"linear rejection verdict changed");
        if (mode==0 || mode==1) {
            bool correct=true;
            for (int i=0;i<C*run.n;++i)
                correct=correct && increment.values[i]==truth(C*run.solver_node.values[i/C]+i%C);
            simple_collective(run.comm,correct,"candidate was not published to ghost nodes");
        }
        bool report=true;
        if (rank==0 && mode!=0) {
            if (mode==4) report=failure.find("no usable candidate")!=std::string::npos;
            else report=message.str().find("[simple-linear]")!=std::string::npos &&
                message.str().find(mode==1?"mars_passed=1":"mars_passed=0")!=std::string::npos &&
                message.str().find(mode==3?"solver_accepted=1":"solver_accepted=0")!=std::string::npos;
        }
        simple_collective(run.comm,report,"missing or incorrect independent rejection diagnostic");
        if (rank==0) report=(message.str().find("[simple-pressure-audit]")!=std::string::npos)
            ==(audit && C==1 && mode!=0 && mode!=4);
        simple_collective(run.comm,report,"pressure audit opt-in or failure-only contract changed");
    }
}

template<int C,class System> void pressure_only_case(Runner& run,System& system,CandidateSolve<C>& solver) {
    Array<double> blocks(static_cast<std::size_t>(C*C)*run.graph.blocks()),rhs(static_cast<std::size_t>(C)*run.n),increment(static_cast<std::size_t>(C)*run.n);
    blocks.zero(); rhs.zero();
    const auto a=run.graph.template view<C>(blocks.data(),rhs.data());
    for (int row=0;row<run.n;++row) for (int k=a.offsets[row];k<a.offsets[row+1];++k)
        if (a.columns[k]==row) for (int c=0;c<C;++c) {
            a.values[C*C*k+C*c+c]=1;
            a.rhs[C*row+c]=truth(C*run.solver_node.values[row]+c);
        }
    for (int mode:{2,3}) {
        solver.mode=mode; bool accepted=true;
        try { run.template solve<C>(system,solver,a,increment,false,{}); }
        catch (const std::exception&) { accepted=false; }
        // This .25 error meets the explicit pressure atol but not momentum's default.
        simple_collective(run.comm,accepted==(C==1 && mode==3),"pressure override changed the wrong verdict");
    }
}

void pressure_configuration(Runner& run,int rank,int ranks) {
    auto reject=[&](bool enabled,double relative,double absolute) {
        bool rejected=false;
        try { run.set_pressure_tolerances(enabled,relative,absolute); }
        catch (const std::exception&) { rejected=true; }
        simple_collective(run.comm,rejected,"invalid pressure configuration accepted");
    };
    reject(true,rank==ranks-1?0:1e-6,0);
    reject(true,1e-6,rank==ranks-1?std::numeric_limits<double>::infinity():0);
    if (ranks>1) {
        reject(rank==ranks-1,1e-6,0);
        reject(true,rank==ranks-1?1e-5:1e-6,0);
        reject(true,1e-6,rank==ranks-1?1e-5:0);
    }
    run.set_pressure_tolerances(false,1e-12,0);
    simple_collective(run.comm,!run.poisson_solve.configured && !run.pressure_tolerance,"default pressure target changed");
    run.set_pressure_tolerances(true,1e-6,.5);
    simple_collective(run.comm,run.poisson_solve.configured && run.poisson_solve.relative==1e-6
        && run.poisson_solve.absolute==.5 && !run.momentum_solve.configured
        && run.tolerance.absolute==1e-13 && run.tolerance.relative==1e-10 && !run.tolerance.maximum,
        "pressure targets not isolated from momentum");
    pressure_only_case<1>(run,run.poisson,run.poisson_solve);
    pressure_only_case<3>(run,run.momentum,run.momentum_solve);
}

template<int C> struct RefiningSolve:CandidateSolve<C> {
    int calls=0,fault=0;
    using CandidateSolve<C>::CandidateSolve;
    template<class System,class Publish>
    PressureRefinementResult refine(System& system,Publish publish,distributed::Tolerance tolerance) {
        static_assert(C==1);
        ++calls;
        struct Operations {
            RefiningSolve& solve; System& system; Publish& publish; distributed::Tolerance tolerance;
            std::vector<double> defect_values,trial;
            auto defect() {
                defect_values.resize(solve.b.size());
                return system.compensated_defect(publish(solve.x),solve.b.data(),defect_values.data(),defect_values.size(),tolerance);
            }
            PressureCorrectionResult correct(int) {
                int rank=0,ranks=1; MPI_Comm_rank(solve.comm,&rank); MPI_Comm_size(solve.comm,&ranks);
                const int local=solve.fault==1 && rank==ranks-1?1:0; int any=0;
                simple_max(solve.comm,&local,&any,1);
                return {!any,1};
            }
            auto candidate() {
                trial.resize(solve.x.size());
                for (std::size_t i=0;i<trial.size();++i) trial[i]=truth(system.first_solver_node()+i);
                if (solve.fault==2) trial=solve.x;
                return system.compensated_defect(publish(trial),solve.b.data(),defect_values.data(),defect_values.size(),tolerance);
            }
            bool verify() { return solve.fault!=3; }
            void keep() { solve.x.swap(trial); }
            void restore() { publish(solve.x); }
        } op{*this,system,publish,tolerance,{}, {}};
        return refine_pressure(op,1,10);
    }
};

template<class Part> void refinement_integration(MPI_Comm comm,const Part& part) {
    DistributedSimpleRunner<HostMatrix,long long,RefiningSolve> run(comm,part.input,part.ownership);
    Array<double> blocks(run.graph.blocks()),rhs(run.n),increment(run.n);
    blocks.zero(); rhs.zero();
    const auto a=run.graph.template view<1>(blocks.data(),rhs.data());
    for (int row=0;row<run.n;++row) for (int k=a.offsets[row];k<a.offsets[row+1];++k) {
        const double value=a.columns[k]==row?4.:-.0625;
        a.values[k]=value; a.rhs[row]+=value*truth(run.solver_node.values[a.columns[k]]);
    }
    for (int fault=0;fault<4;++fault) {
        run.poisson_solve.mode=2; run.poisson_solve.fault=fault;
        bool accepted=true;
        try { run.template solve<1>(run.poisson,run.poisson_solve,a,increment,false,{}); }
        catch (const std::runtime_error&) { accepted=false; }
        simple_collective(comm,accepted==(fault==0),"pressure refinement integration verdict incorrect");
        simple_collective(comm,run.poisson_solve.calls==fault+1,"refinement was not collective");
        if (accepted) {
            bool exact=true;
            for (int i=0;i<run.n;++i) exact=exact && increment.values[i]==truth(run.solver_node.values[i]);
            simple_collective(comm,exact,"refinement left incorrect owned or ghost values");
        }
    }
    const int before=run.poisson_solve.calls;
    run.poisson_solve.mode=0;
    run.template solve<1>(run.poisson,run.poisson_solve,a,increment,false,{});
    simple_collective(comm,run.poisson_solve.calls==before,"successful solve entered refinement");
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        int rank=0,ranks=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        const auto mesh=dsimple_gate::channel(8,2,2);
        const auto part=dsimple_gate::extract(mesh,dsimple_gate::partition(mesh,ranks),rank);
        Runner run(MPI_COMM_WORLD,part.input,part.ownership);
        unsetenv("MARS_SIMPLE_PRESSURE_AUDIT");
        cases<3>(run,run.momentum,run.momentum_solve);
        cases<1>(run,run.poisson,run.poisson_solve);
        // A request on any rank enables the same failure collectives on all ranks.
        if (rank==ranks-1) setenv("MARS_SIMPLE_PRESSURE_AUDIT","1",1);
        cases<3>(run,run.momentum,run.momentum_solve,true);
        cases<1>(run,run.poisson,run.poisson_solve,true);
        pressure_configuration(run,rank,ranks);
        refinement_integration(MPI_COMM_WORLD,part);
        if (!rank) std::cout<<"PASS: accepted/rejected, wrong, missing and nonfinite candidates; original CSR and halo residual; ranks="<<ranks<<'\n';
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); return 1;
    }
    MPI_Finalize();
}
