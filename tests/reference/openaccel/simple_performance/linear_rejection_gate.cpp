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
    std::vector<double> b,x;
    explicit CandidateSolve(MPI_Comm c):comm(c) {}
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

template<int C,class System> void cases(Runner& run,System& system,CandidateSolve<C>& solver) {
    Array<double> blocks(std::size_t(C*C)*run.graph.blocks()),rhs(std::size_t(C)*run.n),increment(std::size_t(C)*run.n);
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
    }
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        int rank=0,ranks=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        const auto mesh=dsimple_gate::channel(8,2,2);
        const auto part=dsimple_gate::extract(mesh,dsimple_gate::partition(mesh,ranks),rank);
        Runner run(MPI_COMM_WORLD,part.input,part.ownership);
        cases<3>(run,run.momentum,run.momentum_solve);
        cases<1>(run,run.poisson,run.poisson_solve);
        if (!rank) std::cout<<"PASS: accepted/rejected, wrong, missing and nonfinite candidates; original CSR and halo residual; ranks="<<ranks<<'\n';
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); return 1;
    }
    MPI_Finalize();
}
