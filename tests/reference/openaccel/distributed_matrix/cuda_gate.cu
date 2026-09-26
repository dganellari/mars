// CUDA/Hypre gate: run_gates() on device buffers, then the real HypreGMRESSolver device-map
// solve through the owned-row adapter with a CUDA-aware ghost exchange before every true
// residual. Options: --probe-empty-rank (Hypre solve with a zero-row rank, policy allow),
// --bench N (timed update/exchange/residual on a larger fixture, N repetitions).
#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_gmres_solver.hpp"
#include "gate_common.hpp"
#include <cuda_runtime.h>
#include <iomanip>
#include <string>
using namespace dmatrix_gate;
using Solver=mars::fem::HypreGMRESSolver<double,int,cstone::GpuTag>;
using Matrix=Solver::Matrix;
using Vector=Solver::Vector;

template<int C> void hypre_gates(MPI_Comm comm,Report& report,Options o,EmptyRanks policy,const std::string& label) {
    int rank=0, ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    const std::string tag="C="+std::to_string(C)+" ranks="+std::to_string(ranks)+" "+label+": ";
    Problem p(C,ranks,o); Local l=extract(p,rank);
    fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,11); });
    Device<C,HYPRE_BigInt> d(l);
    OwnedRowSystem<C,Matrix,HYPRE_BigInt> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes(),policy);
    Vector b, x; b.resize(std::size_t(s.rows())); x.resize(std::size_t(s.rows()));
    Solver solver(comm,2000,1e-12,Solver::BOOMERAMG,100); solver.setVerbose(false); solver.setPointBlock(C);
    GhostExchange exchange(p,l,rank);
    for (unsigned round:{1u,2u}) {   // round 2: new values, RHS and solution; same structure
        if (round==2) { fill(p,l,rank,2,[&](int g,int c) { return p.product(g,c,2,12); }); overwrite(d.blocks,l.blocks); overwrite(d.rhs,l.rhs); }
        s.update(d.view(),b.data(),b.size());
        if (s.rows()) cudaMemset(x.data(),0,std::size_t(s.rows())*sizeof(double));
        const bool solved=solve_owned(solver,s,b,x);
        Buffer<double> local(C*std::size_t(l.nodes()),std::numeric_limits<double>::quiet_NaN());
        s.unpack(x.data(),x.size(),raw(local),local.size());
        const auto before=s.residual(halo_complete(raw(local),local.size()),b.data());   // ghosts still NaN
        exchange.run(raw(local),comm);
        const auto norms=s.residual(halo_complete(raw(local),local.size()),b.data());
        const auto h=download(raw(local),local.size()); double error=0;
        for (int v=0;v<l.nodes();++v) for (int c=0;c<C;++c)
            error=std::max(error,std::abs(h[std::size_t(C)*v+c]-p.solution(l.global[v],c,10+round)));
        double worst=0; MPI_Allreduce(&error,&worst,1,MPI_DOUBLE,MPI_MAX,comm);
        const long long ghosts=sum_all((long long)exchange.recv_nodes.size(),comm);
        std::ostringstream detail;
        detail<<"round "<<round<<" iterations="<<solver.getLastIterations()<<" relative="<<norms.relative()
              <<" max|x-x*| incl. ghosts="<<worst<<" unexchanged_finite="<<before.finite;
        report.result(tag+"Hypre GMRES device-map solve, exchanged true residual, source-keyed solution",
            all_true(solved,comm) && norms.passed && worst<=1e-8 && (ghosts==0 || !before.passed),detail.str());
    }
}

template<int C> void bench(MPI_Comm comm,Report& report,int repetitions) {
    int rank=0, ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    Options o; o.nx=48; o.ny=48; o.nz=24; o.far_ghosts=false;
    Problem p(C,ranks,o); Local l=extract(p,rank);
    fill(p,l,rank,1,[&](int g,int c) { return p.rhs(g,c,1); });
    Device<C,HYPRE_BigInt> d(l);
    OwnedRowSystem<C,Matrix,HYPRE_BigInt> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
    Buffer<double> b(std::size_t(s.rows())), x(C*std::size_t(l.nodes()),1.0);
    GhostExchange exchange(p,l,rank);
    cudaEvent_t start, stop; cudaEventCreate(&start); cudaEventCreate(&stop);
    auto timed=[&](auto f) {
        f(); MPI_Barrier(comm); cudaEventRecord(start);
        for (int i=0;i<repetitions;++i) f();
        cudaEventRecord(stop); cudaEventSynchronize(stop); float ms=0; cudaEventElapsedTime(&ms,start,stop);
        double local=ms/repetitions, worst=0; MPI_Allreduce(&local,&worst,1,MPI_DOUBLE,MPI_MAX,comm); return worst;
    };
    const double update=timed([&] { s.update(d.view(),raw(b),b.size()); });
    const double halo=timed([&] { exchange.run(raw(x),comm); });
    const double residual=timed([&] { s.residual(halo_complete(raw(x),x.size()),raw(b)); });
    const double nnz=double(sum_all((long long)s.nnz(),comm)), rows=double(sum_all((long long)s.rows(),comm));
    std::ostringstream detail; detail<<std::setprecision(4)<<"global nnz="<<nnz<<" rows="<<rows<<" slowest rank ms: update="<<update
        <<" exchange="<<halo<<" residual="<<residual<<" (update moves ~"<<(20*nnz+20*rows)/(update*1e6)/ranks<<" GB/s per rank)";
    report.result("C="+std::to_string(C)+" ranks="+std::to_string(ranks)+" bench (information only)",true,detail.str());
    cudaEventDestroy(start); cudaEventDestroy(stop);
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    Report report; MPI_Comm_rank(MPI_COMM_WORLD,&report.rank);
    MPI_Comm node; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,report.rank,MPI_INFO_NULL,&node);
    int node_rank=0, devices=0; MPI_Comm_rank(node,&node_rank); cudaGetDeviceCount(&devices);
    if (devices>0) cudaSetDevice(node_rank%devices);
    bool probe=false; int repetitions=0;
    for (int i=1;i<argc;++i) {
        const std::string a=argv[i];
        if (a=="--probe-empty-rank") probe=true;
        else if (a=="--bench" && i+1<argc) repetitions=std::stoi(argv[++i]);
    }
    try {
        if (probe) {
            int ranks=1; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
            if (ranks<2) report.skip("zero-row Hypre probe","needs >= 2 ranks");
            else {
                Options o; o.empty_rank=true;
                hypre_gates<1>(MPI_COMM_WORLD,report,o,EmptyRanks::allow,"zero-row rank probe");
                hypre_gates<3>(MPI_COMM_WORLD,report,o,EmptyRanks::allow,"zero-row rank probe");
            }
        } else {
            run_gates<1,Matrix,HYPRE_BigInt>(MPI_COMM_WORLD,report);
            run_gates<3,Matrix,HYPRE_BigInt>(MPI_COMM_WORLD,report);
            hypre_gates<1>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"uneven");
            hypre_gates<3>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"uneven");
            if (repetitions>0) { bench<1>(MPI_COMM_WORLD,report,repetitions); bench<3>(MPI_COMM_WORLD,report,repetitions); }
        }
    } catch (const std::exception& e) {
        std::cerr<<"rank "<<report.rank<<" unexpected exception: "<<e.what()<<std::endl; ++report.failures;
    }
    int failures=0; MPI_Allreduce(&report.failures,&failures,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (report.rank==0) std::cout<<(failures?"FAIL: ":"PASS: ")<<report.passes<<" passed, "<<report.failures<<" failed, "<<report.skips<<" skipped"<<std::endl;
    MPI_Comm_free(&node);
    MPI_Finalize();
    return failures?1:0;
}
