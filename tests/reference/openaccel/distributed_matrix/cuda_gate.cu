// CUDA/Hypre gate: run_gates() on device buffers, then the real HypreGMRESSolver device-map
// solve through the owned-row adapter with a CUDA-aware ghost exchange before every true
// residual. Options: --probe-empty-rank (Hypre solve with a zero-row rank, policy allow),
// --bench N (timed update/exchange/residual on a larger fixture, N repetitions).
#ifndef MARS_REPLAY_CUDA
#define MARS_REPLAY_CUDA
#endif
#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_gmres_solver.hpp"
#include "gate_common.hpp"
#include "mars_segregated_simple_distributed.hpp"
#include "../simple_performance/pressure_profile_fixture.hpp"
#include <cuda_runtime.h>
#include <iomanip>
#include <string>
using namespace dmatrix_gate;
using Solver=mars::fem::HypreGMRESSolver<double,int,cstone::execution::Gpu>;
using Matrix=Solver::Matrix;
using Vector=Solver::Vector;

template<int C> void hypre_gates(MPI_Comm comm,Report& report,Options o,EmptyRanks policy,const std::string& label,int coarse_relax=-1,bool explicit_target=false,int profile_mode=0) {
    int rank=0, ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    const std::string tag="C="+std::to_string(C)+" ranks="+std::to_string(ranks)+" "+label+": ";
    if (profile_mode) o.diffusion=true;
    Problem p(C,ranks,o); Local l=extract(p,rank);
    fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,11); });
    Device<C,HYPRE_BigInt> d(l);
    OwnedRowSystem<C,Matrix,HYPRE_BigInt> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes(),policy);
    Vector b, x; b.resize(std::size_t(s.rows())); x.resize(std::size_t(s.rows()));
    Solver solver(comm,2000,1e-12,Solver::BOOMERAMG,100); solver.setVerbose(false); solver.setPointBlock(C);
    solver.setAMGCoarseRelaxType(coarse_relax);
    const Tolerance acceptance=explicit_target?Tolerance{1e-13,1e-10,true}:Tolerance{};
    if (explicit_target) {
        solver.set_stopping_tolerances(acceptance.relative,acceptance.absolute);
        solver.enable_true_residual_check(acceptance.absolute,acceptance.relative,true);
        solver.enable_fixed_graph_updates();
    }
    const auto profile=pressure_profile_fixture(profile_mode==3,profile_mode==1);
    if (profile_mode) mars::fem::pressure_settings::configure_gpu_profile(solver,profile);
    GhostExchange exchange(p,l,rank);
    for (unsigned round:{1u,2u}) {   // round 2: new values, RHS and solution; same structure
        if (round==2) { fill(p,l,rank,2,[&](int g,int c) { return p.product(g,c,2,12); }); overwrite(d.blocks,l.blocks); overwrite(d.rhs,l.rhs); }
        s.update(d.view(),b.data(),b.size());
        if (s.rows()) cudaMemset(x.data(),0,std::size_t(s.rows())*sizeof(double));
        const bool solved=solve_owned(solver,s,b,x);
        if (profile_mode) {
            bool controls_match=false,levels_match=false;
            std::ostringstream control_detail,level_detail;
            solver.inspect_prepared([&](auto krylov,auto amg,bool flexible) {
                const auto actual=mars::fem::pressure_settings::snapshot(krylov,amg,flexible);
                controls_match=pressure_profile_controls_match(profile,actual,control_detail);
                levels_match=pressure_profile_levels_match(actual,profile_mode==1);
                level_detail<<"round="<<round<<" expected="<<(profile_mode==1?"1":">1")
                            <<" actual="<<actual.at("effective_levels");
            });
            // Synthetic settings are public; print the failing rank's values as well.
            if (!controls_match) std::cerr<<"rank "<<rank<<" "<<tag<<control_detail.str()<<'\n';
            if (!levels_match) std::cerr<<"rank "<<rank<<" "<<tag<<level_detail.str()<<'\n';
            report.result(tag+"actual pressure controls",all_true(controls_match,comm));
            report.result(tag+"hierarchy depth",all_true(levels_match,comm),level_detail.str());
        }
        Buffer<double> local(C*std::size_t(l.nodes()),std::numeric_limits<double>::quiet_NaN());
        s.unpack(x.data(),x.size(),raw(local),local.size());
        const auto before=s.residual(halo_complete(raw(local),local.size()),b.data());   // ghosts still NaN
        exchange.run(raw(local),comm);
        const auto norms=s.residual(halo_complete(raw(local),local.size()),b.data(),acceptance);
        const auto h=download(raw(local),local.size()); double error=0;
        for (int v=0;v<l.nodes();++v) for (int c=0;c<C;++c)
            error=std::max(error,std::abs(h[std::size_t(C)*v+c]-p.solution(l.global[v],c,10+round)));
        double worst=0; MPI_Allreduce(&error,&worst,1,MPI_DOUBLE,MPI_MAX,comm);
        const long long ghosts=sum_all((long long)exchange.recv_nodes.size(),comm);
        std::ostringstream detail;
        detail<<"round "<<round<<" iterations="<<solver.getLastIterations()<<" relative="<<norms.relative()
              <<" max|x-x*| incl. ghosts="<<worst<<" unexchanged_finite="<<before.finite;
        report.result(tag+"Hypre GMRES device-map solve, exchanged true residual, source-keyed solution",
            all_true(solved,comm) && norms.passed && worst<=1e-8 && (ghosts==0 || !before.passed)
            && (!explicit_target || (solver.get_graph_build_count()==1 && solver.get_numeric_update_count()==int(round))),detail.str());
    }
}

// Force an inaccurate initial stop, then exercise the production GPU correction
// and owner-to-ghost publication with an independently known diagonal solution.
void pressure_refinement_gate(MPI_Comm comm,Report& report) {
    using mars::segregated::runtime::HypreSimpleSolve;
    int rank=0,ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    setenv("MARS_HYPRE_MINITER","0",1);
    for (bool cached:{false,true}) {
        Problem p(1,ranks,{}); Local l=extract(p,rank);
        fill(p,l,rank,1,[&](int g,int) { return 2*p.solution(g,0,11); });
        for (int i=0;i<l.nodes();++i) for (int k=l.offsets[i];k<l.offsets[i+1];++k)
            l.blocks[k]=l.columns[k]==i?2.:0.;
        Device<1,HYPRE_BigInt> d(l);
        OwnedRowSystem<1,Matrix,HYPRE_BigInt> system(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
        HypreSimpleSolve<1> solve(comm);
        if (cached) solve.solver.enable_fixed_graph_updates();
        solve.solver.set_stopping_tolerances(1e-12,1e6);
        const Tolerance target{0,1e-10,true};
        solve.solver.enable_true_residual_check(target.absolute,target.relative,true);
        system.update(d.view(),solve.rhs(system.rows()),system.rows());
        const bool initial=solve(system);
        const auto controls=solve.solver.prepared_controls();
        const int setups=solve.solver.get_setup_count(),graphs=solve.solver.get_graph_build_count();
        Buffer<double> local(l.nodes(),std::numeric_limits<double>::quiet_NaN());
        GhostExchange exchange(p,l,rank);
        auto publish=[&](const auto& values) {
            system.unpack(values.data(),values.size(),raw(local),local.size());
            exchange.run(raw(local),comm);
            return halo_complete(raw(local),local.size());
        };
        const auto result=solve.refine(system,publish,target);
        const auto norms=system.residual(halo_complete(raw(local),local.size()),solve.rhs(),target);
        const auto values=download(raw(local),local.size());
        double error=0;
        for (int i=0;i<l.nodes();++i) error=std::max(error,std::abs(values[i]-p.solution(l.global[i],0,11)));
        const auto restored=solve.solver.prepared_controls();
        const bool ok=!initial && result.accepted && result.rounds>0 && result.rounds<=3
            && result.iterations>0 && result.iterations<=solve.solver.get_max_iterations()
            && norms.passed && error<1e-10 && solve.solver.get_setup_count()==setups
            && solve.solver.get_graph_build_count()==graphs
            && restored.relative==controls.relative && restored.absolute==controls.absolute
            && restored.minimum==controls.minimum && restored.maximum==controls.maximum;
        report.result("GPU pressure correction: original target, solution/ghosts, setup reuse, restored controls; cache="+std::to_string(cached),all_true(ok,comm));
    }
    unsetenv("MARS_HYPRE_MINITER");
}

void pressure_expansion_gate(MPI_Comm comm,Report& report) {
    using mars::segregated::runtime::HypreSimpleSolve;
    int rank=0,ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    setenv("MARS_HYPRE_MINITER","0",1);
    for(bool cached:{false,true}) {
        Problem p(1,ranks,{}); Local l=extract(p,rank);
        const double rhs=std::nextafter(3.,INFINITY);
        fill(p,l,rank,1,[&](int,int) { return rhs; });
        for(int i=0;i<l.nodes();++i) for(int k=l.offsets[i];k<l.offsets[i+1];++k) l.blocks[k]=l.columns[k]==i?3.:0.;
        Device<1,HYPRE_BigInt> d(l);
        OwnedRowSystem<1,Matrix,HYPRE_BigInt> system(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
        HypreSimpleSolve<1> solve(comm);
        solve.solver.enable_gpu_aware_mpi();
        if(cached) solve.solver.enable_fixed_graph_updates();
        solve.solver.set_stopping_tolerances(1e-20,1e6);
        const Tolerance target{1e-17,1e-20,true}; solve.solver.enable_true_residual_check(1e-17,1e-20,true);
        system.update(d.view(),solve.rhs(system.rows()),system.rows());
        const bool initial=solve(system);
        // Force an initial rejection without depending on a library-specific stagnation path.
        solve.solver.set_prepared_controls({1e-20,1e-17,0,2000});
        const auto controls=solve.solver.prepared_controls();
        const int setups=solve.solver.get_setup_count(),graphs=solve.solver.get_graph_build_count();
        const auto recovered=solve.expand();
        Buffer<double> high(l.nodes(),NAN),low(l.nodes(),NAN),defect(system.rows());
        system.unpack(solve.solution(),solve.size(),raw(high),high.size());
        system.unpack(solve.expansion_low.data(),solve.expansion_low.size(),raw(low),low.size());
        GhostExchange halo(p,l,rank); halo.run(raw(high),comm); halo.run(raw(low),comm);
        const auto norms=system.expansion_residual(halo_complete(raw(high),high.size()),halo_complete(raw(low),low.size()),
            solve.rhs(),raw(defect),defect.size(),target);
        const auto h=download(raw(high),high.size()),lo=download(raw(low),low.size());
        const auto ordinary=system.residual(halo_complete(raw(high),high.size()),solve.rhs(),target);
        bool retained=true;
        for(int i=0;i<l.nodes();++i) retained=retained && std::isfinite(h[i]) && std::isfinite(lo[i])
            && lo[i]!=0 && std::abs(h[i]-1.)<1e-14;
        const auto restored=solve.solver.prepared_controls();
        const bool ok=!initial && recovered.accepted && recovered.rounds>0 && recovered.rounds<=4
            && recovered.iterations>0 && recovered.iterations<=4*controls.maximum && norms.passed && !ordinary.passed && retained
            && setups==solve.solver.get_setup_count() && graphs==solve.solver.get_graph_build_count()
            && restored.relative==controls.relative && restored.absolute==controls.absolute
            && restored.minimum==controls.minimum && restored.maximum==controls.maximum;
        report.result("GPU retained pressure: original target, halo pair, setup reuse and restored controls; cache="+std::to_string(cached),all_true(ok,comm));
    }
    unsetenv("MARS_HYPRE_MINITER");
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
            // Exercise expansion initialization before any other Hypre solve.
            pressure_expansion_gate(MPI_COMM_WORLD,report);
            run_gates<1,Matrix,HYPRE_BigInt>(MPI_COMM_WORLD,report);
            run_gates<3,Matrix,HYPRE_BigInt>(MPI_COMM_WORLD,report);
            hypre_gates<1>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"uneven");
            hypre_gates<3>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"uneven");
            Options small; small.nx=9; small.ny=3; small.nz=3;
            hypre_gates<1>(MPI_COMM_WORLD,report,small,EmptyRanks::reject,"81-row l1 coarse relaxation",18);
            hypre_gates<3>(MPI_COMM_WORLD,report,small,EmptyRanks::reject,"l1 coarse block relaxation",18);
            hypre_gates<1>(MPI_COMM_WORLD,report,small,EmptyRanks::reject,"explicit pressure target, cached refresh",18,true);
            hypre_gates<1>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"GPU pressure profile one level",18,true,1);
            hypre_gates<1>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"GPU pressure profile multilevel",18,true,2);
            hypre_gates<1>(MPI_COMM_WORLD,report,{},EmptyRanks::reject,"GPU pressure profile FlexGMRES",18,true,3);
            pressure_refinement_gate(MPI_COMM_WORLD,report);
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
