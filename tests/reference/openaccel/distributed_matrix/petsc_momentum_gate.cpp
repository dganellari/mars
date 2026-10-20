#include "gate_common.hpp"
#include "mars_segregated_petsc_momentum.hpp"
#include <sstream>
#ifndef MARS_REPLAY_CUDA
#include "../distributed_simple/channel.hpp"
#include "../simple_duct/duct_host_solver.hpp"
#endif
using namespace dmatrix_gate;
using mars::segregated::runtime::PetscSimpleSolve;

struct PetscGateCsr {
    Buffer<int> offsets,columns; Buffer<double> values;
    void allocate(int rows,int,int nnz) { offsets.resize(rows+1); columns.resize(nnz); values.resize(nnz); }
    int* rowOffsetsPtr() { return raw(offsets); } const int* rowOffsetsPtr() const { return raw(offsets); }
    int* colIndicesPtr() { return raw(columns); } const int* colIndicesPtr() const { return raw(columns); }
    double* valuesPtr() { return raw(values); } const double* valuesPtr() const { return raw(values); }
};

#ifndef MARS_REPLAY_CUDA
void mixed_simple_gate(MPI_Comm comm,Report& report) {
    using namespace mars::segregated::runtime;
    int rank=0,ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    const auto mesh=dsimple_gate::channel(4,2,2);
    const auto part=dsimple_gate::build(mesh,dsimple_gate::partition(mesh,ranks),rank,-1,-1,comm);
    DistributedSimpleRunner<PetscGateCsr,long long,duct::GatheredSolve,PetscSimpleSolve> run(comm,part.input,part.ownership);
    DistributedSimpleRunner<PetscGateCsr,long long,duct::GatheredSolve> reference(comm,part.input,part.ownership);
    run.set_momentum_tolerances(true,1e-12,0);
    run.momentum_solve.solver.configure({{"petsc_fgmres_bjacobi",1},{"rtol",1e-12},{"atol",0},{"maxiter",200}});
    run.momentum_solve.solver.enable_fixed_graph_updates();
    run.poison_unexchanged=true; reference.poison_unexchanged=true;
    for (int step=0;step<2;++step) { run.step(); reference.step(); }
    double error=0;
    for (std::size_t i=0;i<run.velocity.values.size();++i)
        error=std::max(error,std::abs(run.velocity.values[i]-reference.velocity.values[i]));
    for (std::size_t i=0;i<run.pressure.values.size();++i)
        error=std::max(error,std::abs(run.pressure.values[i]-reference.pressure.values[i]));
    report.result("SIMPLE with PETSc momentum and an independent test pressure solver: two steps",
        all_true(error<1e-8 && run.completed==2 && run.momentum_solve.solver.get_numeric_update_count()==2,comm));
}
#endif

void petsc_gates(MPI_Comm comm,Report& report) {
    int rank=0,ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    Problem p(3,ranks,{}); Local l=extract(p,rank);
    fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,11); });
    Device<3,long long> d(l);
    System<3,PetscGateCsr,long long> system(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
    GhostExchange exchange(p,l,rank);
    for (bool cached:{false,true}) {
        PetscSimpleSolve<3> solve(comm);
        solve.set_tolerances(1e-10,0);
        solve.solver.configure({{"petsc_fgmres_bjacobi",1},{"rtol",1e-10},{"atol",0},{"maxiter",100}});
        solve.solver.enable_fixed_graph_updates(cached);
        for (int round=1;round<=3;++round) {
            fill(p,l,rank,round,[&](int g,int c) { return p.product(g,c,round,10+round); });
            const double scale=round==1?1e-9:round==2?1.:1e9;
            for (auto& value:l.blocks) value*=scale;
            for (auto& value:l.rhs) value*=scale;
            overwrite(d.blocks,l.blocks); overwrite(d.rhs,l.rhs);
            system.update(d.view(),solve.rhs(system.rows()),system.rows());
            const bool accepted=solve(system);
            Buffer<double> local(3*l.nodes(),std::numeric_limits<double>::quiet_NaN());
            system.unpack(solve.solution(),solve.size(),raw(local),local.size());
            const auto stale=system.residual(halo_complete(raw(local),local.size()),solve.rhs());
            exchange.run(raw(local),comm);
            const auto residual=system.residual(halo_complete(raw(local),local.size()),solve.rhs(),{0,1e-10,true});
            const auto host=download(raw(local),local.size());
            double error=0;
            for (int i=0;i<l.nodes();++i) for (int c=0;c<3;++c)
                error=std::max(error,std::abs(host[3*i+c]-p.solution(l.global[i],c,10+round)));
            std::ostringstream settings; solve.solver.write_settings(settings);
            std::istringstream input(settings.str()); std::map<std::string,double> actual; std::string key; double value;
            while (input>>key>>value) actual[key]=value;
            report.result("PETSc FGMRES/BJacobi coupled solve, scales, ghosts and numeric refresh; cache="+std::to_string(cached),
                all_true(accepted && residual.passed && error<1e-8 && (ranks==1 || !stale.passed)
                    && solve.solver.get_graph_build_count()==(cached?1:round) && solve.solver.get_numeric_update_count()==round
                    && solve.solver.get_setup_count()==round && actual.at("rtol")==1e-10 && actual.at("maxiter")==100
                    && actual.at("local_blocks")==1 && actual.at("sub_preonly")==1 && actual.at("sub_ilu_levels")==0
                    && actual.at("right_preconditioning")==1 && actual.at("unpreconditioned_norm")==1,comm));
        }
        std::fill(l.rhs.begin(),l.rhs.end(),0); overwrite(d.rhs,l.rhs);
        system.update(d.view(),solve.rhs(system.rows()),system.rows());
        const bool accepted=solve(system);
        report.result("PETSc zero RHS resets the previous solution",all_true(accepted && solve.solver.residual_norm()==0,comm));
    }
    fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,11); });
    overwrite(d.blocks,l.blocks); overwrite(d.rhs,l.rhs);
    PetscSimpleSolve<3> capped(comm);
    capped.set_tolerances(1e-14,0);
    capped.solver.configure({{"petsc_fgmres_bjacobi",1},{"rtol",1e-14},{"atol",0},{"maxiter",1}});
    system.update(d.view(),capped.rhs(system.rows()),system.rows());
    report.result("PETSc rejects an inaccurate capped solve",all_true(!capped(system),comm));
    PetscSimpleSolve<3> invalid(comm); invalid.set_tolerances(1e-10,0);
    bool rejected=false;
    try { invalid.solver.configure({{"petsc_fgmres_bjacobi",1},{"rtol",1e-10},{"atol",0},{"maxiter",rank==ranks-1?0:100}}); }
    catch (const std::exception&) { rejected=true; }
    report.result("PETSc rejects an invalid profile collectively",all_true(rejected,comm));
    if (ranks>1) {
        rejected=false;
        try { invalid.solver.configure({{"petsc_fgmres_bjacobi",1},{"rtol",1e-10},{"atol",0},{"maxiter",rank==ranks-1?101:100}}); }
        catch (const std::exception&) { rejected=true; }
        report.result("PETSc rejects differing valid controls collectively",all_true(rejected,comm));
    }
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
#ifdef MARS_REPLAY_CUDA
    MPI_Comm node; int rank=0,local=0,devices=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
    MPI_Comm_rank(node,&local); MPI_Comm_free(&node);
    if (cudaGetDeviceCount(&devices)!=cudaSuccess || devices<1 || cudaSetDevice(local%devices)!=cudaSuccess) MPI_Abort(MPI_COMM_WORLD,1);
#endif
    char name[]="petsc_momentum_gate",skip[]="-skip_petscrc";
    char* args[]={name,skip,nullptr}; char** argv_petsc=args; int argc_petsc=2;
    if (PetscInitialize(&argc_petsc,&argv_petsc,nullptr,nullptr) || PetscOptionsClear(nullptr)
        || PetscPushErrorHandler(PetscReturnErrorHandler,nullptr)) MPI_Abort(MPI_COMM_WORLD,1);
    Report report; MPI_Comm_rank(MPI_COMM_WORLD,&report.rank);
    try {
        petsc_gates(MPI_COMM_WORLD,report);
#ifndef MARS_REPLAY_CUDA
        mixed_simple_gate(MPI_COMM_WORLD,report);
#endif
    }
    catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    const int failed=report.failures;
    if (!report.rank) std::cout<<(failed?"FAIL":"PASS")<<": "<<report.passes<<" passed, "<<failed<<" failed\n";
    PetscFinalize(); MPI_Finalize(); return failed?1:0;
}
