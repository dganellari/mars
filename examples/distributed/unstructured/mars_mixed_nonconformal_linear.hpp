#pragma once
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#if defined(__CUDACC__) || defined(__HIPCC__)
#include <_hypre_utilities.h>
#include <thrust/equal.h>
#endif
namespace mixedflow {
inline void hypre_check(int code) {check(code==0,"mixed Hypre API");}
class Linear {
public:
    Linear(Graph& graph):graph_(graph) {
        auto rows=host<int64_t>(graph.rowIds),cols=host<int64_t>(graph.columns);auto off=host<int>(graph.offsets);
        check(!rows.empty(),"mixed experiment needs an owned row on every rank");
        std::vector<HYPRE_BigInt> rr(rows.begin(),rows.end()),cc(cols.begin(),cols.end());
        std::vector<HYPRE_Int> counts;for(size_t i=0;i<rows.size();++i) counts.push_back(off[i+1]-off[i]);
        rows_=Buffer<HYPRE_BigInt>(rr.begin(),rr.end());cols_=Buffer<HYPRE_BigInt>(cc.begin(),cc.end());counts_=Buffer<HYPRE_Int>(counts.begin(),counts.end());
        hypre_check(HYPRE_IJMatrixCreate(MPI_COMM_WORLD,rr.front(),rr.back(),rr.front(),rr.back(),&a_));
        hypre_check(HYPRE_IJMatrixSetObjectType(a_,HYPRE_PARCSR));
        for(HYPRE_IJVector* vector:{&b_,&x_}) {hypre_check(HYPRE_IJVectorCreate(MPI_COMM_WORLD,rr.front(),rr.back(),vector));hypre_check(HYPRE_IJVectorSetObjectType(*vector,HYPRE_PARCSR));}
    }
    ~Linear() {HYPRE_IJMatrixDestroy(a_);HYPRE_IJVectorDestroy(b_);HYPRE_IJVectorDestroy(x_);}
    int solve(ChannelOps& ops,bool momentum) {
        ncflow::finish();ops.zero(ops.solution());ncflow::finish();
        hypre_check(HYPRE_IJMatrixInitialize(a_));hypre_check(HYPRE_IJMatrixSetValues(a_,ops.n,pointer(counts_),pointer(rows_),pointer(cols_),pointer(graph_.values)));hypre_check(HYPRE_IJMatrixAssemble(a_));
        hypre_check(HYPRE_IJVectorInitialize(b_));hypre_check(HYPRE_IJVectorInitialize(x_));
        hypre_check(HYPRE_IJVectorSetValues(b_,ops.n,pointer(rows_),ops.rhs()));hypre_check(HYPRE_IJVectorSetValues(x_,ops.n,pointer(rows_),ops.solution()));
        hypre_check(HYPRE_IJVectorAssemble(b_));hypre_check(HYPRE_IJVectorAssemble(x_));
        HYPRE_ParCSRMatrix a;HYPRE_ParVector b,x;hypre_check(HYPRE_IJMatrixGetObject(a_,reinterpret_cast<void**>(&a)));
        hypre_check(HYPRE_IJVectorGetObject(b_,reinterpret_cast<void**>(&b)));hypre_check(HYPRE_IJVectorGetObject(x_,reinterpret_cast<void**>(&x)));
        HYPRE_Solver solver,amg;hypre_check(HYPRE_BoomerAMGCreate(&amg));
        hypre_check(HYPRE_BoomerAMGSetCoarsenType(amg,8));hypre_check(HYPRE_BoomerAMGSetInterpType(amg,6));
        hypre_check(HYPRE_BoomerAMGSetRelaxType(amg,18));hypre_check(HYPRE_BoomerAMGSetStrongThreshold(amg,.5));
        hypre_check(HYPRE_BoomerAMGSetNumSweeps(amg,3));hypre_check(HYPRE_BoomerAMGSetMaxLevels(amg,25));
        hypre_check(HYPRE_BoomerAMGSetTol(amg,0));hypre_check(HYPRE_BoomerAMGSetMaxIter(amg,1));
        if(momentum) hypre_check(HYPRE_BoomerAMGSetNumFunctions(amg,3));
        hypre_check(HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD,&solver));hypre_check(HYPRE_ParCSRGMRESSetKDim(solver,50));
        hypre_check(HYPRE_ParCSRGMRESSetMaxIter(solver,500));hypre_check(HYPRE_ParCSRGMRESSetTol(solver,1e-9));
        hypre_check(HYPRE_ParCSRGMRESSetAbsoluteTol(solver,1e-14));
        hypre_check(HYPRE_ParCSRGMRESSetPrecond(solver,HYPRE_BoomerAMGSolve,HYPRE_BoomerAMGSetup,amg));
        hypre_check(HYPRE_ParCSRGMRESSetup(solver,a,b,x));hypre_check(HYPRE_ParCSRGMRESSolve(solver,a,b,x));
        HYPRE_Int iterations;hypre_check(HYPRE_ParCSRGMRESGetNumIterations(solver,&iterations));
        hypre_check(HYPRE_IJVectorGetValues(x_,ops.n,pointer(rows_),ops.solution()));
#if defined(__CUDACC__) || defined(__HIPCC__)
        // Hypre may use a different stream; its completed correction feeds MARS.
        hypre_check(hypre_ForceSyncComputeStream());
#endif
        ops.apply(ops.solution(),ops.work());ops.axpy(-1,ops.rhs(),ops.work());
        const double residual=ops.norm(ops.work()),scale=ops.norm(ops.rhs());
        check(std::isfinite(residual) && residual<=std::max(1e-13,2e-9*scale),"mixed true linear residual");
        hypre_check(HYPRE_ParCSRGMRESDestroy(solver));hypre_check(HYPRE_BoomerAMGDestroy(amg));return iterations;
    }
private:
    Graph& graph_;Buffer<HYPRE_BigInt> rows_,cols_;Buffer<HYPRE_Int> counts_;HYPRE_IJMatrix a_{};HYPRE_IJVector b_{},x_{};
};
inline bool Channel::solve(int budget,double tolerance,const std::string& output) {
    ChannelOps u(momentum_,input_.mesh,50,3,false),p(pressure_,input_.mesh,50,1,false);Linear velocity(momentum_),pressure(pressure_);
    Buffer<int> previous_closed(closed_.size());
    if(rankId==0) std::puts("Linear solves: Hypre GMRES/PMIS-AMG for both systems, true residual checked; differs from reference PETSc momentum and Falgout pressure.");
    for(int iteration=0;iteration<budget;++iteration) {
        view_.completed=iteration;exchange();
        each<device>(n_,Gradient{view_,0,iteration>0});each<device>(n_,Gradient{view_,1,iteration>0});exchange();
        each<device>(n_,Blend{view_});exchange();assemble(true,u);
        auto mu=reduce(n_,Statistics{view_,u.rhs(),pointer(u.diagonal),3});
        double velocity_scale=mu.maxima[1]<1e-9?1.:mu.maxima[1];double urms=0;
        for(int j=0;j<3;++j) urms=std::max(urms,sqrt(mu.sums[7+j]/input_.mesh.globalNodes)/velocity_scale);
        exchange();int ui=velocity.solve(u,true);each<device>(n_,AddIncrement{view_,u.solution(),true});exchange();
        assemble(false,p);
        auto pr=reduce(n_,Statistics{view_,p.rhs(),pointer(p.diagonal),1});double pressure_scale=std::max(pr.maxima[0]+pr.maxima[3],.5*view_.controls.rho*velocity_scale*velocity_scale);
        double prms=sqrt(pr.sums[7]/input_.mesh.globalNodes)/pressure_scale;
        int pi=pressure.solve(p,false);each<device>(n_,AddIncrement{view_,p.solution(),false});exchange();
        each<device>(n_,Gradient{view_,2,false});each<device>(int(cells_.size()),UpdateCellMass{view_});
        previous_closed=closed_;
        each<device>(int(faces_.size()),UpdateBoundaryMass{view_});update_interface_mass();
        each<device>(n_,Divergence{view_});each<device>(n_,Correct{view_});ncflow::finish();
        const auto balance=reduce(n_,Statistics{view_});
        const double mass_scale=fmax(fabs(balance.sums[0])+fabs(balance.sums[1])+fabs(balance.sums[2]),1e-30);
        const double continuity_relative=balance.sums[5]/mass_scale;
        const double boundary_relative=fabs(balance.sums[0]+balance.sums[1]+balance.sums[2])/mass_scale;
#if defined(__CUDACC__) || defined(__HIPCC__)
        int local_stable=thrust::equal(previous_closed.begin(),previous_closed.end(),closed_.begin());
#else
        int local_stable=std::equal(previous_closed.begin(),previous_closed.end(),closed_.begin());
#endif
        int boundary_stable=0;MPI_Allreduce(&local_stable,&boundary_stable,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
        // A changed outlet state needs another assembly with that boundary condition.
        bool converged=boundary_stable && urms<tolerance && prms<tolerance && continuity_relative<tolerance &&
                       (!view_.conservative || boundary_relative<tolerance);
        if(iteration%10==0 || converged || iteration+1==budget) {
            auto s=reduce(n_,Statistics{view_});check(s.sums[10]==0,"mixed solution has nonfinite values");
            if(rankId==0) {std::printf("iteration %d momentum_rms %.9e pressure_rms %.9e linear %d/%d inlet %.9e outlet %.9e continuity_l1 %.9e interface_in %.9e interface_out %.9e closed_area %.9e\n",iteration+1,urms,prms,ui,pi,s.sums[0],s.sums[1],s.sums[5],s.sums[11],s.sums[12],s.sums[13]);std::fflush(stdout);}
        }
        if(converged) {
            auto s=reduce(n_,Statistics{view_});
            const double boundary=fabs(s.sums[0]+s.sums[1]+s.sums[2]);
            const double scale=fmax(fabs(s.sums[0])+fabs(s.sums[1])+fabs(s.sums[2]),1e-30);
            const double interface=fabs(s.sums[11]+s.sums[12]);
            if(rankId==0) std::printf("Conservation: inlet %.17g outlet %.17g wall %.17g boundary_relative %.17g interface_net %.17g continuity_l1 %.17g\n",s.sums[0],s.sums[1],s.sums[2],boundary/scale,interface,s.sums[5]);
            write(output);return true;
        }
    }
    write(output);return false;
}
}
