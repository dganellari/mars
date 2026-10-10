#pragma once
#include <memory>
#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_compensated_dot.hpp"
#ifdef __CUDACC__
// The C header defines the stream macro; this header declares its GPU accessor.
#include <_hypre_utilities.hpp>
#endif

namespace recovery {
using mars::segregated::CompensatedDot;
#ifdef __CUDACC__
#define REPLAY_HD __host__ __device__
template<class F> __global__ void kernel(int n,F f) { const int i=int(blockIdx.x*blockDim.x+threadIdx.x); if(i<n) f(i); }
template<class F> void each(int n,F f) {
    if(n) kernel<<<(n+255)/256,256,0,hypre_HandleComputeStream(hypre_handle())>>>(n,f);
    frozen::require(cudaGetLastError()==cudaSuccess);
}
#else
#define REPLAY_HD
template<class F> void each(int n,F f) { for(int i=0;i<n;++i) f(i); }
#endif
inline double* data(hypre_ParVector* v) { return hypre_VectorData(hypre_ParVectorLocalVector(v)); }
struct Vector {
    hypre_ParVector* p;
    Vector(hypre_ParVector* like,HYPRE_MemoryLocation memory):p(hypre_ParVectorCloneDeep_v2(like,memory)) {
        frozen::require(p); checked(HYPRE_GetError());
    }
    ~Vector() { hypre_ParVectorDestroy(p); }
    Vector(const Vector&)=delete;
};
struct Gather {
    const int* map; const double* x; double* send;
    REPLAY_HD void operator()(int i) const { send[i]=x[map[i]]; }
};
struct Add {
    const double *x,*delta; double* trial;
    REPLAY_HD void operator()(int i) const { trial[i]=CompensatedDot::add(x[i],delta[i]); }
};
struct AddWithRemainder {
    const double *x,*delta; double *trial,*low;
    REPLAY_HD void operator()(int i) const {
        const double sum=CompensatedDot::add(x[i],delta[i]);
        const double z=CompensatedDot::add(sum,-x[i]);
        trial[i]=sum;
        // TwoSum retains the part discarded by the ordinary solution update.
        low[i]=CompensatedDot::add(CompensatedDot::add(x[i],-CompensatedDot::add(sum,-z)),
                                  CompensatedDot::add(delta[i],-z));
    }
};
struct Defect {
    const int *di,*dj,*oi,*oj;
    const double *da,*oa,*x,*ghost,*b;
    double *r,*error;
    const double *low=nullptr,*ghost_low=nullptr;
    REPLAY_HD void operator()(int row) const {
        CompensatedDot dot;
        for(int k=di[row];k<di[row+1];++k) {
            dot.product(da[k],x[dj[k]]);
            if(low) dot.product(da[k],low[dj[k]]);
        }
        for(int k=oi[row];k<oi[row+1];++k) {
            dot.product(oa[k],ghost[oj[k]]);
            if(low) dot.product(oa[k],ghost_low[oj[k]]);
        }
        dot.product(-1.,b[row]);
        r[row]=-dot.value(); error[row]=fmax(dot.error_bound(),std::sqrt(std::numeric_limits<double>::min()));
        if (r[row]!=0 && r[row]*r[row]<std::numeric_limits<double>::min())
            r[row]=std::numeric_limits<double>::quiet_NaN();
    }
};
struct Norm {
    double lower,upper;
    bool finite() const { return std::isfinite(lower) && std::isfinite(upper); }
};
class Residual {
    hypre_ParCSRMatrix* a_;
    hypre_ParCSRCommPkg* pkg_;
    HYPRE_MemoryLocation memory_;
    double *send_,*ghost_,*ghost_low_=nullptr;
    Vector error_;
    double margin_;
public:
    Residual(hypre_ParCSRMatrix* a,hypre_ParVector* x,HYPRE_MemoryLocation memory)
        :a_(a),memory_(memory),error_(x,memory) {
        if(!hypre_ParCSRMatrixCommPkg(a)) checked(hypre_MatvecCommPkgCreate(a));
        pkg_=hypre_ParCSRMatrixCommPkg(a);
        frozen::require(pkg_);
#ifdef __CUDACC__
        hypre_ParCSRCommPkgCopySendMapElmtsToDevice(pkg_); checked(HYPRE_GetError());
        // This diagnostic requires direct device MPI; no hidden host halo staging.
        frozen::require(hypre_ParCSRCommPkgNumSends(pkg_)+hypre_ParCSRCommPkgNumRecvs(pkg_)==0 || hypre_GetGpuAwareMPI()!=0);
#endif
        const int sends=hypre_ParCSRCommPkgSendMapStart(pkg_,hypre_ParCSRCommPkgNumSends(pkg_));
        const int receives=hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(a));
        frozen::require(receives==hypre_ParCSRCommPkgRecvVecStart(pkg_,hypre_ParCSRCommPkgNumRecvs(pkg_)));
        send_=hypre_TAlloc(double,std::max(1,sends),memory_);
        ghost_=hypre_TAlloc(double,std::max(1,receives),memory_);
        frozen::require(send_ && ghost_);
        // Positive FP64 norm reductions contain at most O(global rows) rounded terms.
        margin_=(8.*double(hypre_ParCSRMatrixGlobalNumRows(a))+128.)*std::numeric_limits<double>::epsilon();
        frozen::require(margin_>0 && margin_<.01);
    }
    ~Residual() { hypre_TFree(send_,memory_); hypre_TFree(ghost_,memory_); hypre_TFree(ghost_low_,memory_); }
    Norm norm(hypre_ParVector* x) const {
        const double square=hypre_ParVectorInnerProd(x,x); checked(HYPRE_GetError());
        if(!(square>=0) || !std::isfinite(square)) return {NAN,INFINITY};
        return {std::nextafter(std::sqrt(square/(1+margin_)),0.),
                std::nextafter(std::sqrt(square/(1-margin_)),INFINITY)};
    }
private:
    void exchange(hypre_ParVector* x,double* ghost) {
        const int n=hypre_ParCSRCommPkgSendMapStart(pkg_,hypre_ParCSRCommPkgNumSends(pkg_));
#ifdef __CUDACC__
        const int* map=hypre_ParCSRCommPkgDeviceSendMapElmts(pkg_);
#else
        const int* map=hypre_ParCSRCommPkgSendMapElmts(pkg_);
#endif
        each(n,Gather{map,data(x),send_});
#ifdef __CUDACC__
        checked(hypre_ForceSyncComputeStream());
#endif
        auto* exchange=hypre_ParCSRCommHandleCreate_v2(1,pkg_,memory_,send_,memory_,ghost);
        frozen::require(exchange); checked(hypre_ParCSRCommHandleDestroy(exchange));
    }
public:
    Norm evaluate(hypre_ParVector* x,hypre_ParVector* b,hypre_ParVector* r,hypre_ParVector* low=nullptr) {
        exchange(x,ghost_);
        if(low) {
            if(!ghost_low_) {
                const int receives=hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(a_));
                ghost_low_=hypre_TAlloc(double,std::max(1,receives),memory_);
                frozen::require(ghost_low_);
            }
            exchange(low,ghost_low_);
        }
        auto* d=hypre_ParCSRMatrixDiag(a_); auto* o=hypre_ParCSRMatrixOffd(a_);
        each(hypre_CSRMatrixNumRows(d),Defect{hypre_CSRMatrixI(d),hypre_CSRMatrixJ(d),
            hypre_CSRMatrixI(o),hypre_CSRMatrixJ(o),hypre_CSRMatrixData(d),hypre_CSRMatrixData(o),
            data(x),ghost_,data(b),data(r),data(error_.p),low?data(low):nullptr,ghost_low_});
        const auto residual=norm(r),error=norm(error_.p);
        if(!residual.finite() || !error.finite()) return {NAN,INFINITY};
        return {std::max(0.,std::nextafter(residual.lower-error.upper,0.)),
                std::nextafter(residual.upper+error.upper,INFINITY)};
    }
};
struct Audit {
    Vector low,work,zero;
    explicit Audit(hypre_ParVector* x,HYPRE_MemoryLocation memory):low(x,memory),work(x,memory),zero(x,memory) {
        checked(hypre_ParVectorSetConstantValues(zero.p,0.));
    }
};
struct SolveResult {
    int error=0,global=0,iterations=0,converged=0;
    double relative=0;
};
inline SolveResult solve(HYPRE_Solver solver,bool flex,HYPRE_ParCSRMatrix a,HYPRE_ParVector b,HYPRE_ParVector x) {
    SolveResult r;
    r.error=(flex?HYPRE_ParCSRFlexGMRESSolve:HYPRE_ParCSRGMRESSolve)(solver,a,b,x);
    r.global=HYPRE_GetError(); HYPRE_ClearAllErrors();
    checked((flex?HYPRE_FlexGMRESGetNumIterations:HYPRE_GMRESGetNumIterations)(solver,&r.iterations));
    checked((flex?HYPRE_FlexGMRESGetConverged:HYPRE_GMRESGetConverged)(solver,&r.converged));
    checked((flex?HYPRE_FlexGMRESGetFinalRelativeResidualNorm:HYPRE_GMRESGetFinalRelativeResidualNorm)(solver,&r.relative));
    return r;
}
inline int maximum(int value,hypre_MPI_Comm comm) {
    int global=0; checked(hypre_MPI_Allreduce(&value,&global,1,HYPRE_MPI_INT,hypre_MPI_MAX,comm)); return global;
}
// Stop codes: target, budget, no certified progress, nonfinite, fatal backend error.
struct Result { int rounds=0,iterations=0,stop=1,audit_steps=0; };
inline Result run(hypre_ParCSRMatrix* a,hypre_ParVector* b,hypre_ParVector* x,
                  HYPRE_Solver solver,bool flex,const settings::Values& controls,HYPRE_MemoryLocation memory,
                  int rounds,const SolveResult& initial,std::ostream& trace,std::ostream* audit_trace=nullptr) {
    const auto comm=hypre_ParCSRMatrixComm(a);
    Result result;
    if(audit_trace) *audit_trace<<"mars-pressure-correction-audit-v1\n"<<std::setprecision(17);
    if(maximum((initial.error|initial.global)&~HYPRE_ERROR_CONV,comm)) { result.stop=4; return result; }
    Residual evaluator(a,x,memory);
    Vector defect(x,memory),delta(x,memory),trial(x,memory);
    std::unique_ptr<Audit> audit;
    if(audit_trace) audit=std::make_unique<Audit>(x,memory);
    const auto rhs=evaluator.norm(b);
    const double limit=std::max(controls.at("atol"),std::nextafter(controls.at("rtol")*rhs.lower,0.));
    const double limit_upper=std::max(controls.at("atol"),std::nextafter(controls.at("rtol")*rhs.upper,INFINITY));
    auto current=evaluator.evaluate(x,b,defect.p);
    trace<<std::setprecision(17)<<"initial "<<current.lower<<' '<<current.upper<<" limit "<<limit<<'\n';
    if(!rhs.finite() || !current.finite() || !std::isfinite(limit)) { result.stop=3; return result; }
    if(current.upper<=limit) { result.stop=0; return result; }
    checked((flex?HYPRE_ParCSRFlexGMRESSetTol:HYPRE_ParCSRGMRESSetTol)(solver,.1));
    checked((flex?HYPRE_ParCSRFlexGMRESSetAbsoluteTol:HYPRE_ParCSRGMRESSetAbsoluteTol)(solver,0.));
    checked((flex?HYPRE_ParCSRFlexGMRESSetMinIter:HYPRE_ParCSRGMRESSetMinIter)(solver,0));
    for(int round=0;round<rounds;++round) {
        const double defect_square=hypre_ParVectorInnerProd(defect.p,defect.p);
        checked(HYPRE_GetError());
        if(defect_square==0) { result.stop=2; break; }
        checked(hypre_ParVectorSetConstantValues(delta.p,0.));
        const auto correction=solve(solver,flex,reinterpret_cast<HYPRE_ParCSRMatrix>(a),
            reinterpret_cast<HYPRE_ParVector>(defect.p),reinterpret_cast<HYPRE_ParVector>(delta.p));
        ++result.rounds;
        const int iterations=maximum(correction.iterations,comm);
        result.iterations+=iterations;
        trace<<"correction "<<result.rounds<<" iterations "<<iterations<<" return "<<correction.error
             <<" global "<<correction.global<<" converged "<<correction.converged<<" reported "<<correction.relative<<'\n';
        if(maximum((correction.error|correction.global)&~HYPRE_ERROR_CONV,comm)) { result.stop=4; break; }
        frozen::require(iterations>=0 && iterations<=controls.at("maxiter"));
        if(iterations==0) { result.stop=2; break; }
        Norm correction_rhs{},correction_residual{},ideal{},rounding{};
        const int local_rows=hypre_CSRMatrixNumRows(hypre_ParCSRMatrixDiag(a));
        if(audit) {
            // Preserve the correction RHS until its independent check is complete.
            correction_rhs=evaluator.norm(defect.p);
            correction_residual=evaluator.evaluate(delta.p,defect.p,audit->work.p);
            each(local_rows,AddWithRemainder{data(x),data(delta.p),data(trial.p),data(audit->low.p)});
            ideal=evaluator.evaluate(trial.p,b,audit->work.p,audit->low.p);
            rounding=evaluator.evaluate(audit->low.p,audit->zero.p,audit->work.p);
        } else each(local_rows,Add{data(x),data(delta.p),data(trial.p)});
        const auto next=evaluator.evaluate(trial.p,b,defect.p);
        trace<<"candidate "<<next.lower<<' '<<next.upper<<'\n';
        if(audit) {
            ++result.audit_steps;
            *audit_trace<<"step "<<result.rounds<<" rhs "<<correction_rhs.lower<<' '<<correction_rhs.upper
                <<" correction "<<correction_residual.lower<<' '<<correction_residual.upper
                <<" ideal "<<ideal.lower<<' '<<ideal.upper<<" rounding "<<rounding.lower<<' '<<rounding.upper
                <<" candidate "<<next.lower<<' '<<next.upper<<" target "<<limit<<' '<<limit_upper<<'\n';
        }
        if(!next.finite()) { result.stop=3; break; }
        if(!(next.upper<current.lower)) { result.stop=2; break; }
        checked(hypre_ParVectorCopy(trial.p,x)); current=next;
        if(current.upper<=limit) { result.stop=0; break; }
    }
    checked((flex?HYPRE_ParCSRFlexGMRESSetTol:HYPRE_ParCSRGMRESSetTol)(solver,controls.at("rtol")));
    checked((flex?HYPRE_ParCSRFlexGMRESSetAbsoluteTol:HYPRE_ParCSRGMRESSetAbsoluteTol)(solver,controls.at("atol")));
    checked((flex?HYPRE_ParCSRFlexGMRESSetMinIter:HYPRE_ParCSRGMRESSetMinIter)(solver,int(controls.at("miniter"))));
    return result;
}
#undef REPLAY_HD
}
