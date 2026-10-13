#pragma once
#include "mars_segregated_distributed_matrix.hpp"
#include <petscksp.h>
#include <petscdevice.h>
#include <map>
#include <iostream>
#include <ostream>
#include <type_traits>

#if PETSC_VERSION_LT(3,24,0)
#error "SIMPLE PETSc momentum requires PETSc 3.24 or newer"
#endif

#if defined(MARS_REPLAY_CUDA) && !defined(PETSC_HAVE_CUDA)
#error "SIMPLE PETSc momentum requires a CUDA-enabled PETSc build"
#endif

namespace mars::segregated::runtime {

// The reference's two flat options leave all other KSP/PC choices at PETSc defaults.
// Keep this adapter restricted to that profile; do not approximate other PETSc PCs.
class PetscMomentumSolver {
public:
    using Values=std::map<std::string,double>;
    struct Timing { double prepare_seconds=0,packing_seconds=0,setup_seconds=0,solve_seconds=0,finish_seconds=0; };
    explicit PetscMomentumSolver(MPI_Comm comm):comm_(comm) {
        static_assert(std::is_same_v<PetscScalar,double>,"SIMPLE requires real FP64 PETSc");
        PetscBool initialized=PETSC_FALSE;
        check(PetscInitialized(&initialized));
        require(initialized,"PETSc is not initialized");
        PetscInt major=0,minor=0,patch=0;
        check(PetscGetVersionNumber(&major,&minor,&patch,nullptr));
        require(major==PETSC_VERSION_MAJOR && minor==PETSC_VERSION_MINOR && patch==PETSC_VERSION_SUBMINOR,
            "PETSc header and runtime versions differ");
    }
    PetscMomentumSolver(const PetscMomentumSolver&)=delete;
    PetscMomentumSolver& operator=(const PetscMomentumSolver&)=delete;
    ~PetscMomentumSolver() { destroy(); }
    void set_stopping_tolerances(double relative,double absolute) {
        require(relative>0 && relative<1 && absolute>=0 && std::isfinite(absolute),"invalid PETSc momentum target");
        relative_=relative; absolute_=absolute;
    }
    void configure(const Values& values) {
        require(values.size()==4 && values.count("petsc_fgmres_bjacobi") && values.count("rtol")
            && values.count("atol") && values.count("maxiter"),"unsupported PETSc momentum profile");
        require(values.at("petsc_fgmres_bjacobi")==1 && values.at("rtol")==relative_ && values.at("atol")==absolute_
            && values.at("maxiter")>=1 && values.at("maxiter")<=std::numeric_limits<int>::max()
            && values.at("maxiter")==std::floor(values.at("maxiter")),"invalid PETSc momentum profile");
        require(!matrix_,"PETSc momentum profile must precede setup");
        const double input[6]={relative_,-relative_,absolute_,-absolute_,values.at("maxiter"),-values.at("maxiter")};
        double global[6];
        if (MPI_Allreduce(input,global,6,MPI_DOUBLE,MPI_MAX,comm_)!=MPI_SUCCESS) MPI_Abort(comm_,1);
        require(global[0]==-global[1] && global[2]==-global[3] && global[4]==-global[5],"PETSc momentum controls differ between ranks");
        maximum_=int(values.at("maxiter")); configured_=true;
    }
    void enable_fixed_graph_updates(bool enabled=true) { cached_=enabled; }
    void enable_timing(bool enabled) { timing_enabled_=enabled; }
    const Timing& get_last_timing() const { return timing_; }
    int getLastIterations() const { return iterations_; }
    int get_graph_build_count() const { return builds_; }
    int get_numeric_update_count() const { return updates_; }
    int get_setup_count() const { return setups_; }
    double residual_norm() const { return residual_; }

    template<class System> bool solve(const System& system,double* rhs,double* solution) {
        require(configured_,"PETSc momentum requires an explicit profile");
        int faults=0; system.synchronize(faults); require(!faults,"PETSc momentum input synchronization failed");
        timing_={}; const double prepare=MPI_Wtime();
        if (!cached_) destroy();
        if (!matrix_) build(system);
        require(system.rows()==rows_ && system.nnz()==nnz_ && system.matrix().rowOffsetsPtr()==offsets_
            && system.matrix().colIndicesPtr()==columns_,"PETSc momentum graph changed");
        check(MatSetValuesCOO(matrix_,system.matrix().valuesPtr(),INSERT_VALUES));
        ++updates_;
        place(b_,rhs); place(x_,solution);
        check(VecSet(x_,0));
        check(KSPSetOperators(ksp_,matrix_,matrix_));
        check(KSPSetTolerances(ksp_,relative_,absolute_,PETSC_DEFAULT,maximum_));
        fence(); if (timing_enabled_) timing_.prepare_seconds=MPI_Wtime()-prepare;
        const double setup=MPI_Wtime();
        check(KSPSetUp(ksp_)); check(KSPSetUpOnBlocks(ksp_)); ++setups_;
        fence(); if (timing_enabled_) timing_.setup_seconds=MPI_Wtime()-setup;
        verify_settings();
        const double start=MPI_Wtime();
        check(KSPSolve(ksp_,b_,x_));
        fence(); if (timing_enabled_) timing_.solve_seconds=MPI_Wtime()-start;
        const double finish=MPI_Wtime();
        PetscInt iterations=0; KSPConvergedReason reason;
        check(KSPGetIterationNumber(ksp_,&iterations)); iterations_=int(iterations);
        check(KSPGetConvergedReason(ksp_,&reason));
        // A backend convergence flag cannot replace the original equation's residual.
        check(MatMult(matrix_,x_,residual_vector_)); check(VecAYPX(residual_vector_,-1,b_));
        double rhs_norm=0;
        check(VecNorm(residual_vector_,NORM_2,&residual_)); check(VecNorm(b_,NORM_2,&rhs_norm));
        fence(); reset(b_); reset(x_);
        if (timing_enabled_) timing_.finish_seconds=MPI_Wtime()-finish;
        const double limit=std::max(absolute_,relative_*rhs_norm);
        return reason>0 && std::isfinite(residual_) && std::isfinite(rhs_norm) && residual_<=limit;
    }

    void write_settings(std::ostream& out) const {
        require(ksp_ && setups_>0,"PETSc momentum settings require a prepared solver");
        PetscInt restart=0,maximum=0,levels=0,count=0; KSP* sub=nullptr;
        PetscReal relative=0,absolute=0,divergence=0;
        PC pc,sub_pc; KSPNormType norm; PCSide side;
        check(KSPGetPC(ksp_,&pc)); check(PCBJacobiGetSubKSP(pc,&count,nullptr,&sub));
        check(KSPGetPC(sub[0],&sub_pc)); check(PCFactorGetLevels(sub_pc,&levels));
        check(KSPGMRESGetRestart(ksp_,&restart)); check(KSPGetTolerances(ksp_,&relative,&absolute,&divergence,&maximum));
        check(KSPGetNormType(ksp_,&norm)); check(KSPGetPCSide(ksp_,&side));
        out.precision(17);
        out<<"petsc_fgmres_bjacobi 1\nrtol "<<relative<<"\natol "<<absolute<<"\nmaxiter "<<maximum
           <<"\npetsc_release "<<10000*PETSC_VERSION_MAJOR+100*PETSC_VERSION_MINOR+PETSC_VERSION_SUBMINOR
           <<"\nrestart "<<restart<<"\ndivergence_tolerance "<<divergence<<"\nlocal_blocks "<<count
           <<"\nsub_preonly 1\nsub_ilu_levels "<<levels<<"\nright_preconditioning "<<(side==PC_RIGHT)
           <<"\nunpreconditioned_norm "<<(norm==KSP_NORM_UNPRECONDITIONED)<<"\ndevice "<<device()<<'\n';
    }
private:
    MPI_Comm comm_;
    Mat matrix_=nullptr; Vec b_=nullptr,x_=nullptr,residual_vector_=nullptr; KSP ksp_=nullptr;
    const int *offsets_=nullptr,*columns_=nullptr;
    int rows_=0,nnz_=0,maximum_=20,iterations_=0,builds_=0,updates_=0,setups_=0;
    double relative_=1e-6,absolute_=1e-16,residual_=0;
    bool configured_=false,cached_=false,timing_enabled_=false;
    Timing timing_;
    static constexpr bool device() {
#ifdef MARS_REPLAY_CUDA
        return true;
#else
        return false;
#endif
    }
    void check(PetscErrorCode error) const {
        if (error) {
            std::cerr<<"ERROR: PETSc momentum API failed, code="<<error<<'\n';
            MPI_Abort(comm_,1); throw std::runtime_error("PETSc momentum API failed");
        }
    }
    void require(bool local,const char* message) const {
        int failed=!local,any=0;
        if (MPI_Allreduce(&failed,&any,1,MPI_INT,MPI_MAX,comm_)!=MPI_SUCCESS) MPI_Abort(comm_,1);
        if (any) throw std::runtime_error(message);
    }
    void fence() const {
#ifdef MARS_REPLAY_CUDA
        PetscDeviceContext context;
        check(PetscDeviceContextGetCurrentContext(&context)); check(PetscDeviceContextSynchronize(context));
#endif
    }
    void place(Vec vector,double* data) const {
#ifdef MARS_REPLAY_CUDA
        check(VecCUDAPlaceArray(vector,data));
#else
        check(VecPlaceArray(vector,data));
#endif
    }
    void reset(Vec vector) const {
#ifdef MARS_REPLAY_CUDA
        check(VecCUDAResetArray(vector));
#else
        check(VecResetArray(vector));
#endif
    }
    void destroy() {
        if (ksp_) check(KSPDestroy(&ksp_));
        if (residual_vector_) check(VecDestroy(&residual_vector_));
        if (x_) check(VecDestroy(&x_));
        if (b_) check(VecDestroy(&b_));
        if (matrix_) check(MatDestroy(&matrix_));
    }
    template<class GlobalId> struct CoordinateRows {
        const int* offsets; const int* columns; const GlobalId* map;
        PetscInt first; PetscInt *row,*column;
#if defined(__CUDACC__)
        __host__ __device__
#endif
        void operator()(int i) const {
            for (int j=offsets[i];j<offsets[i+1];++j) { row[j]=first+i; column[j]=PetscInt(map[columns[j]]); }
        }
    };
    template<class System> void build(const System& system) {
        const auto range=system.hypre_rows();
        require(range.end-range.begin==system.rows() && system.rows()>0 && range.column_begin==0
            && range.column_end<=std::numeric_limits<PetscInt>::max(),"invalid PETSc momentum ownership");
        rows_=system.rows(); nnz_=system.nnz();
        offsets_=system.matrix().rowOffsetsPtr(); columns_=system.matrix().colIndicesPtr();
        check(MatCreate(comm_,&matrix_)); check(MatSetSizes(matrix_,rows_,rows_,range.column_end,range.column_end));
        check(MatSetType(matrix_,device()?MATAIJCUSPARSE:MATAIJ));
        check(MatSetOption(matrix_,MAT_NO_OFF_PROC_ENTRIES,PETSC_TRUE));
        distributed::Buffer<PetscInt> row(nnz_),column(nnz_);
        const auto* map=distributed::raw(system.solver_dof_map());
        using GlobalId=std::remove_cv_t<std::remove_pointer_t<decltype(map)>>;
        int faults=0;
        distributed::apply(rows_,CoordinateRows<GlobalId>{offsets_,columns_,map,PetscInt(range.begin),
            distributed::raw(row),distributed::raw(column)},0,faults);
        system.synchronize(faults); require(!faults,"PETSc momentum graph conversion failed");
        // PETSc stages these indices on host internally during graph setup.
        // Numeric updates below use its device COO map without that staging.
        check(MatSetPreallocationCOO(matrix_,nnz_,distributed::raw(row),distributed::raw(column)));
        fence();
        PetscInt first=0,last=0; check(MatGetOwnershipRange(matrix_,&first,&last));
        require(first==range.begin && last==range.end,"PETSc momentum ownership differs from adapter");
        check(MatCreateVecs(matrix_,&x_,&b_)); check(VecDuplicate(b_,&residual_vector_));
        check(KSPCreate(comm_,&ksp_)); check(KSPSetType(ksp_,KSPFGMRES));
        PC pc; check(KSPGetPC(ksp_,&pc)); check(PCSetType(pc,PCBJACOBI));
        // No SetFromOptions: only the two supported reference choices enter this context.
        check(KSPSetInitialGuessNonzero(ksp_,PETSC_FALSE));
        ++builds_;
    }
    void verify_settings() const {
        PC pc,sub_pc; KSP* sub=nullptr; PetscInt count=0,levels=0;
        KSPNormType norm; PCSide side; KSPType method,sub_method; PCType precond,local_precond;
        check(KSPGetType(ksp_,&method)); check(KSPGetPC(ksp_,&pc)); check(PCGetType(pc,&precond));
        check(KSPGetNormType(ksp_,&norm)); check(KSPGetPCSide(ksp_,&side));
        check(PCBJacobiGetSubKSP(pc,&count,nullptr,&sub));
        require(count==1,"unsupported PETSc block partition");
        check(KSPGetType(sub[0],&sub_method)); check(KSPGetPC(sub[0],&sub_pc)); check(PCGetType(sub_pc,&local_precond));
        require(std::string(method)==KSPFGMRES && std::string(precond)==PCBJACOBI && std::string(sub_method)==KSPPREONLY
            && std::string(local_precond)==PCILU && side==PC_RIGHT && norm==KSP_NORM_UNPRECONDITIONED,
            "PETSc momentum defaults differ from the supported profile");
        check(PCFactorGetLevels(sub_pc,&levels)); require(levels==0,"unsupported PETSc local factorization");
        MatType matrix_type; VecType vector_type; Mat factor; MatSolverType factor_type;
        check(MatGetType(matrix_,&matrix_type)); check(VecGetType(x_,&vector_type));
        check(PCFactorGetMatrix(sub_pc,&factor)); check(MatFactorGetSolverType(factor,&factor_type));
        require(device()?(std::string(matrix_type)==MATSEQAIJCUSPARSE || std::string(matrix_type)==MATMPIAIJCUSPARSE)
                && (std::string(vector_type)==VECSEQCUDA || std::string(vector_type)==VECMPICUDA)
                && factor_type && std::string(factor_type)==MATSOLVERCUSPARSE
            : factor_type && std::string(factor_type)==MATSOLVERPETSC,"PETSc momentum matrix, vector or factor backend changed");
    }
};

template<int C> struct PetscSimpleSolve {
    static_assert(C==3,"PETSc is supported for SIMPLE momentum only");
    PetscMomentumSolver solver;
    distributed::Buffer<double> b,x;
    explicit PetscSimpleSolve(MPI_Comm comm):solver(comm) {}
    void set_tolerances(double relative,double absolute) { solver.set_stopping_tolerances(relative,absolute); }
    double* rhs(std::size_t rows) { b.resize(rows); x.resize(rows); return distributed::raw(b); }
    const double* rhs() const { return distributed::raw(b); }
    const double* solution() const { return distributed::raw(x); }
    std::size_t size() const { return x.size(); }
    template<class System> bool operator()(const System& system) { return solver.solve(system,distributed::raw(b),distributed::raw(x)); }
};
}
