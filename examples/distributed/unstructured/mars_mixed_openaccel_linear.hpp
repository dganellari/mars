#pragma once
#include <petscksp.h>
#include <memory>
#include <type_traits>

namespace mixedflow {
inline void petsc_check(PetscErrorCode code) {
    if(code) {std::fprintf(stderr,"rank %d PETSc error %d\n",rankId,int(code));MPI_Abort(MPI_COMM_WORLD,1);}
}

inline void reference_mpi_check(int code) {
    if(code!=MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,code);
}

// The reference allocates zero entries too: they affect ILU(0) and AMG setup.
inline std::vector<std::set<int64_t>> reference_pattern(const Channel& channel) {
    const auto& mesh=channel.input_.mesh;
    const auto rows=host<int64_t>(channel.pressure_.rowIds);
    std::vector<std::set<int64_t>> pattern(rows.size());
    auto insert=[&](int64_t row,int64_t column) {
        if(row<rows.front() || row>rows.back()) throw std::runtime_error("reference graph row ownership");
        pattern[row-rows.front()].insert(column);
    };
    for(const auto& cell:channel.input_.cells) for(int k=0;k<cell.count;++k) {
        const auto row=mesh.dofs[cell.node[k]];
        if(owner(row,mesh.globalNodes)!=rankId) continue;
        for(int j=0;j<cell.count;++j) {
            int different=0;for(int a=0;a<3;++a) different+=sign(k,a)!=sign(j,a);
            if(cell.count==4 || different<=1) insert(row,mesh.dofs[cell.node[j]]);
        }
    }
    // Send selected-parent graph entries to their row owners. No global gather.
    const auto interfaces=host<Interface>(channel.interfaces_);
    const auto& peers=mesh.candidatePeers;
    std::map<int,std::set<std::pair<int64_t,int64_t>>> outgoing;
    for(int peer:peers) outgoing[peer];
    for(const auto& ip:interfaces) for(int k=0;k<ip.count[0];++k) {
        const auto row=mesh.dofs[ip.node[k]];const int target=owner(row,mesh.globalNodes);
        for(int j=0;j<ip.count[0]+ip.count[1];++j) {
            const auto column=mesh.dofs[ip.node[j]];
            if(target==rankId) insert(row,column);
            else {
                if(!outgoing.count(target)) throw std::runtime_error("reference graph missing row-owner peer");
                outgoing[target].insert({row,column});
            }
        }
    }
    std::vector<std::vector<int64_t>> send(peers.size()),receive(peers.size());
    std::vector<int> send_counts(peers.size()),receive_counts(peers.size());
    bool bounded=true;
    for(size_t i=0;i<peers.size();++i) {
        bounded=bounded && outgoing[peers[i]].size()<=size_t(INT_MAX/2);
        for(auto pair:outgoing[peers[i]]) {send[i].push_back(pair.first);send[i].push_back(pair.second);}
    }
    check(bounded,"reference graph MPI count range");
    std::vector<MPI_Request> requests(2*peers.size());
    for(size_t i=0;i<peers.size();++i) {
        send_counts[i]=int(send[i].size());
        reference_mpi_check(MPI_Irecv(&receive_counts[i],1,MPI_INT,peers[i],310,MPI_COMM_WORLD,&requests[2*i]));
        reference_mpi_check(MPI_Isend(&send_counts[i],1,MPI_INT,peers[i],310,MPI_COMM_WORLD,&requests[2*i+1]));
    }
    reference_mpi_check(MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE));
    for(size_t i=0;i<peers.size();++i) {
        if(receive_counts[i]<0 || receive_counts[i]%2) throw std::runtime_error("invalid reference graph message");
        receive[i].resize(receive_counts[i]);
        reference_mpi_check(MPI_Irecv(receive[i].data(),receive_counts[i],MPI_INT64_T,peers[i],311,MPI_COMM_WORLD,&requests[2*i]));
        reference_mpi_check(MPI_Isend(send[i].data(),send_counts[i],MPI_INT64_T,peers[i],311,MPI_COMM_WORLD,&requests[2*i+1]));
    }
    reference_mpi_check(MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE));
    for(const auto& message:receive) for(size_t i=0;i<message.size();i+=2) insert(message[i],message[i+1]);
    return pattern;
}

// This opt-in profile preserves the CPU reference preconditioners. Staging is
// explicit and reuses buffers; the default GPU solver does not enter this path.
inline void reference_copy(double* destination,const double* source,size_t count,bool to_host) {
#if defined(__CUDACC__) || defined(__HIPCC__)
    check(cudaMemcpy(destination,source,count*sizeof(double),to_host?cudaMemcpyDeviceToHost:cudaMemcpyHostToDevice)==cudaSuccess,
          "reference solver field staging");
#else
    (void)to_host;std::copy_n(source,count,destination);
#endif
}
class OpenAccelLinear {
public:
    OpenAccelLinear(Graph& graph,const std::vector<std::set<int64_t>>& pattern,const Mesh& mesh,bool momentum)
        : graph_(graph),momentum_(momentum) {
        static_assert(std::is_same_v<PetscScalar,double>,"OpenAccel profile requires real FP64 PETSc");
        const auto rows=host<int64_t>(graph.rowIds),columns=host<int64_t>(graph.columns);
        const auto offsets=host<int>(graph.offsets);
        const int components=momentum?3:1;
        check(!rows.empty() && int64_t(components)*mesh.globalNodes<=std::numeric_limits<PetscInt>::max(),
              "reference solver index range");
        std::vector<bool> used(columns.size());
        offsets_.push_back(0);
        for(const auto& row:pattern) for(int a=0;a<components;++a) {
            const int r=int(offsets_.size())-1;rows_.push_back(rows[r]);
            for(auto column:row) for(int b=0;b<components;++b) {
                const int64_t global=components*column+b;
                columns_.push_back(global);
                auto begin=columns.begin()+offsets[r],end=columns.begin()+offsets[r+1];
                const auto found=std::lower_bound(begin,end,global);
                int source=-1;
                if(found!=end && *found==global) {source=int(found-columns.begin());used[source]=true;}
                source_.push_back(source);
            }
            if(columns_.size()>size_t(INT_MAX)) throw std::runtime_error("reference CSR exceeds local index capacity");
            counts_.push_back(int(columns_.size())-offsets_.back());offsets_.push_back(int(columns_.size()));
        }
        for(size_t i=0;i<used.size();++i) if(!used[i]) omitted_.push_back(i);
        assembled_.resize(columns.size());values_.resize(columns_.size());rhs_.resize(rows.size());solution_.resize(rows.size());
        if(momentum_) {
            std::vector<PetscInt> diagonal(rows.size()),off_diagonal(rows.size());
            for(size_t r=0;r<rows.size();++r) for(int p=offsets_[r];p<offsets_[r+1];++p)
                ++(columns_[p]>=rows.front() && columns_[p]<=rows.back()?diagonal[r]:off_diagonal[r]);
            petsc_check(MatCreateAIJ(MPI_COMM_WORLD,rows.size(),rows.size(),components*mesh.globalNodes,components*mesh.globalNodes,
                                    0,diagonal.data(),0,off_diagonal.data(),&matrix_));
            petsc_check(MatSetOption(matrix_,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE));
            PetscInt first,last;petsc_check(MatGetOwnershipRange(matrix_,&first,&last));
            check(first==rows.front() && last==rows.back()+1,"PETSc reference ownership");
            petsc_check(MatCreateVecs(matrix_,&x_,&b_));petsc_check(KSPCreate(MPI_COMM_WORLD,&ksp_));
            petsc_check(KSPSetType(ksp_,KSPFGMRES));PC pc;petsc_check(KSPGetPC(ksp_,&pc));petsc_check(PCSetType(pc,PCBJACOBI));
            petsc_check(KSPSetTolerances(ksp_,1e-6,1e-12,PETSC_DEFAULT,20));
            petsc_columns_.assign(columns_.begin(),columns_.end());
        } else {
            hypre_check(HYPRE_IJMatrixCreate(MPI_COMM_WORLD,rows.front(),rows.back(),rows.front(),rows.back(),&a_));
            hypre_check(HYPRE_IJMatrixSetObjectType(a_,HYPRE_PARCSR));
            for(auto vector:{&hb_,&hx_}) {
                hypre_check(HYPRE_IJVectorCreate(MPI_COMM_WORLD,rows.front(),rows.back(),vector));
                hypre_check(HYPRE_IJVectorSetObjectType(*vector,HYPRE_PARCSR));
            }
        }
    }
    OpenAccelLinear(const OpenAccelLinear&)=delete;
    OpenAccelLinear& operator=(const OpenAccelLinear&)=delete;
    ~OpenAccelLinear() {
        if(momentum_) {petsc_check(KSPDestroy(&ksp_));petsc_check(VecDestroy(&b_));petsc_check(VecDestroy(&x_));petsc_check(MatDestroy(&matrix_));}
        else {hypre_check(HYPRE_IJMatrixDestroy(a_));hypre_check(HYPRE_IJVectorDestroy(hb_));hypre_check(HYPRE_IJVectorDestroy(hx_));}
    }
    int solve(ChannelOps& ops) {
        ncflow::finish();
        reference_copy(assembled_.data(),pointer(graph_.values),assembled_.size(),true);
        reference_copy(rhs_.data(),ops.rhs(),rhs_.size(),true);
        bool zeros=true;for(int p:omitted_) zeros=zeros && assembled_[p]==0;
        check(zeros,"reference storage omits a nonzero coefficient");
        for(size_t p=0;p<source_.size();++p) values_[p]=source_[p]<0?0:assembled_[source_[p]];
        int iterations=momentum_?solve_momentum():solve_pressure();
        reference_copy(ops.solution(),solution_.data(),solution_.size(),false);
        ops.apply(ops.solution(),ops.work());ops.axpy(-1,ops.rhs(),ops.work());
        const double residual=ops.norm(ops.work()),scale=ops.norm(ops.rhs());
        check(std::isfinite(residual) && residual<=std::max(momentum_?1.01e-12:1.01e-10,(momentum_?1.01e-6:1.01e-5)*scale),
              "OpenAccel profile true linear residual");
        return iterations;
    }
private:
    int solve_momentum() {
        for(size_t r=0;r<rows_.size();++r) {
            PetscInt row=rows_[r];
            petsc_check(MatSetValues(matrix_,1,&row,counts_[r],petsc_columns_.data()+offsets_[r],values_.data()+offsets_[r],INSERT_VALUES));
        }
        petsc_check(MatAssemblyBegin(matrix_,MAT_FINAL_ASSEMBLY));petsc_check(MatAssemblyEnd(matrix_,MAT_FINAL_ASSEMBLY));
        PetscScalar* b;petsc_check(VecGetArray(b_,&b));std::copy(rhs_.begin(),rhs_.end(),b);petsc_check(VecRestoreArray(b_,&b));
        petsc_check(VecSet(x_,0));petsc_check(KSPSetOperators(ksp_,matrix_,matrix_));petsc_check(KSPSolve(ksp_,b_,x_));
        PetscInt iterations;KSPConvergedReason reason;petsc_check(KSPGetIterationNumber(ksp_,&iterations));petsc_check(KSPGetConvergedReason(ksp_,&reason));
        check(reason>0,"OpenAccel profile PETSc convergence");
        const PetscScalar* x;petsc_check(VecGetArrayRead(x_,&x));std::copy_n(x,solution_.size(),solution_.data());petsc_check(VecRestoreArrayRead(x_,&x));
        if(first_) {petsc_check(KSPView(ksp_,PETSC_VIEWER_STDOUT_WORLD));first_=false;}
        return int(iterations);
    }
    int solve_pressure() {
        std::fill(solution_.begin(),solution_.end(),0);
        hypre_check(HYPRE_IJMatrixInitialize(a_));
        hypre_check(HYPRE_IJMatrixSetValues(a_,rows_.size(),counts_.data(),rows_.data(),columns_.data(),values_.data()));hypre_check(HYPRE_IJMatrixAssemble(a_));
        hypre_check(HYPRE_IJVectorInitialize(hb_));hypre_check(HYPRE_IJVectorInitialize(hx_));
        hypre_check(HYPRE_IJVectorSetValues(hb_,rows_.size(),rows_.data(),rhs_.data()));hypre_check(HYPRE_IJVectorSetValues(hx_,rows_.size(),rows_.data(),solution_.data()));
        hypre_check(HYPRE_IJVectorAssemble(hb_));hypre_check(HYPRE_IJVectorAssemble(hx_));
        HYPRE_ParCSRMatrix a;HYPRE_ParVector b,x;
        hypre_check(HYPRE_IJMatrixGetObject(a_,reinterpret_cast<void**>(&a)));
        hypre_check(HYPRE_IJVectorGetObject(hb_,reinterpret_cast<void**>(&b)));hypre_check(HYPRE_IJVectorGetObject(hx_,reinterpret_cast<void**>(&x)));
        HYPRE_Solver solver,amg;hypre_check(HYPRE_BoomerAMGCreate(&amg));
        hypre_check(HYPRE_BoomerAMGSetCoarsenType(amg,6));hypre_check(HYPRE_BoomerAMGSetInterpType(amg,6));
        hypre_check(HYPRE_BoomerAMGSetRelaxType(amg,18));hypre_check(HYPRE_BoomerAMGSetStrongThreshold(amg,.5));
        hypre_check(HYPRE_BoomerAMGSetNumSweeps(amg,3));hypre_check(HYPRE_BoomerAMGSetMaxLevels(amg,25));
        hypre_check(HYPRE_BoomerAMGSetTruncFactor(amg,.1));hypre_check(HYPRE_BoomerAMGSetTol(amg,0));hypre_check(HYPRE_BoomerAMGSetMaxIter(amg,1));
        hypre_check(HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD,&solver));hypre_check(HYPRE_ParCSRGMRESSetMaxIter(solver,100));
        hypre_check(HYPRE_ParCSRGMRESSetTol(solver,1e-5));hypre_check(HYPRE_ParCSRGMRESSetAbsoluteTol(solver,1e-10));
        hypre_check(HYPRE_ParCSRGMRESSetPrecond(solver,HYPRE_BoomerAMGSolve,HYPRE_BoomerAMGSetup,amg));
        hypre_check(HYPRE_ParCSRGMRESSetup(solver,a,b,x));hypre_check(HYPRE_ParCSRGMRESSolve(solver,a,b,x));
        HYPRE_Int iterations;hypre_check(HYPRE_ParCSRGMRESGetNumIterations(solver,&iterations));
        hypre_check(HYPRE_IJVectorGetValues(hx_,rows_.size(),rows_.data(),solution_.data()));
        hypre_check(HYPRE_ParCSRGMRESDestroy(solver));hypre_check(HYPRE_BoomerAMGDestroy(amg));return iterations;
    }
    Graph& graph_;bool momentum_,first_=true;
    std::vector<HYPRE_BigInt> rows_,columns_;std::vector<HYPRE_Int> counts_;
    std::vector<PetscInt> petsc_columns_;
    std::vector<int> offsets_,source_,omitted_;
    std::vector<double> assembled_,values_,rhs_,solution_;
    Mat matrix_{};Vec x_{},b_{};KSP ksp_{};
    HYPRE_IJMatrix a_{};HYPRE_IJVector hb_{},hx_{};
};
}
