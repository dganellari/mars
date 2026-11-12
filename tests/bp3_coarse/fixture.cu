#define BP3_CONSOLIDATE_TEST
#include "backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu"
#include <HYPRE_IJ_mv.h>
#include <cmath>
#include <string>

// Failure diagnostics only: these tiny fixture matrices may be copied to the
// host to separate redistribution errors from a device matvec or MPI error.
template<class T> std::vector<T> read_values(const T* data,int count)
{
    std::vector<T> values(count);
#ifdef BP3_CONSOLIDATE_CUDA
    if (count) bp3::gpu(cudaMemcpyAsync(values.data(),data,count*sizeof(T),cudaMemcpyDeviceToHost,bp3::stream()));
    bp3::sync();
#else
    if (count) std::copy(data,data+count,values.begin());
#endif
    return values;
}

std::vector<double> read_matrix(hypre_ParCSRMatrix* a)
{
    const int n=int(hypre_ParCSRMatrixGlobalNumRows(a));
    const int first=int(hypre_ParCSRMatrixFirstRowIndex(a));
    bp3::check(n>0 && n<=128,"diagnostic matrix exceeds fixture limit");
    std::vector<double> dense(n*n,0);
    for (int offd=0;offd<2;++offd) {
        auto* block=offd ? hypre_ParCSRMatrixOffd(a) : hypre_ParCSRMatrixDiag(a);
        const int rows=hypre_CSRMatrixNumRows(block),nnz=hypre_CSRMatrixNumNonzeros(block);
        auto i=read_values(hypre_CSRMatrixI(block),rows+1);
        auto j=read_values(hypre_CSRMatrixJ(block),nnz);
        auto values=read_values(hypre_CSRMatrixData(block),nnz);
        bp3::check(i.front()==0 && i.back()==nnz,"diagnostic CSR endpoints differ");
        for (int row=0;row<rows;++row) {
            bp3::check(i[row]>=0 && i[row]<=i[row+1] && i[row+1]<=nnz,"diagnostic CSR offsets invalid");
            for (int k=i[row];k<i[row+1];++k) {
                bp3::check(j[k]>=0 && j[k]<hypre_CSRMatrixNumCols(block),"diagnostic local column invalid");
                const int col=offd ? int(hypre_ParCSRMatrixColMapOffd(a)[j[k]])
                                   : int(hypre_ParCSRMatrixFirstColDiag(a))+j[k];
                bp3::check(first+row<n && col>=0 && col<n,"diagnostic global index invalid");
                dense[(first+row)*n+col]+=values[k];
            }
        }
    }
    return dense;
}

void diagnose_distributed(bp3::State& s,hypre_ParCSRMatrix* a,hypre_ParVector* v,
                          hypre_ParVector* y,const std::vector<double>& dense)
{
    auto value=[](int col) {return 0.25+(col%13)*0.0625;};
    auto difference=[](const std::vector<double>& actual,const std::vector<double>& expected) {
        double error=0;
        bp3::check(actual.size()==expected.size(),"diagnostic vector size differs");
        for (size_t k=0;k<actual.size();++k) {
            if (!std::isfinite(actual[k])) return HUGE_VAL;
            error=std::max(error,std::abs(actual[k]-expected[k]));
        }
        return error;
    };
    auto* x=hypre_ParVectorLocalVector(v);
    const auto input=read_values(hypre_VectorData(x),s.local);
    const auto original=read_values(hypre_VectorData(hypre_ParVectorLocalVector(y)),s.local);
    std::vector<double> input_expected(s.local),expected(s.local,0),diag_expected(s.local,0),offd_expected(s.local,0);
    for (int row=0;row<s.local;++row) {
        input_expected[row]=value(s.first+row);
        for (int col=0;col<s.rows;++col) {
            const double term=dense[(s.first+row)*s.rows+col]*value(col);
            expected[row]+=term;
            if (col>=s.first && col<s.first+s.local) diag_expected[row]+=term;
            else offd_expected[row]+=term;
        }
    }
    auto* offd=hypre_ParCSRMatrixOffd(a);
    const int ghosts=hypre_CSRMatrixNumCols(offd);
    std::vector<double> ghost_expected(ghosts);
    for (int k=0;k<ghosts;++k) ghost_expected[k]=value(int(hypre_ParCSRMatrixColMapOffd(a)[k]));

#ifdef BP3_CONSOLIDATE_CUDA
    // Inspect the first failed call's buffers before another matvec reuses them.
    auto* pkg=hypre_ParCSRMatrixCommPkg(a);
    bp3::check(pkg!=nullptr,"diagnostic communication package missing");
    const int sends=hypre_ParCSRCommPkgSendMapStart(pkg,hypre_ParCSRCommPkgNumSends(pkg));
    const auto send_map=read_values(hypre_ParCSRCommPkgDeviceSendMapElmts(pkg),sends);
    const auto send_values=read_values(hypre_ParCSRCommPkgBufData(pkg),sends);
    const auto recv_values=read_values(hypre_ParCSRCommPkgTmpData(pkg),ghosts);
    std::vector<double> send_expected(sends);
    int map_errors=0;
    for (int k=0;k<sends;++k) {
        const int col=hypre_ParCSRCommPkgSendMapElmts(pkg)[k];
        map_errors+=send_map[k]!=col;
        send_expected[k]=value(s.first+col);
    }
    const int gpu_mpi=hypre_GetGpuAwareMPI();
    int device=-1;
    bp3::gpu(cudaGetDevice(&device));
#endif

    auto* scratch=hypre_SeqVectorCreate(s.local);
    auto* halo=hypre_SeqVectorCreate(ghosts);
    bp3::hypre(hypre_SeqVectorInitialize_v2(scratch,bp3::location));
    bp3::hypre(hypre_SeqVectorInitialize_v2(halo,bp3::location));
#ifdef BP3_CONSOLIDATE_CUDA
    if (ghosts) bp3::gpu(cudaMemcpyAsync(hypre_VectorData(halo),ghost_expected.data(),ghosts*sizeof(double),
                                       cudaMemcpyHostToDevice,bp3::stream()));
#else
    if (ghosts) std::copy(ghost_expected.begin(),ghost_expected.end(),hypre_VectorData(halo));
#endif
    bp3::hypre(hypre_CSRMatrixMatvec(1,hypre_ParCSRMatrixDiag(a),x,0,scratch));
    const auto diag_y=read_values(hypre_VectorData(scratch),s.local);
    bp3::hypre(hypre_SeqVectorSetConstantValues(scratch,0));
    if (ghosts) bp3::hypre(hypre_CSRMatrixMatvec(1,offd,halo,0,scratch));
    const auto offd_y=read_values(hypre_VectorData(scratch),s.local);
    // Use the same beta=1 accumulation as ParCSR, but with a known halo and a
    // completed diagonal product, to separate SpMV accumulation from transport.
    bp3::hypre(hypre_CSRMatrixMatvec(1,hypre_ParCSRMatrixDiag(a),x,0,scratch));
    bp3::sync();
    if (ghosts) bp3::hypre(hypre_CSRMatrixMatvec(1,offd,halo,1,scratch));
    const double accumulated_error=difference(read_values(hypre_VectorData(scratch),s.local),expected);
    bp3::hypre(hypre_SeqVectorDestroy(halo));
    bp3::hypre(hypre_SeqVectorDestroy(scratch));

    HYPRE_BigInt starts[2]={s.first,s.first+s.local};
    auto* replay=hypre_ParVectorCreate(s.comm,s.rows,starts);
    bp3::hypre(hypre_ParVectorInitialize_v2(replay,bp3::location));
    auto run=[&]() {
        bp3::sync();
        bp3::hypre(hypre_ParCSRMatrixMatvec(1,a,v,0,replay));
        return difference(read_values(hypre_VectorData(hypre_ParVectorLocalVector(replay)),s.local),expected);
    };
    const double repeat_error=run();
#ifdef BP3_CONSOLIDATE_CUDA
    // These controls cannot turn the original failed gate into a pass.
    bp3::hypre(HYPRE_SetGpuAwareMPI(0));
    const double staged_error=run();
    bp3::hypre(HYPRE_SetGpuAwareMPI(gpu_mpi));
    bp3::hypre(HYPRE_SetSpMVUseVendor(0));
    const double own_error=run();
    bp3::hypre(HYPRE_SetSpMVUseVendor(1));
    const double vendor_error=run();
#endif
    bp3::hypre(hypre_ParVectorDestroy(replay));
    int worst=0;
    for (int row=1;row<s.local;++row)
        if (std::abs(original[row]-expected[row])>std::abs(original[worst]-expected[worst])) worst=row;
    for (int rank=0;rank<s.size;++rank) {
        bp3::mpi(MPI_Barrier(s.comm));
        if (s.rank!=rank) continue;
        std::fprintf(stderr,"DIAGNOSTIC rank=%d first=%d vector_first=%lld local=%d input_max_abs=%.17g diag_max_abs=%.17g offd_known_halo_max_abs=%.17g repeat_max_abs=%.17g\n",
                     s.rank,s.first,static_cast<long long>(hypre_ParVectorFirstIndex(v)),s.local,
                     difference(input,input_expected),difference(diag_y,diag_expected),difference(offd_y,offd_expected),repeat_error);
        if (s.local) std::fprintf(stderr,"DIAGNOSTIC rank=%d worst_row=%d original=%.17g expected=%.17g diag_expected=%.17g offd_expected=%.17g\n",
                                 s.rank,s.first+worst,original[worst],expected[worst],diag_expected[worst],offd_expected[worst]);
        std::fprintf(stderr,"DIAGNOSTIC rank=%d offd_rows=%d offd_cols=%d offd_nnz=%d known_halo_beta1_max_abs=%.17g\n",
                     s.rank,hypre_CSRMatrixNumRows(offd),ghosts,hypre_CSRMatrixNumNonzeros(offd),accumulated_error);
#ifdef BP3_CONSOLIDATE_CUDA
        std::fprintf(stderr,"DIAGNOSTIC rank=%d device=%d gpu_mpi=%d send_map_errors=%d send_max_abs=%.17g halo_max_abs=%.17g staged_mpi_max_abs=%.17g own_spmv_max_abs=%.17g vendor_spmv_max_abs=%.17g\n",
                     s.rank,device,gpu_mpi,map_errors,difference(send_values,send_expected),
                     difference(recv_values,ghost_expected),staged_error,own_error,vendor_error);
#endif
        std::fflush(stderr);
    }
}

void diagnose_matrix(bp3::State& s,hypre_ParCSRMatrix* a,hypre_ParVector* v,hypre_ParVector* y)
{
    auto local=read_matrix(a);
    std::vector<double> expected(s.rank==0 ? s.rows*s.rows : 0);
    bp3::mpi(MPI_Reduce(local.data(),expected.data(),s.rows*s.rows,MPI_DOUBLE,MPI_SUM,0,s.comm));
    auto local_y=read_values(hypre_VectorData(hypre_ParVectorLocalVector(y)),s.local);
    std::vector<double> host_y(s.rank==0 ? s.rows : 0);
    bp3::mpi(MPI_Gatherv(local_y.data(),s.local,MPI_DOUBLE,host_y.data(),s.counts.data(),s.offsets.data(),MPI_DOUBLE,0,s.comm));
    // Repeat only the vector gather, comparing device transport to host transport.
    auto* rhs=s.rank==0 ? hypre_VectorData(hypre_ParVectorLocalVector(s.rhs)) : nullptr;
    bp3::mpi(MPI_Gatherv(hypre_VectorData(hypre_ParVectorLocalVector(y)),s.local,MPI_DOUBLE,
                        rhs,s.counts.data(),s.offsets.data(),MPI_DOUBLE,0,s.comm));
    if (s.rank==0) {
        auto actual=read_matrix(s.matrix);
        auto device_y=read_values(rhs,s.rows);
        bp3::hypre(hypre_ParCSRMatrixMatvec(1,s.matrix,s.solution,0,s.rhs));
        auto root_y=read_values(rhs,s.rows);
        double matrix_error=0,matrix_scale=0,action_error=0,action_scale=0,gather_error=0,root_error=0;
        int worst=0;
        for (int k=0;k<s.rows*s.rows;++k) {
            bp3::check(std::isfinite(actual[k]) && std::isfinite(expected[k]),"diagnostic matrix is nonfinite");
            const double delta=std::abs(actual[k]-expected[k]);
            if (delta>matrix_error) {matrix_error=delta;worst=k;}
            matrix_scale=std::max(matrix_scale,std::abs(expected[k]));
        }
        for (int row=0;row<s.rows;++row) {
            double value=0;
            for (int col=0;col<s.rows;++col) value+=expected[row*s.rows+col]*(0.25+(col%13)*0.0625);
            bp3::check(std::isfinite(host_y[row]) && std::isfinite(device_y[row]) && std::isfinite(root_y[row]),"diagnostic action is nonfinite");
            action_error=std::max(action_error,std::abs(host_y[row]-value));
            action_scale=std::max(action_scale,std::abs(value));
            gather_error=std::max(gather_error,std::abs(device_y[row]-host_y[row]));
            root_error=std::max(root_error,std::abs(root_y[row]-value));
        }
        std::fprintf(stderr,"DIAGNOSTIC entries max_abs=%.17g relative=%.17g row=%d col=%d expected=%.17g actual=%.17g\n",
                     matrix_error,matrix_error/(matrix_scale+1e-300),worst/s.rows,worst%s.rows,expected[worst],actual[worst]);
        std::fprintf(stderr,"DIAGNOSTIC distributed_matvec_vs_host relative=%.17g root_matvec_vs_host relative=%.17g repeated_device_gather_vs_host max_abs=%.17g\n",
                     action_error/(action_scale+1e-300),root_error/(action_scale+1e-300),gather_error);
        std::fflush(stderr);
    }
    diagnose_distributed(s,a,v,y,local);
    bp3::mpi(MPI_Barrier(s.comm));
}

#ifdef BP3_CONSOLIDATE_CUDA
__global__ void fill_kernel(double* values,int n,int first,int global,int variant)
{
    int k=int(blockIdx.x*blockDim.x+threadIdx.x);
    if (k<n) values[k]=variant ? sin(3.14159265358979323846*(first+k+1)/(global+1)) : 0.25+((first+k)%13)*0.0625;
}
#endif
void fill(hypre_ParVector* v,int variant=0)
{
    auto* local=hypre_ParVectorLocalVector(v);
    const int n=hypre_VectorSize(local), first=int(hypre_ParVectorFirstIndex(v));
    auto* data=hypre_VectorData(local);
#ifdef BP3_CONSOLIDATE_CUDA
    if (n) fill_kernel<<<(n+127)/128,128,0,bp3::stream()>>>(data,n,first,int(hypre_ParVectorGlobalSize(v)),variant);
    bp3::gpu(cudaGetLastError());
#else
    for (int k=0;k<n;++k) data[k]=variant ? std::sin(3.14159265358979323846*(first+k+1)/(hypre_ParVectorGlobalSize(v)+1)) : 0.25+((first+k)%13)*0.0625;
#endif
}

// Compare the gathered matrix action to the original distributed matrix action.
// This checks off-rank columns as well as values; convergence alone could hide a
// defective but still usable preconditioner.
void check_matrix(HYPRE_Solver solver)
{
    auto* parent=reinterpret_cast<hypre_ParAMGData*>(solver);
    auto found=bp3::states.find(parent);
    if (found==bp3::states.end()) return;
    auto& s=*found->second;
    auto* a=s.source;
    HYPRE_BigInt starts[2]={s.first,s.first+s.local};
    auto* v=hypre_ParVectorCreate(s.comm,s.rows,starts);
    auto* y=hypre_ParVectorCreate(s.comm,s.rows,starts);
    bp3::hypre(hypre_ParVectorInitialize_v2(v,bp3::location));
    bp3::hypre(hypre_ParVectorInitialize_v2(y,bp3::location));
    fill(v);
    bp3::hypre(hypre_ParCSRMatrixMatvec(1,a,v,0,y));
    bp3::sync();
    auto* root_rhs=s.rank==0 ? hypre_VectorData(hypre_ParVectorLocalVector(s.rhs)) : nullptr;
    bp3::mpi(MPI_Gatherv(hypre_VectorData(hypre_ParVectorLocalVector(y)),s.local,MPI_DOUBLE,
                        root_rhs,s.counts.data(),s.offsets.data(),MPI_DOUBLE,0,s.comm));
    double error=0;
    if (s.rank==0) {
        fill(s.solution);
        const double scale=hypre_ParVectorInnerProd(s.rhs,s.rhs);
        bp3::hypre(hypre_ParCSRMatrixMatvec(1,s.matrix,s.solution,-1,s.rhs));
        error=std::sqrt(hypre_ParVectorInnerProd(s.rhs,s.rhs)/(scale+1e-300));
    }
    bp3::mpi(MPI_Bcast(&error,1,MPI_DOUBLE,0,s.comm));
    if (s.rank==0) {std::fprintf(stderr,"gathered matrix relative action error %.17g\n",error);std::fflush(stderr);}
    if (!std::isfinite(error) || error>=1e-12) diagnose_matrix(s,a,v,y);
    bp3::check(std::isfinite(error) && error<1e-12,"gathered matrix action differs");
    bp3::hypre(hypre_ParVectorDestroy(v)); bp3::hypre(hypre_ParVectorDestroy(y));
}

HYPRE_Int fixture_setup(HYPRE_Solver solver,HYPRE_ParCSRMatrix a,HYPRE_ParVector f,HYPRE_ParVector u)
{
    return mars_bp3_amg_setup(solver,a,f,u,bp3::threshold());
}

int main(int argc,char** argv)
{
    MPI_Init(&argc,&argv);
    int rank=0,size=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    const std::string mode=argc>1 ? argv[1] : "balanced";
#ifdef BP3_CONSOLIDATE_CUDA
    MPI_Comm local_comm; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local_comm);
    int local_rank=0,devices=0; MPI_Comm_rank(local_comm,&local_rank);
    bp3::gpu(cudaGetDeviceCount(&devices)); bp3::check(devices>0,"fixture has no GPU");
    bp3::gpu(cudaSetDevice(devices==1 ? 0 : local_rank%devices));
    MPI_Comm_free(&local_comm);
#endif
    bp3::hypre(HYPRE_Initialize());
    const char* profile=std::getenv("BP3_FIXTURE_PROFILE");
    if (!profile) profile="mars";
    const bool mars=std::strcmp(profile,"mars")==0;
    const bool historical=std::strcmp(profile,"historical")==0;
    bp3::check(mars || historical || std::strcmp(profile,"own-spmv")==0,"invalid fixture profile");
#ifdef BP3_CONSOLIDATE_CUDA
    bp3::hypre(HYPRE_SetSpMVUseVendor(historical ? 1 : 0));
    bp3::hypre(HYPRE_SetSpGemmUseVendor(mars ? 0 : 1));
    // Match the measured MARS initializer before the first device allocation.
    // The other profiles retain the old fixture setup for causal comparisons.
    if (mars) bp3::hypre(HYPRE_SetGPUMemoryPoolSize(2,1,27,size_t(2048)<<20));
#endif
    if (rank==0) std::printf("FIXTURE profile=%s threshold=%d\n",profile,bp3::threshold());
    bp3::hypre(HYPRE_SetGpuAwareMPI(1));
    bp3::hypre(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST));
    bp3::hypre(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST));
    const int n=argc>2 ? std::stoi(argv[2]) : 511;
    bp3::check(n>=511 && n<=32767,"fixture row count out of range");
    int active=size,slot=rank;
    if (mode=="empty-root" && size>1) {active=size-1;slot=rank-1;}
    else if (mode=="empty-tail" && size>1) {active=1;slot=rank;}
    else bp3::check(mode=="balanced" || size==1,"unknown fixture partition");
    const int first=slot<0 ? 0 : (slot>=active ? n : n*slot/active);
    const int end=slot<0 ? 0 : (slot>=active ? n : n*(slot+1)/active);
    const int count=end-first;
    HYPRE_IJMatrix ij;
    bp3::hypre(HYPRE_IJMatrixCreate(MPI_COMM_WORLD,first,end-1,first,end-1,&ij));
    bp3::hypre(HYPRE_IJMatrixSetObjectType(ij,HYPRE_PARCSR));
    bp3::hypre(HYPRE_IJMatrixInitialize(ij));
    for (HYPRE_BigInt row=first;row<end;++row) {
        HYPRE_BigInt cols[3]={row}; double values[3]={n==511 ? 2.0 : 2.01}; int len=1;
        if (row>0) {cols[len]=row-1;values[len++]=-1;}
        if (row<n-1) {cols[len]=row+1;values[len++]=-1;}
        bp3::hypre(HYPRE_IJMatrixSetValues(ij,1,&len,&row,cols,values));
    }
    bp3::hypre(HYPRE_IJMatrixAssemble(ij));
    hypre_ParCSRMatrix* a=nullptr;
    bp3::hypre(HYPRE_IJMatrixGetObject(ij,reinterpret_cast<void**>(&a)));
    bp3::hypre(hypre_ParCSRMatrixMigrate(a,bp3::location));
    HYPRE_BigInt starts[2]={first,end};
    auto create=[&]() {
        auto* v=hypre_ParVectorCreate(MPI_COMM_WORLD,n,starts);
        bp3::hypre(hypre_ParVectorInitialize_v2(v,bp3::location)); return v;
    };
    auto* exact=create();auto* rhs=create();auto* x=create();auto* residual=create();
#ifdef BP3_CONSOLIDATE_CUDA
    bp3::hypre(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE));
    bp3::hypre(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE));
#endif
    fill(exact);
    bp3::hypre(hypre_ParCSRMatrixMatvec(1,a,exact,0,rhs));
    double worst=0;
    HYPRE_Solver amg;
    bp3::hypre(HYPRE_BoomerAMGCreate(&amg));
    for (int repeat=0;repeat<2;++repeat) {
        HYPRE_Solver pcg;
        bp3::hypre(HYPRE_BoomerAMGSetCoarsenType(amg,8));
        bp3::hypre(HYPRE_BoomerAMGSetInterpType(amg,6));
        bp3::hypre(HYPRE_BoomerAMGSetPMaxElmts(amg,4));
        bp3::hypre(HYPRE_BoomerAMGSetStrongThreshold(amg,.25));
        bp3::hypre(HYPRE_BoomerAMGSetMaxCoarseSize(amg,9));
        bp3::hypre(HYPRE_BoomerAMGSetRelaxType(amg,18));
        bp3::hypre(HYPRE_BoomerAMGSetRelaxOrder(amg,0));
        bp3::hypre(HYPRE_BoomerAMGSetNumSweeps(amg,1));
        bp3::hypre(HYPRE_BoomerAMGSetCycleRelaxType(amg,18,3));
        bp3::hypre(HYPRE_BoomerAMGSetCycleNumSweeps(amg,4,3));
        bp3::hypre(HYPRE_BoomerAMGSetTol(amg,0));
        bp3::hypre(HYPRE_BoomerAMGSetMaxIter(amg,bp3::threshold()==-2 ? 3 : 1));
        bp3::hypre(HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,&pcg));
        bp3::hypre(HYPRE_ParCSRPCGSetMaxIter(pcg,200));
        bp3::hypre(HYPRE_ParCSRPCGSetTol(pcg,1e-11));
        bp3::hypre(HYPRE_ParCSRPCGSetTwoNorm(pcg,1));
        bp3::hypre(HYPRE_ParCSRPCGSetPrecond(pcg,mars_bp3_amg_solve,fixture_setup,amg));
        auto pa=reinterpret_cast<HYPRE_ParCSRMatrix>(a);
        auto pb=reinterpret_cast<HYPRE_ParVector>(rhs),px=reinterpret_cast<HYPRE_ParVector>(x);
        bp3::hypre(HYPRE_ParCSRPCGSetup(pcg,pa,pb,px));
        check_matrix(amg);
        // A fixed symmetric correction must not depend on a previous right-hand side.
        fill(exact,0);
        bp3::hypre(hypre_ParVectorSetConstantValues(x,0));
        bp3::hypre(mars_bp3_amg_solve(amg,pa,reinterpret_cast<HYPRE_ParVector>(exact),px));
        bp3::hypre(hypre_ParVectorCopy(x,residual));
        const double positive0=hypre_ParVectorInnerProd(exact,residual);
        const double correction_norm=hypre_ParVectorInnerProd(residual,residual);
        fill(rhs,1);
        bp3::hypre(hypre_ParVectorSetConstantValues(x,0));
        bp3::hypre(mars_bp3_amg_solve(amg,pa,pb,px));
        const double positive1=hypre_ParVectorInnerProd(rhs,x);
        const double xy=hypre_ParVectorInnerProd(exact,x),yx=hypre_ParVectorInnerProd(rhs,residual);
        bp3::check(std::isfinite(positive0) && positive0>0 && std::isfinite(positive1) && positive1>0,
                   "preconditioner has nonpositive curvature");
        bp3::check(std::isfinite(xy) && std::isfinite(yx) && std::abs(xy-yx)<=1e-10*(std::abs(xy)+std::abs(yx)),
                   "preconditioner symmetry check failed");
        bp3::hypre(hypre_ParVectorSetConstantValues(x,0));
        bp3::hypre(mars_bp3_amg_solve(amg,pa,reinterpret_cast<HYPRE_ParVector>(exact),px));
        bp3::hypre(hypre_ParVectorAxpy(-1,residual,x));
        const double repeated=std::sqrt(hypre_ParVectorInnerProd(x,x)/correction_norm);
        bp3::check(std::isfinite(repeated) && repeated<=1e-12,"preconditioner depends on earlier calls");
        if (rank==0) std::printf("preconditioner symmetry %.3e repeat %.3e positive PASS\n",
                                std::abs(xy-yx)/(std::abs(xy)+std::abs(yx)),repeated);
        for (int solve=0;solve<2;++solve) {
            fill(exact,solve);
            bp3::hypre(hypre_ParCSRMatrixMatvec(1,a,exact,0,rhs));
            bp3::hypre(hypre_ParVectorSetConstantValues(x,0));
            bp3::hypre(HYPRE_ParCSRPCGSolve(pcg,pa,pb,px));
            bp3::hypre(hypre_ParVectorCopy(rhs,residual));
            bp3::hypre(hypre_ParCSRMatrixMatvec(-1,a,x,1,residual));
            double relative=std::sqrt(hypre_ParVectorInnerProd(residual,residual)/hypre_ParVectorInnerProd(rhs,rhs));
            bp3::check(std::isfinite(relative) && relative<1e-10,"true residual check failed");
            worst=std::max(worst,relative);
            bp3::hypre(hypre_ParVectorCopy(x,residual));
            bp3::hypre(hypre_ParVectorAxpy(-1,exact,residual));
            const double error=std::sqrt(hypre_ParVectorInnerProd(residual,residual)/hypre_ParVectorInnerProd(exact,exact));
            bp3::check(std::isfinite(error) && error<1e-7,"manufactured solution differs");
        }
        bp3::hypre(HYPRE_ParCSRPCGDestroy(pcg));
        bp3::end(amg);
    }
    bp3::hypre(mars_bp3_amg_destroy(amg));
    bp3::check(bp3::audit.setups==bp3::audit.destroyed,"coarse state leaked");
    if (size > 1 && bp3::threshold() == -2)
        bp3::check(bp3::audit.level == 0 && bp3::audit.rows == n && bp3::audit.root_calls > 0 &&
                   bp3::audit.cycles == 3 * bp3::audit.root_calls, "root cycle grouping differs");
    if (size>1 && bp3::threshold()!=0) {
        const bool stalled=n==32767 && bp3::threshold()==100;
        if (stalled) bp3::check(bp3::audit.skipped_setups==2 && bp3::audit.cycles==0,"stalled hierarchy fallback was not exercised");
        else bp3::check(bp3::audit.setups==2 && bp3::audit.cycles>0,"consolidation did not activate");
    }
    if (mode=="empty-root" && size>1 && rank==0)
        bp3::check(count==0 && bp3::audit.empty_setups==(bp3::threshold()!=0 ? 2 : 0),"empty root was not exercised");
    if (rank==0) std::printf("PASS %s ranks=%d true=%.3e setups=%lld cycles=%lld\n",mode.c_str(),size,worst,bp3::audit.setups,bp3::audit.cycles);
    for (auto* v:{exact,rhs,x,residual}) bp3::hypre(hypre_ParVectorDestroy(v));
    bp3::hypre(HYPRE_IJMatrixDestroy(ij));
    bp3::hypre(HYPRE_Finalize());MPI_Finalize();return 0;
}
