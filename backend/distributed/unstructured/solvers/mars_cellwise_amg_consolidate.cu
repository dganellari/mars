#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_parcsr_ls.h>
#include <mpi.h>

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <type_traits>
#include <vector>
#ifdef BP3_CONSOLIDATE_CUDA
#include <cuda_runtime.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
// Hypre 2.33 declares this accessor in its C++ utilities header alongside unrelated
// vendor-solver headers. Keep the exact declaration without requiring those SDKs.
extern "C++" cudaStream_t hypre_DeviceDataComputeStream(hypre_DeviceData* data);
#endif

// Hypre private layouts are supported only for this checked version.
static_assert(HYPRE_RELEASE_NUMBER == 23300, "Hypre 2.33.0 required");
static_assert(std::is_same<HYPRE_Complex, double>::value, "Real FP64 Hypre required");
static_assert(std::is_same<HYPRE_Int, int>::value, "32-bit local indices required");

namespace bp3 {

constexpr int root_row_limit = 300000;
constexpr int root_entry_limit = 30000000;

[[noreturn]] void fail(const char* message)
{
    std::fprintf(stderr, "BP3 consolidation: %s\n", message);
    MPI_Abort(MPI_COMM_WORLD, 91);
    std::abort();
}
void check(bool condition, const char* message) { if (!condition) fail(message); }
void mpi(int code) { check(code == MPI_SUCCESS, "MPI operation failed"); }
void hypre(int code) { check(code == 0, "Hypre operation failed"); }

#ifdef BP3_CONSOLIDATE_CUDA
constexpr HYPRE_MemoryLocation location = HYPRE_MEMORY_DEVICE;
cudaStream_t stream() { return hypre_HandleComputeStream(hypre_handle()); }
void gpu(cudaError_t code) { check(code == cudaSuccess, cudaGetErrorString(code)); }
void sync() { gpu(cudaStreamSynchronize(stream())); }
void copy(void* to, const void* from, size_t bytes)
{
    if (bytes) gpu(cudaMemcpyAsync(to, from, bytes, cudaMemcpyDeviceToDevice, stream()));
}
#else
constexpr HYPRE_MemoryLocation location = HYPRE_MEMORY_HOST;
void sync() {}
void copy(void* to, const void* from, size_t bytes) { if (bytes) std::memcpy(to, from, bytes); }
#endif

template<class T> struct Buffer {
    T* data = nullptr;
    explicit Buffer(size_t count) {
        if (count) { data = hypre_TAlloc(T, count, location); check(data!=nullptr,"scratch allocation failed"); }
    }
    ~Buffer() { hypre_TFree(data, location); }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
};

#ifdef BP3_CONSOLIDATE_CUDA
__host__ __device__
#endif
void pack_row(int row, HYPRE_BigInt first, const int* di, const int* dj, const double* da,
              const int* oi, const int* oj, const double* oa, const HYPRE_BigInt* map,
              int* lengths, int* cols, double* vals)
{
    lengths[row] = di[row+1]-di[row]+oi[row+1]-oi[row];
    int k = di[row]+oi[row];
    for (int j = di[row]; j < di[row+1]; ++j) {
        cols[k] = int(first+dj[j]); vals[k++] = da[j];
    }
    for (int j = oi[row]; j < oi[row+1]; ++j) {
        cols[k] = int(map[oj[j]]); vals[k++] = oa[j];
    }
}

#ifdef BP3_CONSOLIDATE_CUDA
__global__ void pack_kernel(int rows, HYPRE_BigInt first,
    const int* di, const int* dj, const double* da,
    const int* oi, const int* oj, const double* oa, const HYPRE_BigInt* map,
    int* lengths, int* cols, double* vals)
{
    int row = int(blockIdx.x*blockDim.x+threadIdx.x);
    if (row < rows) pack_row(row,first,di,dj,da,oi,oj,oa,map,lengths,cols,vals);
}
#endif

void pack(hypre_ParCSRMatrix* a, int* lengths, int* cols, double* vals)
{
    auto* d = hypre_ParCSRMatrixDiag(a);
    auto* o = hypre_ParCSRMatrixOffd(a);
    const int rows = hypre_CSRMatrixNumRows(d);
    if (!rows) return;
#ifdef BP3_CONSOLIDATE_CUDA
    hypre_ParCSRMatrixCopyColMapOffdToDevice(a);
    const auto* map = hypre_ParCSRMatrixDeviceColMapOffd(a);
    pack_kernel<<<(rows+127)/128,128,0,stream()>>>(rows,
#else
    const auto* map = hypre_ParCSRMatrixColMapOffd(a);
    for (int row = 0; row < rows; ++row) pack_row(row,
#endif
        hypre_ParCSRMatrixFirstColDiag(a),
        hypre_CSRMatrixI(d),hypre_CSRMatrixJ(d),hypre_CSRMatrixData(d),
        hypre_CSRMatrixI(o),hypre_CSRMatrixJ(o),hypre_CSRMatrixData(o),map,lengths,cols,vals);
#ifdef BP3_CONSOLIDATE_CUDA
    gpu(cudaGetLastError());
#endif
}

using Setup = HYPRE_Int (*)(HYPRE_Solver,HYPRE_ParCSRMatrix,HYPRE_ParVector,HYPRE_ParVector);
using Destroy = HYPRE_Int (*)(HYPRE_Solver);
Setup original_setup() { return &HYPRE_BoomerAMGSetup; }
Destroy original_destroy() { return &HYPRE_BoomerAMGDestroy; }

// Collapse when the distributed rows become too few per rank; cap root memory.
int selected_threshold(int requested, int ranks)
{
    check(ranks>0 && (requested==-2 || requested==-1 || requested==0 || (requested>9 && requested<=50000)),
          "consolidation must be off, adaptive, whole-root, or between 10 and 50000 rows");
    if (ranks==1 || requested==0) return 0;
    if (requested==-2) return root_row_limit;
    if (requested>0) return requested;
    return ranks>=50 ? 50000 : (ranks<5 ? 5000 : ranks*1000);
}

#ifdef BP3_CONSOLIDATE_TEST
int threshold()
{
    const char* value=std::getenv("BP3_CONSOLIDATE_ROWS");
    if (!value) return 0;
    char* end=nullptr;
    const long result=std::strtol(value,&end,10);
    check(end!=value && *end=='\0' && (result==-2 || (result>=0 && result<=50000)),"invalid fixture cutoff");
    return int(result);
}
#endif

struct Audit {
    long long setups = 0, destroyed = 0, cycles = 0, empty_setups = 0, skipped_setups = 0;
    int level = -1, rows = 0, nnz = 0, threshold_rows = 0;
    long long root_calls = 0;
    void write() const {
        const char* path = std::getenv("BP3_CONSOLIDATION_AUDIT");
        if (!path) return;
        std::ofstream out(path);
        check(bool(out),"cannot open audit file");
        out << "{\"version\":1,\"setups\":" << setups << ",\"destroyed\":" << destroyed
            << ",\"cycles\":" << cycles << ",\"empty_setups\":" << empty_setups
            << ",\"skipped_setups\":" << skipped_setups
            << ",\"root_calls\":" << root_calls
            << ",\"threshold_rows\":" << threshold_rows << ",\"level\":" << level << ",\"rows\":" << rows << ",\"nnz\":" << nnz << "}\n";
        out.close(); check(bool(out),"cannot write audit file");
    }
} audit;

// A full-communicator duplication gives the persistent point-to-point transfers
// their own tags. Ranks with zero rows still enter the coarse callback.
struct State {
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0, size = 0, level = 0, rows = 0, first = 0, local = 0;
    std::vector<int> counts, offsets;
    std::vector<MPI_Request> gather, scatter;
    double* bound_f = nullptr;
    double* bound_u = nullptr;
    bool bound = false;
    const bool whole;
    const int cycles_per_call;
    hypre_ParCSRMatrix* source = nullptr; // Borrowed until the parent is destroyed.
    hypre_ParCSRMatrix* matrix = nullptr;
    hypre_ParVector* rhs = nullptr;
    hypre_ParVector* solution = nullptr;
    HYPRE_Solver solver = nullptr;

    State(hypre_ParAMGData* parent, int selected, int limit, hypre_ParCSRMatrix* input = nullptr)
        : level(selected), whole(input != nullptr),
          cycles_per_call(whole ? hypre_ParAMGDataMaxIter(parent) : 1)
    {
        auto* a = input ? input : hypre_ParAMGDataAArray(parent)[level];
        source = a;
        check(cycles_per_call > 0, "root correction needs a positive fixed cycle count");
        check(hypre_ParCSRMatrixMemoryLocation(a)==location,"unexpected matrix memory location");
        check(hypre_ParCSRMatrixGlobalNumRows(a)==hypre_ParCSRMatrixGlobalNumCols(a),"matrix is not square");
        check(hypre_ParCSRMatrixFirstRowIndex(a)==hypre_ParCSRMatrixFirstColDiag(a),"row/column ownership differs");
        check(hypre_ParCSRMatrixGlobalNumRows(a)>0 && hypre_ParCSRMatrixGlobalNumRows(a)<=limit && limit<=(whole ? root_row_limit : 50000),
              "selected matrix exceeds the configured coarse limit");
        mpi(MPI_Comm_dup(hypre_ParCSRMatrixComm(a),&comm));
        mpi(MPI_Comm_rank(comm,&rank)); mpi(MPI_Comm_size(comm,&size));
        rows = int(hypre_ParCSRMatrixGlobalNumRows(a));
        first = int(hypre_ParCSRMatrixFirstRowIndex(a));
        local = hypre_CSRMatrixNumRows(hypre_ParCSRMatrixDiag(a));
        const long long local_nnz = static_cast<long long>(hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixDiag(a)))+
                                   hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixOffd(a));
        check(local_nnz>=0 && local_nnz<=INT_MAX,"local coarse matrix count overflows");
        const int nnz = int(local_nnz);
        int metadata[3] = {first,local,nnz};
        std::vector<int> all(rank==0 ? size*3 : 0);
        mpi(MPI_Gather(metadata,3,MPI_INT,all.data(),3,MPI_INT,0,comm));
        int total_nnz = 0;
        std::vector<int> nnz_counts,nnz_offsets;
        if (rank==0) {
            counts.resize(size); offsets.resize(size);
            nnz_counts.resize(size); nnz_offsets.resize(size);
            int row = 0;
            for (int peer=0; peer<size; ++peer) {
                check(all[peer*3]==row,"coarse row partitions are not contiguous");
                counts[peer]=all[peer*3+1]; offsets[peer]=row;
                nnz_counts[peer]=all[peer*3+2]; nnz_offsets[peer]=total_nnz;
                check(counts[peer]>=0 && nnz_counts[peer]>=0 && nnz_counts[peer]<=INT_MAX-total_nnz,
                      "invalid coarse matrix counts");
                row+=counts[peer]; total_nnz+=nnz_counts[peer];
            }
            check(row==rows,"coarse rows do not cover the matrix");
            check(!whole || total_nnz <= root_entry_limit, "root matrix exceeds the entry limit");
            HYPRE_BigInt starts[2]={0,rows};
            matrix=hypre_ParCSRMatrixCreate(MPI_COMM_SELF,rows,rows,starts,starts,0,total_nnz,0);
            hypre(hypre_ParCSRMatrixInitialize_v2(matrix,location));
        }
        mpi(MPI_Bcast(&total_nnz,1,MPI_INT,0,comm));
        Buffer<int> lengths(local), columns(nnz);
        Buffer<double> values(nnz);
        pack(a,lengths.data,columns.data,values.data);
        sync(); // MPI cannot consume the packed rows until their producing stream finishes.
        std::vector<MPI_Request> requests;
        auto transfer=[&](void* data,int n,MPI_Datatype type,int peer,int tag,bool receive) {
            if (!n) return;
            MPI_Request request;
            if (receive) mpi(MPI_Irecv(data,n,type,peer,tag,comm,&request));
            else mpi(MPI_Isend(data,n,type,peer,tag,comm,&request));
            requests.push_back(request);
        };
        if (rank==0) {
            auto* diag=hypre_ParCSRMatrixDiag(matrix);
            for (int peer=1;peer<size;++peer) {
                if (counts[peer]) transfer(hypre_CSRMatrixI(diag)+offsets[peer],counts[peer],MPI_INT,peer,10,true);
                if (nnz_counts[peer]) {
                    transfer(hypre_CSRMatrixJ(diag)+nnz_offsets[peer],nnz_counts[peer],MPI_INT,peer,11,true);
                    transfer(hypre_CSRMatrixData(diag)+nnz_offsets[peer],nnz_counts[peer],MPI_DOUBLE,peer,12,true);
                }
            }
            copy(hypre_CSRMatrixI(diag),lengths.data,local*sizeof(int));
            copy(hypre_CSRMatrixJ(diag),columns.data,nnz*sizeof(int));
            copy(hypre_CSRMatrixData(diag),values.data,nnz*sizeof(double));
        } else {
            transfer(lengths.data,local,MPI_INT,0,10,false);
            transfer(columns.data,nnz,MPI_INT,0,11,false);
            transfer(values.data,nnz,MPI_DOUBLE,0,12,false);
        }
        if (!requests.empty()) mpi(MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE));
        if (rank==0) {
            auto* diag=hypre_ParCSRMatrixDiag(matrix);
            auto* row_ptr=hypre_CSRMatrixI(diag);
#ifdef BP3_CONSOLIDATE_CUDA
            gpu(cudaMemsetAsync(row_ptr+rows,0,sizeof(int),stream()));
            thrust::exclusive_scan(thrust::cuda::par.on(stream()),row_ptr,row_ptr+rows+1,row_ptr);
#else
            int count=0;
            for (int row=0;row<rows;++row) {int n=row_ptr[row]; row_ptr[row]=count; count+=n;}
            row_ptr[rows]=count;
            check(count==total_nnz,"packed CSR length mismatch");
#endif
            sync();
            hypre(hypre_ParCSRMatrixSetNumNonzeros(matrix));
            hypre(hypre_ParCSRMatrixSetDNumNonzeros(matrix));
            HYPRE_BigInt starts[2]={0,rows};
            rhs=hypre_ParVectorCreate(MPI_COMM_SELF,rows,starts);
            solution=hypre_ParVectorCreate(MPI_COMM_SELF,rows,starts);
            hypre(hypre_ParVectorInitialize_v2(rhs,location));
            hypre(hypre_ParVectorInitialize_v2(solution,location));
            hypre(HYPRE_BoomerAMGCreate(&solver));
            hypre(HYPRE_BoomerAMGSetCoarsenType(solver,8));
            hypre(HYPRE_BoomerAMGSetInterpType(solver,6));
            hypre(HYPRE_BoomerAMGSetPMaxElmts(solver,4));
            hypre(HYPRE_BoomerAMGSetStrongThreshold(solver,hypre_ParAMGDataStrongThreshold(parent)));
            hypre(HYPRE_BoomerAMGSetMaxCoarseSize(solver,hypre_ParAMGDataMaxCoarseSize(parent)));
            hypre(HYPRE_BoomerAMGSetRelaxType(solver,18));
            hypre(HYPRE_BoomerAMGSetRelaxOrder(solver,0));
            hypre(HYPRE_BoomerAMGSetNumSweeps(solver,1));
            hypre(HYPRE_BoomerAMGSetCycleRelaxType(solver,18,3));
            hypre(HYPRE_BoomerAMGSetCycleNumSweeps(solver,4,3));
            hypre(HYPRE_BoomerAMGSetMaxLevels(solver,25-level));
            hypre(HYPRE_BoomerAMGSetTol(solver,0));
            hypre(HYPRE_BoomerAMGSetMaxIter(solver,cycles_per_call));
            hypre(original_setup()(solver,reinterpret_cast<HYPRE_ParCSRMatrix>(matrix),
                                  reinterpret_cast<HYPRE_ParVector>(rhs),reinterpret_cast<HYPRE_ParVector>(solution)));
        }
        ++audit.setups; audit.empty_setups+=(local==0);
        audit.level=level; audit.rows=rows; audit.nnz=total_nnz;
        if (rank==0) {
            if (whole) std::printf("BP3 root correction: %d rows, %d entries, %d ranks -> 1, %d cycles per exchange\n",
                                   rows,total_nnz,size,cycles_per_call);
            std::printf("BP3 consolidation: level %d, %d rows, %d entries, %d ranks -> 1\n",level,rows,total_nnz,size);
            auto* child=reinterpret_cast<hypre_ParAMGData*>(solver);
            std::printf("BP3 consolidated hierarchy rows:");
            for (int k=0;k<hypre_ParAMGDataNumLevels(child);++k)
                std::printf(" %lld",static_cast<long long>(hypre_ParCSRMatrixGlobalNumRows(hypre_ParAMGDataAArray(child)[k])));
            std::printf("\n");
        }
    }

    void bind(double* f,double* u)
    {
        bound=true; bound_f=f; bound_u=u;
        auto persistent=[&](double* data,int n,int peer,int tag,bool receive,std::vector<MPI_Request>& list) {
            MPI_Request request;
            if (receive) mpi(MPI_Recv_init(data,n,MPI_DOUBLE,peer,tag,comm,&request));
            else mpi(MPI_Send_init(data,n,MPI_DOUBLE,peer,tag,comm,&request));
            list.push_back(request);
        };
        if (rank==0) {
            auto* f_root=hypre_VectorData(hypre_ParVectorLocalVector(rhs));
            auto* u_root=hypre_VectorData(hypre_ParVectorLocalVector(solution));
            for (int peer=1;peer<size;++peer) if (counts[peer]) {
                persistent(f_root+offsets[peer],counts[peer],peer,20,true,gather);
                persistent(u_root+offsets[peer],counts[peer],peer,21,false,scatter);
            }
        } else if (local) {
            persistent(f,local,0,20,false,gather);
            persistent(u,local,0,21,true,scatter);
        }
    }
    void apply(hypre_ParVector* f,hypre_ParVector* u)
    {
        check(hypre_ParVectorFirstIndex(f)==first && hypre_ParVectorFirstIndex(u)==first &&
              hypre_ParVectorGlobalSize(f)==rows && hypre_ParVectorGlobalSize(u)==rows,
              "coarse vector ownership changed");
        check(hypre_ParVectorNumVectors(f)==1 && hypre_ParVectorNumVectors(u)==1,"multiple vectors are unsupported");
        check(hypre_VectorSize(hypre_ParVectorLocalVector(f))==local &&
              hypre_VectorSize(hypre_ParVectorLocalVector(u))==local,"coarse vector size changed");
        auto* fd=hypre_VectorData(hypre_ParVectorLocalVector(f));
        auto* ud=hypre_VectorData(hypre_ParVectorLocalVector(u));
        // The assembled p=1 interface may wrap different Krylov vectors on each call.
        if (whole && bound && (fd != bound_f || ud != bound_u)) {
            for (auto& request : gather) mpi(MPI_Request_free(&request));
            for (auto& request : scatter) mpi(MPI_Request_free(&request));
            gather.clear(); scatter.clear(); bound = false;
        }
        if (!bound) bind(fd,ud);
        check(fd==bound_f && ud==bound_u,"coarse vector storage changed after request binding");
        auto start=[](std::vector<MPI_Request>& r) {if (!r.empty()) mpi(MPI_Startall(int(r.size()),r.data()));};
        auto wait=[](std::vector<MPI_Request>& r) {if (!r.empty()) mpi(MPI_Waitall(int(r.size()),r.data(),MPI_STATUSES_IGNORE));};
        sync();
        start(gather);
        if (rank!=0) start(scatter);
        if (rank==0 && local) copy(hypre_VectorData(hypre_ParVectorLocalVector(rhs))+first,fd,local*sizeof(double));
        wait(gather);
        if (rank==0) {
            // A complete root solve keeps all its cycles between one gather and scatter.
            hypre(hypre_ParVectorSetConstantValues(solution,0));
            hypre(HYPRE_BoomerAMGSolve(solver,reinterpret_cast<HYPRE_ParCSRMatrix>(matrix),
                                     reinterpret_cast<HYPRE_ParVector>(rhs),reinterpret_cast<HYPRE_ParVector>(solution)));
            sync();
            start(scatter);
            if (local) copy(ud,hypre_VectorData(hypre_ParVectorLocalVector(solution))+first,local*sizeof(double));
        }
        wait(scatter);
        audit.cycles += cycles_per_call;
        if (whole) ++audit.root_calls;
    }
    ~State()
    {
        sync();
        for (auto& r:gather) mpi(MPI_Request_free(&r));
        for (auto& r:scatter) mpi(MPI_Request_free(&r));
        if (solver) hypre(original_destroy()(solver));
        hypre(hypre_ParVectorDestroy(rhs)); hypre(hypre_ParVectorDestroy(solution));
        hypre(hypre_ParCSRMatrixDestroy(matrix));
        if (comm!=MPI_COMM_NULL) mpi(MPI_Comm_free(&comm));
        ++audit.destroyed;
    }
};

std::map<hypre_ParAMGData*,std::unique_ptr<State>> states;
std::map<hypre_ParAMGData*,int> limits;

void end(HYPRE_Solver solver)
{
    auto* data=reinterpret_cast<hypre_ParAMGData*>(solver);
    if (states.erase(data)) hypre_ParAMGDataParticipate(data)=0;
    limits.erase(data);
    audit.write();
}
void begin(HYPRE_Solver solver, int limit)
{
    end(solver);
    auto* data=reinterpret_cast<hypre_ParAMGData*>(solver);
    check(limit==0 || (limit>9 && limit<=50000) || limit==root_row_limit,"invalid selected coarse limit");
    audit.threshold_rows=limit;
    limits[data]=limit;
    hypre(HYPRE_BoomerAMGSetSeqThreshold(solver,limit));
    if (!limit) return;
    check(limit>hypre_ParAMGDataMaxCoarseSize(data),"threshold must exceed the coarsest-row limit");
    check(hypre_ParAMGDataCoarsenType(data)==8 && hypre_ParAMGDataInterpType(data)==6 &&
          hypre_ParAMGDataPMaxElmts(data)==4 && hypre_ParAMGDataAggNumLevels(data)==0 &&
          hypre_ParAMGDataNumFunctions(data)==1 && hypre_ParAMGDataCycleType(data)==1 &&
          hypre_ParAMGDataUserRelaxType(data)==18 && hypre_ParAMGDataSmoothNumLevels(data)==0 &&
          hypre_ParAMGDataRelaxOrder(data)==0 && hypre_ParAMGDataMaxLevels(data)==25 &&
          hypre_ParAMGDataTruncFactor(data)==0 && hypre_ParAMGDataRestriction(data)==0,
          "unsupported AMG settings for the bounded experiment");
    auto* types=hypre_ParAMGDataGridRelaxType(data);
    auto* sweeps=hypre_ParAMGDataNumGridSweeps(data);
    check(types[1]==18 && types[2]==18 && types[3]==18 && sweeps[1]==1 && sweeps[2]==1 && sweeps[3]==4,
          "AMG smoothing differs from the validated BP3 recipe");
}
} // namespace bp3

extern "C" HYPRE_Int hypre_seqAMGSetup(hypre_ParAMGData* data,HYPRE_Int level,HYPRE_Int)
{
    try {
    auto configured=bp3::limits.find(data);
    bp3::check(configured!=bp3::limits.end() && configured->second>0,"coarse setup without explicit opt-in");
    bp3::check(bp3::states.count(data)==0,"duplicate coarse setup");
    auto* matrix=hypre_ParAMGDataAArray(data)[level];
    if (hypre_ParCSRMatrixGlobalNumRows(matrix)>configured->second) {
        // Hypre also calls this hook when coarsening stalls above the cutoff.
        // Retain its distributed coarse relaxation rather than exceed root memory.
        hypre_ParAMGDataParticipate(data)=0;
        ++bp3::audit.skipped_setups;
        int rank=0; bp3::mpi(MPI_Comm_rank(hypre_ParCSRMatrixComm(matrix),&rank));
        if (rank==0) std::printf("BP3 consolidation skipped: coarsening stopped above %d rows\n",configured->second);
        return 0;
    }
    bp3::states.emplace(data,std::unique_ptr<bp3::State>(new bp3::State(data,level,configured->second)));
    hypre_ParAMGDataParticipate(data)=1;
    return 0;
    } catch (...) { bp3::fail("exception in coarse setup"); }
}
extern "C" HYPRE_Int hypre_seqAMGCycle(hypre_ParAMGData* data,HYPRE_Int level,
                                      hypre_ParVector** f,hypre_ParVector** u)
{
    try {
    auto found=bp3::states.find(data);
    bp3::check(found!=bp3::states.end() && found->second->level==level,"missing coarse state");
    found->second->apply(f[level],u[level]);
    return 0;
    } catch (...) { bp3::fail("exception in coarse cycle"); }
}

extern "C" HYPRE_Int mars_bp3_amg_setup(HYPRE_Solver solver,HYPRE_ParCSRMatrix a,
                                        HYPRE_ParVector f,HYPRE_ParVector u,int requested)
{
    try {
    int ranks=0;
    bp3::mpi(MPI_Comm_size(hypre_ParCSRMatrixComm(reinterpret_cast<hypre_ParCSRMatrix*>(a)),&ranks));
    bp3::begin(solver,bp3::selected_threshold(requested,ranks));
    if (requested == -2 && ranks > 1) {
        auto* data = reinterpret_cast<hypre_ParAMGData*>(solver);
        auto* matrix = reinterpret_cast<hypre_ParCSRMatrix*>(a);
        // Assemble the same matrix on one device before constructing any AMG level.
        bp3::states.emplace(data,std::make_unique<bp3::State>(data,0,bp3::root_row_limit,matrix));
        return 0;
    }
    return bp3::original_setup()(solver,a,f,u);
    } catch (...) { bp3::fail("exception in explicit coarse setup"); }
}
extern "C" HYPRE_Int mars_bp3_amg_solve(HYPRE_Solver solver,HYPRE_ParCSRMatrix a,
                                        HYPRE_ParVector f,HYPRE_ParVector u)
{
    try {
    auto* data = reinterpret_cast<hypre_ParAMGData*>(solver);
    auto found = bp3::states.find(data);
    if (found != bp3::states.end() && found->second->whole) {
        bp3::check(found->second->source == reinterpret_cast<hypre_ParCSRMatrix*>(a),
                   "root matrix changed without setup");
        found->second->apply(reinterpret_cast<hypre_ParVector*>(f),reinterpret_cast<hypre_ParVector*>(u));
        return 0;
    }
    bp3::hypre(HYPRE_BoomerAMGSolve(solver,a,f,u));
    return 0;
    } catch (...) { bp3::fail("exception in explicit coarse solve"); }
}

extern "C" HYPRE_Int mars_bp3_amg_destroy(HYPRE_Solver solver)
{
    try {
    bp3::end(solver);
    return bp3::original_destroy()(solver);
    } catch (...) { bp3::fail("exception in explicit coarse destruction"); }
}
