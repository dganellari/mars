#pragma once
// Owned-row adapter: a rank's assembled block CSR -> the scalar CSR, RHS and
// column map that HypreGMRESSolver's device-map solve overload consumes.
//
// Identities (never interchange them):
//   local node      index into this rank's runtime arrays (ElementDomain SFC-local
//                   or runner order); owned and ghost nodes may be interleaved.
//   solver node     contiguous global id; rank r owns [first_r, first_r+owned_r),
//                   first_r = exclusive prefix sum of owned counts in rank order.
//   solver DOF      C*solver_node + component (node-major, matches setPointBlock(C)).
//   local column    C*local_node + component; A.colIndices holds these, and the
//                   separate solver_dof_map() translates them to solver DOFs.
//   source id       Exodus/public-file node id; I/O and comparisons only, never here.
//
// Contract: every owned block row must already hold its complete assembled
// contributions (complete element stars, see the halo-completion work). This adapter
// drops ghost rows and never repairs missing elements or couplings.
//
// Device-first: structure, maps and scratch are built once on the device and reused;
// update() is one fused gather, residual() one fused row-residual pass plus a fixed-order
// reduction. Host transfers are scalars only (see the README next to the gates).
// Every validation ends in one collective, so all ranks throw together.
#include "mars_segregated_assembly.hpp"
#include <mpi.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__CUDACC__)
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#define MARS_DMATRIX_HD __host__ __device__
#else
#include <numeric>
#define MARS_DMATRIX_HD
#endif

namespace mars::segregated::distributed {

enum Fault : int {
    owned_node_out_of_range = 1<<0, duplicate_owned_node = 1<<1, noncontiguous_owned_rows = 1<<2,
    missing_column_id = 1<<3, column_out_of_range = 1<<4, ghost_id_in_owned_range = 1<<5,
    solver_id_out_of_range = 1<<6, missing_diagonal = 1<<7, capacity = 1<<8, overflow = 1<<9,
    empty_rank_rejected = 1<<10, nonfinite = 1<<11, structure_changed = 1<<12,
    caller_error = 1<<13, device_error = 1<<14, values_not_updated = 1<<15
};
inline std::string describe(int faults) {
    static const char* names[] = {"owned node out of range","duplicate owned node",
        "owned rows not contiguous in solver order","owned row references a node without solver id",
        "column outside local nodes","ghost solver id inside own owned range","solver id out of range",
        "owned row without diagonal block","capacity","index overflow for the Hypre wrapper",
        "rank without owned rows (policy reject)","nonfinite owned value or RHS",
        "block graph changed since build","caller-reported error","CUDA error","residual before any update"};
    std::string text;
    for (int bit=0;bit<16;++bit) if (faults&(1<<bit)) text+=(text.empty()?"":"; ")+std::string(names[bit]);
    return text.empty()?"none":text;
}

// Zero-owned-row ranks: HypreGMRESSolver has no explicit empty-range handling, so the
// default rejects them collectively until the Daint probe proves the wrapper path.
enum class EmptyRanks { reject, allow };
struct Tolerance { double absolute=1e-13, relative=1e-10; };
struct ResidualNorms {
    double residual2=0, rhs2=0;
    bool finite=false, passed=false;
    double absolute() const { return std::sqrt(residual2); }
    double relative() const { return absolute()/std::max(std::sqrt(rhs2),1e-300); }
};
// The caller certifies that every ghost entry was refreshed from its owner after the
// last change of owned entries (e.g. ElementDomain::exchangeNodeHaloBlock(x, C)).
struct HaloComplete { const double* values; std::size_t size; };
inline HaloComplete halo_complete(const double* values,std::size_t size) { return {values,size}; }
struct HypreRows { int begin=0, end=0, column_begin=0, column_end=0; };
struct alignas(16) SquareSums { double residual2, rhs2; };

MARS_DMATRIX_HD inline void raise_fault(int* status,int fault) {
#if defined(__CUDA_ARCH__)
    atomicOr(status,fault);
#else
    *status|=fault;
#endif
}
MARS_DMATRIX_HD inline bool finite_value(double v) { return v==v && v<=1.7976931348623157e308 && v>=-1.7976931348623157e308; }

#if defined(__CUDACC__)
using Stream=cudaStream_t;
template<class T> using Buffer=thrust::device_vector<T>;
inline bool cuda_ok(cudaError_t error) { return error==cudaSuccess; }
template<class F> __global__ void dmatrix_apply(int n,F f) { const int i=blockIdx.x*blockDim.x+threadIdx.x; if (i<n) f(i); }
template<class F> void apply(int n,F f,Stream s,int& faults) {
    if (n<=0) return;
    dmatrix_apply<<<(n+255)/256,256,0,s>>>(n,f);
    if (!cuda_ok(cudaGetLastError())) faults|=device_error;
}
template<class T> T* raw(Buffer<T>& b) { return thrust::raw_pointer_cast(b.data()); }
template<class T> const T* raw(const Buffer<T>& b) { return thrust::raw_pointer_cast(b.data()); }
template<class T> T fetch(const T* p,Stream s,int& faults) {
    T value{};
    if (!cuda_ok(cudaMemcpyAsync(&value,p,sizeof(T),cudaMemcpyDeviceToHost,s)) || !cuda_ok(cudaStreamSynchronize(s)))
        faults|=device_error;
    return value;
}
inline void clear_status(int* p,Stream s,int& faults) { if (!cuda_ok(cudaMemsetAsync(p,0,sizeof(int),s))) faults|=device_error; }
template<class T> void copy_in(Buffer<T>& b,const T* p,std::size_t n) {
    b.assign(thrust::device_pointer_cast(p),thrust::device_pointer_cast(p)+n);
}
#else
using Stream=int;
template<class T> using Buffer=std::vector<T>;
template<class F> void apply(int n,F f,Stream,int&) { for (int i=0;i<n;++i) f(i); }
template<class T> T* raw(Buffer<T>& b) { return b.data(); }
template<class T> const T* raw(const Buffer<T>& b) { return b.data(); }
template<class T> T fetch(const T* p,Stream,int&) { return *p; }
inline void clear_status(int* p,Stream,int&) { *p=0; }
template<class T> void copy_in(Buffer<T>& b,const T* p,std::size_t n) { b.assign(p,p+n); }
#endif

inline bool checked_product(long long a,long long b,long long limit,long long& out) {
    if (a<0 || b<0 || (a && b>limit/a)) return false;
    out=a*b; return out<=limit;
}

namespace kernels {
struct MarkOwned {
    const int* owned; int nodes; int* slot; int* status;
    MARS_DMATRIX_HD void operator()(int k) const {
        const int v=owned[k];
        if (v<0 || v>=nodes) { raise_fault(status,owned_node_out_of_range); return; }
#if defined(__CUDA_ARCH__)
        if (atomicCAS(slot+v,0,k+1)!=0) raise_fault(status,duplicate_owned_node);
#else
        if (slot[v]) raise_fault(status,duplicate_owned_node); else slot[v]=k+1;
#endif
    }
};
template<class GlobalId> struct CheckIds {
    const GlobalId* solver_node; const int* slot; long long first, owned, total; int* status;
    MARS_DMATRIX_HD void operator()(int v) const {
        const long long g=static_cast<long long>(solver_node[v]);
        if (g<-1 || g>=total) raise_fault(status,solver_id_out_of_range);
        if (slot[v]) { if (g!=first+slot[v]-1) raise_fault(status,noncontiguous_owned_rows); }
        else if (g>=first && g<first+owned) raise_fault(status,ghost_id_in_owned_range);
    }
};
template<int C,class GlobalId> struct RowLengths {
    BlockCsrView<C> graph; const int* owned; const GlobalId* solver_node; long long* length; int* status;
    MARS_DMATRIX_HD void operator()(int k) const {
        const int v=owned[k];
        if (v<0 || v>=graph.nodes) return;
        const int begin=graph.offsets[v], end=graph.offsets[v+1];
        bool diagonal=false;
        for (int b=begin;b<end;++b) {
            const int u=graph.columns[b];
            if (u<0 || u>=graph.nodes) { raise_fault(status,column_out_of_range); continue; }
            if (solver_node[u]<0) raise_fault(status,missing_column_id);
            diagonal=diagonal || u==v;
        }
        if (!diagonal || end<begin) raise_fault(status,missing_diagonal);
        for (int c=0;c<C;++c) length[C*k+c]=end>begin?C*(long long)(end-begin):0;
    }
};
struct NarrowOffsets {
    const long long* wide; int* offsets;
    MARS_DMATRIX_HD void operator()(int i) const { offsets[i]=int(wide[i]); }
};
// Scalar row C*k+c keeps block order and expands components in place: local column
// C*u+j, value block[C*C*b+C*c+j]. Written once; update() only gathers values.
template<int C> struct FillRows {
    BlockCsrView<C> graph; const int* owned; const int* offsets; int* columns; int* source;
    MARS_DMATRIX_HD void operator()(int row) const {
        const int k=row/C, c=row%C, v=owned[k];
        int out=offsets[row];
        for (int b=graph.offsets[v];b<graph.offsets[v+1];++b)
            for (int j=0;j<C;++j) { columns[out]=C*graph.columns[b]+j; source[out]=C*C*b+C*c+j; ++out; }
    }
};
template<int C,class GlobalId> struct DofMap {
    const GlobalId* solver_node; GlobalId* map;
    MARS_DMATRIX_HD void operator()(int v) const {
        const GlobalId g=solver_node[v];
        for (int j=0;j<C;++j) map[C*v+j]=g<0?GlobalId(-1):GlobalId(C*g+j);
    }
};
// One fused pass: owned values (count entries) then packed RHS (rows entries).
template<int C> struct Gather {
    const double *blocks,*local_rhs; const int *source,*owned; double *values,*rhs; int count; int* status;
    MARS_DMATRIX_HD void operator()(int i) const {
        double v;
        if (i<count) { v=blocks[source[i]]; values[i]=v; }
        else { const int row=i-count; v=local_rhs[C*owned[row/C]+row%C]; rhs[row]=v; }
        if (!finite_value(v)) raise_fault(status,nonfinite);
    }
};
template<int C> struct Unpack {
    const double* x; const int* owned; double* local;
    MARS_DMATRIX_HD void operator()(int row) const { local[C*owned[row/C]+row%C]=x[row]; }
};

#if defined(__CUDACC__)
// W lanes per row (coalesced column/value reads); fixed grid and a fixed-order final
// pass make the reported norms reproducible run to run.
constexpr int residual_threads=256, residual_blocks=1024;
template<int W>
__global__ void owned_residual_partials(int rows,const int* offsets,const int* columns,const double* values,
    const double* x,const double* rhs,SquareSums* partial)
{
    double r2=0, b2=0;
    const int lane=threadIdx.x&(W-1);
    for (long long start=(long long)blockIdx.x*blockDim.x;start<(long long)rows*W;start+=(long long)gridDim.x*blockDim.x) {
        const long long t=start+threadIdx.x; const int row=int(t/W); const bool active=row<rows;
        double sum=0;
        if (active) for (int k=offsets[row]+lane;k<offsets[row+1];k+=W) sum+=values[k]*x[columns[k]];
        for (int o=W/2;o>0;o/=2) sum+=__shfl_down_sync(0xffffffffu,sum,o,W);
        if (active && lane==0) { const double r=sum-rhs[row]; r2+=r*r; b2+=rhs[row]*rhs[row]; }
    }
    __shared__ double s_r[residual_threads/32], s_b[residual_threads/32];
    for (int o=16;o>0;o/=2) { r2+=__shfl_down_sync(0xffffffffu,r2,o); b2+=__shfl_down_sync(0xffffffffu,b2,o); }
    if ((threadIdx.x&31)==0) { s_r[threadIdx.x/32]=r2; s_b[threadIdx.x/32]=b2; }
    __syncthreads();
    if (threadIdx.x==0) {
        SquareSums total{0,0};
        for (int w=0;w<residual_threads/32;++w) { total.residual2+=s_r[w]; total.rhs2+=s_b[w]; }
        partial[blockIdx.x]=total;
    }
}
template<int Threads>
__global__ void owned_residual_finish(int count,const SquareSums* partial,SquareSums* result)
{
    __shared__ double s_r[Threads], s_b[Threads];
    const int t=int(threadIdx.x);
    double r2=0, b2=0;
    for (int i=t;i<count;i+=Threads) { r2+=partial[i].residual2; b2+=partial[i].rhs2; }
    s_r[t]=r2; s_b[t]=b2; __syncthreads();
    for (int o=Threads/2;o>0;o/=2) {
        if (t<o) { s_r[t]+=s_r[t+o]; s_b[t]+=s_b[t+o]; }
        __syncthreads();
    }
    if (t==0) *result={s_r[0],s_b[0]};
}
#endif
} // namespace kernels

// Matrix: mars::fem::SparseMatrix<int,double,cstone::GpuTag> in production (any type with
// allocate(rows,cols,nnz), rowOffsetsPtr(), colIndicesPtr(), valuesPtr()). GlobalId must be
// HYPRE_BigInt for the wrapper's device-map overload. Pointers passed in are device pointers
// in a CUDA build and host pointers otherwise.
template<int C,class Matrix,class GlobalId=std::int64_t>
class OwnedRowSystem {
    static_assert(C==1 || C==3,"segregated SIMPLE solves scalar pressure or three velocity components");
public:
    // graph: the persistent block graph (values/rhs unused here). owned: owned local nodes in
    // solver order. solver_node: solver node per local node, -1 only for nodes no owned row uses.
    OwnedRowSystem(MPI_Comm comm,BlockCsrView<C> graph,const int* owned,int owned_count,
        const GlobalId* solver_node,std::size_t solver_node_size,EmptyRanks empty=EmptyRanks::reject,
        long long row_limit=std::numeric_limits<int>::max(),Stream stream={})
        : comm_(comm), stream_(stream), offsets_(graph.offsets), columns_(graph.columns), nodes_(graph.nodes),
          owned_count_(owned_count), status_(1,0)
    {
        int local=0;
        long long limit=std::min<long long>(row_limit,std::numeric_limits<int>::max());
        limit=std::min<long long>(limit,(long long)std::numeric_limits<GlobalId>::max());
        long long scalar_columns=0, block_values=0;
        if (nodes_<0 || owned_count_<0 || owned_count_>nodes_ || (nodes_>0 && (!offsets_ || !columns_))
            || (owned_count_>0 && !owned) || (nodes_>0 && !solver_node)) local|=capacity;
        if (solver_node_size<std::size_t(std::max(nodes_,0))) local|=capacity;
        if (!checked_product(C,std::max(nodes_,0),std::numeric_limits<int>::max(),scalar_columns)) local|=overflow;
        if (!(local&capacity) && nodes_>0 && !checked_product(C*C,fetch(offsets_+nodes_,stream_,local),
                                                              std::numeric_limits<int>::max(),block_values)) local|=overflow;
        const long long mine=owned_count_<0?0:owned_count_;
        MPI_Exscan(&mine,&first_,1,MPI_LONG_LONG,MPI_SUM,comm_);
        int rank=0; MPI_Comm_rank(comm_,&rank); if (rank==0) first_=0;
        MPI_Allreduce(&mine,&total_,1,MPI_LONG_LONG,MPI_SUM,comm_);
        long long rows_end=0;
        if (!checked_product(C,total_,limit,rows_end)) local|=overflow;
        if (owned_count_==0 && empty==EmptyRanks::reject) local|=empty_rank_rejected;
        rows_=local&(capacity|overflow)?0:C*owned_count_;
        Buffer<long long> length(std::size_t(rows_)+1,0);
        if (!(local&(capacity|overflow))) {
            Buffer<int> slot(nodes_,0);
            apply(owned_count_,kernels::MarkOwned{owned,nodes_,raw(slot),raw(status_)},stream_,local);
            apply(nodes_,kernels::CheckIds<GlobalId>{solver_node,raw(slot),first_,mine,total_,raw(status_)},stream_,local);
            apply(owned_count_,kernels::RowLengths<C,GlobalId>{graph,owned,solver_node,raw(length),raw(status_)},stream_,local);
#if defined(__CUDACC__)
            thrust::exclusive_scan(thrust::cuda::par.on(stream_),length.begin(),length.end(),length.begin());
#else
            std::exclusive_scan(length.begin(),length.end(),length.begin(),0LL);
#endif
            const long long entries=fetch(raw(length)+rows_,stream_,local);
            // The fused gather indexes values and RHS in one int range.
            if (entries+rows_>std::numeric_limits<int>::max()) local|=overflow; else nnz_=int(entries);
            local|=fetch(raw(status_),stream_,local);
        }
        collective(local,"distributed matrix build");
        copy_in(owned_,owned,std::size_t(owned_count_));
        matrix_.allocate(rows_,int(scalar_columns),nnz_);
        source_.resize(nnz_); map_.resize(std::size_t(scalar_columns)); partial_.resize(partial_count); result_.resize(1);
        apply(rows_+1,kernels::NarrowOffsets{raw(length),matrix_.rowOffsetsPtr()},stream_,local);
        apply(rows_,kernels::FillRows<C>{graph,raw(owned_),matrix_.rowOffsetsPtr(),matrix_.colIndicesPtr(),raw(source_)},stream_,local);
        apply(nodes_,kernels::DofMap<C,GlobalId>{solver_node,raw(map_)},stream_,local);
        synchronize(local);
        collective(local,"distributed matrix structure");
    }
    OwnedRowSystem(const OwnedRowSystem&)=delete;
    OwnedRowSystem& operator=(const OwnedRowSystem&)=delete;

    // Gather fresh owned values and RHS from the same graph. No accumulation: every value
    // is overwritten. One status word crosses to the host, then one integer allreduce.
    void update(BlockCsrView<C> assembled,double* rhs,std::size_t rhs_capacity,bool caller_failed=false) {
        int local=caller_failed?caller_error:0;
        if (assembled.offsets!=offsets_ || assembled.columns!=columns_ || assembled.nodes!=nodes_) local|=structure_changed;
        if (rhs_capacity<std::size_t(rows_) || (rows_>0 && (!rhs || !assembled.values || !assembled.rhs))) local|=capacity;
        if (!(local&(structure_changed|capacity))) {
            clear_status(raw(status_),stream_,local);
            apply(nnz_+rows_,kernels::Gather<C>{assembled.values,assembled.rhs,raw(source_),raw(owned_),
                matrix_.valuesPtr(),rhs,nnz_,raw(status_)},stream_,local);
            local|=fetch(raw(status_),stream_,local);
        }
        collective(local,"distributed matrix update");
        updated_=true;
    }
    // Owned solution entries back to local nodes; ghosts are untouched (exchange next).
    // Errors are recorded and reported by the next residual() collective.
    void unpack(const double* x,std::size_t x_size,double* local,std::size_t local_size) {
        if (x_size<std::size_t(rows_) || local_size<std::size_t(C)*std::size_t(nodes_) || (rows_>0 && (!x || !local)))
            deferred_|=capacity;
        else apply(rows_,kernels::Unpack<C>{x,raw(owned_),local},stream_,deferred_);
    }
    // True residual over owned rows: global sums of r^2 and b^2, then one decision on all ranks.
    ResidualNorms residual(HaloComplete x,const double* rhs,Tolerance tolerance={}) {
        int local=deferred_|(updated_?0:values_not_updated); deferred_=0;
        if (x.size<std::size_t(C)*std::size_t(nodes_) || (rows_>0 && (!x.values || !rhs))) local|=capacity;
        SquareSums sums{0,0};
        if (!(local&capacity) && rows_>0) {
#if defined(__CUDACC__)
            constexpr int lanes=C==1?4:8;
            const long long threads=(long long)rows_*lanes;
            const int blocks=int(std::min<long long>(partial_count,(threads+kernels::residual_threads-1)/kernels::residual_threads));
            kernels::owned_residual_partials<lanes><<<blocks,kernels::residual_threads,0,stream_>>>(rows_,
                matrix_.rowOffsetsPtr(),matrix_.colIndicesPtr(),matrix_.valuesPtr(),x.values,rhs,raw(partial_));
            kernels::owned_residual_finish<kernels::residual_threads><<<1,kernels::residual_threads,0,stream_>>>(blocks,raw(partial_),raw(result_));
            if (!cuda_ok(cudaGetLastError())) local|=device_error;
            sums=fetch(raw(result_),stream_,local);
#else
            const int* offsets=matrix_.rowOffsetsPtr(); const int* columns=matrix_.colIndicesPtr(); const double* values=matrix_.valuesPtr();
            for (int row=0;row<rows_;++row) {
                double sum=0; for (int k=offsets[row];k<offsets[row+1];++k) sum+=values[k]*x.values[columns[k]];
                const double r=sum-rhs[row]; sums.residual2+=r*r; sums.rhs2+=rhs[row]*rhs[row];
            }
#endif
        }
        double reduced[3]={sums.residual2,sums.rhs2,local?1.:0.}, global[3]={};
        MPI_Allreduce(reduced,global,3,MPI_DOUBLE,MPI_SUM,comm_);
        if (global[2]>0) collective(local,"distributed residual"); // failure path only: name the faults
        ResidualNorms norms; norms.residual2=global[0]; norms.rhs2=global[1];
        norms.finite=std::isfinite(global[0]) && std::isfinite(global[1]);
        norms.passed=norms.finite && norms.absolute()<=tolerance.absolute+tolerance.relative*std::sqrt(global[1]);
        return norms;
    }
    // Arguments of HypreGMRESSolver::solve(A,b,x,begin,end,column_begin,column_end,map).
    HypreRows hypre_rows() const { return {int(C*first_),int(C*(first_+owned_count_)),0,int(C*total_)}; }
    const Matrix& matrix() const { return matrix_; }
    const Buffer<GlobalId>& solver_dof_map() const { return map_; }
    int rows() const { return rows_; }
    int nodes() const { return nodes_; }
    int nnz() const { return nnz_; }
    long long first_solver_node() const { return first_; }
    long long solver_nodes() const { return total_; }
    void synchronize(int& faults) const {
#if defined(__CUDACC__)
        if (!cuda_ok(cudaStreamSynchronize(stream_)) || !cuda_ok(cudaGetLastError())) faults|=device_error;
#else
        (void)faults;
#endif
    }
private:
    static constexpr int partial_count=1024;
    void collective(int local,const char* stage) const {
        int global=0; MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,comm_);
        if (global) fail(local,global,stage);
    }
    [[noreturn]] void fail(int local,int global,const char* stage) const {
        int rank=0; MPI_Comm_rank(comm_,&rank);
        throw std::runtime_error(std::string(stage)+" rejected on all ranks (global: "+describe(global)
            +"; rank "+std::to_string(rank)+": "+describe(local)+")");
    }
    MPI_Comm comm_; Stream stream_;
    const int *offsets_, *columns_;
    int nodes_, owned_count_, rows_=0, nnz_=0, deferred_=0;
    bool updated_=false;
    long long first_=0, total_=0;
    Matrix matrix_;
    Buffer<int> owned_, source_, status_;
    Buffer<GlobalId> map_;
    Buffer<SquareSums> partial_, result_;
};

// The real wrapper call; Solver is mars::fem::HypreGMRESSolver<double,int,cstone::GpuTag>.
template<class Solver,class System,class Vector>
bool solve_owned(Solver& solver,const System& system,const Vector& rhs,Vector& x) {
    const auto r=system.hypre_rows();
    return solver.solve(system.matrix(),rhs,x,r.begin,r.end,r.column_begin,r.column_end,system.solver_dof_map());
}
} // namespace mars::segregated::distributed
#undef MARS_DMATRIX_HD
