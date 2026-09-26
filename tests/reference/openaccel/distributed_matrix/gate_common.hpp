#pragma once
// Shared gates for the owned-row adapter. The same code runs on host buffers (CPU MPI)
// and device buffers (CUDA-aware MPI); only Buffer/apply differ. Fixtures are generated,
// replicated on every rank (test-only), and keyed by nontrivial source ids.
#include "mars_segregated_distributed_matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#if defined(__CUDACC__)
#define GATE_HD __host__ __device__
#else
#define GATE_HD
#endif

namespace dmatrix_gate {
using namespace mars::segregated::distributed;
using mars::segregated::BlockCsrView;

inline std::uint64_t mix(std::uint64_t x) {
    x+=0x9e3779b97f4a7c15ull; x=(x^(x>>30))*0xbf58476d1ce4e5b9ull; x=(x^(x>>27))*0x94d049bb133111ebull;
    return x^(x>>31);
}
inline std::uint64_t key(std::uint64_t a,std::uint64_t b,std::uint64_t c,std::uint64_t d) { return mix(mix(mix(mix(a)^b)^c)^d); }
inline double symmetric(std::uint64_t k) { return double(mix(k)>>11)*0x1.0p-52-1.0; } // [-1,1)
template<class T> void shuffle(std::vector<T>& v,std::uint64_t seed) {
    for (std::size_t i=v.size();i>1;--i) std::swap(v[i-1],v[mix(seed^(i*0x100000001b3ull))%i]);
}

struct Options {
    int nx=6, ny=5, nz=4;          // 120 nodes, 27-point couplings
    bool empty_rank=false;         // rank 1 owns no rows (needs >= 2 ranks)
    bool empty_rank_holds_nothing=false;
    bool far_ghosts=true;          // partial second ghost layer, never used by owned rows
    bool far_ghosts_without_id=false;
    bool omit_coupling=false;      // last rank drops its largest ghost coupling from one owned row
};

// One global nonsymmetric, diagonally dominant block system and its partition.
struct Problem {
    int C, ranks, nodes;
    Options options;
    std::vector<std::vector<int>> neighbors, owned_order;
    std::vector<int> owner, weight, solver_to_global;
    std::vector<long long> first, solver_node, source_id;
    Problem(int components,int rank_count,Options o):C(components),ranks(rank_count),nodes(o.nx*o.ny*o.nz),options(o) {
        neighbors.resize(nodes);
        for (int z=0;z<o.nz;++z) for (int y=0;y<o.ny;++y) for (int x=0;x<o.nx;++x)
            for (int dz=-1;dz<=1;++dz) for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
                const int X=x+dx, Y=y+dy, Z=z+dz;
                if (X>=0 && Y>=0 && Z>=0 && X<o.nx && Y<o.ny && Z<o.nz) neighbors[x+o.nx*(y+o.ny*z)].push_back(X+o.nx*(Y+o.ny*Z));
            }
        for (auto& n:neighbors) std::sort(n.begin(),n.end());
        std::vector<int> permutation(nodes); for (int g=0;g<nodes;++g) permutation[g]=g;
        shuffle(permutation,7); source_id.resize(nodes);
        for (int g=0;g<nodes;++g) source_id[g]=1000003+17LL*permutation[g];
        // Uneven weights, then irregular seams: every fifth node moves to the next owning rank.
        weight.resize(ranks); int total=0;
        for (int r=0;r<ranks;++r) total+=weight[r]=(o.empty_rank && r==1)?0:r+1;
        owner.resize(nodes);
        for (int g=0,r=0,acc=0;g<nodes;++g) {
            while (r<ranks-1 && (long long)g*total>=(long long)(acc+weight[r])*nodes) acc+=weight[r++];
            while (weight[r]==0) r=(r+1)%ranks;
            owner[g]=r;
            if (ranks>1 && mix(31ull*g+7)%5==0) { int q=(r+1)%ranks; while (weight[q]==0) q=(q+1)%ranks; owner[g]=q; }
        }
        owned_order.resize(ranks); first.assign(ranks+1,0); solver_node.resize(nodes); solver_to_global.resize(nodes);
        for (int g=0;g<nodes;++g) owned_order[owner[g]].push_back(g);
        for (int r=0;r<ranks;++r) {
            shuffle(owned_order[r],1000+r); first[r+1]=first[r]+(long long)owned_order[r].size();
            for (std::size_t k=0;k<owned_order[r].size();++k) {
                solver_node[owned_order[r][k]]=first[r]+(long long)k; solver_to_global[first[r]+k]=owned_order[r][k];
            }
        }
    }
    // Entries are functions of source ids, so any identity mix-up changes the value read.
    double value(int gi,int gj,int c,int j,unsigned seed) const {
        const double v=symmetric(key(std::uint64_t(source_id[gi]),std::uint64_t(source_id[gj]),std::uint64_t(8*c+j),seed));
        return gi==gj && c==j?40.0*C+1.0+0.5*v:0.5*v;
    }
    double rhs(int g,int c,unsigned seed) const { return symmetric(key(std::uint64_t(source_id[g]),77,std::uint64_t(c),seed)); }
    double solution(int g,int c,unsigned seed) const { return symmetric(key(std::uint64_t(source_id[g]),99,std::uint64_t(c),seed)); }
    double product(int g,int c,unsigned value_seed,unsigned x_seed) const {
        long double s=0;
        for (int u:neighbors[g]) for (int j=0;j<C;++j) s+=(long double)value(g,u,c,j,value_seed)*solution(u,j,x_seed);
        return double(s);
    }
    bool coupled(int g,int u) const { return std::binary_search(neighbors[g].begin(),neighbors[g].end(),u); }
    // Nodes present on rank r: owned, every coupled neighbour (complete owned rows), and far ghosts.
    std::vector<int> holds(int r,std::set<int>* far=nullptr) const {
        std::set<int> present(owned_order[r].begin(),owned_order[r].end()), ghosts, extra;
        for (int g:owned_order[r]) for (int u:neighbors[g]) if (!present.count(u)) ghosts.insert(u);
        if (options.empty_rank && weight[r]==0 && !options.empty_rank_holds_nothing)
            for (int g=0;g<nodes;++g) if (mix(std::uint64_t(g))%11==0) ghosts.insert(g);
        if (options.far_ghosts)
            for (int g:ghosts) for (int u:neighbors[g])
                if (!present.count(u) && !ghosts.count(u) && mix(std::uint64_t(u)+5)%3==0) extra.insert(u);
        present.insert(ghosts.begin(),ghosts.end()); present.insert(extra.begin(),extra.end());
        if (far) *far=extra;
        return {present.begin(),present.end()};
    }
};

// Rank-local runtime view: shuffled local order (owned and ghosts interleaved), block graph over
// every local node, owned rows complete, ghost rows partial and NaN-poisoned.
struct Local {
    std::vector<int> global, local_of, offsets, columns, owned;
    std::vector<char> referenced;
    std::vector<long long> solver_node;
    std::vector<double> blocks, rhs;
    int omitted_row=-1, omitted_column=-1;
    int nodes() const { return int(global.size()); }
};
inline Local extract(const Problem& p,int r) {
    Local l; std::set<int> far;
    l.global=p.holds(r,&far); shuffle(l.global,2000+r);
    l.local_of.assign(p.nodes,-1);
    for (int v=0;v<l.nodes();++v) l.local_of[l.global[v]]=v;
    if (p.options.omit_coupling && r==p.ranks-1)
        for (int g:p.owned_order[r]) {
            double best=-1;
            for (int u:p.neighbors[g]) if (u!=g && (p.ranks==1 || p.owner[u]!=r) && std::abs(p.value(g,u,0,0,1))>best)
                { best=std::abs(p.value(g,u,0,0,1)); l.omitted_row=g; l.omitted_column=u; }
            if (best>=0) break;
        }
    l.offsets.push_back(0); l.referenced.assign(l.nodes(),0);
    for (int v=0;v<l.nodes();++v) {
        const int g=l.global[v]; std::vector<int> row;
        for (int u:p.neighbors[g]) if (l.local_of[u]>=0 && !(g==l.omitted_row && u==l.omitted_column)) row.push_back(l.local_of[u]);
        std::sort(row.begin(),row.end());
        if (p.owner[g]==r) for (int u:row) l.referenced[u]=1;
        l.columns.insert(l.columns.end(),row.begin(),row.end()); l.offsets.push_back(int(l.columns.size()));
    }
    for (int g:p.owned_order[r]) l.owned.push_back(l.local_of[g]);
    l.solver_node.resize(l.nodes());
    for (int v=0;v<l.nodes();++v)
        l.solver_node[v]=p.options.far_ghosts_without_id && far.count(l.global[v])?-1:p.solver_node[l.global[v]];
    return l;
}
// Owned rows get the seed's values; ghost rows and ghost RHS stay NaN so any read of them is caught.
template<class Rhs> void fill(const Problem& p,Local& l,int r,unsigned seed,Rhs rhs) {
    const int C=p.C; const double poison=std::numeric_limits<double>::quiet_NaN();
    l.blocks.assign(std::size_t(C*C)*l.columns.size(),poison); l.rhs.assign(std::size_t(C)*l.nodes(),poison);
    for (int v=0;v<l.nodes();++v) {
        const int g=l.global[v]; if (p.owner[g]!=r) continue;
        for (int c=0;c<C;++c) l.rhs[C*v+c]=rhs(g,c);
        for (int b=l.offsets[v];b<l.offsets[v+1];++b)
            for (int c=0;c<C;++c) for (int j=0;j<C;++j) l.blocks[std::size_t(C*C)*b+C*c+j]=p.value(g,l.global[l.columns[b]],c,j,seed);
    }
}

#if defined(__CUDACC__)
template<class T> Buffer<T> upload(const std::vector<T>& h) { return Buffer<T>(h.begin(),h.end()); }
template<class T> void overwrite(Buffer<T>& d,const std::vector<T>& h) { thrust::copy(h.begin(),h.end(),d.begin()); }
template<class T> std::vector<T> download(const T* p,std::size_t n) {
    std::vector<T> h(n); if (n && cudaMemcpy(h.data(),p,n*sizeof(T),cudaMemcpyDeviceToHost)!=cudaSuccess) throw std::runtime_error("download");
    return h;
}
inline void device_sync() { if (cudaDeviceSynchronize()!=cudaSuccess) throw std::runtime_error("device synchronize"); }
#else
template<class T> Buffer<T> upload(const std::vector<T>& h) { return h; }
template<class T> void overwrite(Buffer<T>& d,const std::vector<T>& h) { std::copy(h.begin(),h.end(),d.begin()); }
template<class T> std::vector<T> download(const T* p,std::size_t n) { return std::vector<T>(p,p+n); }
inline void device_sync() {}
#endif

template<int C,class GlobalId> struct Device {
    int nodes;
    Buffer<int> offsets, columns, owned;
    Buffer<double> blocks, rhs;
    Buffer<GlobalId> solver_node;
    explicit Device(const Local& l):nodes(l.nodes()),offsets(upload(l.offsets)),columns(upload(l.columns)),owned(upload(l.owned)),
        blocks(upload(l.blocks)),rhs(upload(l.rhs)),solver_node(upload(std::vector<GlobalId>(l.solver_node.begin(),l.solver_node.end()))) {}
    BlockCsrView<C> view() { return {nodes,raw(offsets),raw(columns),raw(blocks),raw(rhs)}; }
};

struct Pack {
    const int* nodes; const double* x; double* buffer; int C;
    GATE_HD void operator()(int i) const { for (int c=0;c<C;++c) buffer[C*i+c]=x[C*nodes[i]+c]; }
};
struct UnpackGhosts {
    const int* nodes; const double* buffer; double* x; int C, stale;
    GATE_HD void operator()(int i) const { if (i!=stale) for (int c=0;c<C;++c) x[C*nodes[i]+c]=buffer[C*i+c]; }
};
// Test-side stand-in for ElementDomain::exchangeNodeHaloBlock(x, C): persistent buffers, one
// message per peer carrying all components, device pointers handed to (CUDA-aware) MPI.
struct GhostExchange {
    int C;
    std::vector<int> peers, send_offsets{0}, recv_offsets{0};
    Buffer<int> send_nodes, recv_nodes;
    Buffer<double> send_buffer, recv_buffer;
    std::vector<char> recv_referenced;
    GhostExchange(const Problem& p,const Local& l,int rank):C(p.C) {
        std::vector<int> send, recv;
        for (int q=0;q<p.ranks;++q) {
            if (q==rank) continue;
            const auto theirs=p.holds(q); const std::set<int> held(theirs.begin(),theirs.end());
            const std::size_t s0=send.size(), r0=recv.size();
            for (int g=0;g<p.nodes;++g) {
                if (l.local_of[g]<0) continue;
                if (p.owner[g]==rank && held.count(g) && p.owner[g]!=q) send.push_back(l.local_of[g]);
                if (p.owner[g]==q) { recv.push_back(l.local_of[g]); recv_referenced.push_back(l.referenced[l.local_of[g]]); }
            }
            if (send.size()==s0 && recv.size()==r0) continue;
            peers.push_back(q); send_offsets.push_back(int(send.size())); recv_offsets.push_back(int(recv.size()));
        }
        send_nodes=upload(send); recv_nodes=upload(recv);
        send_buffer.resize(std::size_t(C)*send.size()); recv_buffer.resize(std::size_t(C)*recv.size());
    }
    int first_referenced_recv() const {
        for (std::size_t i=0;i<recv_referenced.size();++i) if (recv_referenced[i]) return int(i);
        return -1;
    }
    void run(double* x,MPI_Comm comm,int stale=-1) {
        int faults=0; std::vector<MPI_Request> requests;
        apply(int(send_nodes.size()),Pack{raw(send_nodes),x,raw(send_buffer),C},Stream{},faults);
        device_sync();
        for (std::size_t p=0;p<peers.size();++p) {
            const int count=C*(recv_offsets[p+1]-recv_offsets[p]);
            if (count) { requests.emplace_back(); MPI_Irecv(raw(recv_buffer)+C*recv_offsets[p],count,MPI_DOUBLE,peers[p],0x4d44,comm,&requests.back()); }
        }
        for (std::size_t p=0;p<peers.size();++p) {
            const int count=C*(send_offsets[p+1]-send_offsets[p]);
            if (count) { requests.emplace_back(); MPI_Isend(raw(send_buffer)+C*send_offsets[p],count,MPI_DOUBLE,peers[p],0x4d44,comm,&requests.back()); }
        }
        MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE);
        apply(int(recv_nodes.size()),UnpackGhosts{raw(recv_nodes),raw(recv_buffer),x,C,stale},Stream{},faults);
        device_sync();
        if (faults) throw std::runtime_error("ghost exchange kernel failed");
    }
};

// Independent oracle: expected entries come from Problem's definition, keyed by global identity.
template<class GlobalId> long long oracle_mismatches(const Problem& p,const Local& l,int r,long long first,
    const std::vector<int>& offsets,const std::vector<int>& columns,const std::vector<double>& values,
    const std::vector<GlobalId>& map,const std::vector<double>& rhs,unsigned seed,
    const std::vector<double>& expected_rhs,std::string& note)
{
    const int C=p.C; long long bad=0;
    auto miss=[&](const std::string& what) { if (note.empty()) note=what; ++bad; };
    const int rows=C*int(l.owned.size());
    if (int(offsets.size())!=rows+1 || offsets.front()!=0 || offsets.back()!=int(values.size()) || values.size()!=columns.size()) { miss("CSR shape"); return bad; }
    if (first!=p.first[r]) miss("first solver node");
    if (map.size()!=std::size_t(C)*l.nodes()) miss("map size");
    else for (int v=0;v<l.nodes();++v) for (int j=0;j<C;++j) {
        const long long expect=l.solver_node[v]<0?-1:C*l.solver_node[v]+j;
        if ((long long)map[C*v+j]!=expect) miss("solver DOF map");
    }
    for (int i=0;i<rows;++i) {
        const long long R=C*first+i; const int g=p.solver_to_global[R/C], c=int(R%C);
        if (p.owner[g]!=r || l.global[l.owned[i/C]]!=g) { miss("row identity"); continue; }
        if (rhs[i]!=expected_rhs[i]) miss("RHS entry");
        std::set<std::pair<int,int>> seen;
        for (int k=offsets[i];k<offsets[i+1];++k) {
            if (columns[k]<0 || columns[k]>=C*l.nodes()) { miss("local column range"); continue; }
            const long long Q=(long long)map[columns[k]];
            if (Q<0) { miss("unmapped column"); continue; }
            const int u=p.solver_to_global[Q/C], j=int(Q%C);
            if (!p.coupled(g,u) || !seen.insert({u,j}).second) { miss("extra or duplicate entry"); continue; }
            if (values[k]!=p.value(g,u,c,j,seed)) {
                std::ostringstream s; s<<"value source("<<p.source_id[g]<<","<<c<<")x("<<p.source_id[u]<<","<<j<<")"; miss(s.str());
            }
        }
        if (seen.size()!=std::size_t(C)*p.neighbors[g].size()) {
            std::ostringstream s; s<<"missing coupling in row source "<<p.source_id[g]; miss(s.str());
            bad+=(long long)(std::size_t(C)*p.neighbors[g].size()-seen.size())-1;
        }
    }
    return bad;
}
// Independent residual: definition times the local (exchanged) vector, never the adapter CSR.
inline long double oracle_residual2(const Problem& p,const Local& l,int r,const std::vector<double>& x,
    unsigned value_seed,const std::vector<double>& expected_rhs,long long& incomplete)
{
    const int C=p.C; long double sum=0;
    for (std::size_t k=0;k<l.owned.size();++k) {
        const int g=p.owned_order[r][k];
        for (int c=0;c<C;++c) {
            long double a=0;
            for (int u:p.neighbors[g]) {
                if (l.local_of[u]<0) { ++incomplete; continue; }
                for (int j=0;j<C;++j) a+=(long double)p.value(g,u,c,j,value_seed)*x[std::size_t(C)*l.local_of[u]+j];
            }
            const long double d=a-expected_rhs[C*k+c]; sum+=d*d;
        }
    }
    return sum;
}

struct Report {
    int rank=0, failures=0, passes=0, skips=0;
    void result(const std::string& name,bool ok,const std::string& detail="") {
        if (rank==0) std::cout<<(ok?"PASS ":"FAIL ")<<name<<(detail.empty()?"":"  ["+detail+"]")<<std::endl;
        ok?++passes:++failures;
    }
    void skip(const std::string& name,const std::string& why) { if (rank==0) std::cout<<"SKIP "<<name<<"  ["<<why<<"]"<<std::endl; ++skips; }
};
inline long long sum_all(long long v,MPI_Comm comm) { long long s=0; MPI_Allreduce(&v,&s,1,MPI_LONG_LONG,MPI_SUM,comm); return s; }
inline long double sum_all(long double v,MPI_Comm comm) { double d=double(v), s=0; MPI_Allreduce(&d,&s,1,MPI_DOUBLE,MPI_SUM,comm); return s; }
inline bool all_true(bool v,MPI_Comm comm) { int in=v?1:0, out=0; MPI_Allreduce(&in,&out,1,MPI_INT,MPI_MIN,comm); return out==1; }
// Every rank must throw at the same point; a stranded rank would deadlock (ctest timeout).
template<class F> bool all_ranks_throw(MPI_Comm comm,F f,const char* expected,bool injecting,std::string& message) {
    bool threw=false;
    try { f(); } catch (const std::exception& e) { threw=true; message=e.what(); }
    const bool named=!injecting || message.find(expected)!=std::string::npos;
    return all_true(threw && named,comm);
}

// Runs one injected fault and reports it with the injecting rank's exception text.
template<class F> void expect_rejection(MPI_Comm comm,Report& report,const std::string& name,const char* fault,int injector,F f) {
    int rank=0, ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    std::string message;
    const bool ok=all_ranks_throw(comm,f,fault,rank==injector,message);
    if (injector!=0 && ranks>1) {
        if (rank==injector) { const int n=int(message.size()); MPI_Send(&n,1,MPI_INT,0,0x4d45,comm); MPI_Send(message.data(),n,MPI_CHAR,0,0x4d46,comm); }
        if (rank==0) { int n=0; MPI_Recv(&n,1,MPI_INT,injector,0x4d45,comm,MPI_STATUS_IGNORE); message.assign(std::size_t(n),' ');
                       MPI_Recv(message.data(),n,MPI_CHAR,injector,0x4d46,comm,MPI_STATUS_IGNORE); }
    }
    report.result(name,ok,"rank "+std::to_string(injector)+": "+message);
}
inline std::string sci(double v) { std::ostringstream s; s<<std::scientific<<std::setprecision(3)<<v; return s.str(); }

template<int C,class Matrix,class GlobalId> using System=OwnedRowSystem<C,Matrix,GlobalId>;
template<int C,class Matrix,class GlobalId> struct Downloaded {
    std::vector<int> offsets, columns; std::vector<double> values; std::vector<GlobalId> map;
    explicit Downloaded(const System<C,Matrix,GlobalId>& s)
        : offsets(download(s.matrix().rowOffsetsPtr(),std::size_t(s.rows())+1)),
          columns(download(s.matrix().colIndicesPtr(),std::size_t(s.nnz()))),
          values(download(s.matrix().valuesPtr(),std::size_t(s.nnz()))),
          map(download(raw(s.solver_dof_map()),s.solver_dof_map().size())) {}
};
inline std::vector<double> owned_rhs(const Problem& p,int r,unsigned seed) {
    std::vector<double> b; for (int g:p.owned_order[r]) for (int c=0;c<p.C;++c) b.push_back(p.rhs(g,c,seed)); return b;
}
inline std::vector<double> owned_product(const Problem& p,int r,unsigned value_seed,unsigned x_seed) {
    std::vector<double> b; for (int g:p.owned_order[r]) for (int c=0;c<p.C;++c) b.push_back(p.product(g,c,value_seed,x_seed)); return b;
}
inline std::vector<double> local_solution(const Problem& p,const Local& l,unsigned seed,double shift) {
    std::vector<double> x(std::size_t(p.C)*l.nodes());
    for (int v=0;v<l.nodes();++v) for (int c=0;c<p.C;++c) x[std::size_t(p.C)*v+c]=p.solution(l.global[v],c,seed)+shift;
    return x;
}

// Structural, update, residual, empty-work and fault-injection gates for one component count.
template<int C,class Matrix,class GlobalId> void run_gates(MPI_Comm comm,Report& report,Options base={}) {
    int rank=0, ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    const std::string tag="C="+std::to_string(C)+" ranks="+std::to_string(ranks)+": ";
    const int last=ranks-1; const bool injecting=rank==last;
    {   // Build, first values, oracle comparison of every entry, RHS and map.
        Problem p(C,ranks,base); Local l=extract(p,rank);
        fill(p,l,rank,1,[&](int g,int c) { return p.rhs(g,c,1); });
        Device<C,GlobalId> d(l); Buffer<double> b(C*l.owned.size());
        System<C,Matrix,GlobalId> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
        s.update(d.view(),raw(b),b.size());
        Downloaded<C,Matrix,GlobalId> first_values(s); std::string note;
        long long bad=oracle_mismatches(p,l,rank,s.first_solver_node(),first_values.offsets,first_values.columns,
            first_values.values,first_values.map,download(raw(b),b.size()),1,owned_rhs(p,rank,1),note);
        const auto rows=s.hypre_rows();
        if (rows.begin!=C*p.first[rank] || rows.end!=C*p.first[rank+1] || rows.column_end!=C*p.nodes || s.solver_nodes()!=p.nodes) ++bad;
        bad=sum_all(bad,comm);
        std::string owned_rows;
        for (int q=0;q<ranks;++q) owned_rows+=(q?"/":"")+std::to_string(p.owned_order[q].size());
        report.result(tag+"every CSR entry, RHS and DOF map match the source-keyed oracle",bad==0,
            std::to_string(sum_all((long long)first_values.values.size(),comm))+" entries, owned nodes per rank "+owned_rows
            +", ghosts "+std::to_string(sum_all((long long)(l.nodes()-int(l.owned.size())),comm))+", mismatches="+std::to_string(bad)+(note.empty()?"":" first: "+note));

        // Second update: new values and RHS, same structure, no accumulation.
        fill(p,l,rank,2,[&](int g,int c) { return p.rhs(g,c,2); }); overwrite(d.blocks,l.blocks); overwrite(d.rhs,l.rhs);
        s.update(d.view(),raw(b),b.size());
        Downloaded<C,Matrix,GlobalId> second(s); note.clear();
        long long changed=(second.offsets!=first_values.offsets)+(second.columns!=first_values.columns)+(second.map!=first_values.map);
        long long bad2=oracle_mismatches(p,l,rank,s.first_solver_node(),second.offsets,second.columns,second.values,second.map,
            download(raw(b),b.size()),2,owned_rhs(p,rank,2),note);
        changed=sum_all(changed,comm); bad2=sum_all(bad2,comm);
        report.result(tag+"second update reuses structure and overwrites values",changed==0 && bad2==0,
            "structure changes="+std::to_string(changed)+", mismatches="+std::to_string(bad2));

        // Known solution: unpack owned entries, exchange ghosts, true residual passes.
        fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,5); }); overwrite(d.blocks,l.blocks); overwrite(d.rhs,l.rhs);
        s.update(d.view(),raw(b),b.size());
        std::vector<double> xs; for (int g:p.owned_order[rank]) for (int c=0;c<C;++c) xs.push_back(p.solution(g,c,5));
        Buffer<double> x_owned=upload(xs), x_local=upload(local_solution(p,l,5,1.0)); // previous iterate: x*+1
        std::vector<double> sentinel=download(raw(x_local),x_local.size());
        s.unpack(raw(x_owned),x_owned.size(),raw(x_local),x_local.size()); device_sync();
        const auto unpacked=download(raw(x_local),x_local.size()); long long unpack_bad=0;
        for (int v=0;v<l.nodes();++v) for (int c=0;c<C;++c) {
            const std::size_t i=std::size_t(C)*v+c;
            const bool mine=p.owner[l.global[v]]==rank;
            unpack_bad+=mine?unpacked[i]!=p.solution(l.global[v],c,5):unpacked[i]!=sentinel[i];
        }
        unpack_bad=sum_all(unpack_bad,comm);
        report.result(tag+"unpack writes owned entries only (interleaved local order)",unpack_bad==0,"mismatches="+std::to_string(unpack_bad));
        GhostExchange exchange(p,l,rank);
        const long long ghosts=sum_all((long long)exchange.recv_nodes.size(),comm);
        const auto expected=owned_product(p,rank,1,5); long long incomplete=0;
        if (ghosts>0) {
            const auto stale=s.residual(halo_complete(raw(x_local),x_local.size()),raw(b));
            const long double oracle=sum_all(oracle_residual2(p,l,rank,download(raw(x_local),x_local.size()),1,expected,incomplete),comm);
            report.result(tag+"stale ghosts (exchange skipped) fail the true residual",!stale.passed
                && std::abs(stale.residual2-double(oracle))<=1e-9*double(oracle),
                "relative="+sci(stale.relative())+" adapter_r2="+sci(stale.residual2)+" oracle_r2="+sci(double(oracle)));
        } else report.skip(tag+"stale ghosts (exchange skipped) fail the true residual","no ghosts on one rank");
        exchange.run(raw(x_local),comm);
        const auto fresh=s.residual(halo_complete(raw(x_local),x_local.size()),raw(b));
        const long double oracle=sum_all(oracle_residual2(p,l,rank,download(raw(x_local),x_local.size()),1,expected,incomplete),comm);
        std::ostringstream detail; detail<<"relative="<<fresh.relative()<<" oracle_r2="<<double(oracle)<<" adapter_r2="<<fresh.residual2;
        report.result(tag+"exchanged known solution passes the true residual",fresh.passed && sum_all(incomplete,comm)==0
            && fresh.relative()<1e-14,detail.str());
        if (ghosts>0) {   // One referenced ghost left stale on the last rank.
            Buffer<double> again=upload(local_solution(p,l,5,1.0));
            s.unpack(raw(x_owned),x_owned.size(),raw(again),again.size());
            const int slot=injecting?exchange.first_referenced_recv():-1;
            exchange.run(raw(again),comm,slot);
            const auto one=s.residual(halo_complete(raw(again),again.size()),raw(b));
            if (all_true(!injecting || slot>=0,comm))
                report.result(tag+"a single stale referenced ghost fails the true residual",!one.passed,"relative="+sci(one.relative()));
            else report.skip(tag+"a single stale referenced ghost fails the true residual","last rank has no referenced ghost");
        }
        // Zero RHS: the absolute tolerance decides.
        Buffer<double> zero_b(b.size(),0.0), zero_x(x_local.size(),0.0);
        const auto zero=s.residual(halo_complete(raw(zero_x),zero_x.size()),raw(zero_b));
        const auto nonzero=s.residual(halo_complete(raw(x_local),x_local.size()),raw(zero_b));
        report.result(tag+"zero RHS: absolute tolerance accepts x=0 and rejects x!=0",zero.passed && zero.residual2==0 && !nonzero.passed);
    }
    {   // An omitted ghost-column coupling must fail both the oracle and the true residual.
        Options o=base; o.omit_coupling=true; Problem p(C,ranks,o); Local l=extract(p,rank);
        fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,5); });
        Device<C,GlobalId> d(l); Buffer<double> b(C*l.owned.size());
        System<C,Matrix,GlobalId> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
        s.update(d.view(),raw(b),b.size());
        Downloaded<C,Matrix,GlobalId> got(s); std::string note;
        const long long bad=sum_all(oracle_mismatches(p,l,rank,s.first_solver_node(),got.offsets,got.columns,got.values,got.map,
            download(raw(b),b.size()),1,owned_product(p,rank,1,5),note),comm);
        Buffer<double> x=upload(local_solution(p,l,5,0.0)); GhostExchange(p,l,rank).run(raw(x),comm);
        const auto norms=s.residual(halo_complete(raw(x),x.size()),raw(b));
        const bool ghost=sum_all((long long)(l.omitted_column>=0 && p.owner[l.omitted_column]!=rank),comm)>0;
        report.result(tag+"omitted "+(ghost?"ghost-column":"off-diagonal")+" coupling fails oracle and residual",bad>0 && !norms.passed,
            "oracle mismatches="+std::to_string(bad)+" relative="+sci(norms.relative()));
    }
    for (bool nothing:{false,true}) {   // Empty work: a rank without owned rows (and optionally without nodes).
        const std::string name=tag+(nothing?"rank without any local nodes":"rank without owned rows");
        if (ranks<2) { report.skip(name,"needs >= 2 ranks"); continue; }
        Options o=base; o.empty_rank=true; o.empty_rank_holds_nothing=nothing; Problem p(C,ranks,o); Local l=extract(p,rank);
        fill(p,l,rank,1,[&](int g,int c) { return p.product(g,c,1,5); });
        Device<C,GlobalId> d(l); Buffer<double> b(C*l.owned.size());
        std::string message;
        const bool rejected=all_ranks_throw(comm,[&] { System<C,Matrix,GlobalId> s(comm,d.view(),raw(d.owned),int(l.owned.size()),
            raw(d.solver_node),l.nodes()); },"policy reject",rank==1,message);
        System<C,Matrix,GlobalId> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes(),EmptyRanks::allow);
        s.update(d.view(),raw(b),b.size());
        Downloaded<C,Matrix,GlobalId> got(s); std::string note;
        const long long bad=sum_all(oracle_mismatches(p,l,rank,s.first_solver_node(),got.offsets,got.columns,got.values,got.map,
            download(raw(b),b.size()),1,owned_product(p,rank,1,5),note),comm);
        std::vector<double> xs; for (int g:p.owned_order[rank]) for (int c=0;c<C;++c) xs.push_back(p.solution(g,c,5));
        Buffer<double> x_owned=upload(xs), x=upload(local_solution(p,l,5,1.0));
        s.unpack(raw(x_owned),x_owned.size(),raw(x),x.size()); GhostExchange(p,l,rank).run(raw(x),comm);
        const auto norms=s.residual(halo_complete(raw(x),x.size()),raw(b));
        int empty_rows=rank==1?s.rows():-1, reported=0; MPI_Allreduce(&empty_rows,&reported,1,MPI_INT,MPI_MAX,comm);
        report.result(name+": reject policy collective, allow policy exact",all_true(rejected && bad==0 && norms.passed && (rank!=1 || s.rows()==0),comm),
            "rank 1 rows="+std::to_string(reported)+" local nodes="+std::to_string(int(sum_all((long long)(rank==1?l.nodes():0),comm)))
            +", oracle mismatches="+std::to_string(bad)+", relative="+sci(norms.relative()));
    }
    {   // Invalid inputs injected on the last rank only: every rank must reject at the same collective.
        Problem p(C,ranks,base); Local l=extract(p,rank);
        fill(p,l,rank,1,[&](int g,int c) { return p.rhs(g,c,1); });
        const int ghost=[&] { for (int v=0;v<l.nodes();++v) if (p.owner[l.global[v]]!=rank && l.referenced[v]) return v; return -1; }();
        {
            Device<C,GlobalId> fresh(l); Buffer<double> x(C*std::size_t(l.nodes()),0.0), b(C*l.owned.size(),0.0);
            System<C,Matrix,GlobalId> s(comm,fresh.view(),raw(fresh.owned),int(l.owned.size()),raw(fresh.solver_node),l.nodes());
            expect_rejection(comm,report,tag+"collective rejection: residual before the first update","before any update",0,[&] {
                s.residual(halo_complete(raw(x),x.size()),raw(b)); });
        }
        auto build=[&](Local bad_local,std::size_t map_size,long long limit) {
            Device<C,GlobalId> d(bad_local);
            System<C,Matrix,GlobalId> s(comm,d.view(),raw(d.owned),int(bad_local.owned.size()),raw(d.solver_node),map_size,EmptyRanks::reject,limit);
        };
        struct Case { const char* name; const char* fault; bool applicable; std::function<void()> run; };
        auto with=[&](auto edit) { Local m=l; if (injecting) edit(m); return m; };
        const bool two=l.owned.size()>=2;
        std::vector<Case> cases={
            {"owned local node out of range","owned node out of range",true,[&] { build(with([](Local& m) { m.owned[0]=m.nodes(); }),l.nodes(),1LL<<40); }},
            {"duplicate owned node","duplicate owned node",two,[&] { build(with([](Local& m) { m.owned[1]=m.owned[0]; }),l.nodes(),1LL<<40); }},
            {"owned rows out of solver order","not contiguous",two,[&] { build(with([](Local& m) { std::swap(m.solver_node[m.owned[0]],m.solver_node[m.owned[1]]); }),l.nodes(),1LL<<40); }},
            {"referenced ghost without solver id","without solver id",true,[&] { build(with([&](Local& m) { m.solver_node[ghost>=0?ghost:m.owned[0]]=-1; }),l.nodes(),1LL<<40); }},
            {"ghost id inside own owned range","inside own owned range",ghost>=0 || !injecting,[&] { build(with([&](Local& m) { m.solver_node[ghost]=p.first[rank]; }),l.nodes(),1LL<<40); }},
            {"solver id beyond global count","solver id out of range",true,[&] { build(with([&](Local& m) { m.solver_node[m.owned[0]]=p.nodes; }),l.nodes(),1LL<<40); }},
            {"map capacity below local nodes","capacity",true,[&] { build(l,injecting?std::size_t(l.nodes()-1):std::size_t(l.nodes()),1LL<<40); }},
            {"global rows beyond the wrapper index range","overflow",true,[&] { build(l,l.nodes(),(long long)C*p.nodes-1); }},
        };
        for (auto& c:cases) {
            const std::string name=tag+"collective rejection: "+c.name;
            if (!all_true(c.applicable,comm)) { report.skip(name,"not applicable to this partition"); continue; }
            expect_rejection(comm,report,name,c.fault,last,c.run);
        }
        Device<C,GlobalId> d(l); Buffer<double> b(C*l.owned.size());
        System<C,Matrix,GlobalId> s(comm,d.view(),raw(d.owned),int(l.owned.size()),raw(d.solver_node),l.nodes());
        expect_rejection(comm,report,tag+"collective rejection: RHS capacity","capacity",last,[&] {
            s.update(d.view(),raw(b),injecting?b.size()-1:b.size()); });
        expect_rejection(comm,report,tag+"collective rejection: caller-reported assembly error","caller-reported",last,[&] {
            s.update(d.view(),raw(b),b.size(),injecting); });
        Buffer<int> other=d.offsets;
        expect_rejection(comm,report,tag+"collective rejection: block graph replaced after build","graph changed",last,[&] {
            auto v=d.view(); if (injecting) v.offsets=raw(other); s.update(v,raw(b),b.size()); });
        auto poisoned=l.blocks; poisoned[std::size_t(C*C)*l.offsets[l.owned[0]]]=std::numeric_limits<double>::infinity();
        if (injecting) overwrite(d.blocks,poisoned);
        expect_rejection(comm,report,tag+"collective rejection: nonfinite owned value (ghost rows are NaN and ignored)","nonfinite",last,[&] {
            s.update(d.view(),raw(b),b.size()); });
        overwrite(d.blocks,l.blocks); s.update(d.view(),raw(b),b.size());
        Buffer<double> x(C*std::size_t(l.nodes()),0.0);
        expect_rejection(comm,report,tag+"collective rejection: undersized solution at unpack, reported by residual","capacity",last,[&] {
            s.unpack(raw(x),injecting?std::size_t(s.rows())-1:std::size_t(s.rows()),raw(x),x.size());
            s.residual(halo_complete(raw(x),x.size()),raw(b)); });
        expect_rejection(comm,report,tag+"collective rejection: residual vector smaller than local columns","capacity",last,[&] {
            s.residual(halo_complete(raw(x),injecting?x.size()-1:x.size()),raw(b)); });
    }
}
} // namespace dmatrix_gate
#undef GATE_HD
