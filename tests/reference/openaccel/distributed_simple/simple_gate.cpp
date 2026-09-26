// Distributed SIMPLE gate. Phase 1 (--write-reference FILE, one rank): the unchanged one-rank
// SimpleRunner on a generated channel writes snapshots keyed by global node/element/face and
// component. Phase 2 (--reference FILE, any rank count): DistributedSimpleRunner on an irregular
// partition writes the same snapshots; rank 0 compares every entry. Host builds solve with a
// test-only gathered dense LU; CUDA builds (MARS_REPLAY_CUDA) use Hypre on both sides.
#include "channel.hpp"
#include "mars_segregated_simple_distributed.hpp"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace mars::segregated;
using namespace mars::segregated::runtime;
using namespace dsimple_gate;

namespace {
using Entries=std::vector<std::pair<long long,double>>;
using Snapshots=std::map<std::string,Entries>;
long long node_key(int g,int c) { return 4LL*g+c; }
long long entry_key(int g,int c,int u,int j) { return node_key(g,c)*(1LL<<24)+node_key(u,j); }
std::string tag(const char* name,int k) { return std::string(name)+"@"+std::to_string(k); }

#ifdef MARS_REPLAY_CUDA
using Matrix=HypreSimpleSolve<1>::Solver::Matrix;
using GlobalId=HYPRE_BigInt;
template<int C> using Solve=HypreSimpleSolve<C>;
void graph_host(Graph& g,std::vector<int>& offsets,std::vector<int>& columns) {
    thrust::host_vector<int> o=g.graph.offsets, c=g.graph.columns; offsets.assign(o.begin(),o.end()); columns.assign(c.begin(),c.end());
}
template<class T> std::vector<T> download(const T* p,std::size_t n) {
    std::vector<T> h(n); if (n) assembly_cuda_check(cudaMemcpy(h.data(),p,n*sizeof(T),cudaMemcpyDeviceToHost)); return h;
}
#else
struct Matrix {
    std::vector<int> offsets, columns; std::vector<double> values;
    void allocate(int rows,int,int nnz) { offsets.assign(std::size_t(rows)+1,0); columns.assign(std::size_t(nnz),0); values.assign(std::size_t(nnz),0); }
    int* rowOffsetsPtr() { return offsets.data(); } const int* rowOffsetsPtr() const { return offsets.data(); }
    int* colIndicesPtr() { return columns.data(); } const int* colIndicesPtr() const { return columns.data(); }
    double* valuesPtr() { return values.data(); } const double* valuesPtr() const { return values.data(); }
};
using GlobalId=long long;
void graph_host(Graph& g,std::vector<int>& offsets,std::vector<int>& columns) { offsets=g.offsets; columns=g.columns; }
template<class T> std::vector<T> download(const T* p,std::size_t n) { return std::vector<T>(p,p+n); }
// Test-only oracle: gather owned rows by solver DOF, dense partial-pivot LU on rank 0, scatter.
template<int C> struct Solve {
    std::vector<double> b, x;
    double* rhs(std::size_t rows) { b.assign(rows,0); return b.data(); }
    const double* rhs() const { return b.data(); }
    const double* solution() const { return x.data(); }
    std::size_t size() const { return x.size(); }
    template<class System> bool operator()(const System& s) {
        int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        const int rows=s.rows(); const long long first=C*s.first_solver_node(); const int total=int(C*s.solver_nodes());
        const int* o=s.matrix().rowOffsetsPtr(); const int* c=s.matrix().colIndicesPtr(); const double* v=s.matrix().valuesPtr();
        const auto& map=s.solver_dof_map();
        std::vector<long long> ij; std::vector<double> a;
        for (int i=0;i<rows;++i) for (int k=o[i];k<o[i+1];++k) { ij.push_back((first+i)*total+map[c[k]]); a.push_back(v[k]); }
        int count=int(a.size()); std::vector<int> counts(ranks), displs(ranks), row_counts(ranks), row_displs(ranks);
        MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gather(&rows,1,MPI_INT,row_counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
        int all=0, all_rows=0;
        for (int q=0;q<ranks;++q) { displs[q]=all; all+=counts[q]; row_displs[q]=all_rows; all_rows+=row_counts[q]; }
        std::vector<long long> gij(std::size_t(rank==0?all:0)); std::vector<double> ga(gij.size()), gb(std::size_t(rank==0?total:0));
        MPI_Gatherv(ij.data(),count,MPI_LONG_LONG,gij.data(),counts.data(),displs.data(),MPI_LONG_LONG,0,MPI_COMM_WORLD);
        MPI_Gatherv(a.data(),count,MPI_DOUBLE,ga.data(),counts.data(),displs.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(b.data(),rows,MPI_DOUBLE,gb.data(),row_counts.data(),row_displs.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        int ok=1; std::vector<double> gx(std::size_t(rank==0?total:0));
        if (rank==0) {
            std::vector<double> m(std::size_t(total)*total,0.0);
            for (std::size_t k=0;k<gij.size();++k) m[std::size_t(gij[k])]+=ga[k];
            for (int k=0;k<total && ok;++k) {
                int pivot=k;
                for (int i=k+1;i<total;++i) if (std::abs(m[std::size_t(i)*total+k])>std::abs(m[std::size_t(pivot)*total+k])) pivot=i;
                if (!(std::abs(m[std::size_t(pivot)*total+k])>1e-20)) { ok=0; break; }
                if (pivot!=k) { for (int j=k;j<total;++j) std::swap(m[std::size_t(k)*total+j],m[std::size_t(pivot)*total+j]); std::swap(gb[k],gb[pivot]); }
                for (int i=k+1;i<total;++i) {
                    const double r=m[std::size_t(i)*total+k]/m[std::size_t(k)*total+k];
                    for (int j=k+1;j<total;++j) m[std::size_t(i)*total+j]-=r*m[std::size_t(k)*total+j];
                    gb[i]-=r*gb[k];
                }
            }
            for (int i=total-1;i>=0 && ok;--i) { double t=gb[i]; for (int j=i+1;j<total;++j) t-=m[std::size_t(i)*total+j]*gx[j]; gx[i]=t/m[std::size_t(i)*total+i]; }
        }
        MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);
        x.assign(std::size_t(rows),0);
        MPI_Scatterv(gx.data(),row_counts.data(),row_displs.data(),MPI_DOUBLE,x.data(),rows,MPI_DOUBLE,0,MPI_COMM_WORLD);
        return ok!=0;
    }
};
#endif

struct Options {
    int nx=16, ny=4, nz=4, iterations=2, converge=0;
    std::string reference, write, fault;
    double backflow=0;   // initial outlet-region velocity, see initial()
    double tolerance=1e-10;
};

// Global-keyed snapshots of one rank's copy (owned rows/nodes; every held element/face).
struct Collector {
    Snapshots s;
    void node_field(const std::string& name,const std::vector<double>& a,int C,const std::vector<int>& global,const std::vector<char>& owned) {
        auto& out=s[name];
        for (std::size_t v=0;v<global.size();++v) if (owned[v]) for (int c=0;c<C;++c) out.push_back({node_key(global[v],c),a[C*v+c]});
    }
    void samples(const std::string& name,const std::vector<double>& a,int per,const std::vector<int>& global) {
        auto& out=s[name];
        for (std::size_t i=0;i<global.size();++i) for (int j=0;j<per;++j) out.push_back({8LL*global[i]+j,a[per*i+j]});
    }
    template<int C> void blocks(const std::string& name,Graph& g,const std::vector<double>& values,const std::vector<double>& rhs,
                                const std::vector<int>& global,const std::vector<char>& owned) {
        std::vector<int> o, c; graph_host(g,o,c); auto& m=s[name+"_matrix"]; auto& r=s[name+"_rhs"];
        for (std::size_t v=0;v<global.size();++v) {
            if (!owned[v]) continue;
            for (int b=o[v];b<o[v+1];++b) for (int i=0;i<C;++i) for (int j=0;j<C;++j)
                m.push_back({entry_key(global[v],i,global[c[b]],j),values[std::size_t(C*C)*b+C*i+j]});
            for (int i=0;i<C;++i) r.push_back({node_key(global[v],i),rhs[C*v+i]});
        }
    }
    // What the adapter hands to Hypre: owned scalar rows, local columns through the solver map.
    template<class System> void adapter(const std::string& name,const System& sys,const Partition& p,int C) {
        const auto o=download(sys.matrix().rowOffsetsPtr(),std::size_t(sys.rows())+1);
        const auto c=download(sys.matrix().colIndicesPtr(),std::size_t(sys.nnz()));
        const auto v=download(sys.matrix().valuesPtr(),std::size_t(sys.nnz()));
        std::vector<GlobalId> map(sys.solver_dof_map().size());
#ifdef MARS_REPLAY_CUDA
        thrust::copy(sys.solver_dof_map().begin(),sys.solver_dof_map().end(),map.begin());
#else
        map=sys.solver_dof_map();
#endif
        auto& m=s[name+"_adapter"];
        for (int i=0;i<sys.rows();++i) {
            const long long R=C*sys.first_solver_node()+i; const int g=p.solver_to_global[R/C];
            for (int k=o[i];k<o[i+1];++k) {
                const long long Q=(long long)map[c[k]];
                m.push_back({entry_key(g,int(R%C),p.solver_to_global[Q/C],int(Q%C)),v[k]});
            }
        }
    }
    void sums(const std::string& name,const SimpleSums& a) {
        auto& out=s[name];
        const double v[]={a.volume,a.momentum2,a.continuity2,a.continuity,a.velocity_change2,a.pressure_change2,a.speed2,
                          a.inlet,a.outlet,a.inlet_area,a.flux_change,double(a.closed),double(a.changed),double(a.invalid)};
        for (int i=0;i<14;++i) out.push_back({i,v[i]});
    }
};

void write_snapshots(const std::string& path,const Snapshots& s,const std::string& header) {
    std::ofstream out(path,std::ios::binary); if (!out) throw std::runtime_error("cannot write "+path);
    auto put=[&](const void* p,std::size_t n) { out.write(static_cast<const char*>(p),std::streamsize(n)); };
    std::size_t n=header.size(); put(&n,sizeof n); put(header.data(),n);
    n=s.size(); put(&n,sizeof n);
    for (const auto& [name,entries]:s) {
        std::size_t k=name.size(); put(&k,sizeof k); put(name.data(),k);
        k=entries.size(); put(&k,sizeof k); put(entries.data(),k*sizeof(entries[0]));
    }
    if (!out) throw std::runtime_error("write failed: "+path);
}
Snapshots read_snapshots(const std::string& path,std::string& header) {
    std::ifstream in(path,std::ios::binary); if (!in) throw std::runtime_error("cannot read "+path);
    auto get=[&](void* p,std::size_t n) { in.read(static_cast<char*>(p),std::streamsize(n)); if (!in) throw std::runtime_error("truncated "+path); };
    std::size_t n=0; get(&n,sizeof n); header.resize(n); get(header.data(),n);
    Snapshots s; get(&n,sizeof n);
    for (std::size_t i=0;i<n;++i) {
        std::size_t k=0; get(&k,sizeof k); std::string name(k,' '); get(name.data(),k);
        get(&k,sizeof k); Entries e(k); get(e.data(),k*sizeof(e[0])); s[name]=std::move(e);
    }
    return s;
}
// Gather every rank's entries per name on rank 0 (names are identical on all ranks).
Snapshots gather(const Snapshots& mine) {
    int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    Snapshots all;
    for (const auto& [name,e]:mine) {
        int count=int(e.size()); std::vector<int> counts(ranks), displs(ranks);
        MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
        int total=0; for (int q=0;q<ranks;++q) { displs[q]=total*int(sizeof(e[0])); total+=counts[q]; counts[q]*=int(sizeof(e[0])); }
        Entries out(std::size_t(rank==0?total:0));
        MPI_Gatherv(e.data(),count*int(sizeof(e[0])),MPI_BYTE,out.data(),counts.data(),displs.data(),MPI_BYTE,0,MPI_COMM_WORLD);
        if (rank==0) all[name]=std::move(out);
    }
    return all;
}

// Every reference key must be present; unique-ownership names exactly once, copies all equal.
struct Verdict { int failures=0; double worst=0; std::string first; };
Verdict compare(const Snapshots& ref,const Snapshots& got,double tolerance) {
    Verdict v;
    for (const auto& [name,r]:ref) {
        // Element/face histories and the allreduced diagnostics exist once per holding rank.
        const bool copies=name.rfind("eflux",0)==0 || name.rfind("bflux",0)==0 || name.rfind("trace",0)==0
            || name.rfind("flags",0)==0 || name.rfind("diagnostics",0)==0;
        auto it=got.find(name);
        if (it==got.end()) { ++v.failures; if (v.first.empty()) v.first=name+": missing"; continue; }
        std::map<long long,double> expect; double scale=0;
        for (const auto& [k,x]:r) { expect[k]+=x; scale=std::max(scale,std::abs(x)); }
        std::map<long long,int> seen; double worst=0; long long bad=0;
        for (const auto& [k,x]:it->second) {
            auto e=expect.find(k);
            if (e==expect.end()) { ++bad; continue; }
            ++seen[k]; const double d=std::abs(x-e->second)/std::max(scale,1e-300);
            if (!(d<=tolerance)) ++bad;
            worst=std::max(worst,std::isfinite(d)?d:1e300);
        }
        for (const auto& [k,x]:expect) { (void)x; const int c=seen[k]; if (c==0 || (!copies && c!=1)) ++bad; }
        v.worst=std::max(v.worst,worst);
        std::cout<<"  "<<std::left<<std::setw(26)<<name<<" entries="<<std::setw(7)<<r.size()<<" received="<<std::setw(7)<<it->second.size()
                 <<" worst_scaled="<<std::scientific<<std::setprecision(2)<<worst<<std::defaultfloat<<(bad?"  MISMATCH "+std::to_string(bad):"")<<'\n';
        if (bad) { ++v.failures; if (v.first.empty()) v.first=name; }
    }
    return v;
}

std::string header_of(const Options& o) {
    std::ostringstream h; h<<"MARS_DSIMPLE_V1 "<<o.nx<<'x'<<o.ny<<'x'<<o.nz<<" iterations="<<o.iterations<<" converge="<<o.converge<<" backflow="<<o.backflow;
    return h.str();
}
// lx is the global channel length: a rank's local extent must not change the initial field.
template<class Runner> void initial(Runner& run,const std::vector<double>& x,const std::vector<double>& y,double lx,double backflow,const SimpleControls& c) {
    if (backflow==0) return;
    // backflow>0: uniform reverse flow near the outlet (closes every outlet face);
    // backflow<0: sheared, forward below y=0.5 and reverse above (closes part of the outlet).
    std::vector<double> u(3*x.size(),0.0);
    for (std::size_t v=0;v<x.size();++v) if (x[v]>0.75*lx)
        u[3*v]=backflow>0?-backflow*c.inlet_speed:(y[v]>0.5?backflow:-backflow)*c.inlet_speed;
    run.velocity.values=u;
}

// Runs a runner (reference or distributed) through the schedule and records snapshots.
template<class Runner,class Record> std::string drive(Runner& run,const Options& o,Record record) {
    try {
        for (int k=0;;++k) {
            run.assemble_momentum();
            const auto sums=run.diagnostics(); const auto m=simple_metrics(sums,run.controls);
            record.assembled(run,k,sums);
            const bool done=o.converge>0 && simple_converged(m,run.completed,sums.changed,1e-6,1e-6,1e-6);
            if (done) return "converged@"+std::to_string(run.completed);
            if (k==(o.converge>0?o.converge:o.iterations)) return o.converge>0?"not converged":"done";
            run.advance(); record.advanced(run,k);
        }
    } catch (const std::exception& e) { return std::string("threw: ")+e.what(); }
}

int execute(const Options& o) {
    int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    const auto mesh=channel(o.nx,o.ny,o.nz); const SimpleControls controls;
    if (!o.write.empty()) {
        if (ranks!=1) throw std::runtime_error("--write-reference runs the one-rank SimpleRunner: use one rank");
        SimpleRunner run(mesh); run.momentum.verbose=run.poisson.verbose=false;
        initial(run,mesh.x,mesh.y,double(o.nx)/o.ny,o.backflow,controls);
        std::vector<int> global(mesh.x.size()); for (std::size_t g=0;g<global.size();++g) global[g]=int(g);
        std::vector<int> el(mesh.nodes[0].size()), fa(mesh.faces.size());
        for (std::size_t i=0;i<el.size();++i) el[i]=int(i);
        for (std::size_t i=0;i<fa.size();++i) fa[i]=int(i);
        const std::vector<char> all(global.size(),1);
        Collector c;
        struct { Collector& c; const std::vector<int> &global,&el,&fa; const std::vector<char>& all; bool full;
            void assembled(SimpleRunner& r,int k,const SimpleSums& s) {
                c.sums(tag("diagnostics",k),s);
                if (full) c.blocks<3>(tag("momentum",k),r.graph,r.momentum.blocks.host(),r.momentum.rhs.host(),global,all);
            }
            void advanced(SimpleRunner& r,int k) {
                if (!full) return;
                c.blocks<1>(tag("pressure",k),r.graph,r.poisson.blocks.host(),r.poisson.rhs.host(),global,all);
                c.node_field(tag("velocity",k),r.velocity.host(),3,global,all); c.node_field(tag("pressure",k),r.pressure.host(),1,global,all);
                c.node_field(tag("influence",k),r.d.host(),3,global,all); c.node_field(tag("divergence",k),r.div.host(),1,global,all);
                c.samples(tag("eflux",k),r.eflux.host(),6,el); c.samples(tag("bflux",k),r.bflux.host(),3,fa); c.samples(tag("trace",k),r.trace.host(),3,fa);
                const auto f=r.flags.host(); c.samples(tag("flags",k),std::vector<double>(f.begin(),f.end()),3,fa);
            }
        } rec{c,global,el,fa,all,o.converge==0};
        const std::string outcome=drive(run,o,rec);
        if (o.converge>0) c.node_field("velocity@final",run.velocity.host(),3,global,all), c.node_field("pressure@final",run.pressure.host(),1,global,all);
        c.s["outcome:"+outcome];
        write_snapshots(o.write,c.s,header_of(o));
        std::cout<<"reference "<<header_of(o)<<" outcome="<<outcome<<" snapshots="<<c.s.size()<<" -> "<<o.write<<'\n';
        return 0;
    }
    std::string header; const auto ref=read_snapshots(o.reference,header);
    if (header!=header_of(o)) throw std::runtime_error("reference header '"+header+"' does not match options '"+header_of(o)+"'");
    const auto p=partition(mesh,ranks);
    int drop=-1, drop_rank=-1;
    if (o.fault=="missing-element") {
        // Lowest-rank ownership leaves the highest rank without halo elements; rank 0 has the most.
        drop_rank=0;
        for (int el:present_elements(mesh,p,drop_rank)) if (p.element_owner[el]!=drop_rank) { drop=el; break; }
        if (drop<0) throw std::runtime_error("missing-element fault: rank 0 holds no halo element (use >= 2 ranks)");
    }
    auto part=extract(mesh,p,rank,drop,drop_rank);
    if (o.fault=="duplicate-face") {
        int injected=0;
        if (rank==0) for (int i:part.ownership.owned_faces) if (part.input.faces[i].kind==0) { part.ownership.owned_faces.push_back(i); injected=1; break; }
        MPI_Bcast(&injected,1,MPI_INT,0,MPI_COMM_WORLD);
        if (!injected) throw std::runtime_error("duplicate-face fault: rank 0 owns no inlet face");
    }
    if (o.fault=="stale-ghost") {
        int ghosts=int(std::count(part.node_owned.begin(),part.node_owned.end(),0)), last=0;
        if (rank==ranks-1) last=ghosts;
        MPI_Bcast(&last,1,MPI_INT,ranks-1,MPI_COMM_WORLD);
        if (!last) throw std::runtime_error("stale-ghost fault: the last rank holds no ghost node (use >= 2 ranks)");
    }
    if (!o.fault.empty() && o.fault!="missing-element" && o.fault!="duplicate-face" && o.fault!="stale-ghost")
        throw std::runtime_error("unknown fault "+o.fault);
    DistributedSimpleRunner<Matrix,GlobalId,Solve> run(MPI_COMM_WORLD,part.input,part.ownership,controls);
    run.poison_unexchanged=true;
    initial(run,part.input.x,part.input.y,double(o.nx)/o.ny,o.backflow,controls);
    Collector c;
    int closed_owned=0;   // outlet faces this rank owns that were ever closed
    struct { Collector& c; Part& part; const Partition& p; const std::string& fault; int rank, ranks; bool full; int& closed;
        void assembled(DistributedSimpleRunner<Matrix,GlobalId,Solve>& r,int k,const SimpleSums& s) {
            c.sums(tag("diagnostics",k),s);
            if (full) c.blocks<3>(tag("momentum",k),r.graph,r.momentum_blocks.host(),r.momentum_rhs.host(),part.node_global,part.node_owned);
        }
        void advanced(DistributedSimpleRunner<Matrix,GlobalId,Solve>& r,int k) {
            if (fault=="stale-ghost" && k==0 && rank==ranks-1) {   // one ghost misses the final publish
                auto u=r.velocity.host(); const auto old=r.old_velocity.host();
                for (std::size_t v=0;v<part.node_owned.size();++v) if (!part.node_owned[v]) { for (int j=0;j<3;++j) u[3*v+j]=old[3*v+j]; break; }
                r.velocity.values=u;
            }
            const auto flags=r.flags.host(); int now=0;
            for (int i:part.ownership.owned_faces) now+=flags[3*i]!=0;
            closed=std::max(closed,now);
            if (!full) return;
            c.blocks<1>(tag("pressure",k),r.graph,r.poisson_blocks.host(),r.poisson_rhs.host(),part.node_global,part.node_owned);
            c.adapter(tag("pressure",k),r.poisson,p,1);
            c.node_field(tag("velocity",k),r.velocity.host(),3,part.node_global,part.node_owned);
            c.node_field(tag("pressure",k),r.pressure.host(),1,part.node_global,part.node_owned);
            c.node_field(tag("influence",k),r.d.host(),3,part.node_global,part.node_owned);
            c.node_field(tag("divergence",k),r.div.host(),1,part.node_global,part.node_owned);
            c.samples(tag("eflux",k),r.eflux.host(),6,part.element_global); c.samples(tag("bflux",k),r.bflux.host(),3,part.face_global);
            c.samples(tag("trace",k),r.trace.host(),3,part.face_global);
            const auto f=r.flags.host(); c.samples(tag("flags",k),std::vector<double>(f.begin(),f.end()),3,part.face_global);
        }
    } rec{c,part,p,o.fault,rank,ranks,o.converge==0,closed_owned};
    // The momentum adapter matrix is what Hypre receives; record it right after assembly.
    struct Wrapped { decltype(rec)& inner; Collector& c; const Partition& p;
        void assembled(DistributedSimpleRunner<Matrix,GlobalId,Solve>& r,int k,const SimpleSums& s) { inner.assembled(r,k,s); }
        void advanced(DistributedSimpleRunner<Matrix,GlobalId,Solve>& r,int k) {
            if (inner.full) c.adapter(tag("momentum",k),r.momentum,p,3);
            inner.advanced(r,k);
        }
    } wrapped{rec,c,p};
    const std::string outcome=drive(run,o,wrapped);
    if (o.converge>0) { c.node_field("velocity@final",run.velocity.host(),3,part.node_global,part.node_owned);
                        c.node_field("pressure@final",run.pressure.host(),1,part.node_global,part.node_owned); }
    // Rename the adapter snapshots onto the reference's block-matrix names for comparison.
    Snapshots mine;
    for (auto& [name,e]:c.s) {
        const auto at=name.find("_adapter");
        if (at!=std::string::npos) mine[name.substr(0,at)+"_matrix"+name.substr(at+8)+"#adapter"]=e; else mine[name]=e;
    }
    const auto all=gather(mine);
    long long outlet_ranks=0, mine_outlet=0;
    for (int i:part.ownership.owned_faces) mine_outlet+=part.input.faces[i].kind==1;
    const long long has=mine_outlet>0; MPI_Allreduce(&has,&outlet_ranks,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    long long closing=closed_owned>0, closing_ranks=0; MPI_Allreduce(&closing,&closing_ranks,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    int failures=0;
    if (rank==0) {
        std::cout<<"distributed "<<header_of(o)<<" ranks="<<ranks<<" fault="<<(o.fault.empty()?"none":o.fault)
                 <<" ranks_with_outlet_faces="<<outlet_ranks<<" ranks_owning_closed_faces="<<closing_ranks
                 <<" exchange_rounds="<<run.exchange.rounds()<<" ghost_values_received="<<run.exchange.received_values()<<'\n';
        Snapshots reference=ref, adapter;
        for (auto& [name,e]:all) if (name.size()>8 && name.compare(name.size()-8,8,"#adapter")==0) adapter[name.substr(0,name.size()-8)]=e;
        Snapshots blocks; for (auto& [name,e]:all) if (!(name.size()>8 && name.compare(name.size()-8,8,"#adapter")==0)) blocks[name]=e;
        std::string ref_outcome, got_outcome=outcome;   // rank 0's outcome; failures are collective
        for (auto& [name,e]:reference) if (name.rfind("outcome:",0)==0) ref_outcome=name.substr(8);
        std::cout<<" outcome reference='"<<ref_outcome<<"' distributed='"<<got_outcome<<"'\n";
        Snapshots ref_values, ref_adapter;
        for (auto& [name,e]:reference) if (name.rfind("outcome:",0)!=0) ref_values[name]=e;
        for (auto& [name,e]:adapter) ref_adapter[name]=reference.count(name)?reference.at(name):Entries{};
        std::cout<<" block rows and fields:\n"; auto v=compare(ref_values,blocks,o.tolerance);
        std::cout<<" adapter (owned scalar rows handed to Hypre):\n"; auto a=compare(ref_adapter,adapter,o.tolerance);
        const bool same_outcome=ref_outcome==got_outcome || (ref_outcome.rfind("threw",0)==0 && got_outcome.rfind("threw",0)==0);
        failures=v.failures+a.failures+(same_outcome?0:1);
        const bool expect_mismatch=!o.fault.empty() && o.fault!="none";
        std::cout<<(expect_mismatch?(failures?"PASS: injected fault detected":"FAIL: injected fault not detected")
                                   :(failures?"FAIL: "+(v.first.empty()?a.first:v.first):"PASS: distributed SIMPLE matches the one-rank reference entry by entry"))<<'\n';
        failures=expect_mismatch?(failures?0:1):failures;
    }
    MPI_Bcast(&failures,1,MPI_INT,0,MPI_COMM_WORLD);
    return failures?1:0;
}
} // namespace

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
#ifdef MARS_REPLAY_CUDA
    {   // one GPU per rank on the node, as the other multi-rank drivers select it
        int rank=0, local=0, devices=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm node; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        MPI_Comm_rank(node,&local); MPI_Comm_free(&node);
        if (cudaGetDeviceCount(&devices)==cudaSuccess && devices>0) cudaSetDevice(local%devices);
    }
#endif
    Options o; int result=1;
#ifdef MARS_REPLAY_CUDA
    o.tolerance=1e-8;   // both sides stop Hypre at rtol 1e-12 on different partitions
#endif
    try {
        for (int i=1;i+1<argc;i+=2) {
            const std::string k=argv[i], v=argv[i+1];
            if (k=="--mesh") { char x; std::istringstream s(v); s>>o.nx>>x>>o.ny>>x>>o.nz; }
            else if (k=="--iterations") o.iterations=std::stoi(v);
            else if (k=="--converge") o.converge=std::stoi(v);
            else if (k=="--reference") o.reference=v;
            else if (k=="--write-reference") o.write=v;
            else if (k=="--fault") o.fault=v;
            else if (k=="--backflow") o.backflow=std::stod(v);
            else if (k=="--tolerance") o.tolerance=std::stod(v);
            else throw std::runtime_error("unknown option "+k);
        }
        if (o.write.empty()==o.reference.empty()) throw std::runtime_error("give exactly one of --write-reference FILE / --reference FILE");
        result=execute(o);
    } catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; result=1; }
    MPI_Finalize();
    return result;
}
