#pragma once
// Distributed SIMPLE: SimpleRunner's iteration, unchanged kernels, on a rank's partition.
//
// Inputs (from native ingestion and element-halo completion; not computed here):
//   local mesh     every element touching an owned node (complete stars), every boundary
//                  face of those elements, and all their nodes (owned and ghost, any order);
//   ownership      owned nodes in solver order, solver node per local node, the elements and
//                  boundary faces this rank owns uniquely (a face belongs to its element's owner);
//   halo lists     NodeHaloTopology-style peers and send/recv node lists.
//
// Schedule per iteration (fields published owner -> ghost, 4 rounds, 14 values per ghost;
// nothing at setup):
//   gradient(p)            owned rows exact from complete stars          -> publish grad p (3)
//   momentum assembly      owned rows only enter the solve; d at owned rows
//   momentum solve         -> publish [du (3), d (3)]; true residual; u += du on every node
//   outlet trace           area moments over owned open outlet faces, one allreduce; every
//                          rank holding a face recomputes the same trace from the same inputs
//   pressure assembly+solve-> publish phi (1); true residual; p += alpha_p*phi on every node
//   flux update            every held element/face recomputes its stored flux and reversal
//                          flags from identical inputs, so all copies of a history agree
//   velocity correction    owned nodes -> publish [u (3), p (1)]
// No reverse additions: owned sums come from complete stars, ghost partial sums are never read
// (validation can poison them). Diagnostics and outlet moments sum uniquely owned nodes, faces
// and element samples, then allreduce. Every failure is collective.
#include "mars_segregated_simple_runtime.hpp"
#include "mars_segregated_distributed_matrix.hpp"
#include "mars_segregated_halo_exchange.hpp"
#include <limits>
#include <memory>
#if defined(__CUDACC__) && !defined(MARS_REPLAY_CUDA)
#error "Distributed SIMPLE: CUDA builds define MARS_REPLAY_CUDA, as the SIMPLE drivers do"
#endif
#if defined(__CUDACC__)
#define MARS_DSIMPLE_HD __host__ __device__
#else
#define MARS_DSIMPLE_HD
#endif

namespace mars::segregated::runtime {

template<class GlobalId> struct SimpleOwnership {
    std::vector<int> owned_nodes, owned_elements, owned_faces;
    std::vector<GlobalId> solver_node;
    std::vector<int> peers, send_offsets, send_nodes, recv_offsets, recv_nodes;
};

// Apply a per-index functor to a list of indices (owned nodes, elements or faces).
template<class F> struct OnList {
    F f; const int* list;
    MARS_DSIMPLE_HD auto operator()(int k) const { return f(list[k]); }
};
// Samples of listed entities: entity list[i/per], sample i%per of a per-entity array.
template<class F> struct OnSamples {
    F f; const int* list; int per;
    MARS_DSIMPLE_HD auto operator()(int i) const { return f(per*list[i/per]+i%per); }
};
struct PoisonGhosts {
    double* values; int components; const unsigned char* owned; double value;
    MARS_DSIMPLE_HD void operator()(int i) const { if (!owned[i/components]) values[i]=value; }
};
struct MarkOwnedNodes {
    const int* list; unsigned char* owned;
    MARS_DSIMPLE_HD void operator()(int k) const { owned[list[k]]=1; }
};

inline void simple_collective(MPI_Comm comm,bool local_ok,const char* message) {
    int bad=local_ok?0:1, any=0; MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);
    if (any) throw std::runtime_error(std::string(message)+" (on "+(local_ok?"another rank":"this rank")+")");
}
inline SimpleSums allreduce_sums(MPI_Comm comm,const SimpleSums& a) {
    double sum[12]={a.volume,a.momentum2,a.continuity2,a.continuity,a.velocity_change2,a.pressure_change2,
                    a.inlet,a.outlet,a.inlet_area,double(a.closed),double(a.changed),double(a.invalid)};
    double max[2]={a.speed2,a.flux_change}, gs[12], gm[2];
    MPI_Allreduce(sum,gs,12,MPI_DOUBLE,MPI_SUM,comm); MPI_Allreduce(max,gm,2,MPI_DOUBLE,MPI_MAX,comm);
    SimpleSums g; g.volume=gs[0]; g.momentum2=gs[1]; g.continuity2=gs[2]; g.continuity=gs[3];
    g.velocity_change2=gs[4]; g.pressure_change2=gs[5]; g.inlet=gs[6]; g.outlet=gs[7]; g.inlet_area=gs[8];
    g.closed=int(gs[9]); g.changed=int(gs[10]); g.invalid=int(gs[11]); g.speed2=gm[0]; g.flux_change=gm[1];
    return g;
}

#ifdef MARS_REPLAY_CUDA
// Production linear solve: the Hypre device-map overload through the owned-row adapter.
template<int C> struct HypreSimpleSolve {
    using Solver=mars::fem::HypreGMRESSolver<double,int,cstone::GpuTag>;
    Solver solver{MPI_COMM_WORLD,2000,1e-12,Solver::BOOMERAMG,100};
    typename Solver::Vector b,x;
    HypreSimpleSolve() { solver.setVerbose(false); solver.setPointBlock(C); }
    double* rhs(std::size_t rows) { if (b.size()!=rows) { b.resize(rows); x.resize(rows); } return b.data(); }
    const double* rhs() const { return b.data(); }
    const double* solution() const { return x.data(); }
    std::size_t size() const { return x.size(); }
    template<class System> bool operator()(const System& s) {
        if (s.rows()) assembly_cuda_check(cudaMemset(x.data(),0,std::size_t(s.rows())*sizeof(double)));
        return distributed::solve_owned(solver,s,b,x);
    }
};
#endif

template<class Matrix,class GlobalId,template<int> class Solve>
struct DistributedSimpleRunner {
    MPI_Comm comm;
    int n,e,b,completed=0,owned_nodes,owned_elements,owned_faces;
    bool poison_unexchanged=false; // validation: NaN in ghost entries that must never be read
    SimpleControls controls;
    Array<double> x,y,z,velocity,pressure,vg,pg,d,volume,div,eflux,bflux,trace,factor,sum,gp,moment,old_velocity,old_pressure,old_eflux,old_bflux;
    Array<int> n0,n1,n2,n3,error,flags,old_flags,owned,elements,boundary;
    Array<unsigned char> owned_mask;
    Array<GlobalId> solver_node;
    Array<SimpleFace> faces; Array<TetGeometry<double>> geometry;
    SimpleMesh mesh; SimpleState state;
    Graph graph;
    Array<double> momentum_blocks,momentum_rhs,poisson_blocks,poisson_rhs,du,phi;
    distributed::OwnedRowSystem<3,Matrix,GlobalId> momentum;
    distributed::OwnedRowSystem<1,Matrix,GlobalId> poisson;
    Solve<3> momentum_solve; Solve<1> poisson_solve;
    distributed::FieldExchange exchange;
    bool assembled=false;
    distributed::Tolerance tolerance{1e-13,1e-10};

    // Id may be wider than GlobalId (e.g. int64 ids for a 32-bit HYPRE_BigInt build); ids
    // that do not fit are rejected collectively instead of narrowed.
    template<class Input,class Id> DistributedSimpleRunner(MPI_Comm c,const Input& f,const SimpleOwnership<Id>& o,SimpleControls ctl={},
        distributed::EmptyRanks empty=distributed::EmptyRanks::reject):
        comm(c),n(int(f.x.size())),e(int(f.nodes[0].size())),b(int(f.faces.size())),
        owned_nodes(int(o.owned_nodes.size())),owned_elements(int(o.owned_elements.size())),owned_faces(int(o.owned_faces.size())),controls(ctl),
        x(f.x),y(f.y),z(f.z),velocity(3*n),pressure(n),vg(9*n),pg(3*n),d(3*n),volume(n),div(n),
        eflux(6*e),bflux(3*b),trace(3*b),factor(n),sum(9*n),gp(3*n),moment(2),old_velocity(3*n),old_pressure(n),old_eflux(6*e),old_bflux(3*b),
        n0(f.nodes[0]),n1(f.nodes[1]),n2(f.nodes[2]),n3(f.nodes[3]),error(1),flags(3*b),old_flags(3*b),
        owned(o.owned_nodes),elements(o.owned_elements),boundary(o.owned_faces),owned_mask(std::size_t(n)),solver_node(narrow(c,o.solver_node)),
        faces(f.faces),geometry(e),
        mesh{n,e,b,{n0.data(),n1.data(),n2.data(),n3.data()},x.data(),y.data(),z.data(),faces.data(),geometry.data()},
        state{velocity.data(),pressure.data(),vg.data(),pg.data(),d.data(),volume.data(),div.data(),
              eflux.data(),bflux.data(),trace.data(),factor.data(),error.data(),flags.data()},
        graph(mesh),
        momentum_blocks(std::size_t(graph.blocks())*9),momentum_rhs(3*n),poisson_blocks(std::size_t(graph.blocks())),poisson_rhs(n),du(3*n),phi(n),
        momentum(c,graph.template view<3>(momentum_blocks.data(),momentum_rhs.data()),owned.data(),owned_nodes,solver_node.data(),solver_node.values.size(),empty),
        poisson(c,graph.template view<1>(poisson_blocks.data(),poisson_rhs.data()),owned.data(),owned_nodes,solver_node.data(),solver_node.values.size(),empty),
        exchange(c,o.peers,o.send_offsets,o.send_nodes,o.recv_offsets,o.recv_nodes,n)
    {
        simple_collective(comm,ctl.density>0 && ctl.viscosity>0 && ctl.pseudo_dt>0 && ctl.inlet_speed>0,"invalid material or pseudo-time");
        for (double alpha:{ctl.alpha_u,ctl.alpha_p,ctl.alpha_mass,ctl.beta})
            simple_collective(comm,alpha>0 && alpha<=1,"relaxation and beta must be in (0,1]");
        bool lists_ok=true;
        for (int v:o.owned_elements) lists_ok=lists_ok && v>=0 && v<e;
        for (int v:o.owned_faces) lists_ok=lists_ok && v>=0 && v<b;
        simple_collective(comm,lists_ok,"owned element or face list outside the local mesh");
        launch(owned_nodes,MarkOwnedNodes{owned.data(),owned_mask.data()});
        launch(e,SimpleGeometry{mesh,state}); check("native geometry failed");
        launch(b,SimpleBoundaryFactor{mesh,factor.data()});
        // Owned volumes and boundary factors are complete. Ghost entries stay partial on purpose:
        // they only feed ghost gradients that are republished or multiplied by a zero blend.
    }
    template<class Id> static std::vector<GlobalId> narrow(MPI_Comm c,const std::vector<Id>& ids) {
        bool fits=true;
        for (Id v:ids) fits=fits && (long double)v>=(long double)std::numeric_limits<GlobalId>::min()
                                  && (long double)v<=(long double)std::numeric_limits<GlobalId>::max();
        simple_collective(c,fits,"solver node id does not fit the solver's global index type");
        return std::vector<GlobalId>(ids.begin(),ids.end());
    }
    DistributedSimpleRunner(const DistributedSimpleRunner&)=delete;
    DistributedSimpleRunner& operator=(const DistributedSimpleRunner&)=delete;
    void check(const char* message) { simple_collective(comm,error.host()[0]==0,message); }
    // NaN where ghost entries are never read; a huge finite value where they are read but
    // multiplied by a zero blend (velocity gradient), so 0*value must stay 0.
    void poison(Array<double>& a,int components,bool never_read=true) {
        if (poison_unexchanged) launch(components*n,PoisonGhosts{a.data(),components,owned_mask.data(),
            never_read?std::numeric_limits<double>::quiet_NaN():1e30});
    }
    void assemble_momentum() {
        gradient<3>(mesh,state,state.velocity,sum,state.velocity_gradient);
        poison(vg,9,false); // only multiplied by velocity_blend=0; never published
        gradient<1>(mesh,state,state.pressure,sum,state.pressure_gradient);
        exchange({{pg.data(),3}});
        momentum_blocks.zero(); momentum_rhs.zero(); auto am=graph.template view<3>(momentum_blocks.data(),momentum_rhs.data());
        launch(e,SimpleInterior<3>{mesh,state,controls,am});
        launch(b,SimpleBoundary<3>{mesh,state,controls,am,completed>0,false});
        launch(owned_nodes,OnList<SimpleMomentumNode>{{state,controls,am},owned.data()});
        check("momentum assembly failed"); assembled=true;
    }
    SimpleSums diagnostics() {
        ensure(assembled,"diagnostics require a fresh momentum assembly"); // same program order on every rank
        auto a=reduce_sums(owned_nodes,OnList<SimpleNodeSums>{{state,momentum_rhs.data(),old_velocity.data(),old_pressure.data()},owned.data()});
        a=SimpleSumCombine{}(a,reduce_sums(owned_faces,OnList<SimpleFaceSums>{{mesh,state,old_flags.data()},boundary.data()}));
        a=SimpleSumCombine{}(a,reduce_sums(6*owned_elements,OnSamples<SimpleFluxChange>{{eflux.data(),old_eflux.data()},elements.data(),6}));
        a=SimpleSumCombine{}(a,reduce_sums(3*owned_faces,OnSamples<SimpleFluxChange>{{bflux.data(),old_bflux.data()},boundary.data(),3}));
        return allreduce_sums(comm,a);
    }
    template<int C,class System,class Solver> void solve(System& system,Solver& solver,BlockCsrView<C> view,Array<double>& increment,
        bool assembly_failed,std::initializer_list<distributed::Field> with) {
        system.update(view,solver.rhs(std::size_t(system.rows())),std::size_t(system.rows()),assembly_failed);
        // Solve policies return a rank-consistent result (Hypre: global norms; the test oracle: a broadcast).
        ensure(solver(system),"linear solve failed");
        system.unpack(solver.solution(),solver.size(),increment.data(),increment.values.size());
        if (with.size()==0) exchange({{increment.data(),C}});
        else { auto it=with.begin(); exchange({{increment.data(),C},*it}); }
        const auto norms=system.residual(distributed::halo_complete(increment.data(),increment.values.size()),solver.rhs(),tolerance);
        ensure(norms.passed,"true linear residual failed"); // from allreduced sums: identical on every rank
    }
    template<class Observer=NoSimpleObserver> void advance(Observer observe={}) {
        ensure(assembled,"advance requires momentum assembly");
        old_velocity.copy_from(velocity); old_pressure.copy_from(pressure); old_flags.copy_from(flags);
        old_eflux.copy_from(eflux); old_bflux.copy_from(bflux);
        poison(momentum_rhs,3);
        solve<3>(momentum,momentum_solve,graph.template view<3>(momentum_blocks.data(),momentum_rhs.data()),du,false,{{d.data(),3}});
        launch(3*n,SimpleAddIncrement{state.velocity,du.data(),1});
        observe("momentum",velocity); observe("influence",d);
        moment.zero(); launch(owned_faces,OnList<SimpleTraceMoment>{{mesh,state,moment.data()},boundary.data()});
        const auto local=moment.host(); std::vector<double> global(2);
        MPI_Allreduce(local.data(),global.data(),2,MPI_DOUBLE,MPI_SUM,comm);
        if (!(std::isfinite(global[0]) && std::isfinite(global[1]) && global[1]>0))
            throw std::runtime_error("all outlet faces closed: no open pressure anchor; cannot solve this prescribed-inflow case");
        moment.values=global;
        launch(b,SimpleTrace{mesh,state,controls,moment.data()}); observe("trace",trace);
        poisson_blocks.zero(); poisson_rhs.zero(); auto ap=graph.template view<1>(poisson_blocks.data(),poisson_rhs.data());
        launch(e,SimpleInterior<1>{mesh,state,controls,ap}); launch(b,SimpleBoundary<1>{mesh,state,controls,ap});
        solve<1>(poisson,poisson_solve,ap,phi,error.host()[0]!=0,{});
        observe("raw_pressure_increment",phi);
        launch(n,SimpleAddIncrement{state.pressure,phi.data(),controls.alpha_p});
        gradient<1>(mesh,state,phi.data(),sum,gp.data());
        poison(gp,3);
        // Match the reference: new p, predicted u, old grad(p), then reversal and velocity correction.
        div.zero(); launch(e,SimpleInterior<1>{mesh,state,controls,ap,true});
        launch(b,SimpleBoundary<1>{mesh,state,controls,ap,completed>0,true}); check("flux update failed");
        poison(div,1);
        launch(owned_nodes,OnList<SimpleCorrectVelocity>{{state,gp.data()},owned.data()});
        exchange({{velocity.data(),3},{pressure.data(),1}});
        observe("pressure",pressure); observe("velocity",velocity); observe("interior_flux",eflux);
        observe("boundary_flux",bflux); observe("mass_divergence",div);
        ++completed; assembled=false;
    }
    template<class Observer=NoSimpleObserver> void step(Observer observe={}) { assemble_momentum(); advance(observe); }
};
} // namespace mars::segregated::runtime
#undef MARS_DSIMPLE_HD
