#pragma once
// Distributed SIMPLE: SimpleRunner's iteration, unchanged kernels, on a rank's partition.
//
// Inputs (from native ingestion and element-halo completion; not computed here):
//   local mesh     every element touching an owned node (complete stars), every boundary
//                  face of those elements, and all their nodes (owned and ghost, any order);
//   ownership      owned nodes in solver order, solver node per local node, the elements and
//                  boundary faces this rank owns uniquely (the owner must hold the face);
//   halo lists     NodeHaloTopology-style peers and send/recv node lists.
//
// Schedule per iteration (fields published owner -> ghost, 4 rounds, 14 values per ghost;
// in addition to setup metadata checks; high-resolution adds 12 values to the first round):
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
// and element samples, then allreduce. Callers must abort the communicator on unexpected
// rank-local exceptions; exchange failures abort before peers can use invalid buffers.
#include "mars_segregated_simple_runtime.hpp"
#include "mars_segregated_distributed_matrix.hpp"
#include "mars_segregated_halo_exchange.hpp"
#include "mars_segregated_simple_reduction.hpp"
#include "mars_segregated_simple_profile.hpp"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#ifdef MARS_REPLAY_CUDA
#include <thrust/logical.h>
#include <thrust/partition.h>
#include <thrust/sequence.h>
#endif
#if defined(__CUDACC__) && !defined(MARS_REPLAY_CUDA)
#error "Distributed SIMPLE: CUDA builds define MARS_REPLAY_CUDA, as the SIMPLE drivers do"
#endif
#if defined(__CUDACC__)
#define MARS_DSIMPLE_HD __host__ __device__
#else
#define MARS_DSIMPLE_HD
#endif

namespace mars::segregated::runtime {

template<class T> using HostStorage=std::vector<T>;
template<class GlobalId,template<class> class Storage=HostStorage> struct SimpleOwnership {
    Storage<int> owned_nodes, owned_elements, owned_faces;
    Storage<GlobalId> solver_node;
    std::vector<int> peers, send_offsets, recv_offsets;
    Storage<int> send_nodes, recv_nodes;
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

struct LocalMomentumEntity {
    SimpleMesh mesh; const unsigned char* owned; bool boundary;
    MARS_DSIMPLE_HD bool operator()(int i) const {
        const int element=boundary?mesh.faces[i].element:i;
        // Boundary reconstruction may use the fourth node, not just the face nodes.
        for (int k=0;k<4;++k) if (!owned[mesh.nodes[k][element]]) return false;
        return true;
    }
};

struct CheckTraceMoment {
    const double* moment; int* error;
    MARS_DSIMPLE_HD void operator()(int) const {
        if (!(distributed::finite_value(moment[0]) && distributed::finite_value(moment[1]) && moment[1]>0))
            distributed::raise_fault(error,1);
    }
};
struct OutsideRange {
    int size;
    MARS_DSIMPLE_HD bool operator()(int i) const { return i<0 || i>=size; }
};
inline bool valid_indices(const std::vector<int>& indices,int size) {
    return std::none_of(indices.begin(),indices.end(),OutsideRange{size});
}
#ifdef MARS_REPLAY_CUDA
inline bool valid_indices(const thrust::device_vector<int>& indices,int size) {
    return !thrust::any_of(indices.begin(),indices.end(),OutsideRange{size});
}
#endif
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
    using Solver=mars::fem::HypreGMRESSolver<double,int,cstone::execution::Gpu>;
    Solver solver;
    typename Solver::Vector b,x;
    explicit HypreSimpleSolve(MPI_Comm comm):solver(comm,2000,1e-12,Solver::BOOMERAMG,100) {
        solver.setVerbose(false); solver.setPointBlock(C);
        const distributed::Tolerance acceptance;
        solver.enable_true_residual_check(acceptance.absolute,acceptance.relative);
        // A one-level hierarchy needs l1 row norms too; the default coarse type omits them.
        // l1-Jacobi also keeps coarse relaxation on the device.
        solver.setAMGCoarseRelaxType(18);
    }
    void set_tolerances(double relative,double absolute) {
        solver.set_stopping_tolerances(relative,absolute);
        solver.enable_true_residual_check(absolute,relative,true);
    }
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
    int n,e,b,completed=0,limiter_iteration=-1,owned_nodes,owned_elements,owned_faces;
    int local_momentum_elements=0,local_momentum_faces=0;
    bool overlap_assembly=true;
    bool poison_unexchanged=false; // validation: NaN in ghost entries that must never be read
    SimpleControls controls;
    Array<double> x,y,z,velocity,pressure,vg,pg,d,volume,div,eflux,bflux,trace,factor,sum,gp,moment,old_velocity,old_pressure,old_eflux,old_bflux;
    Array<double> blend,blend_lower,blend_upper,blend_candidate;
    Array<int> n0,n1,n2,n3,error,flags,old_flags,owned,elements,boundary,momentum_elements,momentum_faces;
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
    SimpleProfile profile;
#ifdef MARS_REPLAY_CUDA
    SimpleDeviceReduction reduction;
    Array<SimpleReport> report{1};
#endif
    bool assembled=false;
    distributed::Tolerance tolerance{1e-13,1e-10};
    std::optional<distributed::Tolerance> pressure_tolerance;

    // All ranks enter, including those without the override, before either solve.
    void set_pressure_tolerances(bool enabled,double relative,double absolute) {
        simple_collective(comm,completed==0 && std::isfinite(relative) && std::isfinite(absolute)
            && relative>0 && relative<1 && absolute>=0,"invalid or late pressure tolerances");
        const double local[6]={double(enabled),-double(enabled),relative,-relative,absolute,-absolute};
        double global[6];
        if (MPI_Allreduce(local,global,6,MPI_DOUBLE,MPI_MAX,comm)!=MPI_SUCCESS) {
            MPI_Abort(comm,1); throw std::runtime_error("pressure tolerance reduction failed");
        }
        simple_collective(comm,global[0]==-global[1] && global[2]==-global[3] && global[4]==-global[5],
            "pressure tolerances differ between ranks");
        if (enabled) {
            poisson_solve.set_tolerances(relative,absolute);
            pressure_tolerance=distributed::Tolerance{absolute,relative,true};
        }
    }

    // Id may be wider than GlobalId (e.g. int64 ids for a 32-bit HYPRE_BigInt build); ids
    // that do not fit are rejected collectively instead of narrowed.
    template<class Input,class Id,template<class> class Storage> DistributedSimpleRunner(MPI_Comm c,const Input& f,const SimpleOwnership<Id,Storage>& o,SimpleControls ctl={},
        distributed::EmptyRanks empty=distributed::EmptyRanks::reject):
        comm(c),n(int(f.x.size())),e(int(f.nodes[0].size())),b(int(f.faces.size())),
        owned_nodes(int(o.owned_nodes.size())),owned_elements(int(o.owned_elements.size())),owned_faces(int(o.owned_faces.size())),controls(ctl),
        x(f.x),y(f.y),z(f.z),velocity(3*n),pressure(n),vg(9*n),pg(3*n),d(3*n),volume(n),div(n),
        eflux(6*e),bflux(3*b),trace(3*b),factor(n),sum(9*n),gp(3*n),moment(2),old_velocity(3*n),old_pressure(n),old_eflux(6*e),old_bflux(3*b),
        blend(ctl.high_resolution?3*n:0),blend_lower(ctl.high_resolution?3*n:0),blend_upper(ctl.high_resolution?3*n:0),blend_candidate(ctl.high_resolution?3*n:0),
        n0(f.nodes[0]),n1(f.nodes[1]),n2(f.nodes[2]),n3(f.nodes[3]),error(1),flags(3*b),old_flags(3*b),
        owned(o.owned_nodes),elements(o.owned_elements),boundary(o.owned_faces),momentum_elements(e),momentum_faces(b),owned_mask(std::size_t(n)),solver_node(narrow(c,o.solver_node)),
        faces(f.faces),geometry(e),
        mesh{n,e,b,{n0.data(),n1.data(),n2.data(),n3.data()},x.data(),y.data(),z.data(),faces.data(),geometry.data()},
        state{velocity.data(),pressure.data(),vg.data(),pg.data(),d.data(),volume.data(),div.data(),
              eflux.data(),bflux.data(),trace.data(),factor.data(),error.data(),flags.data(),blend.data()},
        graph(mesh),
        momentum_blocks(std::size_t(graph.blocks())*9),momentum_rhs(3*n),poisson_blocks(std::size_t(graph.blocks())),poisson_rhs(n),du(3*n),phi(n),
        momentum(c,graph.template view<3>(momentum_blocks.data(),momentum_rhs.data()),owned.data(),owned_nodes,solver_node.data(),solver_node.values.size(),empty),
        poisson(c,graph.template view<1>(poisson_blocks.data(),poisson_rhs.data()),owned.data(),owned_nodes,solver_node.data(),solver_node.values.size(),empty),
        momentum_solve(c),poisson_solve(c),
        exchange(c,o.peers,o.send_offsets,o.send_nodes,o.recv_offsets,o.recv_nodes,n,ctl.high_resolution?15:8)
    {
        simple_collective(comm,valid_simple_controls(ctl),"invalid SIMPLE controls");
        simple_collective(comm,valid_indices(o.owned_elements,e) && valid_indices(o.owned_faces,b),
                          "owned element or face list outside the local mesh");
        launch(owned_nodes,MarkOwnedNodes{owned.data(),owned_mask.data()});
        local_momentum_elements=partition_momentum(momentum_elements,false);
        local_momentum_faces=partition_momentum(momentum_faces,true);
        launch(e,SimpleGeometry{mesh,state}); check("native geometry failed");
        launch(b,SimpleBoundaryFactor{mesh,factor.data()});
        // Owned volumes and boundary factors are complete. Ghost entries stay partial on purpose:
        // only owned gradients and limiter bounds are used before publication.
    }
    int partition_momentum(Array<int>& list,bool boundary) {
        const LocalMomentumEntity local{mesh,owned_mask.data(),boundary};
#ifdef MARS_REPLAY_CUDA
        thrust::sequence(list.values.begin(),list.values.end());
        return int(thrust::stable_partition(list.values.begin(),list.values.end(),local)-list.values.begin());
#else
        std::iota(list.values.begin(),list.values.end(),0);
        return int(std::stable_partition(list.values.begin(),list.values.end(),local)-list.values.begin());
#endif
    }
    template<class Id> static std::vector<GlobalId> narrow(MPI_Comm c,const std::vector<Id>& ids) {
        bool fits=true;
        for (Id v:ids) fits=fits && (long double)v>=(long double)std::numeric_limits<GlobalId>::min()
                                  && (long double)v<=(long double)std::numeric_limits<GlobalId>::max();
        simple_collective(c,fits,"solver node id does not fit the solver's global index type");
        return std::vector<GlobalId>(ids.begin(),ids.end());
    }
#ifdef MARS_REPLAY_CUDA
    template<class Id> static thrust::device_vector<GlobalId> narrow(MPI_Comm c,const thrust::device_vector<Id>& ids) {
        // Solver ids are nonnegative, and the wrapper only supports int-sized global row ranges.
        simple_collective(c,!thrust::any_of(ids.begin(),ids.end(),IdTooWide<Id>{}),
                          "solver node id does not fit the solver's global index type");
        return thrust::device_vector<GlobalId>(ids.begin(),ids.end());
    }
    template<class Id> struct IdTooWide {
        __host__ __device__ bool operator()(Id id) const {
            return id<0 || static_cast<unsigned long long>(id)>static_cast<unsigned long long>(std::numeric_limits<GlobalId>::max());
        }
    };
#endif
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
        auto timing=profile.scope(SimpleProfile::assembly);
        gradient<3>(mesh,state,state.velocity,sum,state.velocity_gradient,controls.velocity_shifted);
        if (controls.high_resolution && limiter_iteration!=completed) {
            const auto a=graph.template view<3>(nullptr,nullptr);
            launch(owned_nodes,OnList<SimpleBlendBounds>{{a,state.velocity,blend_lower.data(),blend_upper.data(),blend_candidate.data(),state.error},owned.data()});
            const SimpleBlendSamples samples{mesh,state,blend_lower.data(),blend_upper.data(),blend_candidate.data(),owned_mask.data()};
            launch(e,SimpleBlendInterior{samples}); launch(b,SimpleBlendBoundary{samples});
            launch(owned_nodes,OnList<SimpleBlendFinish>{{blend_candidate.data(),blend.data()},owned.data()});
            limiter_iteration=completed;
        }
        gradient<1>(mesh,state,state.pressure,sum,state.pressure_gradient);
        momentum_blocks.zero(); momentum_rhs.zero(); auto am=graph.template view<3>(momentum_blocks.data(),momentum_rhs.data());
        poison(pg,3);
        if (controls.high_resolution) {
            poison(vg,9); poison(blend,3);
            // Owners have complete stars: compute locally, then batch both fields with grad(p).
            exchange.begin({{pg.data(),3},{vg.data(),9},{blend.data(),3}});
        } else {
            poison(vg,9,false); // Upwind multiplies the unused gradient by zero.
            exchange.begin({{pg.data(),3}});
        }
        if (overlap_assembly) {
            launch(local_momentum_elements,OnList<SimpleInterior<3>>{{mesh,state,controls,am},momentum_elements.data()});
            launch(local_momentum_faces,OnList<SimpleBoundary<3>>{{mesh,state,controls,am,completed>0,false},momentum_faces.data()});
            exchange.end();
            if (e>local_momentum_elements)
                launch(e-local_momentum_elements,OnList<SimpleInterior<3>>{{mesh,state,controls,am},momentum_elements.data()+local_momentum_elements});
            if (b>local_momentum_faces)
                launch(b-local_momentum_faces,OnList<SimpleBoundary<3>>{{mesh,state,controls,am,completed>0,false},momentum_faces.data()+local_momentum_faces});
        } else {
            exchange.end();
            launch(e,SimpleInterior<3>{mesh,state,controls,am});
            launch(b,SimpleBoundary<3>{mesh,state,controls,am,completed>0,false});
        }
        launch(owned_nodes,OnList<SimpleMomentumNode>{{state,controls,am},owned.data()});
        check("momentum assembly failed"); assembled=true;
    }
#ifdef MARS_REPLAY_CUDA
    void reduce_diagnostics() {
        ensure(assembled,"diagnostics require a fresh momentum assembly");
        reduction.reduce(owned_nodes,OnList<SimpleNodeSums>{{state,momentum_rhs.data(),old_velocity.data(),old_pressure.data()},owned.data()},0);
        reduction.reduce(owned_faces,OnList<SimpleFaceSums>{{mesh,state,old_flags.data()},boundary.data()},1);
        reduction.reduce(6*owned_elements,OnSamples<SimpleFluxChange>{{eflux.data(),old_eflux.data()},elements.data(),6},2);
        reduction.reduce(3*owned_faces,OnSamples<SimpleFluxChange>{{bflux.data(),old_bflux.data()},boundary.data(),3},3);
        reduction.finish(comm);
    }
    SimpleSums diagnostics() {
        auto timing=profile.scope(SimpleProfile::diagnostics);
        reduce_diagnostics(); return reduction.result.host()[0];
    }
    SimpleReport diagnostic_report(double residual,double mass,double change) {
        auto timing=profile.scope(SimpleProfile::diagnostics);
        reduce_diagnostics();
        launch(1,FinishSimpleReport{reduction.result.data(),controls,completed,residual,mass,change,report.data()});
        return report.host()[0]; // Fixed-size logging/convergence report; no field or reduction staging.
    }
#else
    SimpleSums diagnostics() {
        auto timing=profile.scope(SimpleProfile::diagnostics);
        ensure(assembled,"diagnostics require a fresh momentum assembly"); // same program order on every rank
        auto a=reduce_sums(owned_nodes,OnList<SimpleNodeSums>{{state,momentum_rhs.data(),old_velocity.data(),old_pressure.data()},owned.data()});
        a=SimpleSumCombine{}(a,reduce_sums(owned_faces,OnList<SimpleFaceSums>{{mesh,state,old_flags.data()},boundary.data()}));
        a=SimpleSumCombine{}(a,reduce_sums(6*owned_elements,OnSamples<SimpleFluxChange>{{eflux.data(),old_eflux.data()},elements.data(),6}));
        a=SimpleSumCombine{}(a,reduce_sums(3*owned_faces,OnSamples<SimpleFluxChange>{{bflux.data(),old_bflux.data()},boundary.data(),3}));
        return allreduce_sums(comm,a);
    }
#endif
    template<int C,class System,class Solver> void solve(System& system,Solver& solver,BlockCsrView<C> view,Array<double>& increment,
        bool assembly_failed,std::initializer_list<distributed::Field> with) {
        auto timing=profile.scope(C==3?SimpleProfile::momentum:SimpleProfile::pressure);
        system.update(view,solver.rhs(std::size_t(system.rows())),std::size_t(system.rows()),assembly_failed);
        // Keep a candidate on nonconvergence so both matrix representations can be checked.
        const bool solved=solver(system);
        const auto context=[&](const char* reason) {
            return std::string(C==3?"momentum":"pressure correction")+" at SIMPLE iteration "+std::to_string(completed+1)+": "+reason;
        };
        const bool candidate=solver.size()==std::size_t(system.rows()) && (!system.rows() || solver.solution());
        int verdict[2]={candidate?1:0,solved?1:0};
        if (MPI_Allreduce(MPI_IN_PLACE,verdict,2,MPI_INT,MPI_MIN,comm)!=MPI_SUCCESS) {
            MPI_Abort(comm,1); throw std::runtime_error("linear verdict reduction failed");
        }
        if (!verdict[0]) throw std::runtime_error(context("linear solver returned no usable candidate"));
        system.unpack(solver.solution(),solver.size(),increment.data(),increment.values.size());
        if (with.size()==0) exchange({{increment.data(),C}});
        else { auto it=with.begin(); exchange({{increment.data(),C},*it}); }
        const auto acceptance=C==1 && pressure_tolerance?*pressure_tolerance:tolerance;
        const auto norms=system.residual(distributed::halo_complete(increment.data(),increment.values.size()),solver.rhs(),acceptance);
        // Check rejected candidates too, before discarding the only independent evidence.
        // Neither a solver rejection nor a failed MARS residual can become a success.
        if (!verdict[1] || !norms.passed) {
            int rank=0; MPI_Comm_rank(comm,&rank);
            if (!rank) std::cerr<<"[simple-linear] stage="<<(C==3?"momentum":"pressure")
                <<" iteration="<<completed+1<<" solver_accepted="<<verdict[1]
                <<" mars_absolute_residual="<<norms.absolute()<<" rhs_norm="<<std::sqrt(norms.rhs2)
                <<" acceptance_limit="<<acceptance.limit(std::sqrt(norms.rhs2))
                <<" mars_passed="<<norms.passed<<'\n';
            if constexpr (C==1) {
                const char* option=std::getenv("MARS_SIMPLE_PRESSURE_AUDIT");
                const int enabled=option && std::string(option)=="1";
                int any=0;
                if (MPI_Allreduce(&enabled,&any,1,MPI_INT,MPI_MAX,comm)!=MPI_SUCCESS) {
                    MPI_Abort(comm,1); throw std::runtime_error("pressure audit selection failed");
                }
                if (any) {
                    const auto audit=system.pressure_audit(distributed::halo_complete(increment.data(),increment.values.size()),solver.rhs(),acceptance);
                    if (!rank) std::cerr<<"[simple-pressure-audit] finite="<<audit.finite
                        <<" zero_row="<<audit.zero_row<<" nonpositive_diagonal="<<audit.nonpositive_diagonal
                        <<" positive_offdiagonal="<<audit.positive_offdiagonal
                        <<" constant_mode_detected="<<audit.constant_mode_detected
                        <<" residual_within_roundoff_bound="<<audit.residual_within_roundoff_bound
                        <<" roundoff_bound_exceeds_limit="<<audit.roundoff_bound_exceeds_limit
                        <<" compensated_residual_finite="<<audit.compensated_residual_finite
                        <<" compensated_residual_passed="<<audit.compensated_residual_passed<<'\n';
                }
            }
            throw std::runtime_error(context(verdict[1]?"true linear residual failed":"linear solve failed"));
        }
    }
    template<class Observer=NoSimpleObserver> void advance(Observer observe={}) {
        ensure(assembled,"advance requires momentum assembly");
        old_velocity.copy_from(velocity); old_pressure.copy_from(pressure); old_flags.copy_from(flags);
        old_eflux.copy_from(eflux); old_bflux.copy_from(bflux);
        poison(momentum_rhs,3);
        solve<3>(momentum,momentum_solve,graph.template view<3>(momentum_blocks.data(),momentum_rhs.data()),du,false,{{d.data(),3}});
        launch(3*n,SimpleAddIncrement{state.velocity,du.data(),1});
        observe("momentum",velocity); observe("influence",d);
        {
            auto timing=profile.scope(SimpleProfile::outlet);
            moment.zero(); launch(owned_faces,OnList<SimpleTraceMoment>{{mesh,state,moment.data()},boundary.data()});
#ifdef MARS_REPLAY_CUDA
            assembly_cuda_check(cudaStreamSynchronize(nullptr));
            ensure(MPI_Allreduce(MPI_IN_PLACE,moment.data(),2,MPI_DOUBLE,MPI_SUM,comm)==MPI_SUCCESS,"outlet moment reduction failed");
            launch(1,CheckTraceMoment{moment.data(),error.data()});
            check("all outlet faces closed: no open pressure anchor or nonfinite outlet moments");
#else
            const auto local=moment.host(); std::vector<double> global(2);
            MPI_Allreduce(local.data(),global.data(),2,MPI_DOUBLE,MPI_SUM,comm);
            if (!(std::isfinite(global[0]) && std::isfinite(global[1]) && global[1]>0))
                throw std::runtime_error("all outlet faces closed: no open pressure anchor; cannot solve this prescribed-inflow case");
            moment.values=global;
#endif
            launch(b,SimpleTrace{mesh,state,controls,moment.data()}); observe("trace",trace);
        }
        auto ap=graph.template view<1>(poisson_blocks.data(),poisson_rhs.data());
        {
            auto timing=profile.scope(SimpleProfile::pressure_assembly);
            poisson_blocks.zero(); poisson_rhs.zero();
            launch(e,SimpleInterior<1>{mesh,state,controls,ap}); launch(b,SimpleBoundary<1>{mesh,state,controls,ap});
        }
        solve<1>(poisson,poisson_solve,ap,phi,error.host()[0]!=0,{});
        observe("raw_pressure_increment",phi);
        {
            auto timing=profile.scope(SimpleProfile::correction);
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
        }
        ++completed; assembled=false;
    }
    template<class Observer=NoSimpleObserver> void step(Observer observe={}) { assemble_momentum(); advance(observe); }
};
} // namespace mars::segregated::runtime
#undef MARS_DSIMPLE_HD
