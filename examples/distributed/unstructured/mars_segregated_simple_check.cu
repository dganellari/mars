#include "mars_segregated_simple_runtime.hpp"
#include "fixture.hpp"
using namespace mars::segregated;
using namespace mars::segregated::runtime;
using namespace simple_check;

void execute(const Fixture& f) {
    SimpleRunner run(f);
    for (int it=0;it<2;++it) {
        const auto& expected=f.expected[it]; std::string label="iteration="+std::to_string(it+1)+" ";
        run.step([&](const char* stage,Array<double>& field) {
            const std::string name=stage;
            const auto* oracle=&expected.mass_divergence;
            if (name=="momentum") oracle=&expected.momentum;
            else if (name=="influence") oracle=&expected.influence;
            else if (name=="trace") oracle=&expected.trace;
            else if (name=="raw_pressure_increment") oracle=&expected.increment;
            else if (name=="pressure") oracle=&expected.pressure;
            else if (name=="velocity") oracle=&expected.velocity;
            else if (name=="interior_flux") oracle=&expected.interior_flux;
            else if (name=="boundary_flux") oracle=&expected.boundary_flux;
            else require(name=="mass_divergence","unknown observed stage");
            compare(field.host(),*oracle,label+name);
        });
        auto hdiv=run.div.host(),hb=run.bflux.host();
        double balance=std::accumulate(hb.begin(),hb.end(),0.),total=std::accumulate(hdiv.begin(),hdiv.end(),0.);
        require(std::abs(balance-total)<=1e-12,"global flux cancellation failed");
        std::cout<<label<<"net_boundary_mass="<<balance<<" sum_continuity="<<total<<'\n';
    }
}

void reversal_check() {
    SimpleInput f; f.x={0,1,0,0}; f.y={0,0,1,0}; f.z={0,0,0,1};
    for (int k=0;k<4;++k) f.nodes[k]={k}; f.faces={{0,1,1}};
    SimpleRunner run(f);
    auto fill=[](Array<double>& a,double value) {
        Array<double> source(std::vector<double>(a.values.size(),value)); a.copy_from(source);
    };
    auto flags=[&](int expected) { for (int value:run.flags.host()) require(value==expected,"reversal flag mismatch"); };
    auto zero=[](Array<double>& a) { for (double value:a.host()) require(value==0,"closed opening has nonzero contribution"); };
    auto flux=[&]() { launch(1,SimpleBoundary<1>{run.mesh,run.state,run.controls,{},true,true}); run.check("reversal kernel failed"); };
    fill(run.velocity,-1); flux(); flags(1); zero(run.bflux);
    auto am=run.graph.view<3>(run.momentum.blocks.data(),run.momentum.rhs.data());
    auto ap=run.graph.view<1>(run.poisson.blocks.data(),run.poisson.rhs.data());
    launch(1,SimpleBoundary<3>{run.mesh,run.state,run.controls,am,true,false});
    launch(1,SimpleBoundary<1>{run.mesh,run.state,run.controls,ap,true,false});
    run.check("closed boundary assembly failed");
    zero(run.momentum.blocks); zero(run.momentum.rhs); zero(run.poisson.blocks); zero(run.poisson.rhs);
    fill(run.trace,2);
    launch(1,SimpleTraceMoment{run.mesh,run.state,run.moment.data()}); zero(run.moment);
    launch(1,SimpleTrace{run.mesh,run.state,run.controls,run.moment.data()});
    for (double value:run.trace.host()) require(value==2,"closed trace changed");
    fill(run.velocity,1); flux(); flags(1); zero(run.bflux);
    fill(run.pressure,3); flux(); flags(0); zero(run.bflux);
    launch(1,SimpleTraceMoment{run.mesh,run.state,run.moment.data()});
    launch(1,SimpleTrace{run.mesh,run.state,run.controls,run.moment.data()});
    for (double value:run.trace.host()) require(std::abs(value)<1e-13,"reopened trace mismatch");
    flux(); for (double q:run.bflux.host()) require(std::abs(q-.375)<1e-13,"reopened flux mismatch");
    fill(run.d,2); launch(1,SimpleBoundary<1>{run.mesh,run.state,run.controls,ap,true,false});
    run.check("reopened boundary assembly failed");
    const auto values=run.poisson.blocks.host();
    require(std::abs(std::accumulate(values.begin(),values.end(),0.)-3)<1e-13,"reopened pressure anchor mismatch");
    // Exercise the device reduction too; the integration driver must not dispatch it on the host.
    auto sums=reduce_sums(1,SimpleFaceSums{run.mesh,run.state,run.old_flags.data()});
    require(sums.invalid==0 && std::abs(sums.outlet-1.125)<1e-13,"boundary reduction mismatch");
    auto changes=reduce_sums(3,SimpleFluxChange{run.bflux.data(),run.old_bflux.data()});
    require(changes.invalid==0 && std::abs(changes.flux_change-.375)<1e-13,"flux-change reduction mismatch");
    fill(run.factor,1); fill(run.momentum.rhs,1.5);
    const auto node_sums=reduce_sums(4,SimpleNodeSums{run.state,run.momentum.rhs.data(),run.old_velocity.data(),run.old_pressure.data()});
    auto near=[](double a,double b) { require(std::abs(a-b)<1e-12*std::max(1.,std::abs(b)),"node reduction mismatch"); };
    require(node_sums.invalid==0,"invalid node reduction"); near(node_sums.volume,1./6);
    near(node_sums.momentum2,1152); near(node_sums.continuity2,10.125); near(node_sums.continuity,1.125);
    near(node_sums.velocity_change2,.5); near(node_sums.pressure_change2,1.5); near(node_sums.speed2,3);
    std::cout<<"PASS: native outlet close/retain/reopen kernels and pressure anchor\n";
}
int main(int argc,char** argv) {
#ifdef MARS_REPLAY_CUDA
    MPI_Init(&argc,&argv); int ranks=0; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
#endif
    int result=0;
    try {
#ifdef MARS_REPLAY_CUDA
        require(ranks==1,"this SIMPLE gate supports exactly one rank");
#endif
        require(argc==2,"usage: mars_segregated_simple_check PUBLIC_SIMPLE_FILE"); execute(load(argv[1])); reversal_check();
#ifdef MARS_REPLAY_CUDA
        const char* backend="CUDA/Hypre";
#else
        const char* backend="host/direct";
#endif
        std::cout<<"PASS: "<<backend<<" two native SIMPLE iterations scalar_checks="<<checks
                 <<"; no injected stage fields; MPI and nonlinear convergence not tested\n";
    } catch (const std::exception& error) { std::cerr<<"FAIL: "<<error.what()<<'\n'; result=1; }
#ifdef MARS_REPLAY_CUDA
    MPI_Finalize();
#endif
    return result;
}
