#include "mars_segregated_simple_runtime.hpp"
#include "mars_segregated_simple_input.hpp"
#include "mars_segregated_native_input.hpp"
#include <filesystem>
#include <iomanip>
#include <limits>
using namespace mars::segregated;
using namespace mars::segregated::runtime;

struct Options {
    std::string mesh,output,format="prepared";
    int iterations=2000,report=10;
    double residual=1e-6,mass=1e-6,change=1e-6;
};
Options options(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;++i) {
        std::string key=argv[i];
        ensure(i+1<argc,"each option requires a value; see --help"); const std::string value=argv[++i];
        if (key=="--mesh") o.mesh=value;
        else if (key=="--mesh-format") { ensure(value=="prepared" || value=="exodus","mesh format must be prepared or exodus"); o.format=value; }
        else if (key=="--output-prefix") o.output=value;
        else {
            std::size_t end=0; const double number=std::stod(value,&end);
            ensure(end==value.size() && std::isfinite(number) && number>0,"options require positive finite numbers");
            if (key=="--iterations" || key=="--report-every") {
                ensure(number<=std::numeric_limits<int>::max() && number==std::floor(number),"iteration counts must be integers");
                (key=="--iterations"?o.iterations:o.report)=int(number);
            } else if (key=="--residual-tol") o.residual=number;
            else if (key=="--mass-tol") o.mass=number;
            else if (key=="--change-tol") o.change=number;
            else throw std::runtime_error("unknown option: "+key);
        }
    }
    ensure(!o.mesh.empty() && !o.output.empty(),"--mesh and --output-prefix are required"); return o;
}
void save_fields(const Options& o,const std::vector<int>& source_nodes,SimpleRunner& run) {
    // Explicit final public-field export; no field downloads inside the iteration loop.
    const auto u=run.velocity.host(),p=run.pressure.host(),x=run.x.host(),y=run.y.host(),z=run.z.host();
    std::ofstream out(o.output+"-fields.csv"); ensure(bool(out),"cannot write fields");
    out<<std::setprecision(17)<<"node,x,y,z,u,v,w,p\n";
    for (int n=0;n<run.n;++n) out<<source_nodes[n]<<','<<x[n]<<','<<y[n]<<','<<z[n]<<','
        <<u[3*n]<<','<<u[3*n+1]<<','<<u[3*n+2]<<','<<p[n]<<'\n';
    out.close(); ensure(bool(out),"field output failed");
}
template<class Input>
int solve(const Options& o,const Input& input,const std::vector<int>& source_nodes) {
    SimpleRunner run(input); run.momentum.verbose=run.poisson.verbose=false;
    std::ofstream csv(o.output+"-metrics.csv"); ensure(bool(csv),"cannot write metrics");
    csv<<std::setprecision(17)<<"iteration,momentum,continuity,mass_balance,du,dp,dflux,cancellation,inlet_kg_s,outlet_kg_s,umax_m_s,closed_faces,changed_faces\n";
    std::cout<<"SIMPLE public Tet4 channel, one rank, upwind, laminar, rho=1 mu=0.1 U=0.1 L=1\n"
             <<"alpha_u=0.3 alpha_p=0.3 alpha_mass=0.75 beta=0.05 pseudo_dt=0.01 (steady, no physical time)\n"
             <<"Norms are dimensionless MARS residuals; not the OpenAccel printed RMS normalization.\n";
    bool converged=false;
    for (;;) {
        // Evaluate the current state before advancing it; no mixed-iteration stopping test.
        run.assemble_momentum(); const auto sums=run.diagnostics(); const auto m=simple_metrics(sums,run.controls);
        ensure(m.finite,"nonfinite nonlinear diagnostics");
        ensure(m.cancellation<=1e-10,"assembled continuity does not match boundary mass flux");
        converged=simple_converged(m,run.completed,sums.changed,o.residual,o.mass,o.change);
        csv<<run.completed<<','<<m.momentum<<','<<m.continuity<<','<<m.flux<<','<<m.velocity_change<<','
           <<m.pressure_change<<','<<m.flux_change<<','<<m.cancellation<<','<<sums.inlet<<','<<sums.outlet<<','<<sqrt(sums.speed2)<<','
           <<sums.closed<<','<<sums.changed<<'\n';
        ensure(bool(csv),"metric output failed");
        if (run.completed%o.report==0 || converged || run.completed==o.iterations)
            std::cout<<"[simple] iteration="<<run.completed<<" momentum="<<m.momentum<<" continuity="<<m.continuity
                     <<" balance="<<m.flux<<" du="<<m.velocity_change<<" dp="<<m.pressure_change
                     <<" dflux="<<m.flux_change<<" umax="<<sqrt(sums.speed2)<<" closed="<<sums.closed<<" changed="<<sums.changed<<std::endl;
        if (converged || run.completed==o.iterations) break;
        run.advance();
    }
    csv.close(); ensure(bool(csv),"metric output failed"); save_fields(o,source_nodes,run);
    std::cout<<(converged?"CONVERGED":"NOT CONVERGED: iteration limit")<<" iterations="<<run.completed<<'\n';
    return converged?0:2;
}
int execute(const Options& o) {
    for (const char* suffix:{"-metrics.csv","-fields.csv"})
        ensure(!std::filesystem::exists(o.output+suffix),"output exists; choose a fresh prefix");
    if (o.format=="exodus") {
#ifdef MARS_REPLAY_CUDA
        NativeSimpleInput input(o.mesh);
        std::cout<<"Native SIMPLE Exodus input: "<<std::filesystem::canonical(o.mesh).string()<<'\n';
        // Exodus storage rows identify final output; the comparator applies node_num_map.
        std::vector<int> ids(input.source_node.size());
        thrust::copy(input.source_node.begin(),input.source_node.end(),ids.begin());
        return solve(o,input,ids);
#else
        throw std::runtime_error("native ElementDomain input requires the CUDA build");
#endif
    }
    const auto input=load_simple_input(o.mesh.c_str());
    std::vector<int> ids(input.x.size()); std::iota(ids.begin(),ids.end(),0);
    return solve(o,input,ids);
}
int main(int argc,char** argv) {
    if (argc==2 && std::string(argv[1])=="--help") {
        std::cout<<"mars_segregated_simple --mesh PUBLIC_MESH --output-prefix PATH [--iterations 2000] "
                 <<"[--mesh-format prepared|exodus] "
                 <<"[--report-every 10] [--residual-tol 1e-6] [--mass-tol 1e-6] [--change-tol 1e-6]\n"; return 0;
    }
#ifdef MARS_REPLAY_CUDA
    MPI_Init(&argc,&argv); int ranks=0; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if (ranks!=1) { std::cerr<<"SIMPLE public driver currently requires one rank\n"; MPI_Finalize(); return 1; }
#endif
    int result=1;
    try { result=execute(options(argc,argv)); }
    catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; }
#ifdef MARS_REPLAY_CUDA
    MPI_Finalize();
#endif
    return result;
}
