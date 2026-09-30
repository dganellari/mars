// Host CPU/MPI rectangular-duct run through the production distributed SIMPLE code:
//   duct lattice (duct_mesh.hpp, identical to duct_mesh.py's Exodus mesh)
//   -> test slab partition in ElementDomain's shape (slab_view)
//   -> production simple_partition (solver ids, star check, exterior faces, face owners)
//   -> production DistributedSimpleRunner (host build), linear solves by duct_host_solver.hpp.
// SIMPLE controls go through the production option parser, so every --rho/--mu/... flag means
// what it means for mars_segregated_simple; the files written (-metrics.csv, -fields.csv) and
// the final CONVERGED line match that executable, so duct_compare.py reads both alike.
//
//   mpirun -n P duct_host_run --cells N [--length 7 --width 2 --height 1 --stretch 2]
//          --output-prefix PREFIX [production SIMPLE options]
//   duct_host_run --cells N --dump-mesh FILE        canonical text (compare with duct_mesh.py --dump)
#include "duct_mesh.hpp"
#include "duct_host_solver.hpp"
#include "mars_segregated_simple_options.hpp"
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace mars::segregated;
using namespace mars::segregated::runtime;
using Runner=DistributedSimpleRunner<duct::HostMatrix,long long,duct::GatheredSolve>;

struct DuctOptions {
    int cells=0; double length=7, width=2, height=1, stretch=2;
    std::string dump;
    SimpleOptions simple;
};
DuctOptions parse(int argc,char** argv) {
    DuctOptions o; std::vector<std::string> rest{"duct_host_run","--mesh=duct-lattice"};
    for (int i=1;i<argc;++i) {
        const std::string key=argv[i];
        auto value=[&]() { ensure(i+1<argc,"missing option value"); return std::string(argv[++i]); };
        if (key=="--cells") o.cells=std::stoi(value());
        else if (key=="--length") o.length=std::stod(value());
        else if (key=="--width") o.width=std::stod(value());
        else if (key=="--height") o.height=std::stod(value());
        else if (key=="--stretch") o.stretch=std::stod(value());
        else if (key=="--dump-mesh") o.dump=value();
        else rest.push_back(key);
    }
    if (o.dump.empty()) {
        std::vector<char*> args; for (auto& word:rest) args.push_back(word.data());
        o.simple=simple_options(int(args.size()),args.data());
        ensure(!o.simple.pressure_tolerances,"pressure linear overrides require the production Hypre driver");
    }
    return o;
}

int execute(const DuctOptions& d) {
    int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    const auto lattice=duct::Lattice::make(d.cells,d.length,d.width,d.height,d.stretch);
    if (!d.dump.empty()) {
        if (rank==0) {
            ensure(!std::filesystem::exists(d.dump),"dump exists; choose a fresh name");
            std::ofstream out(d.dump,std::ios::binary); out<<duct::canonical(lattice); ensure(bool(out),"dump failed");
        }
        return 0;
    }
    const auto& o=d.simple;
    if (rank==0) for (const char* suffix:{"-metrics.csv","-fields.csv"})
        ensure(!std::filesystem::exists(o.output+suffix),"output exists; choose a fresh prefix");

    auto view=duct::slab_view(lattice,rank,ranks);
    const std::vector<long long> key=view.key;
    auto kind=[&](const int* face) { const long long g[3]={key[face[0]],key[face[1]],key[face[2]]}; return lattice.kind(g); };
    auto part=simple_partition<long long>(MPI_COMM_WORLD,view,kind);

    // Global coverage, as the production driver checks it: every node and duct face owned once.
    long long mine[4]={(long long)part.ownership.owned_nodes.size(),0,0,0}, total[4]={};
    for (int i:part.ownership.owned_faces) ++mine[1+part.input.faces[i].kind];
    MPI_Allreduce(mine,total,4,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    const auto faces=lattice.boundary_faces();
    simple_collective(MPI_COMM_WORLD,total[0]==lattice.nodes() && total[1]==faces[0] && total[2]==faces[1] && total[3]==faces[2],
                      "ownership does not cover every node and duct boundary face exactly once");

    Runner run(MPI_COMM_WORLD,part.input,part.ownership,o.controls);
    std::ofstream csv;
    const auto& c=o.controls;
    if (rank==0) {
        csv.open(o.output+"-metrics.csv"); ensure(bool(csv),"cannot write metrics");
        csv<<std::setprecision(17)<<"iteration,momentum,continuity,mass_balance,du,dp,dflux,cancellation,inlet_kg_s,outlet_kg_s,umax_m_s,closed_faces,changed_faces\n";
        std::cout<<std::setprecision(17)<<"SIMPLE Tet4 duct (host CPU/MPI oracle), "<<ranks<<" ranks, "
                 <<(c.high_resolution?"high-resolution":"upwind")<<", laminar\n"
                 <<"lattice "<<lattice.nx<<"x"<<lattice.ny<<"x"<<lattice.nz<<" L="<<lattice.length<<" W="<<lattice.width<<" H="<<lattice.height
                 <<" nodes="<<lattice.nodes()<<" tets="<<lattice.elements()<<'\n'
                 <<"rho="<<c.density<<" mu="<<c.viscosity<<" inlet_speed="<<c.inlet_speed<<" outlet_pressure="<<c.pressure_reference
                 <<" reference_length="<<c.reference_length<<" alpha_u="<<c.alpha_u<<" alpha_p="<<c.alpha_p
                 <<" alpha_mass="<<c.alpha_mass<<" beta="<<c.beta<<" pseudo_dt="<<c.pseudo_dt<<'\n';
    }
    const auto start=std::chrono::steady_clock::now();
    bool converged=false;
    for (;;) {
        run.assemble_momentum();
        const auto sums=run.diagnostics(); const auto m=simple_metrics(sums,c);
        ensure(m.finite,"nonfinite nonlinear diagnostics");
        ensure(m.cancellation<=1e-10,"assembled continuity does not match boundary mass flux");
        converged=simple_converged(m,run.completed,sums.changed,o.residual,o.mass,o.change);
        if (rank==0) {
            csv<<run.completed<<','<<m.momentum<<','<<m.continuity<<','<<m.flux<<','<<m.velocity_change<<','
               <<m.pressure_change<<','<<m.flux_change<<','<<m.cancellation<<','<<sums.inlet<<','<<sums.outlet<<','<<std::sqrt(sums.speed2)<<','
               <<sums.closed<<','<<sums.changed<<'\n';
            if (run.completed%o.report==0 || converged || run.completed==o.iterations)
                std::cout<<"[simple] iteration="<<run.completed<<" momentum="<<m.momentum<<" continuity="<<m.continuity
                         <<" balance="<<m.flux<<" du="<<m.velocity_change<<" dp="<<m.pressure_change<<" dflux="<<m.flux_change
                         <<" umax="<<std::sqrt(sums.speed2)<<" closed="<<sums.closed<<" changed="<<sums.changed
                         <<" gmres(u,p)="<<run.momentum_solve.stats.last<<','<<run.poisson_solve.stats.last<<std::endl;
        }
        if (converged || run.completed==o.iterations) break;
        run.advance();
    }
    if (rank==0) { csv.close(); ensure(bool(csv),"metric output failed"); }

    // Owned nodes to rank 0, one row per lattice node (the production driver's CSV layout).
    const auto u=run.velocity.host(), p=run.pressure.host(), x=run.x.host(), y=run.y.host(), z=run.z.host();
    std::vector<double> rows;
    for (int v:part.ownership.owned_nodes) {
        rows.push_back(double(key[v])); rows.push_back(x[v]); rows.push_back(y[v]); rows.push_back(z[v]);
        for (int j=0;j<3;++j) rows.push_back(u[3*v+j]);
        rows.push_back(p[v]);
    }
    const int count=int(rows.size()); std::vector<int> counts(ranks), displs(ranks);
    MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    long long all=0; if (rank==0) for (int q=0;q<ranks;++q) { displs[q]=int(all); all+=counts[q]; }
    std::vector<double> gathered(std::size_t(rank==0?all:0));
    MPI_Gatherv(rows.data(),count,MPI_DOUBLE,gathered.data(),counts.data(),displs.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (rank==0) {
        ensure(gathered.size()==8*std::size_t(lattice.nodes()),"field gather does not cover every node");
        std::vector<const double*> order(std::size_t(lattice.nodes()),nullptr);
        for (std::size_t i=0;i<gathered.size();i+=8) {
            const long long g=(long long)gathered[i];
            ensure(g>=0 && g<lattice.nodes() && !order[g],"field gather returned duplicate or missing nodes");
            order[g]=&gathered[i];
        }
        std::ofstream out(o.output+"-fields.csv"); ensure(bool(out),"cannot write fields");
        out<<std::setprecision(17)<<"node,x,y,z,u,v,w,p\n";
        for (const double* r:order) { for (int j=0;j<8;++j) out<<(j?",":"")<<r[j]; out<<'\n'; }
        ensure(bool(out),"field output failed");
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"linear solves: momentum "<<run.momentum_solve.stats.solves<<" (mean GMRES "
                 <<double(run.momentum_solve.stats.iterations)/std::max(1LL,run.momentum_solve.stats.solves)<<"), pressure "
                 <<run.poisson_solve.stats.solves<<" (mean GMRES "
                 <<double(run.poisson_solve.stats.iterations)/std::max(1LL,run.poisson_solve.stats.solves)<<"); "<<seconds<<" s\n";
        std::cout<<(converged?"CONVERGED":"NOT CONVERGED: iteration limit")<<" iterations="<<run.completed<<" ranks="<<ranks
                 <<" exchange_rounds="<<run.exchange.rounds()<<'\n';
    }
    return converged?0:2;
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int result=1;
    try { result=execute(parse(argc,argv)); }
    catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Bcast(&result,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Finalize();
    return result;
}
