// Multi-rank public-channel SIMPLE driver on ElementDomain (cstone): each rank passes a slice of
// the elements, cstone's sync distributes them by SFC and builds the halos; domain_view +
// simple_partition turn the domain into DistributedSimpleRunner input. Same options, metrics and
// CSV formats as mars_segregated_simple, so a one-rank and a P-rank run compare line by line.
// Validation driver: every rank reads the (small) public mesh to restore exact coordinates and
// boundary tags by SFC key; production ingestion supplies those without replication.
#include "mars_segregated_simple_partition.hpp"
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
using namespace mars;
using namespace mars::segregated;
using namespace mars::segregated::runtime;
using Domain=ElementDomain<TetTag,double,uint64_t,cstone::GpuTag>;
using Runner=DistributedSimpleRunner<HypreSimpleSolve<1>::Solver::Matrix,HYPRE_BigInt,HypreSimpleSolve>;

struct Options { std::string mesh,output; int iterations=2000,report=10; double residual=1e-6,mass=1e-6,change=1e-6; };
Options options(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;i+=2) {
        ensure(i+1<argc,"each option requires a value");
        const std::string k=argv[i], v=argv[i+1];
        if (k=="--mesh") o.mesh=v; else if (k=="--output-prefix") o.output=v;
        else if (k=="--iterations") o.iterations=std::stoi(v); else if (k=="--report-every") o.report=std::stoi(v);
        else if (k=="--residual-tol") o.residual=std::stod(v); else if (k=="--mass-tol") o.mass=std::stod(v);
        else if (k=="--change-tol") o.change=std::stod(v); else throw std::runtime_error("unknown option: "+k);
    }
    ensure(!o.mesh.empty() && !o.output.empty() && o.iterations>0 && o.report>0,"--mesh and --output-prefix are required");
    return o;
}

int execute(const Options& o) {
    int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if (rank==0) for (const char* suffix:{"-metrics.csv","-fields.csv"})
        ensure(!std::filesystem::exists(o.output+suffix),"output exists; choose a fresh prefix");
    const auto input=load_simple_input(o.mesh.c_str());
    const int nodes=int(input.x.size()), elements=int(input.nodes[0].size());

    // Every rank passes all coordinates and its slice of elements; cstone distributes by SFC.
    std::vector<double> hx(input.x), hy(input.y), hz(input.z);
    const int begin=int((long long)elements*rank/ranks), end=int((long long)elements*(rank+1)/ranks);
    std::array<std::vector<uint64_t>,4> conn;
    for (int k=0;k<4;++k) for (int el=begin;el<end;++el) conn[k].push_back(uint64_t(input.nodes[k][el]));
    Domain domain(std::make_tuple(hx,hy,hz),std::make_tuple(conn[0],conn[1],conn[2],conn[3]),rank,ranks);

    // Exact coordinates and source node ids by SFC key (Tet4 domains store decoded coordinates).
    auto view=domain_view(domain);
    std::vector<std::array<double,3>> coords(nodes);
    for (int g=0;g<nodes;++g) coords[g]={input.x[g],input.y[g],input.z[g]};
    const auto local=domain.resolveSideSetNodesToLocalKeepMisses(coords);
    std::vector<int> source(view.x.size(),-1);
    for (int g=0;g<nodes;++g) if (local[g]>=0) { source[local[g]]=g; view.x[local[g]]=input.x[g]; view.y[local[g]]=input.y[g]; view.z[local[g]]=input.z[g]; }
    simple_collective(MPI_COMM_WORLD,std::find(source.begin(),source.end(),-1)==source.end(),"a local node did not match a public mesh node by SFC key");

    // Boundary tags from the public face list, keyed by sorted source node ids.
    std::map<std::array<int,3>,int> tags;
    for (const auto& f:input.faces) {
        std::array<int,3> k3; for (int j=0;j<3;++j) k3[j]=input.nodes[tet_face_node(f.ordinal,j)][f.element];
        std::sort(k3.begin(),k3.end()); tags[k3]=f.kind;
    }
    auto kind=[&](const int* face) {
        std::array<int,3> k3{source[face[0]],source[face[1]],source[face[2]]}; std::sort(k3.begin(),k3.end());
        auto it=tags.find(k3); return it==tags.end()?-1:it->second;
    };
    auto part=simple_partition<HYPRE_BigInt>(MPI_COMM_WORLD,view,kind);

    // Global coverage: every node owned once, every public boundary face owned once.
    long long mine[4]={(long long)part.ownership.owned_nodes.size(),0,0,0}, total[4]={};
    for (int i:part.ownership.owned_faces) ++mine[1+part.input.faces[i].kind];
    MPI_Allreduce(mine,total,4,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    long long kinds[3]={0,0,0}; for (const auto& f:input.faces) ++kinds[f.kind];
    simple_collective(MPI_COMM_WORLD,total[0]==nodes && total[1]==kinds[0] && total[2]==kinds[1] && total[3]==kinds[2],
                      "ownership does not cover every node and boundary face exactly once");

    Runner run(MPI_COMM_WORLD,part.input,part.ownership);
    std::ofstream csv;
    if (rank==0) {
        csv.open(o.output+"-metrics.csv"); ensure(bool(csv),"cannot write metrics");
        csv<<std::setprecision(17)<<"iteration,momentum,continuity,mass_balance,du,dp,dflux,cancellation,inlet_kg_s,outlet_kg_s,umax_m_s,closed_faces,changed_faces\n";
        std::cout<<"SIMPLE public Tet4 channel, "<<ranks<<" ranks (ElementDomain/cstone), upwind, laminar, rho=1 mu=0.1 U=0.1 L=1\n";
    }
    bool converged=false;
    for (;;) {
        run.assemble_momentum(); const auto sums=run.diagnostics(); const auto m=simple_metrics(sums,run.controls);
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
                         <<" umax="<<std::sqrt(sums.speed2)<<" closed="<<sums.closed<<" changed="<<sums.changed<<std::endl;
        }
        if (converged || run.completed==o.iterations) break;
        run.advance();
    }
    // Owned nodes to rank 0, written in public node order (same file as the one-rank driver).
    const auto u=run.velocity.host(), p=run.pressure.host();
    std::vector<double> rows;
    for (int v:part.ownership.owned_nodes) { rows.push_back(source[v]); for (int j=0;j<3;++j) rows.push_back(u[3*v+j]); rows.push_back(p[v]); }
    int count=int(rows.size()); std::vector<int> counts(ranks), displs(ranks);
    MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    int all=0; if (rank==0) for (int q=0;q<ranks;++q) { displs[q]=all; all+=counts[q]; }
    std::vector<double> gathered(std::size_t(rank==0?all:0));
    MPI_Gatherv(rows.data(),count,MPI_DOUBLE,gathered.data(),counts.data(),displs.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (rank==0) {
        std::vector<std::array<double,4>> field(nodes); std::vector<int> seen(nodes,0);
        for (std::size_t i=0;i+4<gathered.size();i+=5) { const int g=int(gathered[i]); ++seen[g]; field[g]={gathered[i+1],gathered[i+2],gathered[i+3],gathered[i+4]}; }
        ensure(std::count(seen.begin(),seen.end(),1)==nodes,"field gather did not return every node exactly once");
        std::ofstream out(o.output+"-fields.csv"); ensure(bool(out),"cannot write fields");
        out<<std::setprecision(17)<<"node,x,y,z,u,v,w,p\n";
        for (int g=0;g<nodes;++g) out<<g<<','<<input.x[g]<<','<<input.y[g]<<','<<input.z[g]<<','<<field[g][0]<<','<<field[g][1]<<','<<field[g][2]<<','<<field[g][3]<<'\n';
        ensure(bool(out),"field output failed");
        std::cout<<(converged?"CONVERGED":"NOT CONVERGED: iteration limit")<<" iterations="<<run.completed<<" ranks="<<ranks
                 <<" exchange_rounds="<<run.exchange.rounds()<<'\n';
    }
    return converged?0:2;
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    {
        int rank=0, local=0, devices=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm node; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        MPI_Comm_rank(node,&local); MPI_Comm_free(&node);
        if (cudaGetDeviceCount(&devices)==cudaSuccess && devices>0) cudaSetDevice(local%devices);
    }
    int result=1;
    try { result=execute(options(argc,argv)); }
    catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Bcast(&result,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Finalize();
    return result;
}
