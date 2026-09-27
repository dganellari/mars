// Native Exodus -> device ElementDomain -> distributed SIMPLE. Host work is file I/O and API control.
#include "mars_segregated_simple_native_mesh.hpp"
#include <filesystem>
#include <iomanip>
#include <limits>

using namespace mars;
using namespace mars::segregated;
using namespace mars::segregated::runtime;
using Runner=DistributedSimpleRunner<HypreSimpleSolve<1>::Solver::Matrix,HYPRE_BigInt,HypreSimpleSolve>;

struct Options { std::string mesh,output; int iterations=2000,report=10; double residual=1e-6,mass=1e-6,change=1e-6; bool setup_only=false; };
Options options(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;i+=2) {
        ensure(i+1<argc,"each option requires a value");
        const std::string k=argv[i], v=argv[i+1];
        if (k=="--mesh") o.mesh=v; else if (k=="--output-prefix") o.output=v;
        else if (k=="--mesh-format") ensure(v=="exodus","distributed native SIMPLE requires --mesh-format exodus");
        else if (k=="--iterations") o.iterations=std::stoi(v); else if (k=="--report-every") o.report=std::stoi(v);
        else if (k=="--residual-tol") o.residual=std::stod(v); else if (k=="--mass-tol") o.mass=std::stod(v);
        else if (k=="--change-tol") o.change=std::stod(v);
        else if (k=="--setup-only") { ensure(v=="0" || v=="1","--setup-only expects 0 or 1"); o.setup_only=v=="1"; }
        else throw std::runtime_error("unknown option: "+k);
    }
    ensure(!o.mesh.empty() && !o.output.empty() && o.iterations>0 && o.report>0,"--mesh and --output-prefix are required");
    return o;
}

struct FieldRow { double values[8]; };
struct FieldRowLess {
    __host__ __device__ bool operator()(const FieldRow& a,const FieldRow& b) const { return a.values[0]<b.values[0]; }
};
struct PackOutput {
    const int *owned,*source; const double *x,*y,*z,*u,*p; FieldRow* rows;
    __device__ void operator()(int i) const {
        const int n=owned[i]; rows[i]={{double(source[n]),x[n],y[n],z[n],u[3*n],u[3*n+1],u[3*n+2],p[n]}};
    }
};
struct CheckOutput {
    const FieldRow* rows; int* error;
    __device__ void operator()(int i) const { if (rows[i].values[0]!=double(i)) atomicExch(error,1); }
};

int execute(const Options& o) {
    int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if (rank==0) for (const char* suffix:{"-metrics.csv","-fields.csv"})
        ensure(!std::filesystem::exists(o.output+suffix),"output exists; choose a fresh prefix");
    int nodes=0;
    Buffer<int> source_node;
    auto make_runner=[&]() {
        const auto input=read_simple_mesh(MPI_COMM_WORLD,o.mesh); nodes=int(input.x.size());
        auto domain=distribute_simple_mesh(MPI_COMM_WORLD,input);
        NativeSimpleMesh<HYPRE_BigInt> native(MPI_COMM_WORLD,*domain,input);
        source_node=std::move(native.source_node);
        return std::make_unique<Runner>(MPI_COMM_WORLD,native.partition.input,native.partition.ownership);
    };
    auto runner=make_runner(); // Release replicated file arrays and setup scratch before iterating.
    auto& run=*runner;
    if (o.setup_only) {
        if (rank==0) std::cout<<"PASS: ElementDomain SIMPLE setup ranks="<<ranks<<"; no iterations run"<<std::endl;
        return 0;
    }
    std::ofstream csv;
    if (rank==0) {
        csv.open(o.output+"-metrics.csv"); ensure(bool(csv),"cannot write metrics");
        csv<<std::setprecision(17)<<"iteration,momentum,continuity,mass_balance,du,dp,dflux,cancellation,inlet_kg_s,outlet_kg_s,umax_m_s,closed_faces,changed_faces\n";
        std::cout<<"SIMPLE public Tet4 channel, "<<ranks<<" ranks (ElementDomain/cstone), upwind, laminar, rho=1 mu=0.1 U=0.1 L=1\n";
    }
    bool converged=false;
    for (;;) {
        run.assemble_momentum(); const auto report=run.diagnostic_report(o.residual,o.mass,o.change);
        const auto& sums=report.sums; const auto& m=report.metrics;
        ensure(m.finite,"nonfinite nonlinear diagnostics");
        ensure(m.cancellation<=1e-10,"assembled continuity does not match boundary mass flux");
        converged=report.converged;
        if (rank==0) {
            csv<<run.completed<<','<<m.momentum<<','<<m.continuity<<','<<m.flux<<','<<m.velocity_change<<','
               <<m.pressure_change<<','<<m.flux_change<<','<<m.cancellation<<','<<sums.inlet<<','<<sums.outlet<<','<<report.speed<<','
               <<sums.closed<<','<<sums.changed<<'\n';
            if (run.completed%o.report==0 || converged || run.completed==o.iterations)
                std::cout<<"[simple] iteration="<<run.completed<<" momentum="<<m.momentum<<" continuity="<<m.continuity
                         <<" balance="<<m.flux<<" du="<<m.velocity_change<<" dp="<<m.pressure_change<<" dflux="<<m.flux_change
                         <<" umax="<<report.speed<<" closed="<<sums.closed<<" changed="<<sums.changed<<std::endl;
        }
        if (converged || run.completed==o.iterations) break;
        run.advance();
    }
    // Output is the only field download. MPI gathers device rows before rank zero writes the CSV.
    Buffer<FieldRow> rows(run.owned_nodes);
    launch(run.owned_nodes,PackOutput{run.owned.data(),raw(source_node),run.x.data(),run.y.data(),run.z.data(),
                                     run.velocity.data(),run.pressure.data(),raw(rows)});
    static_assert(sizeof(FieldRow)==8*sizeof(double));
    simple_collective(MPI_COMM_WORLD,run.owned_nodes<=INT_MAX/8,"field output exceeds MPI count capacity");
    const int count=run.owned_nodes*8; std::vector<int> counts(ranks),displacements(ranks);
    ensure(MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS,"field counts failed");
    long long all=0;
    if (!rank) for (int q=0;q<ranks;++q) { displacements[q]=int(all); all+=counts[q]; ensure(all<=INT_MAX,"field output exceeds MPI count capacity"); }
    simple_collective(MPI_COMM_WORLD,rank!=0 || all==8LL*nodes,"field gather does not cover every source node");
    Buffer<FieldRow> gathered(size_t(rank==0?nodes:0));
    mesh_mpi_ready();
    ensure(MPI_Gatherv(raw(rows),count,MPI_DOUBLE,raw(gathered),counts.data(),displacements.data(),MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS,"device field gather failed");
    if (rank==0) {
        mesh_sort(gathered,FieldRowLess{});
        Array<int> output_error(1);
        launch(nodes,CheckOutput{raw(gathered),output_error.data()});
        ensure(output_error.host()[0]==0,"field gather returned duplicate or missing source nodes");
        std::vector<FieldRow> field(gathered.size()); thrust::copy(gathered.begin(),gathered.end(),field.begin());
        std::ofstream out(o.output+"-fields.csv"); ensure(bool(out),"cannot write fields");
        out<<std::setprecision(17)<<"node,x,y,z,u,v,w,p\n";
        for (const auto& row:field) { for (int j=0;j<8;++j) out<<(j?",":"")<<row.values[j]; out<<'\n'; }
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
        if (cudaGetDeviceCount(&devices)!=cudaSuccess || devices<=0 || cudaSetDevice(local%devices)!=cudaSuccess) {
            std::cerr<<"ERROR: cannot select a CUDA device"<<std::endl;
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    int result=1;
    try { result=execute(options(argc,argv)); }
    catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Bcast(&result,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Finalize();
    return result;
}
