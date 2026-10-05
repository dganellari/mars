// Native Exodus -> device ElementDomain -> distributed SIMPLE. Host work is file I/O and API control.
#include "mars_segregated_simple_native_mesh.hpp"
#include "mars_segregated_simple_output.hpp"
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace mars;
using namespace mars::segregated;
using namespace mars::segregated::runtime;
using Runner=DistributedSimpleRunner<HypreSimpleSolve<1>::Solver::Matrix,HYPRE_BigInt,HypreSimpleSolve>;

int execute(const SimpleOptions& o) {
    int rank=0, ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    simple_output_preflight(MPI_COMM_WORLD,o.output,o.field_output);
    const double setup_start=MPI_Wtime();
    int nodes=0;
    Buffer<int> source_node;
    auto make_runner=[&]() {
        const auto input=read_simple_mesh_root(MPI_COMM_WORLD,o.mesh,o.boundaries); nodes=input.global_nodes;
        auto domain=distribute_simple_mesh(MPI_COMM_WORLD,input);
        NativeSimpleMesh<HYPRE_BigInt> native(MPI_COMM_WORLD,*domain,input);
        source_node=std::move(native.source_node);
        return std::make_unique<Runner>(MPI_COMM_WORLD,native.partition.input,native.partition.ownership,o.controls);
    };
    auto runner=make_runner(); // Release root file arrays and setup scratch before iterating.
    auto& run=*runner;
    run.set_pressure_tolerances(o.pressure_tolerances,o.pressure_rtol,o.pressure_atol);
    run.profile.configure(o.profile,o.profile_warmup);
    run.exchange.enable_profiling(o.profile);
    run.overlap_assembly=o.halo_overlap;
    if (o.linear_cache) { run.momentum_solve.solver.enable_fixed_graph_updates(); run.poisson_solve.solver.enable_fixed_graph_updates(); }
    run.momentum_solve.solver.enable_timing(o.profile); run.poisson_solve.solver.enable_timing(o.profile);
    const double setup_seconds=MPI_Wtime()-setup_start;
    if (o.setup_only) {
        if (rank==0) std::cout<<"PASS: ElementDomain SIMPLE setup ranks="<<ranks<<"; no iterations run"<<std::endl;
        return 0;
    }
    std::ofstream csv;
    if (rank==0) {
        csv.open(o.output+"-metrics.csv"); ensure(bool(csv),"cannot write metrics");
        csv<<std::setprecision(17)<<"iteration,momentum,continuity,mass_balance,du,dp,dflux,cancellation,inlet_kg_s,outlet_kg_s,umax_m_s,closed_faces,changed_faces\n";
        const auto& c=o.controls;
        std::cout<<std::setprecision(17)<<"SIMPLE Tet4, "<<ranks<<" ranks (ElementDomain/cstone), "
                 <<(c.high_resolution?"high-resolution":"upwind")<<", laminar\n"
                 <<"velocity_interpolation="<<(c.velocity_shifted?"linear-linear":"trilinear")<<'\n'
                 <<"rho="<<c.density<<" mu="<<c.viscosity<<" nu="<<c.viscosity/c.density
                 <<" inlet_speed="<<c.inlet_speed<<" (inward normal) outlet_pressure="<<c.pressure_reference
                 <<" reference_length="<<c.reference_length<<'\n'
                 <<"alpha_u="<<c.alpha_u<<" alpha_p="<<c.alpha_p<<" alpha_mass="<<c.alpha_mass
                 <<" beta="<<c.beta<<" pseudo_dt="<<c.pseudo_dt<<" (steady, no physical time)\n"
                 <<"linear_cache="<<o.linear_cache<<" halo_overlap="<<o.halo_overlap
                 <<" field_output="<<o.field_output<<" profile="<<o.profile<<'\n'
                 <<"Norms are dimensionless MARS residuals; not OpenAccel printed RMS.\n";
        if (o.pressure_tolerances)
            std::cout<<"pressure_linear_rtol="<<o.pressure_rtol<<" pressure_linear_atol="<<o.pressure_atol
                     <<" pressure_acceptance=max(atol,rtol*rhs_norm); momentum targets unchanged\n";
    }
    if (o.profile) {
        char host[MPI_MAX_PROCESSOR_NAME]; int length=0,device=0;
        ensure(MPI_Get_processor_name(host,&length)==MPI_SUCCESS,"cannot read profile host name");
        assembly_cuda_check(cudaGetDevice(&device));
        cudaDeviceProp properties{}; assembly_cuda_check(cudaGetDeviceProperties(&properties,device));
        std::ostringstream line;
        line<<"[simple-profile] rank="<<rank<<" host="<<std::string(host,length)
            <<" device="<<device<<" model="<<properties.name<<'\n';
        std::cout<<line.str();
    }
    Buffer<FieldRow> snapshot_rows(o.snapshot_iterations?run.owned_nodes:0);
    double snapshot_seconds=0;
    const double iteration_start=MPI_Wtime();
    distributed::FieldExchange::Profile halo_baseline;
    bool converged=false;
    for (;;) {
        run.assemble_momentum(); const auto report=run.diagnostic_report(o.residual,o.mass,o.change);
        run.profile.collect(run.completed);
        if (o.profile && run.completed<=o.profile_warmup) halo_baseline=run.exchange.profile();
        const auto& sums=report.sums; const auto& m=report.metrics;
        ensure(m.finite,"nonfinite nonlinear diagnostics");
        ensure(m.cancellation<=1e-10,"assembled continuity does not match boundary mass flux");
        converged=report.converged;
        if (rank==0) {
            csv<<run.completed<<','<<m.momentum<<','<<m.continuity<<','<<m.flux<<','<<m.velocity_change<<','
               <<m.pressure_change<<','<<m.flux_change<<','<<m.cancellation<<','<<sums.inlet<<','<<sums.outlet<<','<<report.speed<<','
               <<sums.closed<<','<<sums.changed<<'\n';
            ensure(bool(csv),"metric output failed");
            if (run.completed%o.report==0 || converged || run.completed==o.iterations)
                std::cout<<"[simple] iteration="<<run.completed<<" momentum="<<m.momentum<<" continuity="<<m.continuity
                         <<" balance="<<m.flux<<" du="<<m.velocity_change<<" dp="<<m.pressure_change<<" dflux="<<m.flux_change
                         <<" umax="<<report.speed<<" closed="<<sums.closed<<" changed="<<sums.changed<<std::endl;
        }
        if (o.snapshot_iterations && run.completed<=o.snapshot_iterations) {
            const double start=MPI_Wtime();
            const auto prefix=simple_snapshot_prefix(o.output,run.completed);
            simple_output_preflight(MPI_COMM_WORLD,prefix,o.field_output);
            launch(run.owned_nodes,PackOutput{run.owned.data(),raw(source_node),run.x.data(),run.y.data(),run.z.data(),
                                             run.velocity.data(),run.pressure.data(),raw(snapshot_rows)});
            write_simple_fields(MPI_COMM_WORLD,prefix,o.field_output,nodes,snapshot_rows);
            snapshot_seconds+=MPI_Wtime()-start;
        }
        if (converged || run.completed==o.iterations) break;
        run.advance();
        run.profile.record_linear(3,run.momentum_solve.solver,run.completed);
        run.profile.record_linear(1,run.poisson_solve.solver,run.completed);
    }
    if (!rank) { csv.close(); ensure(bool(csv),"metric output failed"); }
    const double iteration_seconds=MPI_Wtime()-iteration_start-snapshot_seconds;
    const double output_start=MPI_Wtime();
    if (o.field_output!="none") {
        Buffer<FieldRow> rows(run.owned_nodes);
        launch(run.owned_nodes,PackOutput{run.owned.data(),raw(source_node),run.x.data(),run.y.data(),run.z.data(),
                                         run.velocity.data(),run.pressure.data(),raw(rows)});
        write_simple_fields(MPI_COMM_WORLD,o.output,o.field_output,nodes,rows);
    }
    const double output_seconds=MPI_Wtime()-output_start+snapshot_seconds;
    if (!rank) std::cout<<(converged?"CONVERGED":"NOT CONVERGED: iteration limit")<<" iterations="<<run.completed<<" ranks="<<ranks
                       <<" exchange_rounds="<<run.exchange.rounds()<<'\n';
    double wall[3]={setup_seconds,iteration_seconds,output_seconds},maximum[3];
    ensure(MPI_Reduce(wall,maximum,3,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD)==MPI_SUCCESS,"run timing reduction failed");
    if (!rank) std::cout<<"[simple-time] scope=rank_max_wall setup_seconds="<<maximum[0]<<" iteration_seconds="<<maximum[1]
                       <<" output_seconds="<<maximum[2]<<" seconds_per_iteration="<<(run.completed?maximum[1]/run.completed:0)<<'\n';
    run.profile.write(MPI_COMM_WORLD,std::cout);
    if (o.profile) {
        const auto& h=run.exchange.profile();
        double local[4]={h.pack_seconds-halo_baseline.pack_seconds,h.pack_ready_wait_seconds-halo_baseline.pack_ready_wait_seconds,
                         h.wait_seconds-halo_baseline.wait_seconds,h.unpack_seconds-halo_baseline.unpack_seconds},global[4];
        long long counts[3]={h.rounds-halo_baseline.rounds,h.sent_bytes-halo_baseline.sent_bytes,h.received_bytes-halo_baseline.received_bytes},totals[3];
        ensure(MPI_Reduce(local,global,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD)==MPI_SUCCESS,"halo timing reduction failed");
        ensure(MPI_Reduce(counts,totals,3,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD)==MPI_SUCCESS,"halo count reduction failed");
        int builds[6]={run.momentum_solve.solver.get_graph_build_count(),run.momentum_solve.solver.get_numeric_update_count(),run.momentum_solve.solver.get_setup_count(),
                       run.poisson_solve.solver.get_graph_build_count(),run.poisson_solve.solver.get_numeric_update_count(),run.poisson_solve.solver.get_setup_count()},max_builds[6];
        ensure(MPI_Reduce(builds,max_builds,6,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD)==MPI_SUCCESS,"linear count reduction failed");
        if (!rank) {
            std::cout<<"[simple-profile] halo_scope=post_warmup_application_only pack_seconds="<<global[0]<<" pack_ready_wait_seconds="<<global[1]
                     <<" mpi_wait_seconds="<<global[2]<<" unpack_seconds="<<global[3]<<" rank_rounds_sum="<<totals[0]
                     <<" sent_bytes_sum="<<totals[1]<<" received_bytes_sum="<<totals[2]<<'\n';
            for (int c=0;c<2;++c) std::cout<<"[simple-profile] linear="<<(c?"pressure":"momentum")<<" scope=lifetime_rank_max graph_builds="<<max_builds[3*c]
                                        <<" numeric_updates="<<max_builds[3*c+1]<<" setups="<<max_builds[3*c+2]<<'\n';
        }
    }
    return converged?0:2;
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    SimpleOptions parsed;
    try { parsed=simple_options(argc,argv); }
    catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    if (parsed.help) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (!rank) std::cout<<simple_help();
        MPI_Finalize(); return 0;
    }
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
    try { result=execute(parsed); }
    catch (const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Bcast(&result,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Finalize();
    return result;
}
