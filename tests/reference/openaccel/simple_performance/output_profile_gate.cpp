#include "mars_segregated_simple_output.hpp"
#include "mars_segregated_simple_profile.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
using namespace mars::segregated::runtime;
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        int rank,ranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
#ifdef MARS_REPLAY_CUDA
        int local,devices=0; MPI_Comm node;
        ensure(MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node)==MPI_SUCCESS,"cannot create node communicator");
        MPI_Comm_rank(node,&local); MPI_Comm_free(&node);
        assembly_cuda_check(cudaGetDeviceCount(&devices)); ensure(devices>0,"no CUDA device");
        assembly_cuda_check(cudaSetDevice(local%devices));
#endif
        ensure(argc==2,"output directory required");
        std::filesystem::create_directories(argv[1]);
        const std::string prefix=std::string(argv[1])+"/field";
        // Interleaved source IDs ensure partition output cannot assume a contiguous range.
        Buffer<FieldRow> rows(2);
        std::vector<FieldRow> file(2);
        for (int j=0;j<2;++j) { const int n=rank+(1-j)*ranks; file[j]={{double(n),double(n),0,0,.1*n,0,0,2.*n}}; }
        rows.assign(file.begin(),file.end());
        write_simple_fields(MPI_COMM_WORLD,prefix,"distributed",2*ranks,rows);
        write_simple_fields(MPI_COMM_WORLD,prefix,"gathered",2*ranks,rows);
        ensure(std::filesystem::exists(simple_part_path(prefix,rank)),"missing local output");
        if (!rank) ensure(std::filesystem::exists(prefix+"-fields.json") && std::filesystem::exists(prefix+"-fields.csv"),"missing output manifest or reference");
        bool rejected=false;
        try { simple_output_preflight(MPI_COMM_WORLD,prefix,"distributed"); } catch (const std::runtime_error&) { rejected=true; }
        ensure(rejected,"existing output was not rejected collectively");
        SimpleProfile profile;
        { auto timing=profile.scope(SimpleProfile::assembly); }
        profile.collect(0); ensure(!profile.samples && !profile.totals[0].calls,"disabled profiler performed work");
        profile.configure(true,1);
        for (int iteration=0;iteration<4;++iteration) {
            { auto timing=profile.scope(SimpleProfile::assembly); }
            { auto timing=profile.scope(SimpleProfile::diagnostics); }
            profile.collect(iteration);
        }
        ensure(profile.samples==2 && profile.totals[0].calls==2 && profile.totals[1].calls==2,"profile warmup or counts failed");
        std::ostringstream report; profile.write(MPI_COMM_WORLD,report);
        if (!rank) ensure(report.str().find("samples=2")!=std::string::npos,"missing rank-max profile");
        if (!rank) std::cout<<"PASS: distributed field output and optional phase profile ranks="<<ranks<<'\n';
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
