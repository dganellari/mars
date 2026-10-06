#include "mars_segregated_simple_output.hpp"
#include "mars_segregated_simple_profile.hpp"
#include "mars_segregated_simple_audit.hpp"
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
        mars::segregated::assembly_cuda_check(cudaGetDeviceCount(&devices)); ensure(devices>0,"no CUDA device");
        mars::segregated::assembly_cuda_check(cudaSetDevice(local%devices));
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
        for (int step=0;step<3;++step) {
            const auto snapshot=simple_snapshot_prefix(prefix,step);
            simple_output_preflight(MPI_COMM_WORLD,snapshot,"distributed");
            for (auto& row:file) row.values[7]+=1.;
            rows.assign(file.begin(),file.end());
            write_simple_fields(MPI_COMM_WORLD,snapshot,"distributed",2*ranks,rows);
            std::ifstream saved(simple_part_path(snapshot,rank));
            std::string line; std::getline(saved,line);
            int count=0;
            while (std::getline(saved,line)) {
                std::istringstream input(line); double value[8];
                for (int j=0;j<8;++j) { input>>value[j]; if (j!=7) ensure(input.get()==',',"invalid snapshot CSV"); }
                ensure(value[7]==2.*value[0]+step+1,"snapshot contains stale state"); ++count;
            }
            ensure(count==2,"snapshot lost owned rows");
        }
        const std::vector<int> h_source{2*ranks-1-rank,rank},h_owned{1,0},h_offsets{0,2,4},h_columns{0,1,0,1};
        Buffer<int> source(h_source.begin(),h_source.end()),owned(h_owned.begin(),h_owned.end()),
                    offsets(h_offsets.begin(),h_offsets.end()),columns(h_columns.begin(),h_columns.end());
        Buffer<double> blocks(36,2.),vector(6,3.),pressure(4,4.),scalar(2,5.);
        SimpleFirstStepAudit audit(MPI_COMM_WORLD,prefix,2,2,4,raw(source),raw(owned),raw(offsets),raw(columns));
        audit.momentum({2,raw(offsets),raw(columns),raw(blocks),raw(vector)},raw(vector),raw(vector),raw(vector));
        audit.pressure({2,raw(offsets),raw(columns),raw(pressure),raw(scalar)},raw(scalar));
        audit.finish(raw(vector));
        std::ostringstream audit_path; audit_path<<prefix<<"-audit-rank"<<std::setw(6)<<std::setfill('0')<<rank<<".bin";
        std::ifstream binary(audit_path.str(),std::ios::binary);
        std::uint64_t header[7]; binary.read(reinterpret_cast<char*>(header),sizeof(header));
        ensure(header[0]==0x4d53415544495431ULL && header[1]==1 && header[2]==std::uint64_t(rank) && header[3]==std::uint64_t(ranks),"audit identity failed");
        ensure(header[4]==2 && header[5]==2 && header[6]==4,"audit counts failed");
        int metadata[11]; binary.read(reinterpret_cast<char*>(metadata),sizeof(metadata));
        ensure(metadata[0]==2*ranks-1-rank && metadata[1]==rank && metadata[2]==1 && metadata[3]==0,"audit node mapping failed");
        double data[74]; binary.read(reinterpret_cast<char*>(data),sizeof(data));
        ensure(bool(binary) && binary.peek()==std::char_traits<char>::eof(),"audit binary size failed");
        for (int i=0;i<74;++i) ensure(data[i]==(i<36?2.:(i<60?3.:(i<64?4.:(i<68?5.:3.)))),"audit stage content failed");
        rejected=false;
        try { SimpleFirstStepAudit duplicate(MPI_COMM_WORLD,prefix,2,2,4,raw(source),raw(owned),raw(offsets),raw(columns)); }
        catch (const std::runtime_error&) { rejected=true; }
        ensure(rejected,"audit overwrite was not rejected collectively");
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
