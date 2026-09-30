#pragma once
#include "mars_segregated_simple_mesh.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#if defined(__CUDACC__)
#define MARS_SOUTPUT_HD __host__ __device__
#else
#define MARS_SOUTPUT_HD
#endif
namespace mars::segregated::runtime {
struct FieldRow { double values[8]; };
struct FieldRowLess {
    MARS_SOUTPUT_HD bool operator()(const FieldRow& a,const FieldRow& b) const { return a.values[0]<b.values[0]; }
};
struct PackOutput {
    const int *owned,*source; const double *x,*y,*z,*u,*p; FieldRow* rows;
    MARS_SOUTPUT_HD void operator()(int i) const {
        const int n=owned[i]; rows[i]={{double(source[n]),x[n],y[n],z[n],u[3*n],u[3*n+1],u[3*n+2],p[n]}};
    }
};
struct CheckOutput {
    const FieldRow* rows; int global_nodes; bool complete; int* error;
    MARS_SOUTPUT_HD void operator()(int i) const {
        const double node=rows[i].values[0];
        if (!(node>=0 && node<global_nodes) || node!=double(int(node)) ||
            (i && !(rows[i-1].values[0]<node)) || (complete && node!=double(i))) distributed::raise_fault(error,1);
        for (double value:rows[i].values) if (!geometry_finite(value)) distributed::raise_fault(error,1);
    }
};
inline std::string simple_part_path(const std::string& prefix,int rank) {
    std::ostringstream name; name<<prefix<<"-fields-rank"<<std::setw(6)<<std::setfill('0')<<rank<<".csv";
    return name.str();
}
inline void simple_output_preflight(MPI_Comm comm,const std::string& prefix,const std::string& mode) {
    int rank; MPI_Comm_rank(comm,&rank);
    bool free=!std::filesystem::exists(simple_part_path(prefix,rank));
    if (!rank) for (const char* suffix:{"-metrics.csv","-fields.csv","-fields.json"})
        free=free && !std::filesystem::exists(prefix+suffix);
    simple_collective(comm,free,"output exists; choose a fresh prefix");
    simple_collective(comm,mode=="gathered" || mode=="distributed" || mode=="none","invalid field output mode");
}
inline std::string simple_json_string(const std::string& value) {
    std::ostringstream out; out<<'"';
    for (unsigned char c:value) {
        if (c=='"' || c=='\\') out<<'\\'<<char(c);
        else if (c<32) out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
        else out<<char(c);
    }
    out<<'"'; return out.str();
}
inline bool write_simple_csv(const std::string& path,const Buffer<FieldRow>& rows) {
    std::vector<FieldRow> fields(rows.size());
#ifdef MARS_REPLAY_CUDA
    thrust::copy(rows.begin(),rows.end(),fields.begin());
#else
    fields=rows;
#endif
    std::ofstream out(path); out<<std::setprecision(17)<<"node,x,y,z,u,v,w,p\n";
    for (const auto& row:fields) { for (int j=0;j<8;++j) out<<(j?",":"")<<row.values[j]; out<<'\n'; }
    out.close(); return bool(out);
}
inline void write_simple_fields(MPI_Comm comm,const std::string& prefix,const std::string& mode,int nodes,Buffer<FieldRow>& rows) {
    if (mode=="none") return;
    int rank,ranks; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    static_assert(sizeof(FieldRow)==8*sizeof(double));
    mesh_sort(rows,FieldRowLess{});
    Buffer<int> error(1,0);
    launch(int(rows.size()),CheckOutput{raw(rows),nodes,false,raw(error)});
    mesh_check(comm,error,"invalid output node identity or nonfinite field");
    long long count=static_cast<long long>(rows.size()),total=0;
    ensure(MPI_Allreduce(&count,&total,1,MPI_LONG_LONG,MPI_SUM,comm)==MPI_SUCCESS,"field count reduction failed");
    simple_collective(comm,total==nodes,"field output does not cover every source node");
    if (mode=="distributed") {
        simple_collective(comm,write_simple_csv(simple_part_path(prefix,rank),rows),"field output failed");
        bool ok=true;
        if (!rank) {
            std::ofstream out(prefix+"-fields.json");
            out<<"{\n  \"format\": \"mars-simple-fields-v1\",\n  \"nodes\": "<<nodes<<",\n  \"parts\": [";
            for (int q=0;q<ranks;++q) out<<(q?", ":"")<<simple_json_string(std::filesystem::path(simple_part_path(prefix,q)).filename().string());
            out<<"]\n}\n"; out.close(); ok=bool(out);
        }
        simple_collective(comm,ok,"field manifest output failed");
        return;
    }
    simple_collective(comm,count<=INT_MAX/8,"field output exceeds MPI count capacity; use distributed output");
    const int values=int(count)*8; std::vector<int> counts(ranks),offsets(ranks);
    ensure(MPI_Gather(&values,1,MPI_INT,counts.data(),1,MPI_INT,0,comm)==MPI_SUCCESS,"field counts failed");
    long long all=0;
    if (!rank) for (int q=0;q<ranks;++q) { offsets[q]=int(std::min(all,static_cast<long long>(INT_MAX))); all+=counts[q]; }
    simple_collective(comm,rank!=0 || all<=INT_MAX,"field gather exceeds MPI count capacity; use distributed output");
    Buffer<FieldRow> gathered(rank?0:std::size_t(nodes));
#ifdef MARS_REPLAY_CUDA
    assembly_cuda_check(cudaStreamSynchronize(nullptr));
#endif
    ensure(MPI_Gatherv(raw(rows),values,MPI_DOUBLE,raw(gathered),counts.data(),offsets.data(),MPI_DOUBLE,0,comm)==MPI_SUCCESS,"device field gather failed");
    if (!rank) {
        mesh_sort(gathered,FieldRowLess{});
        launch(nodes,CheckOutput{raw(gathered),nodes,true,raw(error)});
    }
    mesh_check(comm,error,"field gather returned duplicate or missing source nodes");
    simple_collective(comm,rank!=0 || write_simple_csv(prefix+"-fields.csv",gathered),"field output failed");
}
}

#undef MARS_SOUTPUT_HD
