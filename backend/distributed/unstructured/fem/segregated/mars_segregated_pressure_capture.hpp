#pragma once
#include <algorithm>
#include <cmath>
#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace mars::segregated::frozen {
inline void require(bool pass) { if (!pass) throw std::runtime_error("invalid private pressure capture"); }
inline std::string part_name(int rank,const char* suffix=".bin") {
    std::ostringstream out; out<<"rank-"<<std::setw(6)<<std::setfill('0')<<rank<<suffix; return out.str();
}
// Only opt-in diagnostic output uses this host staging. Solver storage stays on device.
template<class T> std::vector<T> copy_input(const T* input,std::size_t size) {
    std::vector<T> host(size);
    if (size) {
#ifdef MARS_REPLAY_CUDA
        require(cudaMemcpy(host.data(),input,size*sizeof(T),cudaMemcpyDeviceToHost)==cudaSuccess);
#else
        std::copy_n(input,size,host.data());
#endif
    }
    return host;
}
class Writer {
    int fd_=-1;
public:
    explicit Writer(const std::filesystem::path& path):fd_(::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600)) { require(fd_>=0); }
    ~Writer() { if (fd_>=0) ::close(fd_); }
    Writer(const Writer&)=delete;
    void bytes(const void* data,std::size_t count) {
        auto p=static_cast<const char*>(data);
        while (count) {
            const auto n=::write(fd_,p,count);
            if (n<0 && errno==EINTR) continue;
            require(n>0); p+=n; count-=std::size_t(n);
        }
    }
    template<class T> void array(const std::vector<T>& data) { bytes(data.data(),data.size()*sizeof(T)); }
    void finish() { require(::fsync(fd_)==0); const int old=fd_; fd_=-1; require(::close(old)==0); }
};
struct Part {
    static constexpr std::uint64_t magic=0x4d53505245535331ULL;
    std::uint64_t rank=0,ranks=1,first=0,last=0,total=0,nodes=0,nnz=0,iteration=0;
    bool solver_passed=false,mars_passed=false,maximum=true;
    double absolute=0,relative=0;
    std::vector<std::int32_t> offsets,columns;
    std::vector<std::int64_t> map;
    std::vector<double> values,rhs,candidate;
    std::size_t rows() const { return std::size_t(last-first); }
    void validate() const {
        require(ranks>0 && rank<ranks && first<last && last<=total && total<=INT32_MAX
            && nodes<=INT32_MAX && nnz<=INT32_MAX && iteration>0 && nnz>=rows());
        require(offsets.size()==rows()+1 && columns.size()==nnz && values.size()==nnz
            && map.size()==nodes && rhs.size()==rows() && candidate.size()==nodes);
        require(offsets.front()==0 && offsets.back()==std::int64_t(nnz));
        require(std::isfinite(absolute) && std::isfinite(relative) && absolute>=0 && relative>=0 && absolute+relative>0);
        for (std::size_t row=0;row<rows();++row) require(offsets[row]>=0 && offsets[row]<offsets[row+1]);
        for (auto id:map) require(id>=-1 && id<std::int64_t(total));
        for (std::size_t k=0;k<nnz;++k) require(columns[k]>=0 && columns[k]<std::int64_t(nodes)
            && map[columns[k]]>=0 && std::isfinite(values[k]));
        for (double value:rhs) require(std::isfinite(value));
        // Unreferenced ghosts may be poisoned. Owners and referenced columns are checked by the reader.
    }
    void write(const std::filesystem::path& path) const {
        validate(); Writer out(path);
        const std::vector<std::uint64_t> header={magic,1,rank,ranks,first,last,total,nodes,nnz,iteration,
            std::uint64_t(solver_passed),std::uint64_t(mars_passed),std::uint64_t(maximum)};
        out.array(header); out.array(std::vector<double>{absolute,relative});
        out.array(offsets); out.array(columns); out.array(map); out.array(values); out.array(rhs); out.array(candidate); out.finish();
    }
    static Part read(const std::filesystem::path& path) {
        std::ifstream in(path,std::ios::binary); require(bool(in));
        auto read=[&](auto& v) { in.read(reinterpret_cast<char*>(v.data()),std::streamsize(v.size()*sizeof(v[0]))); require(bool(in)); };
        std::vector<std::uint64_t> h(13); read(h);
        require(h[0]==magic && h[1]==1 && h[10]<=1 && h[11]<=1 && h[12]<=1);
        Part p; p.rank=h[2]; p.ranks=h[3]; p.first=h[4]; p.last=h[5]; p.total=h[6]; p.nodes=h[7]; p.nnz=h[8]; p.iteration=h[9];
        p.solver_passed=h[10]; p.mars_passed=h[11]; p.maximum=h[12];
        require(p.first<p.last && p.last<=p.total && p.total<=INT32_MAX && p.nodes<=INT32_MAX && p.nnz<=INT32_MAX);
        const std::uint64_t bytes=13*8+16+4*(p.rows()+1)+12*p.nnz+16*p.nodes+8*p.rows();
        require(std::filesystem::file_size(path)==bytes);
        std::vector<double> tol(2); read(tol); p.absolute=tol[0]; p.relative=tol[1];
        p.offsets.resize(p.rows()+1); p.columns.resize(p.nnz); p.map.resize(p.nodes);
        p.values.resize(p.nnz); p.rhs.resize(p.rows()); p.candidate.resize(p.nodes);
        read(p.offsets); read(p.columns); read(p.map); read(p.values); read(p.rhs); read(p.candidate); p.validate(); return p;
    }
};
template<class System,class Candidate,class Id> Part capture_part(const System& system,const Id* ids_input,const double* rhs,Candidate candidate,
    int rank,int ranks,int iteration,bool solved,bool passed,double absolute,double relative,bool maximum) {
    Part p; const auto range=system.hypre_rows(); p.rank=rank; p.ranks=ranks;
    p.first=range.begin; p.last=range.end; p.total=range.column_end; p.nodes=system.nodes(); p.nnz=system.nnz(); p.iteration=iteration;
    p.solver_passed=solved; p.mars_passed=passed; p.absolute=absolute; p.relative=relative; p.maximum=maximum;
    const auto& a=system.matrix();
    auto offsets=copy_input(a.rowOffsetsPtr(),p.rows()+1),columns=copy_input(a.colIndicesPtr(),p.nnz);
    p.offsets.assign(offsets.begin(),offsets.end()); p.columns.assign(columns.begin(),columns.end());
    const auto ids=copy_input(ids_input,p.nodes); p.map.assign(ids.begin(),ids.end());
    p.values=copy_input(a.valuesPtr(),p.nnz); p.rhs=copy_input(rhs,p.rows());
    require(candidate.size>=p.nodes); p.candidate=copy_input(candidate.values,p.nodes); return p;
}
}
