#pragma once
#include "mars_segregated_simple_output.hpp"
#include <cstdint>

namespace mars::segregated::runtime {
// Private, opt-in file I/O. No copies or extra synchronization occur in ordinary runs.
class SimpleFirstStepAudit {
    MPI_Comm comm_;
    std::ofstream file_;
    int nodes_,blocks_,stage_=0;
    template<class T> void write(const T* values,std::size_t count) {
        static_assert(sizeof(int)==4 && sizeof(double)==8);
        std::vector<T> host(count);
        if (count) {
#ifdef MARS_REPLAY_CUDA
            assembly_cuda_check(cudaMemcpy(host.data(),values,count*sizeof(T),cudaMemcpyDeviceToHost));
#else
            std::copy_n(values,count,host.data());
#endif
            file_.write(reinterpret_cast<const char*>(host.data()),std::streamsize(count*sizeof(T)));
        }
    }
    void checked() { file_.flush(); simple_collective(comm_,bool(file_),"private first-step audit output failed"); }
public:
    SimpleFirstStepAudit(MPI_Comm comm,const std::string& prefix,int nodes,int owned,int blocks,
                        const int* source,const int* owned_nodes,const int* offsets,const int* columns)
        :comm_(comm),nodes_(nodes),blocks_(blocks) {
        int rank,ranks; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
        std::ostringstream name; name<<prefix<<"-audit-rank"<<std::setw(6)<<std::setfill('0')<<rank<<".bin";
        simple_collective(comm,!std::filesystem::exists(name.str()),"private audit exists; choose a fresh prefix");
        file_.open(name.str(),std::ios::binary);
        // Fixed-width native-endian header; the reader checks its byte-order marker.
        const std::uint64_t header[]={0x4d53415544495431ULL,1,std::uint64_t(rank),std::uint64_t(ranks),
                                     std::uint64_t(nodes),std::uint64_t(owned),std::uint64_t(blocks)};
        file_.write(reinterpret_cast<const char*>(header),sizeof(header));
        write(source,nodes); write(owned_nodes,owned); write(offsets,std::size_t(nodes)+1); write(columns,blocks);
        checked();
    }
    void momentum(BlockCsrView<3> a,const double* increment,const double* predictor,const double* influence) {
        ensure(stage_++==0,"invalid first-step audit sequence");
        write(a.values,9*std::size_t(blocks_)); write(a.rhs,3*std::size_t(nodes_));
        write(increment,3*std::size_t(nodes_)); write(predictor,3*std::size_t(nodes_)); write(influence,3*std::size_t(nodes_));
        checked();
    }
    void pressure(BlockCsrView<1> a,const double* increment) {
        ensure(stage_++==1,"invalid first-step audit sequence");
        write(a.values,blocks_); write(a.rhs,nodes_); write(increment,nodes_); checked();
    }
    void finish(const double* gradient) {
        ensure(stage_++==2,"invalid first-step audit sequence");
        write(gradient,3*std::size_t(nodes_)); file_.close();
        simple_collective(comm_,bool(file_),"private first-step audit output failed");
    }
};
}
