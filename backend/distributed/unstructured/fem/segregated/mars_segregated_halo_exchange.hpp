#pragma once
// Fused node-field halo exchange from explicit per-peer lists, in NodeHaloTopology's layout
// (peers, CSR send/recv offsets, local node ids). Owners' values of every listed field go to
// each ghost copy in one message per peer: node-major, fields concatenated per node.
// Persistent buffers; device pointers go straight to (CUDA-aware) MPI in a CUDA build.
// The node lists themselves come from the halo/ownership work; nothing here decides ownership.
#include "mars_segregated_distributed_matrix.hpp"
#include <initializer_list>
#include <string>
#include <vector>
#if defined(__CUDACC__)
#define MARS_HALO_HD __host__ __device__
#else
#define MARS_HALO_HD
#endif

namespace mars::segregated::distributed {
struct Field { double* values; int components; };
constexpr int max_fields=4;
struct FieldSet { double* values[max_fields]; int components[max_fields]; int count, stride; };

namespace kernels {
struct PackFields {
    FieldSet f; const int* nodes; double* buffer;
    MARS_HALO_HD void operator()(int i) const {
        const int node=nodes[i]; int out=f.stride*i;
        for (int k=0;k<f.count;++k) for (int c=0;c<f.components[k];++c) buffer[out++]=f.values[k][f.components[k]*node+c];
    }
};
struct UnpackFields {
    FieldSet f; const int* nodes; const double* buffer;
    MARS_HALO_HD void operator()(int i) const {
        const int node=nodes[i]; int in=f.stride*i;
        for (int k=0;k<f.count;++k) for (int c=0;c<f.components[k];++c) f.values[k][f.components[k]*node+c]=buffer[in++];
    }
};
} // namespace kernels

class FieldExchange {
public:
    // Lists are host copies (ElementDomain: getNodeHaloTopology().peers_, sendOffsets_, recvOffsets_,
    // and downloads of sendNodeIds_/recvNodeIds_). Construction checks, collectively, that every
    // peer sends exactly as many nodes as the receiver expects and that ids address local nodes.
    FieldExchange(MPI_Comm comm,const std::vector<int>& peers,const std::vector<int>& send_offsets,
        const std::vector<int>& send_nodes,const std::vector<int>& recv_offsets,const std::vector<int>& recv_nodes,
        int nodes,int max_stride=8,Stream stream={})
        : comm_(comm), stream_(stream), peers_(peers), send_offsets_(send_offsets), recv_offsets_(recv_offsets),
          max_stride_(max_stride)
    {
        int local=0;
        const std::size_t p=peers.size();
        if (send_offsets.size()!=p+1 || recv_offsets.size()!=p+1 || send_offsets.front()!=0 || recv_offsets.front()!=0
            || std::size_t(send_offsets.back())!=send_nodes.size() || std::size_t(recv_offsets.back())!=recv_nodes.size())
            local|=capacity;
        for (int v:send_nodes) if (v<0 || v>=nodes) local|=capacity;
        for (int v:recv_nodes) if (v<0 || v>=nodes) local|=capacity;
        long long buffer=0;
        if (!checked_product(max_stride,(long long)std::max(send_nodes.size(),recv_nodes.size()),std::numeric_limits<int>::max(),buffer))
            local|=overflow;
        // Counts, not a peer handshake: an asymmetric peer list must fail, not hang. q's
        // receive count from me must equal my send count to q (build time, one int per rank).
        int ranks=1, rank=0; MPI_Comm_size(comm_,&ranks); MPI_Comm_rank(comm_,&rank);
        std::vector<int> receive(ranks,0), expected_by(ranks,0);
        if (!(local&capacity)) for (std::size_t i=0;i<p;++i) {
            if (peers[i]<0 || peers[i]>=ranks || peers[i]==rank) { local|=capacity; continue; }
            receive[peers[i]]+=recv_offsets[i+1]-recv_offsets[i];
        }
        MPI_Alltoall(receive.data(),1,MPI_INT,expected_by.data(),1,MPI_INT,comm_);
        if (!(local&capacity)) {
            std::vector<int> send(ranks,0);
            for (std::size_t i=0;i<p;++i) send[peers[i]]+=send_offsets[i+1]-send_offsets[i];
            for (int q=0;q<ranks;++q) if (send[q]!=expected_by[q]) local|=capacity;
        }
        int global=0; MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,comm_);
        if (global) throw std::runtime_error("halo exchange lists rejected on all ranks (global: "+describe(global)+"; this rank: "+describe(local)+")");
        copy_in(send_nodes_,send_nodes.data(),send_nodes.size()); copy_in(recv_nodes_,recv_nodes.data(),recv_nodes.size());
        send_buffer_.resize(std::size_t(max_stride)*send_nodes.size());
        recv_buffer_.resize(std::size_t(max_stride)*recv_nodes.size());
    }
    FieldExchange(const FieldExchange&)=delete;
    FieldExchange& operator=(const FieldExchange&)=delete;

    // Refresh the ghost entries of every field from its owner, in one round. Owned entries
    // must be final; ghost entries are overwritten. Field arrays hold components*nodes values.
    void operator()(std::initializer_list<Field> fields) {
        FieldSet set{}; set.count=int(fields.size()); set.stride=0;
        if (set.count<1 || set.count>max_fields) throw std::runtime_error("halo exchange: 1..4 fields per round");
        int k=0;
        for (const Field& f:fields) { set.values[k]=f.values; set.components[k]=f.components; set.stride+=f.components; ++k; }
        if (set.stride>max_stride_) throw std::runtime_error("halo exchange: fields exceed the buffer stride");
        int faults=0;
        apply(int(send_nodes_.size()),kernels::PackFields{set,raw(send_nodes_),raw(send_buffer_)},stream_,faults);
        sync(faults);
        std::vector<MPI_Request> requests; requests.reserve(2*peers_.size());
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=set.stride*(recv_offsets_[i+1]-recv_offsets_[i]);
            if (count) { requests.emplace_back(); MPI_Irecv(raw(recv_buffer_)+std::size_t(set.stride)*recv_offsets_[i],count,MPI_DOUBLE,peers_[i],tag_values,comm_,&requests.back()); }
        }
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=set.stride*(send_offsets_[i+1]-send_offsets_[i]);
            if (count) { requests.emplace_back(); MPI_Isend(raw(send_buffer_)+std::size_t(set.stride)*send_offsets_[i],count,MPI_DOUBLE,peers_[i],tag_values,comm_,&requests.back()); }
        }
        MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE);
        apply(int(recv_nodes_.size()),kernels::UnpackFields{set,raw(recv_nodes_),raw(recv_buffer_)},stream_,faults);
        sync(faults);
        if (faults) throw std::runtime_error("halo exchange: CUDA error");
        ++rounds_; values_+=(long long)set.stride*(long long)recv_nodes_.size();
    }
    std::size_t ghosts() const { return recv_nodes_.size(); }
    long long rounds() const { return rounds_; }
    long long received_values() const { return values_; }
private:
    static constexpr int tag_values=0x4d48;
    void sync(int& faults) const {
#if defined(__CUDACC__)
        if (!cuda_ok(cudaStreamSynchronize(stream_))) faults|=device_error;
#else
        (void)faults;
#endif
    }
    MPI_Comm comm_; Stream stream_;
    std::vector<int> peers_, send_offsets_, recv_offsets_;
    int max_stride_;
    Buffer<int> send_nodes_, recv_nodes_;
    Buffer<double> send_buffer_, recv_buffer_;
    long long rounds_=0, values_=0;
};
} // namespace mars::segregated::distributed
#undef MARS_HALO_HD
