#pragma once
// Fused node-field halo exchange from explicit per-peer lists, in NodeHaloTopology's layout
// (peers, CSR send/recv offsets, local node ids). Owners' values of every listed field go to
// each ghost copy in one message per peer: node-major, fields concatenated per node.
// Persistent buffers; device pointers go straight to (CUDA-aware) MPI in a CUDA build.
// The node lists themselves come from the halo/ownership work; nothing here decides ownership.
#include "mars_segregated_distributed_matrix.hpp"
#include <initializer_list>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <type_traits>
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
// Reverse direction: a ghost's partial value is added into its owner's entry.
struct AddFields {
    FieldSet f; const int* nodes; const double* buffer;
    MARS_HALO_HD void operator()(int i) const {
        const int node=nodes[i]; int in=f.stride*i;
        for (int k=0;k<f.count;++k) for (int c=0;c<f.components[k];++c) assembly_add(f.values[k]+f.components[k]*node+c,buffer[in++]);
    }
};
template<class T> struct PackInteger {
    const T* values; const int* nodes; T* buffer;
    MARS_HALO_HD void operator()(int i) const { buffer[i]=values[nodes[i]]; }
};
template<class T> struct UnpackInteger {
    T* values; const int* nodes; const T* buffer;
    MARS_HALO_HD void operator()(int i) const { values[nodes[i]]=buffer[i]; }
};
struct CheckHaloNode {
    const int* nodes; int count; int* seen; int* status;
    MARS_HALO_HD void operator()(int i) const {
        const int v=nodes[i];
        if (v<0 || v>=count) { raise_fault(status,capacity); return; }
        if (seen) {
#if defined(__CUDA_ARCH__)
            if (atomicExch(seen+v,1)) raise_fault(status,capacity);
#else
            if (seen[v]) raise_fault(status,capacity);
            seen[v]=1;
#endif
        }
    }
};
} // namespace kernels

class FieldExchange {
public:
    // CPU fixtures upload their lists once. Native callers pass device lists directly.
    FieldExchange(MPI_Comm comm,const std::vector<int>& peers,const std::vector<int>& send_offsets,
        const std::vector<int>& send_nodes,const std::vector<int>& recv_offsets,const std::vector<int>& recv_nodes,
        int nodes,int max_stride=8,Stream stream={})
        :comm_(comm),stream_(stream),peers_(peers),send_offsets_(send_offsets),recv_offsets_(recv_offsets),
         max_stride_(max_stride),nodes_(nodes),send_nodes_(send_nodes.begin(),send_nodes.end()),recv_nodes_(recv_nodes.begin(),recv_nodes.end())
    { initialize(); }

    struct NodeLists { const int* send; std::size_t send_size; const int* recv; std::size_t recv_size; };
    FieldExchange(MPI_Comm comm,const std::vector<int>& peers,const std::vector<int>& send_offsets,
        const std::vector<int>& recv_offsets,NodeLists lists,int nodes,int max_stride=8,Stream stream={})
        :comm_(comm),stream_(stream),peers_(peers),send_offsets_(send_offsets),recv_offsets_(recv_offsets),
         max_stride_(max_stride),nodes_(nodes)
    {
        reject_lists((lists.send_size && !lists.send) || (lists.recv_size && !lists.recv)?capacity:0);
        if (lists.send_size) copy_in(send_nodes_,lists.send,lists.send_size);
        if (lists.recv_size) copy_in(recv_nodes_,lists.recv,lists.recv_size);
        initialize();
    }
#if defined(__CUDACC__)
    FieldExchange(MPI_Comm comm,const std::vector<int>& peers,const std::vector<int>& send_offsets,
        const Buffer<int>& send_nodes,const std::vector<int>& recv_offsets,const Buffer<int>& recv_nodes,
        int nodes,int max_stride=8,Stream stream={})
        :FieldExchange(comm,peers,send_offsets,recv_offsets,
                       NodeLists{raw(send_nodes),send_nodes.size(),raw(recv_nodes),recv_nodes.size()},nodes,max_stride,stream) {}
#endif
    FieldExchange(const FieldExchange&)=delete;
    FieldExchange& operator=(const FieldExchange&)=delete;

    // Refresh the ghost entries of every field from its owner, in one round. Owned entries
    // must be final; ghost entries are overwritten. Field arrays hold components*nodes values.
    void operator()(std::initializer_list<Field> fields) {
        const FieldSet set=field_set(fields);
        int faults=0;
        apply(int(send_nodes_.size()),kernels::PackFields{set,raw(send_nodes_),raw(send_buffer_)},stream_,faults);
        finish_device(faults);
        auto& requests=requests_; requests.clear();
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=set.stride*(recv_offsets_[i+1]-recv_offsets_[i]);
            if (count) { requests.emplace_back(); check_mpi(MPI_Irecv(raw(recv_buffer_)+std::size_t(set.stride)*recv_offsets_[i],count,MPI_DOUBLE,peers_[i],tag_values,comm_,&requests.back())); }
        }
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=set.stride*(send_offsets_[i+1]-send_offsets_[i]);
            if (count) { requests.emplace_back(); check_mpi(MPI_Isend(raw(send_buffer_)+std::size_t(set.stride)*send_offsets_[i],count,MPI_DOUBLE,peers_[i],tag_values,comm_,&requests.back())); }
        }
        check_mpi(MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE));
        apply(int(recv_nodes_.size()),kernels::UnpackFields{set,raw(recv_nodes_),raw(recv_buffer_)},stream_,faults);
        finish_device(faults);
        ++rounds_; values_+=(long long)set.stride*(long long)recv_nodes_.size();
    }
    // Transpose of the publish: ghost entries are added into their owners' entries (atomically:
    // several peers may hold the same node). Ghost entries are left as they were. Used for setup
    // checks such as star completeness, not in the SIMPLE iteration.
    void reverse_add(Field field) {
        const FieldSet set=field_set({field});
        int faults=0;
        apply(int(recv_nodes_.size()),kernels::PackFields{set,raw(recv_nodes_),raw(recv_buffer_)},stream_,faults);
        finish_device(faults);
        auto& requests=requests_; requests.clear();
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=set.stride*(send_offsets_[i+1]-send_offsets_[i]);
            if (count) { requests.emplace_back(); check_mpi(MPI_Irecv(raw(send_buffer_)+std::size_t(set.stride)*send_offsets_[i],count,MPI_DOUBLE,peers_[i],tag_reverse,comm_,&requests.back())); }
        }
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=set.stride*(recv_offsets_[i+1]-recv_offsets_[i]);
            if (count) { requests.emplace_back(); check_mpi(MPI_Isend(raw(recv_buffer_)+std::size_t(set.stride)*recv_offsets_[i],count,MPI_DOUBLE,peers_[i],tag_reverse,comm_,&requests.back())); }
        }
        check_mpi(MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE));
        apply(int(send_nodes_.size()),kernels::AddFields{set,raw(send_nodes_),raw(send_buffer_)},stream_,faults);
        finish_device(faults);
    }
    // Integer identities travel as integers directly from device memory, never through double.
    template<class T> void publish(T* values,std::size_t size) {
        static_assert(std::is_integral_v<T> && (sizeof(T)==4 || sizeof(T)==8));
        if (size!=std::size_t(nodes_) || (size && !values)) fatal("metadata does not cover local nodes");
        MPI_Datatype type;
        if constexpr (sizeof(T)==8) {
            if constexpr (std::is_signed<T>::value) type=MPI_INT64_T;
            else type=MPI_UINT64_T;
        } else {
            if constexpr (std::is_signed<T>::value) type=MPI_INT32_T;
            else type=MPI_UINT32_T;
        }
        Buffer<T> send(send_nodes_.size()),recv(recv_nodes_.size());
        int faults=0;
        apply(int(send.size()),kernels::PackInteger<T>{values,raw(send_nodes_),raw(send)},stream_,faults);
        finish_device(faults);
        auto& requests=requests_; requests.clear();
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=recv_offsets_[i+1]-recv_offsets_[i];
            if (count) { requests.emplace_back(); check_mpi(MPI_Irecv(raw(recv)+recv_offsets_[i],count,type,peers_[i],tag_metadata,comm_,&requests.back())); }
        }
        for (std::size_t i=0;i<peers_.size();++i) {
            const int count=send_offsets_[i+1]-send_offsets_[i];
            if (count) { requests.emplace_back(); check_mpi(MPI_Isend(raw(send)+send_offsets_[i],count,type,peers_[i],tag_metadata,comm_,&requests.back())); }
        }
        check_mpi(MPI_Waitall(int(requests.size()),requests.data(),MPI_STATUSES_IGNORE));
        apply(int(recv.size()),kernels::UnpackInteger<T>{values,raw(recv_nodes_),raw(recv)},stream_,faults);
        finish_device(faults);
    }
    template<class T> void publish_host(std::vector<T>& values) {
        Buffer<T> device(values.begin(),values.end());
        publish(raw(device),device.size());
#if defined(__CUDACC__)
        thrust::copy(device.begin(),device.end(),values.begin());
#else
        values=device;
#endif
    }
    std::size_t ghosts() const { return recv_nodes_.size(); }
    long long rounds() const { return rounds_; }
    long long received_values() const { return values_; }
private:
    void initialize() {
        int local=0;
        const std::size_t p=peers_.size();
        if (send_offsets_.size()!=p+1 || recv_offsets_.size()!=p+1 || send_offsets_.front()!=0 || recv_offsets_.front()!=0
            || std::size_t(send_offsets_.back())!=send_nodes_.size() || std::size_t(recv_offsets_.back())!=recv_nodes_.size())
            local|=capacity;
        if (nodes_<0 || max_stride_<1 || !std::is_sorted(send_offsets_.begin(),send_offsets_.end())
            || !std::is_sorted(recv_offsets_.begin(),recv_offsets_.end())) local|=capacity;
        long long buffer=0;
        if (!checked_product(max_stride_,(long long)std::max(send_nodes_.size(),recv_nodes_.size()),std::numeric_limits<int>::max(),buffer)
            || !checked_product(max_stride_,nodes_,std::numeric_limits<int>::max(),buffer)) local|=overflow;
        int ranks=1,rank=0; MPI_Comm_size(comm_,&ranks); MPI_Comm_rank(comm_,&rank);
        for (int peer:peers_) if (peer<0 || peer>=ranks || peer==rank) local|=capacity;
        if (std::set<int>(peers_.begin(),peers_.end()).size()!=p || p>std::size_t(std::numeric_limits<int>::max()/2)) local|=capacity;
        reject_lists(local);
        Buffer<int> status(1,0),seen(std::size_t(nodes_),0);
        apply(int(send_nodes_.size()),kernels::CheckHaloNode{raw(send_nodes_),nodes_,nullptr,raw(status)},stream_,local);
        apply(int(recv_nodes_.size()),kernels::CheckHaloNode{raw(recv_nodes_),nodes_,raw(seen),raw(status)},stream_,local);
        local|=fetch(raw(status),stream_,local);
        reject_lists(local);
        validate_peer_counts();
        try {
            send_buffer_.resize(std::size_t(max_stride_)*send_nodes_.size());
            recv_buffer_.resize(std::size_t(max_stride_)*recv_nodes_.size());
            requests_.reserve(2*p);
        } catch (const std::exception&) { fatal("cannot allocate halo buffers"); }
    }
    static constexpr int tag_counts=0x4d47, tag_values=0x4d48, tag_reverse=0x4d49, tag_metadata=0x4d4a;
    [[noreturn]] void fatal(const char* message) const {
        std::fprintf(stderr,"ERROR: node-field halo: %s\n",message);
        MPI_Abort(comm_,1);
        std::abort();
    }
    void check_mpi(int code) const { if (code!=MPI_SUCCESS) fatal("MPI exchange failed"); }
    void reject_lists(int local) const {
        int global=0; check_mpi(MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,comm_));
        if (global) throw std::runtime_error("halo exchange lists rejected on all ranks (global: "+describe(global)+"; this rank: "+describe(local)+")");
    }
    // Sparse NBX handshake: synchronous sends complete only when received. Keep probing
    // until all ranks finish their sends, including ranks with no peers or one-sided lists.
    void validate_peer_counts() {
        std::vector<std::array<int,2>> send(peers_.size());
        std::vector<MPI_Request> pending(peers_.size());
        std::map<int,std::array<int,2>> received;
        for (std::size_t i=0;i<peers_.size();++i) {
            send[i]={send_offsets_[i+1]-send_offsets_[i],recv_offsets_[i+1]-recv_offsets_[i]};
            check_mpi(MPI_Issend(send[i].data(),2,MPI_INT,peers_[i],tag_counts,comm_,&pending[i]));
        }
        MPI_Request barrier=MPI_REQUEST_NULL;
        bool started=false; int done=0, local=0;
        while (!done) {
            int ready=0; MPI_Status status;
            check_mpi(MPI_Iprobe(MPI_ANY_SOURCE,tag_counts,comm_,&ready,&status));
            if (ready) {
                std::array<int,2> counts;
                check_mpi(MPI_Recv(counts.data(),2,MPI_INT,status.MPI_SOURCE,tag_counts,comm_,MPI_STATUS_IGNORE));
                if (!received.emplace(status.MPI_SOURCE,counts).second) local|=capacity;
            }
            if (started) check_mpi(MPI_Test(&barrier,&done,MPI_STATUS_IGNORE));
            else {
                int sent=0; check_mpi(MPI_Testall(int(pending.size()),pending.data(),&sent,MPI_STATUSES_IGNORE));
                if (sent) { check_mpi(MPI_Ibarrier(comm_,&barrier)); started=true; }
            }
        }
        if (received.size()!=peers_.size()) local|=capacity;
        for (std::size_t i=0;i<peers_.size();++i) {
            const auto it=received.find(peers_[i]);
            if (it==received.end() || it->second[0]!=send[i][1] || it->second[1]!=send[i][0]) local|=capacity;
        }
        reject_lists(local);
    }
    FieldSet field_set(std::initializer_list<Field> fields) const {
        FieldSet set{};
        if (fields.size()<1 || fields.size()>max_fields) fatal("1..4 fields required per round");
        for (const Field& f:fields) {
            if (f.components<1 || f.components>max_stride_-set.stride || (nodes_ && !f.values))
                fatal("invalid field pointer or buffer stride");
            set.values[set.count]=f.values; set.components[set.count++]=f.components; set.stride+=f.components;
        }
        return set;
    }
    void finish_device(int faults) const {
#if defined(__CUDACC__)
        if (!cuda_ok(cudaStreamSynchronize(stream_))) faults|=device_error;
#endif
        if (faults) fatal("CUDA pack or unpack failed");
    }
    MPI_Comm comm_; Stream stream_;
    std::vector<int> peers_, send_offsets_, recv_offsets_;
    int max_stride_, nodes_;
    std::vector<MPI_Request> requests_;
    Buffer<int> send_nodes_, recv_nodes_;
    Buffer<double> send_buffer_, recv_buffer_;
    long long rounds_=0, values_=0;
};
} // namespace mars::segregated::distributed
#undef MARS_HALO_HD
