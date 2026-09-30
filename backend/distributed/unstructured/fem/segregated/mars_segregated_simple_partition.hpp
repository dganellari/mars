#pragma once
// Rank-local SIMPLE input and ownership from ElementDomain-shaped state: node keys (sorted
// local order), the cstone local element range, node ownership flags and NodeHaloTopology
// lists. Setup only, host-side; the iteration itself stays on the device.
//
// Boundary faces: a face of a held element is exterior when no other held element shares it.
// That is exact for every face touching an owned node, because complete stars hold the
// neighbour across such a face. Faces touching no owned node may be misjudged; they then only
// feed ghost rows and ghost sums, which are never read. A face belongs to the owner of its
// smallest-key node, so each exterior face is counted by exactly one rank.
// Checks (all collective): list sizes, peer key identities, every ghost received from an owner, and
// complete element stars (global incidence via a reverse add must equal held incidence).
#include "mars_segregated_simple_distributed.hpp"
#include "mars_segregated_simple_input.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace mars::segregated::runtime {

template<class Key> struct DomainView {
    std::vector<double> x,y,z;             // per local node
    std::vector<Key> key;                  // node identity, strictly increasing (ElementDomain's SFC map)
    std::array<std::vector<int>,4> nodes;  // local node ids per held element, source node order
    int element_begin=0, element_end=0;    // this rank's own elements; the others are halo
    std::vector<unsigned char> owned;      // 1 where this rank owns the node
    std::vector<int> peers, send_offsets{0}, send_nodes, recv_offsets{0}, recv_nodes;
};
template<class GlobalId> struct SimplePartition {
    SimpleInput input;
    SimpleOwnership<GlobalId> ownership;
    long long incomplete_stars=0, unreceived_ghosts=0;
};

// kind(local[3]) returns 0 inlet, 1 outlet, 2 wall, or -1 for an exterior face SIMPLE ignores.
template<class GlobalId,class Key,class Kind>
SimplePartition<GlobalId> simple_partition(MPI_Comm comm,const DomainView<Key>& v,Kind kind) {
    const int n=int(v.x.size()), e=int(v.nodes[0].size());
    bool ok=v.y.size()==std::size_t(n) && v.z.size()==std::size_t(n) && v.key.size()==std::size_t(n) && v.owned.size()==std::size_t(n)
        && v.element_begin>=0 && v.element_begin<=v.element_end && v.element_end<=e;
    for (int k=1;k<4 && ok;++k) ok=v.nodes[k].size()==std::size_t(e);
    for (int k=0;k<4 && ok;++k) for (int u:v.nodes[k]) ok=ok && u>=0 && u<n;
    for (int i=1;i<n && ok;++i) ok=v.key[i-1]<v.key[i];   // multi-block (duplicate keys) is not supported
    simple_collective(comm,ok,"SIMPLE partition: inconsistent local mesh, or duplicate node keys (multi-block)");

    SimplePartition<GlobalId> part;
    auto& o=part.ownership;
    o.peers=v.peers; o.send_offsets=v.send_offsets; o.send_nodes=v.send_nodes; o.recv_offsets=v.recv_offsets; o.recv_nodes=v.recv_nodes;
    distributed::FieldExchange exchange(comm,v.peers,v.send_offsets,v.send_nodes,v.recv_offsets,v.recv_nodes,n);

    // Every ghost must be refreshed by some peer, or it would keep a stale value forever.
    std::vector<int> received(n,0);
    for (int u:v.recv_nodes) ++received[u];
    for (int u=0;u<n;++u) part.unreceived_ghosts+=!v.owned[u] && received[u]!=1;
    for (int u=0;u<n;++u) ok=ok && !(v.owned[u] && received[u]);
    for (int u:v.send_nodes) ok=ok && v.owned[u];
    simple_collective(comm,ok && part.unreceived_ghosts==0,"SIMPLE partition: halo must send owned nodes and receive each ghost exactly once");

    auto received_keys=v.key;
    exchange.publish_host(received_keys);
    simple_collective(comm,received_keys==v.key,"SIMPLE partition: halo peer node keys do not match");

    // Solver ids: owned nodes in local order, contiguous per rank in rank order.
    for (int u=0;u<n;++u) if (v.owned[u]) o.owned_nodes.push_back(u);
    long long mine=(long long)o.owned_nodes.size(), first=0, total=0; int rank=0;
    MPI_Exscan(&mine,&first,1,MPI_LONG_LONG,MPI_SUM,comm); MPI_Comm_rank(comm,&rank); if (rank==0) first=0;
    MPI_Allreduce(&mine,&total,1,MPI_LONG_LONG,MPI_SUM,comm);
    static_assert(std::is_integral_v<GlobalId>);
    simple_collective(comm,total>0 && (long double)(total-1)<=(long double)std::numeric_limits<GlobalId>::max(),
                      "SIMPLE partition: solver node ids exceed the global index type");
    std::vector<long long> id(n,-1);
    for (std::size_t k=0;k<o.owned_nodes.size();++k) id[o.owned_nodes[k]]=first+(long long)k;
    exchange.publish_host(id);
    for (long long value:id) ok=ok && value>=0 && value<total;
    simple_collective(comm,ok,"SIMPLE partition: a node received no valid solver id");
    o.solver_node.resize(n);
    for (int u=0;u<n;++u) o.solver_node[u]=GlobalId(id[u]);

    // Complete stars: global incidence (own elements, reverse-added to owners) == held incidence.
    std::vector<double> own(n,0.0), held(n,0.0);
    for (int el=0;el<e;++el) for (int k=0;k<4;++k) { held[v.nodes[k][el]]+=1; if (el>=v.element_begin && el<v.element_end) own[v.nodes[k][el]]+=1; }
    Array<double> incidence(own); exchange.reverse_add({incidence.data(),1}); own=incidence.host();
    for (int u=0;u<n;++u) part.incomplete_stars+=v.owned[u] && own[u]!=held[u];
    simple_collective(comm,part.incomplete_stars==0,"SIMPLE partition: incomplete element star at an owned node (halo completion)");

    // Exterior faces of held elements, kind by the caller, owner = owner of the smallest key.
    std::map<std::array<Key,3>,int> count;
    auto face_key=[&](int el,int ordinal) {
        std::array<Key,3> k3; for (int j=0;j<3;++j) k3[j]=v.key[v.nodes[tet_face_node(ordinal,j)][el]];
        std::sort(k3.begin(),k3.end()); return k3;
    };
    for (int el=0;el<e;++el) for (int f=0;f<4;++f) ++count[face_key(el,f)];
    auto& in=part.input;
    in.x=v.x; in.y=v.y; in.z=v.z; in.nodes=v.nodes;
    for (int el=0;el<e;++el) for (int f=0;f<4;++f) {
        if (count[face_key(el,f)]!=1) continue;
        int local[3]; for (int j=0;j<3;++j) local[j]=v.nodes[tet_face_node(f,j)][el];
        const int k=kind(local);
        if (k<0) continue;
        int smallest=local[0]; for (int j=1;j<3;++j) if (v.key[local[j]]<v.key[smallest]) smallest=local[j];
        if (v.owned[smallest]) o.owned_faces.push_back(int(in.faces.size()));
        in.faces.push_back({el,f,k});
    }
    for (int el=v.element_begin;el<v.element_end;++el) o.owned_elements.push_back(el);
    return part;
}

// Host copy of a tet ElementDomain's state for simple_partition. Coordinates are the domain's
// (SFC-decoded for Tet4 unless the caller overwrites x/y/z with exact values afterwards).
template<class KeyType,class Domain,class Copy>
DomainView<KeyType> domain_view_host(const Domain& d,Copy copy) {
    // Building the lazy key map replaces the input node count with this rank's held-node count.
    const auto& keys=d.getLocalToGlobalSfcMap();
    DomainView<KeyType> v; const std::size_t n=keys.size(), e=d.getElementCount();
    auto read=[&](auto& host,const auto& source,std::size_t count,const char* field) {
        if (source.size()!=count)
            throw std::runtime_error(std::string("SIMPLE domain snapshot: inconsistent ")+field+" size");
        copy(host,source,count,field);
    };
    read(v.key,keys,n,"node keys");
    read(v.x,d.getNodeX(),n,"x coordinates"); read(v.y,d.getNodeY(),n,"y coordinates"); read(v.z,d.getNodeZ(),n,"z coordinates");
    const auto& c=d.getElementToNodeConnectivity();
    read(v.nodes[0],std::get<0>(c),e,"corner 0"); read(v.nodes[1],std::get<1>(c),e,"corner 1");
    read(v.nodes[2],std::get<2>(c),e,"corner 2"); read(v.nodes[3],std::get<3>(c),e,"corner 3");
    v.element_begin=int(d.startIndex()); v.element_end=int(d.endIndex());
    read(v.owned,d.getNodeOwnershipMap(),n,"ownership");
    if (d.numRanks()>1) {   // one rank has no NodeHaloTopology object
        const auto& t=d.getNodeHaloTopology();
        v.peers=t.peers_; v.send_offsets=t.sendOffsets_; v.recv_offsets=t.recvOffsets_;
        if (v.send_offsets.empty() || v.recv_offsets.empty() || v.send_offsets.back()<0 || v.recv_offsets.back()<0)
            throw std::runtime_error("SIMPLE domain snapshot: invalid halo offsets");
        read(v.send_nodes,t.sendNodeIds_,std::size_t(v.send_offsets.back()),"halo send nodes");
        read(v.recv_nodes,t.recvNodeIds_,std::size_t(v.recv_offsets.back()),"halo receive nodes");
    }
    return v;
}

#ifdef MARS_REPLAY_CUDA
template<class RealType,class KeyType>
DomainView<KeyType> domain_view(const ElementDomain<TetTag,RealType,KeyType,cstone::execution::Gpu>& d) {
    return domain_view_host<KeyType>(d,[](auto& host,const auto& device,std::size_t count,const char* field) {
        using S=std::remove_cv_t<std::remove_pointer_t<decltype(thrust::raw_pointer_cast(device.data()))>>;
        std::vector<S> raw_values(count);
        if (count) {
            const auto error=cudaMemcpy(raw_values.data(),thrust::raw_pointer_cast(device.data()),count*sizeof(S),cudaMemcpyDeviceToHost);
            if (error!=cudaSuccess)
                throw std::runtime_error(std::string("SIMPLE domain snapshot: ")+field+": "+cudaGetErrorString(error));
        }
        host.assign(raw_values.begin(),raw_values.end());
    });
}
#endif
} // namespace mars::segregated::runtime
