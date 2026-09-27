#include "mars_segregated_simple_partition.hpp"
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>

using mars::segregated::runtime::domain_view_host;

namespace {
void require(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}

struct Topology {
    std::vector<int> peers_,sendOffsets_{0},recvOffsets_{0},sendNodeIds_,recvNodeIds_;
};

// The mesh input is replicated, but its lazy local map may shrink or grow after redistribution.
struct LazyDomain {
    mutable std::size_t node_count;
    mutable bool built=false;
    int ranks;
    std::vector<std::uint64_t> keys;
    std::vector<float> x,y,z;
    std::array<std::vector<std::uint64_t>,4> connectivity;
    std::vector<unsigned char> owned;
    Topology topology;

    LazyDomain(std::size_t input_count,std::size_t held_count,int rank_count)
        :node_count(input_count),ranks(rank_count),keys(held_count),x(held_count),y(held_count),z(held_count),owned(held_count,1) {
        for (std::size_t i=0;i<held_count;++i) {
            keys[i]=(std::uint64_t(1)<<54)+i;
            x[i]=float(i); y[i]=float(i)+0.25f; z[i]=-float(i);
        }
        if (held_count>=4) {
            for (std::size_t k=0;k<4;++k) connectivity[k]={k};
            if (ranks>1) {
                owned.back()=0;
                topology={{1},{0,1},{0,1},{0},{int(held_count-1)}};
            }
        }
    }
    std::size_t getNodeCount() const { return node_count; }
    std::size_t getElementCount() const { return connectivity[0].size(); }
    const auto& getLocalToGlobalSfcMap() const { built=true; node_count=keys.size(); return keys; }
    const auto& getNodeX() const { require(built,"coordinates requested before the local map"); return x; }
    const auto& getNodeY() const { return y; }
    const auto& getNodeZ() const { return z; }
    const auto& getElementToNodeConnectivity() const { return connectivity; }
    const auto& getNodeOwnershipMap() const { return owned; }
    const auto& getNodeHaloTopology() const { require(ranks>1,"one-rank topology requested"); return topology; }
    std::size_t startIndex() const { return 0; }
    std::size_t endIndex() const { return getElementCount(); }
    int numRanks() const { return ranks; }
};

struct HostCopy {
    std::vector<std::string> fields;
    template<class H,class D> void operator()(H& host,const D& source,std::size_t count,const char* field) {
        require(count==source.size(),"invalid extent reached the transfer");
        fields.emplace_back(field);
        host.assign(source.begin(),source.end());
    }
};

void snapshot(std::size_t input_count,std::size_t held_count,int ranks) {
    LazyDomain domain(input_count,held_count,ranks);
    HostCopy copy;
    const auto view=domain_view_host<std::uint64_t>(domain,std::ref(copy));
    require(domain.built && domain.node_count==held_count,"local map not initialized");
    require(view.key==domain.keys && view.x.size()==held_count,"snapshot uses input count");
    require(view.owned==domain.owned,"ownership differs");
    for (std::size_t i=0;i<held_count;++i)
        require(view.x[i]==double(domain.x[i]) && view.y[i]==double(domain.y[i]) && view.z[i]==double(domain.z[i]),
                "coordinate conversion differs");
    for (std::size_t k=0;k<4;++k)
        require(view.nodes[k]==std::vector<int>(domain.connectivity[k].begin(),domain.connectivity[k].end()),
                "connectivity conversion differs");
    require(view.element_begin==0 && view.element_end==int(domain.getElementCount()),"element range differs");
    require(view.peers==domain.topology.peers_ && view.send_nodes==domain.topology.sendNodeIds_
            && view.recv_nodes==domain.topology.recvNodeIds_,"halo lists differ");
    require(view.send_offsets==domain.topology.sendOffsets_ && view.recv_offsets==domain.topology.recvOffsets_,
            "halo offsets differ");
}

template<class Change> void invalid(Change change,const char* field) {
    LazyDomain domain(425,4,2);
    change(domain);
    HostCopy copy;
    bool rejected=false;
    try { (void)domain_view_host<std::uint64_t>(domain,std::ref(copy)); }
    catch (const std::runtime_error& e) {
        rejected=std::string(e.what()).find(field)!=std::string::npos;
    }
    require(rejected,"invalid extent was not rejected with its field name");
    require(std::find(copy.fields.begin(),copy.fields.end(),field)==copy.fields.end(),"invalid field was transferred");
}
} // namespace

int main() {
    try {
        snapshot(425,4,2);
        snapshot(1,6,2);
        snapshot(425,0,2);
        snapshot(4,4,1);
        invalid([](auto& d) { d.x.pop_back(); },"x coordinates");
        invalid([](auto& d) { d.connectivity[2].clear(); },"corner 2");
        invalid([](auto& d) { d.owned.push_back(1); },"ownership");
        invalid([](auto& d) { d.topology.sendNodeIds_.clear(); },"halo send nodes");
        invalid([](auto& d) { d.topology.recvNodeIds_.clear(); },"halo receive nodes");
        invalid([](auto& d) { d.topology.sendOffsets_.clear(); },"halo offsets");
        invalid([](auto& d) { d.topology.recvOffsets_.back()=-1; },"halo offsets");
        std::cout<<"PASS: 11 lazy domain snapshot and transfer-bound cases\n";
    } catch (const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
