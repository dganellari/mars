#include "mars_segregated_simple_partition.hpp"
#include <iostream>
#include <string>

using namespace mars::segregated::runtime;
using mars::segregated::distributed::FieldExchange;

void require(MPI_Comm comm,bool ok,const char* message) { simple_collective(comm,ok,message); }

void ring(MPI_Comm comm,bool isolate_last,const std::string& fault) {
    int rank=0, ranks=0; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    const int active=ranks-(isolate_last && ranks>2?1:0);
    std::vector<int> peers, so{0}, ro{0}, sn, rn;
    if (rank<active && active>1) for (int q=0;q<active;++q) {
        const bool sends=q==(rank+1)%active, receives=q==(rank+active-1)%active;
        if (!sends && !receives) continue;
        peers.push_back(q);
        if (sends) sn.push_back(0);
        if (receives) rn.push_back(1);
        so.push_back(int(sn.size())); ro.push_back(int(rn.size()));
    }
    bool reject=false;
    if (rank==0 && ranks>1) {
        if (fault=="asymmetric") { peers.clear(); so={0}; ro={0}; sn.clear(); rn.clear(); reject=true; }
        if (fault=="counts") { sn.push_back(0); ++so.back(); reject=true; }
        if (fault=="offsets") { so.front()=1; reject=true; }
        if (fault=="duplicate-peer") { peers.push_back(peers.back()); so.push_back(so.back()); ro.push_back(ro.back()); reject=true; }
        if (fault=="duplicate-recv") { rn.push_back(1); ++ro.back(); reject=true; }
        if (fault=="node-overflow") reject=true;
    }
    int any_reject=reject; MPI_Allreduce(MPI_IN_PLACE,&any_reject,1,MPI_INT,MPI_MAX,comm);
    int caught=0;
    try {
        FieldExchange exchange(comm,peers,so,sn,ro,rn,rank==0 && fault=="node-overflow"?std::numeric_limits<int>::max():2);
        if (any_reject) throw std::logic_error("bad lists accepted");
        Array<double> scalar(std::vector<double>{rank+0.25,-9});
        Array<double> vector(std::vector<double>{double(rank),double(rank)+1,double(rank)+2,-9,-9,-9});
        if (fault=="stride" || fault=="reverse-stride") {
            if (fault=="stride") exchange({{scalar.data(),rank==0?9:1}});
            else exchange.reverse_add({scalar.data(),rank==0?0:1});
            throw std::logic_error("invalid field did not abort");
        }
        exchange({{scalar.data(),1},{vector.data(),3}});
        auto s=scalar.host(), v=vector.host();
        const int previous=(rank+active-1)%active;
        require(comm,s[0]==rank+0.25 && v[0]==rank && (rn.empty() ||
                (s[1]==previous+0.25 && v[3]==previous && v[4]==previous+1 && v[5]==previous+2)),"fused publish mismatch");
        exchange.reverse_add({scalar.data(),1});
        s=scalar.host();
        require(comm,s[0]==(sn.empty()?1:2)*(rank+0.25),"reverse sum mismatch");
        std::vector<std::uint64_t> keys{(UINT64_C(1)<<62)+std::uint64_t(rank),0};
        exchange.publish_host(keys);
        require(comm,rn.empty() || keys[1]==(UINT64_C(1)<<62)+std::uint64_t(previous),"integer metadata lost bits");
    } catch (const std::runtime_error& e) {
        if (any_reject && std::string(e.what()).find("halo exchange lists rejected")==0) caught=1;
        else throw;
    }
    require(comm,caught==any_reject,"list rejection was not collective");
}

void partition_identity(MPI_Comm comm,bool swap_slots) {
    int rank=0, ranks=0; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    require(comm,ranks<=4 && (!swap_slots || ranks==2),"identity fixture needs 1/2/4 ranks (swap: 2)");
    DomainView<std::uint64_t> v;
    v.x={0,1,0,0}; v.y={0,0,1,0}; v.z={0,0,0,1};
    for (int i=0;i<4;++i) { v.key.push_back((UINT64_C(1)<<62)+i); v.nodes[i]={i}; v.owned.push_back(i%ranks==rank); }
    v.element_begin=rank==0?0:1; v.element_end=1;
    for (int q=0;q<ranks;++q) if (q!=rank) {
        v.peers.push_back(q);
        for (int i=0;i<4;++i) {
            if (i%ranks==rank) v.send_nodes.push_back(i);
            if (i%ranks==q) v.recv_nodes.push_back(i);
        }
        v.send_offsets.push_back(int(v.send_nodes.size())); v.recv_offsets.push_back(int(v.recv_nodes.size()));
    }
    if (swap_slots && rank==1) std::swap(v.recv_nodes[0],v.recv_nodes[1]);
    int caught=0;
    try {
        const auto part=simple_partition<long long>(comm,v,[](const int*) { return 2; });
        for (int i=0;i<4;++i) {
            int expected=0;
            for (int j=0;j<4;++j) expected+=(j%ranks<i%ranks) || (j%ranks==i%ranks && j<i);
            require(comm,part.ownership.solver_node[i]==expected,"wrong solver identity");
        }
    } catch (const std::runtime_error& e) {
        if (swap_slots && std::string(e.what()).find("halo peer node keys do not match")!=std::string::npos) caught=1;
        else throw;
    }
    require(comm,caught==int(swap_slots),"swapped halo slots were not rejected collectively");
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank=0, ranks=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
#ifdef MARS_REPLAY_CUDA
    MPI_Comm node; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
    int local=0, devices=0; MPI_Comm_rank(node,&local); MPI_Comm_free(&node);
    mars::segregated::assembly_cuda_check(cudaGetDeviceCount(&devices));
    if (!devices) MPI_Abort(MPI_COMM_WORLD,1);
    mars::segregated::assembly_cuda_check(cudaSetDevice(local%devices));
#endif
        const std::string fault=argc>1?argv[1]:"";
        if (fault=="identity") partition_identity(MPI_COMM_WORLD,true);
        else ring(MPI_COMM_WORLD,false,fault);
        if (fault.empty()) {
            ring(MPI_COMM_WORLD,true,"");
            partition_identity(MPI_COMM_WORLD,false);
        }
        if (rank==0) std::cout<<"PASS: node-field exchange ranks="<<ranks<<" fault="<<(fault.empty()?"none":fault)<<std::endl;
    } catch (const std::exception& e) { std::cerr<<e.what()<<std::endl; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
