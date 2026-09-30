#include "channel.hpp"
#include "mars_segregated_simple_mesh.hpp"
#include <iostream>
using namespace mars::segregated;
using namespace mars::segregated::runtime;
#ifdef MARS_REPLAY_CUDA
#define GATE_HD __host__ __device__
#else
#define GATE_HD
#endif
struct PlaneTags {
    const double *x,*y,*z; double length; bool missing;
    GATE_HD int operator()(const int* f) const {
        bool inlet=true,outlet=true,y0=true,y1=true,z0=true,z1=true;
        for (int j=0;j<3;++j) {
            const int i=f[j]; inlet=inlet && x[i]==0; outlet=outlet && x[i]==length;
            y0=y0 && y[i]==0; y1=y1 && y[i]==1; z0=z0 && z[i]==0; z1=z1 && z[i]==1;
        }
        return inlet?(missing?-1:0):outlet?1:y0 || y1 || z0 || z1?2:-1;
    }
};
template<class T> std::vector<T> host(const Buffer<T>& values) {
#ifdef MARS_REPLAY_CUDA
    std::vector<T> out(values.size()); thrust::copy(values.begin(),values.end(),out.begin()); return out;
#else
    return values;
#endif
}
template<class T> void upload(Buffer<T>& out,const std::vector<T>& in) { out.assign(in.begin(),in.end()); }
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv); int rank,ranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
#ifdef MARS_REPLAY_CUDA
        int devices=0,local; MPI_Comm shared; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&shared);
        MPI_Comm_rank(shared,&local); MPI_Comm_free(&shared);
        assembly_cuda_check(cudaGetDeviceCount(&devices)); ensure(devices>0,"no CUDA device"); assembly_cuda_check(cudaSetDevice(local%devices));
#endif
        const std::string fault=argc>1?argv[1]:"none";
        const auto source=dsimple_gate::channel(8,2,2); const auto ownership=dsimple_gate::partition(source,ranks);
        auto reference=dsimple_gate::build(source,ownership,rank);
        SimpleMeshData<long long> v;
        upload(v.x,reference.input.x); upload(v.y,reference.input.y); upload(v.z,reference.input.z);
        std::vector<long long> keys; for (int id:reference.node_global) keys.push_back((1LL<<54)+id);
        upload(v.key,keys);
        std::vector<unsigned char> mask(reference.node_owned.begin(),reference.node_owned.end()); upload(v.owned,mask);
        auto& o=reference.ownership;
        v.element_begin=int(std::count_if(reference.element_global.begin(),reference.element_global.end(),[&](int e) { return ownership.element_owner[e]<rank; }));
        v.element_end=v.element_begin+int(o.owned_elements.size());
        bool altered=false;
        if (fault=="star" && rank==0 && ranks>1) {
            for (int e=v.element_end;e<int(reference.input.nodes[0].size());++e) {
                bool touches=false; for (int j=0;j<4;++j) touches=touches || mask[reference.input.nodes[j][e]];
                if (!touches) continue;
                for (auto& column:reference.input.nodes) column.erase(column.begin()+e);
                altered=true; break;
            }
        }
        for (int j=0;j<4;++j) upload(v.nodes[j],reference.input.nodes[j]);
        v.peers=o.peers; v.send_offsets=o.send_offsets; v.recv_offsets=o.recv_offsets;
        upload(v.send_nodes,o.send_nodes); upload(v.recv_nodes,o.recv_nodes);
        if (fault=="identity" && rank==0 && o.recv_nodes.size()>1) {
            auto wrong=o.recv_nodes; std::swap(wrong[0],wrong[1]); upload(v.recv_nodes,wrong); altered=true;
        }
        if (fault=="tag" && rank==0) altered=true;
        int mutation=altered?1:0; MPI_Allreduce(MPI_IN_PLACE,&mutation,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        if (fault!="none") ensure(mutation>0,"fault injection was not exercised");
        bool caught=false;
        try {
            auto result=build_simple_partition<long long>(MPI_COMM_WORLD,v,PlaneTags{raw(v.x),raw(v.y),raw(v.z),4.,fault=="tag" && rank==0});
            bool equal=host(result.ownership.owned_nodes)==o.owned_nodes && host(result.ownership.solver_node)==o.solver_node
                && host(result.ownership.owned_elements)==o.owned_elements && host(result.ownership.owned_faces)==o.owned_faces;
            const auto faces=host(result.input.faces);
            equal=equal && faces.size()==reference.input.faces.size();
            for (size_t i=0;i<faces.size() && equal;++i) {
                const auto a=faces[i],b=reference.input.faces[i]; equal=a.element==b.element && a.ordinal==b.ordinal && a.kind==b.kind;
            }
            simple_collective(MPI_COMM_WORLD,equal,"device mesh builder differs from independent host builder");
        } catch (const std::runtime_error& e) {
            const std::string message=e.what();
            caught=(fault=="star" && message.find("incomplete element star")!=std::string::npos)
                || (fault=="identity" && message.find("peer node keys")!=std::string::npos)
                || (fault=="tag" && message.find("exterior boundary coverage")!=std::string::npos);
            if (!caught) throw;
        }
        simple_collective(MPI_COMM_WORLD,caught==(fault!="none"),"fault was not rejected collectively");
        if (!rank) std::cout<<"PASS: native partition topology, ownership, integer identities and complete stars ranks="<<ranks<<" fault="<<fault<<'\n';
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
