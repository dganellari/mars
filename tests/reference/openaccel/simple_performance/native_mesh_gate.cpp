#include "mars_segregated_simple_native_mesh.hpp"
#include <iostream>
using namespace mars::segregated::runtime;

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv); int rank,ranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
        Buffer<int> error(1,0);
        Buffer<SourceKey> sorted; Buffer<double> x,y,z;
        if (!rank) { sorted={{10,1},{20,0}}; x={1.25,-2.5}; y={3.75,4.5}; z={-5.25,6.5}; }
        // Non-root ranks deliberately cross a chunk boundary; rank zero can have no requests.
        Buffer<uint64_t> requests(rank || ranks==1?native_mesh_chunk+17:0);
        for (size_t i=0;i<requests.size();++i) requests[i]=i%2?20:10;
        Buffer<NativeCoordinate> response;
        auto coordinates=[&](const uint64_t* query,NativeCoordinate* reply) {
            return NativeCoordinateLookup{query,reply,raw(sorted),int(sorted.size()),raw(x),raw(y),raw(z),raw(error)};
        };
        native_root_lookup(MPI_COMM_WORLD,requests,response,coordinates);
        mesh_check(MPI_COMM_WORLD,error,"coordinate lookup failed");
        bool correct=response.size()==requests.size();
        for (size_t i=0;i<response.size();++i) {
            const int id=i%2?0:1; const auto c=response[i];
            correct=correct && c.source==id && c.x==(id?-2.5:1.25) && c.y==(id?4.5:3.75) && c.z==(id?6.5:-5.25);
        }
        simple_collective(MPI_COMM_WORLD,correct,"exact coordinate or source ID changed");
        Buffer<SourceTag> tags;
        if (!rank) tags={{{{10,20,30}},2,19}};
        Buffer<SimpleFaceKey<uint64_t>> faces={{{10,20,30}},{{20,30,40}}}; Buffer<SourceTag> tagged;
        native_root_lookup(MPI_COMM_WORLD,faces,tagged,[&](const SimpleFaceKey<uint64_t>* query,SourceTag* reply) {
            return NativeTagLookup{query,reply,raw(tags),int(tags.size())};
        });
        simple_collective(MPI_COMM_WORLD,tagged[0].kind==2 && tagged[0].source==19 && tagged[1].kind==-1 && tagged[1].source==-1,
                          "boundary tag lookup changed");
        Buffer<int> coverage(rank?0:ranks,0),ids={rank}; Buffer<unsigned char> ack;
        auto count=[&](const int* query,unsigned char* reply) {
            return NativeCoverageLookup{query,reply,raw(coverage),ranks,raw(error)};
        };
        native_root_lookup(MPI_COMM_WORLD,ids,ack,count);
        if (!rank) launch(ranks,CoverageCheck{raw(coverage),raw(error)});
        mesh_check(MPI_COMM_WORLD,error,"coverage failed");
        native_root_lookup(MPI_COMM_WORLD,ids,ack,count);
        if (!rank) launch(ranks,CoverageCheck{raw(coverage),raw(error)});
        bool rejected=false; try { mesh_check(MPI_COMM_WORLD,error,"duplicate coverage"); } catch (...) { rejected=true; }
        simple_collective(MPI_COMM_WORLD,rejected,"duplicate coverage was accepted");
        error[0]=0; requests={999}; native_root_lookup(MPI_COMM_WORLD,requests,response,coordinates);
        rejected=false; try { mesh_check(MPI_COMM_WORLD,error,"unknown SFC key"); } catch (...) { rejected=true; }
        simple_collective(MPI_COMM_WORLD,rejected,"unknown source key was accepted");
        if (argc==2) {
            const auto source=read_simple_mesh_root(MPI_COMM_WORLD,argv[1]);
            simple_collective(MPI_COMM_WORLD,source.global_nodes==425 && source.global_elements==1536 && source.global_faces==576,
                              "public fixture dimensions changed");
            simple_collective(MPI_COMM_WORLD,rank==0 || (source.x.empty() && source.y.empty() && source.z.empty()
                && source.nodes[0].empty() && source.nodes[1].empty() && source.nodes[2].empty() && source.nodes[3].empty() && source.faces.empty()),
                "source payload was replicated");
            const auto local=native_initial_partition(MPI_COMM_WORLD,source);
            const auto reference=read_simple_mesh(MPI_COMM_WORLD,argv[1]);
            const int first=int(static_cast<long long>(source.global_elements)*rank/ranks);
            const int last=int(static_cast<long long>(source.global_elements)*(rank+1)/ranks);
            std::set<int> used;
            for (int e=first;e<last;++e) for (int j=0;j<4;++j) used.insert(reference.nodes[j][e]);
            correct=local.nodes[0].size()==size_t(last-first) && local.x.size()==used.size();
            for (int e=0;e<last-first;++e) for (int j=0;j<4;++j) {
                const int a=local.nodes[j][e],b=reference.nodes[j][first+e];
                correct=correct && a==int(std::distance(used.begin(),used.find(b))) && local.x[a]==reference.x[b] && local.y[a]==reference.y[b] && local.z[a]==reference.z[b];
            }
            simple_collective(MPI_COMM_WORLD,correct,"initial partition changed coordinates or element order");
            rejected=false; try { read_simple_mesh_root(MPI_COMM_WORLD,std::string(argv[1])+".missing"); } catch (...) { rejected=true; }
            simple_collective(MPI_COMM_WORLD,rejected,"root reader fault was not collective");
        }
        if (!rank) std::cout<<"PASS native root routing, exact identities, tags, chunk boundary and collective faults ranks="<<ranks<<'\n';
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
