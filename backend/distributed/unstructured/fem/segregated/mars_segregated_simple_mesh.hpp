#pragma once
// Device partition construction. The host branch is an executable oracle for the same kernels.
#include "mars_segregated_simple_distributed.hpp"
#include "mars_segregated_native_mapping.hpp"
#include <climits>
#ifdef MARS_REPLAY_CUDA
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#endif

#if defined(__CUDACC__)
#define MARS_SMESH_HD __host__ __device__
#else
#define MARS_SMESH_HD
#endif
namespace mars::segregated::runtime {
using distributed::Buffer;
using distributed::raw;

template<class KeyType> struct SimpleMeshData {
    Buffer<double> x,y,z;
    Buffer<KeyType> key;
    std::array<Buffer<int>,4> nodes;
    Buffer<unsigned char> owned;
    int element_begin=0,element_end=0;
    std::vector<int> peers,send_offsets{0},recv_offsets{0};
    Buffer<int> send_nodes,recv_nodes;
};
struct SimpleDeviceInput {
    Buffer<double> x,y,z;
    std::array<Buffer<int>,4> nodes;
    Buffer<SimpleFace> faces;
};
template<class GlobalId> struct NativeSimplePartition {
    SimpleDeviceInput input;
    SimpleOwnership<GlobalId,Buffer> ownership;
};

// Only error words return to the CPU. Topology and identities never do.
inline void mesh_check(MPI_Comm comm,Buffer<int>& status,const char* message) {
    int faults=0;
    const int local=distributed::fetch(raw(status),{},faults);
    simple_collective(comm,local==0 && faults==0,message);
}
template<class T,class F> void mesh_sort(Buffer<T>& values,F less) {
#ifdef MARS_REPLAY_CUDA
    thrust::sort(values.begin(),values.end(),less);
#else
    std::sort(values.begin(),values.end(),less);
#endif
}
template<class T,class Predicate> void mesh_compact(Buffer<T>& values,Predicate predicate) {
    Buffer<T> out(values.size());
#ifdef MARS_REPLAY_CUDA
    auto end=thrust::copy_if(values.begin(),values.end(),out.begin(),predicate);
#else
    auto end=std::copy_if(values.begin(),values.end(),out.begin(),predicate);
#endif
    out.resize(std::size_t(end-out.begin())); values=std::move(out);
}
inline Buffer<int> mesh_sequence(int first,int last) {
    Buffer<int> result(std::size_t(last-first));
#ifdef MARS_REPLAY_CUDA
    thrust::sequence(result.begin(),result.end(),first);
#else
    std::iota(result.begin(),result.end(),first);
#endif
    return result;
}
template<class KeyType> struct MeshNodeCheck {
    const KeyType* keys; const unsigned char* owned; int* status;
    MARS_SMESH_HD void operator()(int i) const {
        if ((i && !(keys[i-1]<keys[i])) || owned[i]>1) distributed::raise_fault(status,1);
    }
};
struct MeshCellCheck {
    const int* nodes[4]; int count; int* status;
    MARS_SMESH_HD void operator()(int e) const {
        for (int i=0;i<4;++i) {
            if (nodes[i][e]<0 || nodes[i][e]>=count) distributed::raise_fault(status,1);
            for (int j=0;j<i;++j) if (nodes[i][e]==nodes[j][e]) distributed::raise_fault(status,1);
        }
    }
};
struct MeshOwned { const unsigned char* owned; MARS_SMESH_HD bool operator()(int i) const { return owned[i]!=0; } };
struct MeshHaloMark {
    const int* nodes; const unsigned char* owned; int* received; int* status; bool receive;
    MARS_SMESH_HD void operator()(int i) const {
        const int n=nodes[i];
        if (receive) {
            if (owned[n]) distributed::raise_fault(status,1);
#if defined(__CUDA_ARCH__)
            atomicAdd(received+n,1);
#else
            ++received[n];
#endif
        } else if (!owned[n]) distributed::raise_fault(status,1);
    }
};
struct MeshGhostCheck {
    const unsigned char* owned; const int* received; int* status;
    MARS_SMESH_HD void operator()(int i) const {
        if (received[i]!=(owned[i]?0:1)) distributed::raise_fault(status,1);
    }
};
template<class T> struct MeshEqual {
    const T *a,*b; int* status;
    MARS_SMESH_HD void operator()(int i) const { if (a[i]!=b[i]) distributed::raise_fault(status,1); }
};
template<class GlobalId> struct MeshSolverId {
    const int* nodes; GlobalId* ids; long long first;
    MARS_SMESH_HD void operator()(int i) const { ids[nodes[i]]=GlobalId(first+i); }
};
template<class GlobalId> struct MeshIdCheck {
    const GlobalId* ids; long long total; int* status;
    MARS_SMESH_HD void operator()(int i) const {
        if (ids[i]<0 || static_cast<unsigned long long>(ids[i])>=static_cast<unsigned long long>(total)) distributed::raise_fault(status,1);
    }
};
struct MeshIncidence {
    const int* nodes[4]; int first,last; double *own,*held;
    MARS_SMESH_HD void operator()(int e) const {
        for (int j=0;j<4;++j) {
            assembly_add(held+nodes[j][e],1.);
            if (e>=first && e<last) assembly_add(own+nodes[j][e],1.);
        }
    }
};
struct MeshStarCheck {
    const unsigned char* owned; const double *own,*held; int* status;
    MARS_SMESH_HD void operator()(int i) const { if (owned[i] && own[i]!=held[i]) distributed::raise_fault(status,1); }
};
template<class KeyType> struct SimpleFaceKey {
    KeyType nodes[3]{};
    MARS_SMESH_HD bool operator<(const SimpleFaceKey& b) const {
        for (int j=0;j<3;++j) { if (nodes[j]<b.nodes[j]) return true; if (b.nodes[j]<nodes[j]) return false; }
        return false;
    }
    MARS_SMESH_HD bool operator==(const SimpleFaceKey& b) const {
        return nodes[0]==b.nodes[0] && nodes[1]==b.nodes[1] && nodes[2]==b.nodes[2];
    }
};
template<class KeyType> MARS_SMESH_HD SimpleFaceKey<KeyType> mesh_face_key(KeyType a,KeyType b,KeyType c) {
    if (b<a) { auto t=a; a=b; b=t; }
    if (c<b) { auto t=b; b=c; c=t; }
    if (b<a) { auto t=a; a=b; b=t; }
    return {{a,b,c}};
}
template<class KeyType> struct MeshFaceRecord {
    SimpleFaceKey<KeyType> key; int element,ordinal;
};
template<class KeyType> struct MeshFaceLess {
    MARS_SMESH_HD bool operator()(const MeshFaceRecord<KeyType>& a,const MeshFaceRecord<KeyType>& b) const { return a.key<b.key; }
};
template<class KeyType> struct MeshFaces {
    const int* nodes[4]; const KeyType* keys; MeshFaceRecord<KeyType>* faces;
    MARS_SMESH_HD void operator()(int i) const {
        const int e=i/4,f=i%4;
        faces[i]={mesh_face_key(keys[nodes[tet_face_node(f,0)][e]],keys[nodes[tet_face_node(f,1)][e]],keys[nodes[tet_face_node(f,2)][e]]),e,f};
    }
};
template<class KeyType,class Kind> struct MeshBoundary {
    const MeshFaceRecord<KeyType>* records; int count; const int* nodes[4]; const unsigned char* owned;
    Kind kind; SimpleFace* faces; int* status;
    MARS_SMESH_HD void operator()(int i) const {
        const auto r=records[i];
        const bool previous=i>0 && r.key==records[i-1].key, next=i+1<count && r.key==records[i+1].key;
        if (previous && next) distributed::raise_fault(status,1); // nonmanifold
        if (previous || next) return;
        int local[3]; bool touches=false;
        for (int j=0;j<3;++j) { local[j]=nodes[tet_face_node(r.ordinal,j)][r.element]; touches=touches || owned[local[j]]; }
        const int tag=kind(local);
        if (tag<0) { if (touches) distributed::raise_fault(status,1); return; }
        if (tag>2) { distributed::raise_fault(status,1); return; }
        faces[4*r.element+r.ordinal]={r.element,r.ordinal,tag};
    }
};
struct MeshTagged { MARS_SMESH_HD bool operator()(SimpleFace f) const { return f.kind>=0; } };
template<class KeyType> struct MeshOwnedFace {
    const SimpleFace* faces; const int* nodes[4]; const KeyType* keys; const unsigned char* owned;
    MARS_SMESH_HD bool operator()(int i) const {
        const auto f=faces[i]; int smallest=nodes[tet_face_node(f.ordinal,0)][f.element];
        for (int j=1;j<3;++j) { const int n=nodes[tet_face_node(f.ordinal,j)][f.element]; if (keys[n]<keys[smallest]) smallest=n; }
        return owned[smallest]!=0;
    }
};

template<class GlobalId,class KeyType,class Kind>
NativeSimplePartition<GlobalId> build_simple_partition(MPI_Comm comm,const SimpleMeshData<KeyType>& v,Kind kind) {
    const std::size_t ns=v.key.size(),es=v.nodes[0].size();
    bool sizes=ns<=std::size_t(INT_MAX/9) && es<=std::size_t(INT_MAX/16) && v.x.size()==ns && v.y.size()==ns && v.z.size()==ns && v.owned.size()==ns;
    for (int k=1;k<4;++k) sizes=sizes && v.nodes[k].size()==es;
    sizes=sizes && v.element_begin>=0 && v.element_end>=v.element_begin && std::size_t(v.element_end)<=es;
    simple_collective(comm,sizes,"SIMPLE mesh: inconsistent local array sizes");
    const int n=int(ns),e=int(es); Buffer<int> status(1,0);
    launch(n,MeshNodeCheck<KeyType>{raw(v.key),raw(v.owned),raw(status)});
    launch(e,MeshCellCheck{{raw(v.nodes[0]),raw(v.nodes[1]),raw(v.nodes[2]),raw(v.nodes[3])},n,raw(status)});
    mesh_check(comm,status,"SIMPLE mesh: invalid keys, ownership mask or connectivity");
    NativeSimplePartition<GlobalId> part; auto& o=part.ownership;
    o.peers=v.peers; o.send_offsets=v.send_offsets; o.recv_offsets=v.recv_offsets; o.send_nodes=v.send_nodes; o.recv_nodes=v.recv_nodes;
    distributed::FieldExchange exchange(comm,o.peers,o.send_offsets,o.recv_offsets,
        {raw(o.send_nodes),o.send_nodes.size(),raw(o.recv_nodes),o.recv_nodes.size()},n);
    Buffer<int> received(ns,0);
    launch(int(o.recv_nodes.size()),MeshHaloMark{raw(o.recv_nodes),raw(v.owned),raw(received),raw(status),true});
    launch(int(o.send_nodes.size()),MeshHaloMark{raw(o.send_nodes),raw(v.owned),raw(received),raw(status),false});
    launch(n,MeshGhostCheck{raw(v.owned),raw(received),raw(status)});
    mesh_check(comm,status,"SIMPLE mesh: halo must send owned nodes and receive every ghost once");
    auto received_keys=v.key; exchange.publish(raw(received_keys),ns);
    launch(n,MeshEqual<KeyType>{raw(v.key),raw(received_keys),raw(status)});
    mesh_check(comm,status,"SIMPLE mesh: halo peer node keys do not match");
    o.owned_nodes=mesh_sequence(0,n); mesh_compact(o.owned_nodes,MeshOwned{raw(v.owned)});
    const long long mine=static_cast<long long>(o.owned_nodes.size()); long long first=0,total=0; int rank;
    MPI_Comm_rank(comm,&rank); MPI_Exscan(&mine,&first,1,MPI_LONG_LONG,MPI_SUM,comm); if (!rank) first=0;
    MPI_Allreduce(&mine,&total,1,MPI_LONG_LONG,MPI_SUM,comm);
    simple_collective(comm,total>0 && static_cast<unsigned long long>(total-1)<=static_cast<unsigned long long>(std::numeric_limits<GlobalId>::max()),
                      "SIMPLE mesh: solver node ids exceed the global index type");
    o.solver_node.resize(ns,GlobalId(-1));
    launch(int(mine),MeshSolverId<GlobalId>{raw(o.owned_nodes),raw(o.solver_node),first});
    exchange.publish(raw(o.solver_node),ns);
    launch(n,MeshIdCheck<GlobalId>{raw(o.solver_node),total,raw(status)});
    mesh_check(comm,status,"SIMPLE mesh: a node received no valid solver id");
    Buffer<double> own(ns,0),held(ns,0);
    launch(e,MeshIncidence{{raw(v.nodes[0]),raw(v.nodes[1]),raw(v.nodes[2]),raw(v.nodes[3])},v.element_begin,v.element_end,raw(own),raw(held)});
    exchange.reverse_add({raw(own),1});
    launch(n,MeshStarCheck{raw(v.owned),raw(own),raw(held),raw(status)});
    mesh_check(comm,status,"SIMPLE mesh: incomplete element star at an owned node");
    Buffer<MeshFaceRecord<KeyType>> records(4*es);
    launch(4*e,MeshFaces<KeyType>{{raw(v.nodes[0]),raw(v.nodes[1]),raw(v.nodes[2]),raw(v.nodes[3])},raw(v.key),raw(records)});
    mesh_sort(records,MeshFaceLess<KeyType>{});
    auto& f=part.input; f.x=v.x; f.y=v.y; f.z=v.z; f.nodes=v.nodes; f.faces.resize(4*es,SimpleFace{0,0,-1});
    launch(4*e,MeshBoundary<KeyType,Kind>{raw(records),4*e,{raw(v.nodes[0]),raw(v.nodes[1]),raw(v.nodes[2]),raw(v.nodes[3])},raw(v.owned),kind,raw(f.faces),raw(status)});
    mesh_check(comm,status,"SIMPLE mesh: invalid exterior boundary coverage");
    mesh_compact(f.faces,MeshTagged{});
    o.owned_faces=mesh_sequence(0,int(f.faces.size()));
    mesh_compact(o.owned_faces,MeshOwnedFace<KeyType>{raw(f.faces),{raw(v.nodes[0]),raw(v.nodes[1]),raw(v.nodes[2]),raw(v.nodes[3])},raw(v.key),raw(v.owned)});
    o.owned_elements=mesh_sequence(v.element_begin,v.element_end);
    return part;
}
} // namespace mars::segregated::runtime
#undef MARS_SMESH_HD
