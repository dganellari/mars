#pragma once
#include "mars_segregated_simple_mesh.hpp"
#include "../../utils/mars_read_exodus_raw.hpp"
#if defined(__CUDACC__)
#define MARS_SINPUT_HD __host__ __device__
#else
#define MARS_SINPUT_HD
#endif
namespace mars::segregated::runtime {
inline void mesh_mpi_ready() {
#ifdef MARS_REPLAY_CUDA
    assembly_cuda_check(cudaStreamSynchronize(nullptr));
#endif
}
template<class T> void broadcast_mesh_array(MPI_Comm comm,const std::vector<T>& file,Buffer<T>& values,MPI_Datatype type) {
    int rank; MPI_Comm_rank(comm,&rank);
    unsigned long long count=file.size();
    ensure(MPI_Bcast(&count,1,MPI_UNSIGNED_LONG_LONG,0,comm)==MPI_SUCCESS,"mesh size broadcast failed");
    ensure(count<=static_cast<unsigned long long>(INT_MAX),"mesh array exceeds MPI count capacity");
    values.resize(size_t(count));
    if (!rank) values.assign(file.begin(),file.end());
    mesh_mpi_ready();
    ensure(MPI_Bcast(raw(values),int(count),type,0,comm)==MPI_SUCCESS,"mesh device broadcast failed");
}
struct SourceCoordinatesCheck {
    const double *x,*y,*z; int* error;
    MARS_SINPUT_HD void operator()(int n) const {
        if (!(geometry_finite(x[n]) && geometry_finite(y[n]) && geometry_finite(z[n]))) distributed::raise_fault(error,1);
    }
};
struct SourceConnectivity {
    const long long* source; int* nodes[4]; int count; int* error;
    MARS_SINPUT_HD void operator()(int e) const {
        for (int j=0;j<4;++j) {
            const long long id=source[4*e+j];
            if (id<1 || id>count) { distributed::raise_fault(error,1); nodes[j][e]=-1; }
            else nodes[j][e]=int(id-1);
        }
    }
};
struct SourceBoundary {
    const long long *elements,*sides; int element_count,kind,offset; SimpleFace* faces; int* tags; int* error;
    MARS_SINPUT_HD void operator()(int i) const {
        const long long element=elements[i],side=sides[i];
        if (element<1 || element>element_count || side<1 || side>4) { distributed::raise_fault(error,1); return; }
        const int e=int(element-1),f=int(side-1);
        faces[offset+i]={e,f,kind};
#if defined(__CUDA_ARCH__)
        if (atomicCAS(tags+4*e+f,-1,kind)!=-1) distributed::raise_fault(error,1);
#else
        if (tags[4*e+f]!=-1) distributed::raise_fault(error,1);
        tags[4*e+f]=kind;
#endif
    }
};
struct SourceExteriorCheck {
    const MeshFaceRecord<int>* faces; const int* tags; int count; int* error;
    MARS_SINPUT_HD void operator()(int i) const {
        const auto f=faces[i];
        const bool previous=i>0 && faces[i-1].key==f.key,next=i+1<count && faces[i+1].key==f.key;
        const int tag=tags[4*f.element+f.ordinal];
        if ((previous && next) || ((tag>=0)!=( !previous && !next))) distributed::raise_fault(error,1);
    }
};

// Rank zero reads bytes. Broadcasts, index conversion and all topology validation use device buffers.
inline SimpleDeviceInput read_simple_mesh(MPI_Comm comm,const std::string& path) {
    int rank; MPI_Comm_rank(comm,&rank);
    ExodusRawTet4 file;
    bool ok=true; std::string message;
    if (!rank) try { file=readExodusTet4Raw(path); } catch (const std::exception& e) { ok=false; message=e.what(); }
    if (!ok) std::cerr<<message<<'\n';
    simple_collective(comm,ok,"native Exodus read failed");
    SimpleDeviceInput source;
    broadcast_mesh_array(comm,file.coordinates[0],source.x,MPI_DOUBLE);
    broadcast_mesh_array(comm,file.coordinates[1],source.y,MPI_DOUBLE);
    broadcast_mesh_array(comm,file.coordinates[2],source.z,MPI_DOUBLE);
    Buffer<long long> connectivity;
    broadcast_mesh_array(comm,file.connectivity,connectivity,MPI_LONG_LONG);
    const size_t n=source.x.size(),e=connectivity.size()/4;
    simple_collective(comm,n>0 && n<=size_t(INT_MAX/9) && source.y.size()==n && source.z.size()==n && e>0 && e<=size_t(INT_MAX/16)
                      && connectivity.size()==4*e,"invalid native mesh sizes");
    Buffer<int> error(1,0);
    launch(int(n),SourceCoordinatesCheck{raw(source.x),raw(source.y),raw(source.z),raw(error)});
    for (auto& c:source.nodes) c.resize(e);
    launch(int(e),SourceConnectivity{raw(connectivity),{raw(source.nodes[0]),raw(source.nodes[1]),raw(source.nodes[2]),raw(source.nodes[3])},int(n),raw(error)});
    mesh_check(comm,error,"invalid source coordinates or connectivity");
    launch(int(e),MeshCellCheck{{raw(source.nodes[0]),raw(source.nodes[1]),raw(source.nodes[2]),raw(source.nodes[3])},int(n),raw(error)});
    mesh_check(comm,error,"degenerate source connectivity");
    // Side-set names are file metadata. All element/face matching is performed below on the device.
    int sets=int(file.side_sets.size()); MPI_Bcast(&sets,1,MPI_INT,0,comm);
    simple_collective(comm,sets==3,"SIMPLE requires inlet, outlet and walls side sets");
    Buffer<int> tags(4*e,-1); int seen=0;
    for (int i=0;i<sets;++i) {
        int kind=-1;
        if (!rank) {
            const auto& name=file.side_sets[size_t(i)].name;
            kind=name=="inlet"?0:name=="outlet"?1:name=="walls"?2:-1;
        }
        MPI_Bcast(&kind,1,MPI_INT,0,comm);
        simple_collective(comm,kind>=0 && !(seen&(1<<kind)),"missing, repeated or unsupported boundary name"); seen|=1<<kind;
        Buffer<long long> elements,sides; const std::vector<long long> empty;
        broadcast_mesh_array(comm,rank?empty:file.side_sets[size_t(i)].elements,elements,MPI_LONG_LONG);
        broadcast_mesh_array(comm,rank?empty:file.side_sets[size_t(i)].sides,sides,MPI_LONG_LONG);
        simple_collective(comm,elements.size()==sides.size() && !elements.empty() && source.faces.size()+elements.size()<=size_t(INT_MAX/3),"invalid side-set lengths");
        const int first=int(source.faces.size()); source.faces.resize(source.faces.size()+elements.size());
        launch(int(elements.size()),SourceBoundary{raw(elements),raw(sides),int(e),kind,first,raw(source.faces),raw(tags),raw(error)});
        mesh_check(comm,error,"invalid or repeated boundary face");
    }
    auto ids=mesh_sequence(0,int(n)); Buffer<MeshFaceRecord<int>> faces(4*e);
    launch(int(4*e),MeshFaces<int>{{raw(source.nodes[0]),raw(source.nodes[1]),raw(source.nodes[2]),raw(source.nodes[3])},raw(ids),raw(faces)});
    mesh_sort(faces,MeshFaceLess<int>{});
    launch(int(4*e),SourceExteriorCheck{raw(faces),raw(tags),int(4*e),raw(error)});
    mesh_check(comm,error,"source side sets do not cover the manifold exterior exactly once");
    return source;
}

#ifdef MARS_REPLAY_CUDA
using NativeSimpleDomain=ElementDomain<TetTag,double,uint64_t,cstone::GpuTag>;
inline std::unique_ptr<NativeSimpleDomain> distribute_simple_mesh(MPI_Comm comm,const SimpleDeviceInput& source) {
    int rank,ranks; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    // ElementDomain currently uses MPI_COMM_WORLD internally.
    int relation; MPI_Comm_compare(comm,MPI_COMM_WORLD,&relation);
    simple_collective(comm,relation==MPI_IDENT || relation==MPI_CONGRUENT,"ElementDomain requires the world communicator");
    const char* mode=std::getenv("MARS_OWNERSHIP");
    simple_collective(comm,!mode || std::string(mode)!="vote","native SIMPLE requires SFC node ownership");
    NativeSimpleDomain::DeviceCoordsTuple coordinates;
    auto copy=[](auto& out,const auto& in) {
        out.resize(in.size()); thrust::copy(in.begin(),in.end(),thrust::device_pointer_cast(out.data()));
    };
    copy(std::get<0>(coordinates),source.x); copy(std::get<1>(coordinates),source.y); copy(std::get<2>(coordinates),source.z);
    NativeSimpleDomain::DeviceConnectivityTuple connectivity;
    const auto e=static_cast<long long>(source.nodes[0].size());
    const size_t first=size_t(e*rank/ranks),last=size_t(e*(rank+1)/ranks);
    auto slice=[&](auto& out,const auto& in) {
        out.resize(last-first); thrust::copy(in.begin()+first,in.begin()+last,thrust::device_pointer_cast(out.data()));
    };
    slice(std::get<0>(connectivity),source.nodes[0]); slice(std::get<1>(connectivity),source.nodes[1]);
    slice(std::get<2>(connectivity),source.nodes[2]); slice(std::get<3>(connectivity),source.nodes[3]);
    return std::make_unique<NativeSimpleDomain>(std::move(coordinates),std::move(connectivity),rank,ranks);
}
struct InvalidNativeNode {
    size_t count;
    template<class T> __host__ __device__ bool operator()(T node) const {
        return static_cast<unsigned long long>(node)>=count;
    }
};
struct SourceKey {
    uint64_t key; int source;
    __host__ __device__ bool operator<(const SourceKey& b) const { return key<b.key; }
};
struct EncodeSourceKeys {
    const double *x,*y,*z; cstone::Box<double> box; SourceKey* records; uint64_t* keys;
    __device__ void operator()(int i) const {
        const uint64_t k=cstone::sfc3D<cstone::HilbertKey<uint64_t>>(x[i],y[i],z[i],box).value(); records[i]={k,i}; keys[i]=k;
    }
};
struct SourceKeyLess { __host__ __device__ bool operator()(SourceKey a,SourceKey b) const { return a<b; } };
struct UniqueSourceKeys {
    const SourceKey* keys; int* error;
    __device__ void operator()(int i) const { if (i && keys[i-1].key==keys[i].key) atomicExch(error,1); }
};
struct RestoreCoordinates {
    const SourceKey* sorted; int count; const uint64_t* keys;
    const double *sx,*sy,*sz; double *x,*y,*z; int *source,*error;
    __device__ void operator()(int i) const {
        const uint64_t key=keys[i]; int lo=0,hi=count;
        while (lo<hi) { const int mid=lo+(hi-lo)/2; if (sorted[mid].key<key) lo=mid+1; else hi=mid; }
        if (lo==count || sorted[lo].key!=key) { atomicExch(error,1); return; }
        const int g=sorted[lo].source; source[i]=g; x[i]=sx[g]; y[i]=sy[g]; z[i]=sz[g];
    }
};
struct SourceTag { SimpleFaceKey<uint64_t> key; int kind,source; };
struct SourceTagLess { __host__ __device__ bool operator()(const SourceTag& a,const SourceTag& b) const { return a.key<b.key; } };
struct EncodeSourceTags {
    const SimpleFace* faces; const int* nodes[4]; const uint64_t* keys; SourceTag* tags;
    __device__ void operator()(int i) const {
        const auto f=faces[i]; tags[i]={mesh_face_key(keys[nodes[tet_face_node(f.ordinal,0)][f.element]],
            keys[nodes[tet_face_node(f.ordinal,1)][f.element]],keys[nodes[tet_face_node(f.ordinal,2)][f.element]]),f.kind,i};
    }
};
struct NativeBoundaryKind {
    const uint64_t* keys; const SourceTag* tags; int count;
    __host__ __device__ int find(const int* face) const {
        const auto key=mesh_face_key(keys[face[0]],keys[face[1]],keys[face[2]]); int lo=0,hi=count;
        while (lo<hi) { const int mid=lo+(hi-lo)/2; if (tags[mid].key<key) lo=mid+1; else hi=mid; }
        return lo<count && tags[lo].key==key?lo:-1;
    }
    __host__ __device__ int operator()(const int* face) const { const int i=find(face); return i<0?-1:tags[i].kind; }
};
struct NativeCoverage {
    const int *owned_nodes,*source; const SimpleFace* faces; const int* owned_faces; const int* nodes[4];
    int node_work,face_work,source_nodes; NativeBoundaryKind lookup; int* counts;
    __device__ void operator()(int i) const {
        if (i<node_work) atomicAdd(counts+source[owned_nodes[i]],1);
        if (i<face_work) {
            const auto f=faces[owned_faces[i]]; int local[3]; for (int j=0;j<3;++j) local[j]=nodes[tet_face_node(f.ordinal,j)][f.element];
            const int match=lookup.find(local);
            if (match>=0) atomicAdd(counts+source_nodes+lookup.tags[match].source,1);
        }
    }
};
struct CoverageCheck { const int* counts; int* error; __device__ void operator()(int i) const { if (counts[i]!=1) atomicExch(error,1); } };

template<class GlobalId> struct NativeSimpleMesh {
    NativeSimplePartition<GlobalId> partition;
    Buffer<int> source_node;
    NativeSimpleMesh(MPI_Comm comm,const NativeSimpleDomain& domain,const SimpleDeviceInput& source) {
        // Force the lazy local map first: the constructor's input node count is not rank-local.
        simple_collective(comm,domain.numRanks()==1 || domain.sfcOwnership(),"native SIMPLE requires SFC node ownership");
        const auto& keys=domain.getLocalToGlobalSfcMap();
        const auto& conn=domain.getElementToNodeConnectivity();
        const size_t n=keys.size(),e=domain.getElementCount();
        simple_collective(comm,n<=size_t(INT_MAX/9) && e<=size_t(INT_MAX/16),"native local mesh exceeds index capacity");
        SimpleMeshData<uint64_t> v; Buffer<int> error(1,0);
        auto copy=[&](auto& out,const auto& in,size_t count) {
            simple_collective(comm,in.size()==count,"native local array size mismatch");
            out.resize(count); if (count) thrust::copy(thrust::device_pointer_cast(in.data()),thrust::device_pointer_cast(in.data())+count,out.begin());
        };
        auto valid_column=[&](const auto& column) {
            return !thrust::any_of(thrust::device_pointer_cast(column.data()),
                                  thrust::device_pointer_cast(column.data())+column.size(),InvalidNativeNode{n});
        };
        simple_collective(comm,valid_column(std::get<0>(conn)) && valid_column(std::get<1>(conn))
                          && valid_column(std::get<2>(conn)) && valid_column(std::get<3>(conn)),"invalid native node index before narrowing");
        copy(v.key,keys,n); copy(v.nodes[0],std::get<0>(conn),e); copy(v.nodes[1],std::get<1>(conn),e);
        copy(v.nodes[2],std::get<2>(conn),e); copy(v.nodes[3],std::get<3>(conn),e);
        copy(v.owned,domain.getNodeOwnershipMap(),n);
        v.element_begin=int(domain.startIndex()); v.element_end=int(domain.endIndex());
        if (domain.numRanks()>1) {
            const auto& h=domain.getNodeHaloTopology(); v.peers=h.peers_; v.send_offsets=h.sendOffsets_; v.recv_offsets=h.recvOffsets_;
            copy(v.send_nodes,h.sendNodeIds_,h.sendNodeIds_.size()); copy(v.recv_nodes,h.recvNodeIds_,h.recvNodeIds_.size());
        }
        v.x.resize(n); v.y.resize(n); v.z.resize(n); source_node.resize(n);
        Buffer<SourceKey> sorted(source.x.size()); Buffer<uint64_t> source_keys(source.x.size());
        launch(int(source.x.size()),EncodeSourceKeys{raw(source.x),raw(source.y),raw(source.z),domain.getBoundingBox(),raw(sorted),raw(source_keys)});
        mesh_sort(sorted,SourceKeyLess{});
        launch(int(sorted.size()),UniqueSourceKeys{raw(sorted),raw(error)});
        mesh_check(comm,error,"source nodes collide under native SFC identity");
        launch(int(n),RestoreCoordinates{raw(sorted),int(sorted.size()),raw(v.key),raw(source.x),raw(source.y),raw(source.z),raw(v.x),raw(v.y),raw(v.z),raw(source_node),raw(error)});
        mesh_check(comm,error,"a local node has no source coordinate");
        Buffer<SourceTag> tags(source.faces.size());
        launch(int(tags.size()),EncodeSourceTags{raw(source.faces),{raw(source.nodes[0]),raw(source.nodes[1]),raw(source.nodes[2]),raw(source.nodes[3])},raw(source_keys),raw(tags)});
        mesh_sort(tags,SourceTagLess{});
        NativeBoundaryKind lookup{raw(v.key),raw(tags),int(tags.size())};
        partition=build_simple_partition<GlobalId>(comm,v,lookup);
        auto& o=partition.ownership; const auto& f=partition.input;
        // Identity coverage is a public validation gate; its counters and collective stay on-device.
        simple_collective(comm,source.x.size()+source.faces.size()<=size_t(INT_MAX),"coverage exceeds MPI count capacity");
        Buffer<int> coverage(source.x.size()+source.faces.size(),0);
        launch(int(std::max(o.owned_nodes.size(),o.owned_faces.size())),NativeCoverage{raw(o.owned_nodes),raw(source_node),raw(f.faces),raw(o.owned_faces),
            {raw(f.nodes[0]),raw(f.nodes[1]),raw(f.nodes[2]),raw(f.nodes[3])},int(o.owned_nodes.size()),int(o.owned_faces.size()),int(source.x.size()),lookup,raw(coverage)});
        mesh_mpi_ready();
        ensure(MPI_Allreduce(MPI_IN_PLACE,raw(coverage),int(coverage.size()),MPI_INT,MPI_SUM,comm)==MPI_SUCCESS,"native coverage reduction failed");
        launch(int(coverage.size()),CoverageCheck{raw(coverage),raw(error)});
        mesh_check(comm,error,"native ownership must cover every source node and boundary face once");
    }
};
#endif
} // namespace mars::segregated::runtime
#undef MARS_SINPUT_HD
