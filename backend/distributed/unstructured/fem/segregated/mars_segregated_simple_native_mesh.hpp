#pragma once
#include "mars_segregated_simple_options.hpp"
#include "mars_segregated_simple_mesh.hpp"
#include "../../utils/mars_read_exodus_raw.hpp"
#include <map>
#include <type_traits>
#ifdef MARS_REPLAY_CUDA
#include <thrust/unique.h>
#endif
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

inline std::vector<int> source_boundary_kinds(const std::vector<ExodusRawSideSet>& sets,const SimpleBoundaryNames& names) {
    ensure(names.valid() && sets.size()==names.walls.size()+2,
           "SIMPLE boundary selection must name every side set exactly once");
    std::map<std::string,size_t> aliases;
    auto add=[&](const std::string& name,size_t index) {
        if (name.empty()) return;
        const auto entry=aliases.emplace(simple_boundary_name(name),index);
        ensure(entry.second || entry.first->second==index,"ambiguous Exodus boundary names or aliases");
    };
    for (size_t i=0;i<sets.size();++i) {
        const auto& set=sets[i];
        auto name=simple_boundary_name(set.name);
        if (set.id>0) {
            const auto id=std::to_string(set.id);
            // Ioss replaces a generated surface name whose embedded ID is stale.
            if (name.rfind("surface_",0)==0) {
                const auto suffix=name.substr(8);
                if (!suffix.empty() && suffix[0]>='1' && suffix[0]<='9' &&
                    suffix.find_first_not_of("0123456789")==std::string::npos) name="surface_"+id;
            }
            add("surface_"+id,i); add("sideset_"+id,i);
        }
        add(name,i);
    }
    std::vector<int> kinds(sets.size(),-1);
    auto select=[&](const std::string& name,int kind) {
        const auto found=aliases.find(simple_boundary_name(name));
        ensure(found!=aliases.end(),"selected boundary name has no Exodus name or ID alias match");
        ensure(kinds[found->second]<0,"multiple boundary selections refer to the same Exodus side set");
        kinds[found->second]=kind;
    };
    select(names.inlet,0); select(names.outlet,1);
    for (const auto& wall:names.walls) select(wall,2);
    return kinds;
}

// Rank zero reads bytes. Broadcasts, index conversion and all topology validation use device buffers.
inline SimpleDeviceInput read_simple_mesh(MPI_Comm comm,const std::string& path,const SimpleBoundaryNames& names={}) {
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
    std::vector<int> kinds;
    if (!rank) try { kinds=source_boundary_kinds(file.side_sets,names); }
    catch (const std::exception& ex) { ok=false; std::cerr<<ex.what()<<'\n'; }
    simple_collective(comm,ok,"native boundary name resolution failed");
    int sets=int(file.side_sets.size()); MPI_Bcast(&sets,1,MPI_INT,0,comm);
    Buffer<int> tags(4*e,-1);
    for (int i=0;i<sets;++i) {
        int kind=rank?-1:kinds[size_t(i)];
        MPI_Bcast(&kind,1,MPI_INT,0,comm);
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

struct NativeSimpleSource : SimpleDeviceInput {
    int global_nodes=0,global_elements=0,global_faces=0;
};
inline NativeSimpleSource read_simple_mesh_root(MPI_Comm comm,const std::string& path,const SimpleBoundaryNames& names={}) {
    int rank; MPI_Comm_rank(comm,&rank);
    NativeSimpleSource source; bool ok=true;
    if (!rank) try {
        static_cast<SimpleDeviceInput&>(source)=read_simple_mesh(MPI_COMM_SELF,path,names);
        source.global_nodes=int(source.x.size()); source.global_elements=int(source.nodes[0].size());
        source.global_faces=int(source.faces.size());
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; ok=false; }
    simple_collective(comm,ok,"root native mesh read or validation failed");
    int counts[3]={source.global_nodes,source.global_elements,source.global_faces};
    ensure(MPI_Bcast(counts,3,MPI_INT,0,comm)==MPI_SUCCESS,"source count broadcast failed");
    source.global_nodes=counts[0]; source.global_elements=counts[1]; source.global_faces=counts[2];
    return source;
}

// Payloads remain on the device. Only counts and fault words reach the host.
// A bounded root service avoids replicating the source or a global routing directory.
inline constexpr int native_mesh_chunk=65536;
template<class Request,class Response,class Kernel>
void native_root_lookup(MPI_Comm comm,const Buffer<Request>& requests,Buffer<Response>& responses,Kernel kernel) try {
    static_assert(std::is_trivially_copyable_v<Request> && std::is_trivially_copyable_v<Response>);
    static_assert(sizeof(Request)<=INT_MAX/native_mesh_chunk && sizeof(Response)<=INT_MAX/native_mesh_chunk);
    int rank,ranks; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    simple_collective(comm,requests.size()<=size_t(INT_MAX),"native lookup count overflow");
    responses.resize(requests.size());
    Buffer<Request> incoming; Buffer<Response> outgoing;
    if (!rank) { incoming.resize(native_mesh_chunk); outgoing.resize(native_mesh_chunk); }
    for (int peer=0;peer<ranks;++peer) {
        int count=rank==peer?int(requests.size()):0;
        ensure(MPI_Bcast(&count,1,MPI_INT,peer,comm)==MPI_SUCCESS,"native lookup count broadcast failed");
        for (int first=0;first<count;) {
            const int work=std::min(native_mesh_chunk,count-first);
            if (rank==peer || !rank) mesh_mpi_ready();
            if (peer && rank==peer)
                ensure(MPI_Send(raw(requests)+first,work*sizeof(Request),MPI_BYTE,0,27101,comm)==MPI_SUCCESS,"native lookup request failed");
            if (!rank) {
                const Request* request=nullptr;
                Response* response=nullptr;
                if (peer) {
                    ensure(MPI_Recv(raw(incoming),work*sizeof(Request),MPI_BYTE,peer,27101,comm,MPI_STATUS_IGNORE)==MPI_SUCCESS,"native lookup receive failed");
                    request=raw(incoming); response=raw(outgoing);
                } else { request=raw(requests)+first; response=raw(responses)+first; }
                launch(work,kernel(request,response)); mesh_mpi_ready();
                if (peer) ensure(MPI_Send(response,work*sizeof(Response),MPI_BYTE,peer,27102,comm)==MPI_SUCCESS,"native lookup reply failed");
            }
            if (peer && rank==peer)
                ensure(MPI_Recv(raw(responses)+first,work*sizeof(Response),MPI_BYTE,0,27102,comm,MPI_STATUS_IGNORE)==MPI_SUCCESS,"native lookup reply receive failed");
            first+=work;
        }
    }
} catch (const std::exception&) {
    // A failed allocation or device launch cannot leave a peer waiting for its reply.
    MPI_Abort(comm,1); throw;
}
struct NativeCell {
    double x[4],y[4],z[4]; int source[4];
};
struct NativeCellLookup {
    const int* requests; NativeCell* result; const int* nodes[4]; const double *x,*y,*z;
    MARS_SINPUT_HD void operator()(int i) const {
        for (int j=0;j<4;++j) { const int n=nodes[j][requests[i]]; result[i].x[j]=x[n]; result[i].y[j]=y[n]; result[i].z[j]=z[n]; result[i].source[j]=n; }
    }
};
struct NativeVertex { int source; double x,y,z; };
struct NativeVertexLess { MARS_SINPUT_HD bool operator()(NativeVertex a,NativeVertex b) const { return a.source<b.source; } };
struct NativeVertexEqual { MARS_SINPUT_HD bool operator()(NativeVertex a,NativeVertex b) const { return a.source==b.source; } };
struct NativeCellVertices {
    const NativeCell* cells; NativeVertex* vertices;
    MARS_SINPUT_HD void operator()(int i) const {
        const auto& cell=cells[i/4]; const int j=i%4;
        vertices[i]={cell.source[j],cell.x[j],cell.y[j],cell.z[j]};
    }
};
struct NativeVertexUnpack {
    const NativeVertex* vertices; double *x,*y,*z;
    MARS_SINPUT_HD void operator()(int i) const { const auto v=vertices[i]; x[i]=v.x; y[i]=v.y; z[i]=v.z; }
};
struct NativeCellUnpack {
    const NativeCell* cells; const NativeVertex* vertices; int count; int* nodes[4];
    MARS_SINPUT_HD void operator()(int e) const {
        for (int j=0;j<4;++j) {
            const int source=cells[e].source[j]; int lo=0,hi=count;
            while (lo<hi) { const int mid=lo+(hi-lo)/2; if (vertices[mid].source<source) lo=mid+1; else hi=mid; }
            nodes[j][e]=lo;
        }
    }
};
inline SimpleDeviceInput native_initial_partition(MPI_Comm comm,const NativeSimpleSource& source) {
    int rank,ranks; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    simple_collective(comm,source.global_elements>=ranks,"native mesh needs at least one source element per rank");
    const int first=int(static_cast<long long>(source.global_elements)*rank/ranks);
    const int last=int(static_cast<long long>(source.global_elements)*(rank+1)/ranks),count=last-first;
    auto ids=mesh_sequence(first,last); Buffer<NativeCell> cells;
    native_root_lookup(comm,ids,cells,[&](const int* requests,NativeCell* response) {
        return NativeCellLookup{requests,response,{raw(source.nodes[0]),raw(source.nodes[1]),raw(source.nodes[2]),raw(source.nodes[3])},
                                raw(source.x),raw(source.y),raw(source.z)};
    });
    // Shared source IDs must remain shared locally: ElementDomain averages incident cell sizes per node.
    Buffer<NativeVertex> vertices(4*count);
    launch(4*count,NativeCellVertices{raw(cells),raw(vertices)}); mesh_sort(vertices,NativeVertexLess{});
#ifdef MARS_REPLAY_CUDA
    vertices.resize(size_t(thrust::unique(vertices.begin(),vertices.end(),NativeVertexEqual{})-vertices.begin()));
#else
    vertices.resize(size_t(std::unique(vertices.begin(),vertices.end(),NativeVertexEqual{})-vertices.begin()));
#endif
    SimpleDeviceInput local; local.x.resize(vertices.size()); local.y.resize(vertices.size()); local.z.resize(vertices.size());
    for (auto& column:local.nodes) column.resize(count);
    launch(int(vertices.size()),NativeVertexUnpack{raw(vertices),raw(local.x),raw(local.y),raw(local.z)});
    launch(count,NativeCellUnpack{raw(cells),raw(vertices),int(vertices.size()),
        {raw(local.nodes[0]),raw(local.nodes[1]),raw(local.nodes[2]),raw(local.nodes[3])}});
    return local;
}

#ifdef MARS_REPLAY_CUDA
using NativeSimpleDomain=ElementDomain<TetTag,double,uint64_t,cstone::execution::Gpu>;
inline std::unique_ptr<NativeSimpleDomain> distribute_simple_mesh(MPI_Comm comm,const NativeSimpleSource& source) {
    int rank,ranks; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
    // ElementDomain currently uses MPI_COMM_WORLD internally.
    int relation; MPI_Comm_compare(comm,MPI_COMM_WORLD,&relation);
    simple_collective(comm,relation==MPI_IDENT || relation==MPI_CONGRUENT,"ElementDomain requires the world communicator");
    const char* mode=std::getenv("MARS_OWNERSHIP");
    simple_collective(comm,!mode || std::string(mode)!="vote","native SIMPLE requires SFC node ownership");
    mars::requestFullNodeHaloExchange(); // not yet checked for halo reads beyond one layer
    auto local=native_initial_partition(comm,source);
    NativeSimpleDomain::DeviceCoordsTuple coordinates;
    auto copy=[](auto& out,const auto& in) {
        out.resize(in.size()); thrust::copy(in.begin(),in.end(),thrust::device_pointer_cast(out.data()));
    };
    copy(std::get<0>(coordinates),local.x); copy(std::get<1>(coordinates),local.y); copy(std::get<2>(coordinates),local.z);
    NativeSimpleDomain::DeviceConnectivityTuple connectivity;
    copy(std::get<0>(connectivity),local.nodes[0]); copy(std::get<1>(connectivity),local.nodes[1]);
    copy(std::get<2>(connectivity),local.nodes[2]); copy(std::get<3>(connectivity),local.nodes[3]);
    return std::make_unique<NativeSimpleDomain>(std::move(coordinates),std::move(connectivity),rank,ranks);
}
#endif
MARS_SINPUT_HD inline void native_count(int* value) {
#if defined(__CUDA_ARCH__)
    atomicAdd(value,1);
#else
    ++*value;
#endif
}
struct InvalidNativeNode {
    size_t count;
    template<class T> MARS_SINPUT_HD bool operator()(T node) const {
        return static_cast<unsigned long long>(node)>=count;
    }
};
struct SourceKey {
    uint64_t key; int source;
    MARS_SINPUT_HD bool operator<(const SourceKey& b) const { return key<b.key; }
};
#ifdef MARS_REPLAY_CUDA
struct EncodeSourceKeys {
    const double *x,*y,*z; cstone::Box<double> box; SourceKey* records; uint64_t* keys;
    MARS_SINPUT_HD void operator()(int i) const {
        const uint64_t k=cstone::sfc3D<cstone::HilbertKey<uint64_t>>(x[i],y[i],z[i],box).value(); records[i]={k,i}; keys[i]=k;
    }
};
#endif
struct SourceKeyLess { MARS_SINPUT_HD bool operator()(SourceKey a,SourceKey b) const { return a<b; } };
struct UniqueSourceKeys {
    const SourceKey* keys; int* error;
    MARS_SINPUT_HD void operator()(int i) const { if (i && keys[i-1].key==keys[i].key) distributed::raise_fault(error,1); }
};
struct NativeCoordinate { double x,y,z; int source; };
struct NativeCoordinateLookup {
    const uint64_t* requests; NativeCoordinate* response; const SourceKey* sorted; int count;
    const double *x,*y,*z; int* error;
    MARS_SINPUT_HD void operator()(int i) const {
        const uint64_t key=requests[i]; int lo=0,hi=count;
        while (lo<hi) { const int mid=lo+(hi-lo)/2; if (sorted[mid].key<key) lo=mid+1; else hi=mid; }
        if (lo==count || sorted[lo].key!=key) { distributed::raise_fault(error,1); response[i]={0,0,0,-1}; return; }
        const int g=sorted[lo].source; response[i]={x[g],y[g],z[g],g};
    }
};
struct NativeCoordinateUnpack {
    const NativeCoordinate* coordinates; double *x,*y,*z; int* source;
    MARS_SINPUT_HD void operator()(int i) const { const auto c=coordinates[i]; x[i]=c.x; y[i]=c.y; z[i]=c.z; source[i]=c.source; }
};
struct SourceTag { SimpleFaceKey<uint64_t> key; int kind,source; };
struct SourceTagLess { MARS_SINPUT_HD bool operator()(const SourceTag& a,const SourceTag& b) const { return a.key<b.key; } };
struct EncodeSourceTags {
    const SimpleFace* faces; const int* nodes[4]; const uint64_t* keys; SourceTag* tags;
    MARS_SINPUT_HD void operator()(int i) const {
        const auto f=faces[i]; tags[i]={mesh_face_key(keys[nodes[tet_face_node(f.ordinal,0)][f.element]],
            keys[nodes[tet_face_node(f.ordinal,1)][f.element]],keys[nodes[tet_face_node(f.ordinal,2)][f.element]]),f.kind,i};
    }
};
struct NativeBoundaryKind {
    const uint64_t* keys; const SourceTag* tags; int count;
    MARS_SINPUT_HD int find(const int* face) const {
        const auto key=mesh_face_key(keys[face[0]],keys[face[1]],keys[face[2]]); int lo=0,hi=count;
        while (lo<hi) { const int mid=lo+(hi-lo)/2; if (tags[mid].key<key) lo=mid+1; else hi=mid; }
        return lo<count && tags[lo].key==key?lo:-1;
    }
    MARS_SINPUT_HD int operator()(const int* face) const { const int i=find(face); return i<0?-1:tags[i].kind; }
};
struct NativeFaceRequest {
    const int* nodes[4]; const uint64_t* keys; SimpleFaceKey<uint64_t>* requests;
    MARS_SINPUT_HD void operator()(int i) const {
        const int e=i/4,f=i%4;
        requests[i]=mesh_face_key(keys[nodes[tet_face_node(f,0)][e]],keys[nodes[tet_face_node(f,1)][e]],keys[nodes[tet_face_node(f,2)][e]]);
    }
};
struct NativeTagLookup {
    const SimpleFaceKey<uint64_t>* requests; SourceTag* response; const SourceTag* tags; int count;
    MARS_SINPUT_HD void operator()(int i) const {
        const auto key=requests[i]; int lo=0,hi=count;
        while (lo<hi) { const int mid=lo+(hi-lo)/2; if (tags[mid].key<key) lo=mid+1; else hi=mid; }
        response[i]=lo<count && tags[lo].key==key?tags[lo]:SourceTag{key,-1,-1};
    }
};
struct NativeTagged { MARS_SINPUT_HD bool operator()(SourceTag tag) const { return tag.source>=0; } };
struct NativeSameTag { MARS_SINPUT_HD bool operator()(SourceTag a,SourceTag b) const { return a.key==b.key; } };
struct NativeCoverageIds {
    const int *owned_nodes,*source; const SimpleFace* faces; const int* owned_faces; const int* nodes[4];
    int node_work,source_nodes; NativeBoundaryKind lookup; int* ids;
    MARS_SINPUT_HD void operator()(int i) const {
        if (i<node_work) ids[i]=source[owned_nodes[i]];
        else {
            const auto f=faces[owned_faces[i-node_work]]; int local[3];
            for (int j=0;j<3;++j) local[j]=nodes[tet_face_node(f.ordinal,j)][f.element];
            const int match=lookup.find(local); ids[i]=match<0?-1:source_nodes+lookup.tags[match].source;
        }
    }
};
struct NativeCoverageLookup {
    const int* ids; unsigned char* response; int* counts; int size; int* error;
    MARS_SINPUT_HD void operator()(int i) const {
        const int id=ids[i]; response[i]=0;
        if (id<0 || id>=size) distributed::raise_fault(error,1); else native_count(counts+id);
    }
};
struct CoverageCheck { const int* counts; int* error; MARS_SINPUT_HD void operator()(int i) const { if (counts[i]!=1) distributed::raise_fault(error,1); } };

#ifdef MARS_REPLAY_CUDA
template<class GlobalId> struct NativeSimpleMesh {
    NativeSimplePartition<GlobalId> partition;
    Buffer<int> source_node;
    NativeSimpleMesh(MPI_Comm comm,const NativeSimpleDomain& domain,const NativeSimpleSource& source) {
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
        int rank; MPI_Comm_rank(comm,&rank);
        v.x.resize(n); v.y.resize(n); v.z.resize(n); source_node.resize(n);
        Buffer<SourceKey> sorted(source.x.size()); Buffer<uint64_t> source_keys(source.x.size());
        if (!rank) {
            launch(int(source.x.size()),EncodeSourceKeys{raw(source.x),raw(source.y),raw(source.z),domain.getBoundingBox(),raw(sorted),raw(source_keys)});
            mesh_sort(sorted,SourceKeyLess{});
            launch(int(sorted.size()),UniqueSourceKeys{raw(sorted),raw(error)});
        }
        mesh_check(comm,error,"source nodes collide under native SFC identity");
        Buffer<NativeCoordinate> coordinates;
        native_root_lookup(comm,v.key,coordinates,[&](const uint64_t* requests,NativeCoordinate* response) {
            return NativeCoordinateLookup{requests,response,raw(sorted),int(sorted.size()),raw(source.x),raw(source.y),raw(source.z),raw(error)};
        });
        mesh_check(comm,error,"a local node has no source coordinate");
        launch(int(n),NativeCoordinateUnpack{raw(coordinates),raw(v.x),raw(v.y),raw(v.z),raw(source_node)});
        Buffer<SourceTag> source_tags(source.faces.size());
        if (!rank) {
            launch(int(source_tags.size()),EncodeSourceTags{raw(source.faces),{raw(source.nodes[0]),raw(source.nodes[1]),raw(source.nodes[2]),raw(source.nodes[3])},raw(source_keys),raw(source_tags)});
            mesh_sort(source_tags,SourceTagLess{});
        }
        Buffer<SimpleFaceKey<uint64_t>> face_keys(4*e);
        launch(int(4*e),NativeFaceRequest{{raw(v.nodes[0]),raw(v.nodes[1]),raw(v.nodes[2]),raw(v.nodes[3])},raw(v.key),raw(face_keys)});
        Buffer<SourceTag> tags;
        native_root_lookup(comm,face_keys,tags,[&](const SimpleFaceKey<uint64_t>* requests,SourceTag* response) {
            return NativeTagLookup{requests,response,raw(source_tags),int(source_tags.size())};
        });
        mesh_compact(tags,NativeTagged{}); mesh_sort(tags,SourceTagLess{});
        tags.resize(size_t(thrust::unique(tags.begin(),tags.end(),NativeSameTag{})-tags.begin()));
        NativeBoundaryKind lookup{raw(v.key),raw(tags),int(tags.size())};
        partition=build_simple_partition<GlobalId>(comm,v,lookup);
        auto& o=partition.ownership; const auto& f=partition.input;
        simple_collective(comm,static_cast<long long>(source.global_nodes)+source.global_faces<=INT_MAX,"coverage exceeds index capacity");
        const int total=source.global_nodes+source.global_faces;
        Buffer<int> coverage(rank?0:size_t(total),0),coverage_ids(o.owned_nodes.size()+o.owned_faces.size());
        launch(int(coverage_ids.size()),NativeCoverageIds{raw(o.owned_nodes),raw(source_node),raw(f.faces),raw(o.owned_faces),
            {raw(f.nodes[0]),raw(f.nodes[1]),raw(f.nodes[2]),raw(f.nodes[3])},int(o.owned_nodes.size()),source.global_nodes,lookup,raw(coverage_ids)});
        Buffer<unsigned char> acknowledgements;
        native_root_lookup(comm,coverage_ids,acknowledgements,[&](const int* requests,unsigned char* response) {
            return NativeCoverageLookup{requests,response,raw(coverage),total,raw(error)};
        });
        if (!rank) launch(total,CoverageCheck{raw(coverage),raw(error)});
        mesh_check(comm,error,"native ownership must cover every source node and boundary face once");
    }
};
#endif
} // namespace mars::segregated::runtime
#undef MARS_SINPUT_HD
