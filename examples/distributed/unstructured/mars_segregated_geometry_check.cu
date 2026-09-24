#include "mars_segregated_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef MARS_REPLAY_CUDA
#include "mars_segregated_geometry_device.hpp"
#define GEOMETRY_CHECK_HD __host__ __device__
#else
#define GEOMETRY_CHECK_HD
#endif
using namespace mars::segregated;
using Geometry = TetGeometry<double>;

void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
template<class T> void read(std::istream& stream,T& value) {
    require(bool(stream >> value) && std::isfinite(double(value)),"invalid/truncated geometry data");
}
template<class T> void read_all(std::istream& stream,T& values) { for (auto& v:values) read(stream,v); }
struct Frame { int components,call; std::vector<double> field,expected; };
struct Fixture {
    std::vector<double> x,y,z,element_expected,face_expected;
    std::array<std::vector<int>,4> nodes;
    std::vector<TetBoundaryFace> faces;
    std::vector<Frame> frames;
};
Fixture load(const char* path) {
    std::ifstream stream(path); std::string magic; int n,e,b,f;
    require(bool(stream >> magic >> n >> e >> b >> f) && magic=="MARS_PUBLIC_GEOMETRY_V1"
            && n==425 && e==1536 && b==576 && f==4,"expected complete public geometry fixture");
    Fixture data; data.x.resize(n); data.y.resize(n); data.z.resize(n);
    unsigned long long id,previous=0;
    for (int i=0;i<n;++i) {
        read(stream,id); require(id>previous,"unordered/duplicate node ID"); previous=id;
        read(stream,data.x[i]); read(stream,data.y[i]); read(stream,data.z[i]);
    }
    for (auto& nodes:data.nodes) nodes.resize(e);
    data.element_expected.resize(e*138); previous=0;
    for (int i=0;i<e;++i) {
        read(stream,id); require(id>previous,"unordered/duplicate element ID"); previous=id;
        for (auto& nodes:data.nodes) { read(stream,nodes[i]); require(nodes[i]>=0 && nodes[i]<n,"invalid node index"); }
        for (int a=0;a<4;++a) for (int c=0;c<a;++c) require(data.nodes[a][i]!=data.nodes[c][i],"repeated element node");
        for (int j=0;j<138;++j) read(stream,data.element_expected[138*i+j]);
    }
    data.faces.resize(b); data.face_expected.resize(b*18);
    std::vector<bool> visited(e*4,false);
    for (int i=0;i<b;++i) {
        auto& face=data.faces[i]; read(stream,face.element); read(stream,face.ordinal);
        require(face.element>=0 && face.element<e && face.ordinal>=0 && face.ordinal<4,"invalid boundary face");
        require(!visited[4*face.element+face.ordinal],"duplicate boundary face");
        visited[4*face.element+face.ordinal]=true;
        for (int j=0;j<18;++j) read(stream,data.face_expected[18*i+j]);
    }
    for (int i=0;i<f;++i) {
        Frame frame; read(stream,frame.components); read(stream,frame.call);
        require(frame.components==(i<2?3:1) && frame.call==i%2+1,"wrong gradient state order");
        frame.field.resize(n*frame.components); frame.expected.resize(n*frame.components*3);
        for (int node=0;node<n;++node) {
            for (int c=0;c<frame.components;++c) read(stream,frame.field[node*frame.components+c]);
            for (int c=0;c<3*frame.components;++c) read(stream,frame.expected[node*frame.components*3+c]);
        }
        data.frames.push_back(std::move(frame));
    }
    require(!(stream>>magic),"trailing geometry data"); return data;
}

GEOMETRY_CHECK_HD void element_values(const Geometry& g,double* out) {
    for (int s=0;s<6;++s) for (int j=0;j<12;++j) out[12*s+j]=g.gradient[j];
    for (int i=0;i<18;++i) out[72+i]=g.area[i];
    for (int s=0;s<6;++s) {
        tet_sample_shape(s,false,out+90+4*s);
        tet_sample_shape(s,false,out+114+4*s);
    }
}
GEOMETRY_CHECK_HD void face_values(const Geometry& g,int face,double* out) {
    for (int s=0;s<3;++s) {
        tet_boundary_area(g,face,out+3*s); tri_sample_shape(s,false,out+9+3*s);
    }
}
std::size_t checks=0;
void compare(const std::vector<double>& actual,const std::vector<double>& expected,const std::string& label) {
    require(actual.size()==expected.size(),"comparison size mismatch"); double worst=0; std::size_t failures=0;
    for (std::size_t i=0;i<actual.size();++i) {
        ++checks; const double error=std::abs(actual[i]-expected[i])/std::max(1.,std::abs(expected[i]));
        worst=std::max(worst,error);
        if (!std::isfinite(actual[i]) || error>1e-12) {
            if (failures++<4) std::cerr << label << " index=" << i << " expected=" << expected[i] << " actual=" << actual[i] << '\n';
        }
    }
    std::cout << label << " scalars=" << actual.size() << " worst_scaled=" << worst << '\n';
    require(failures==0,"native geometry/reconstruction mismatch");
}

#ifdef MARS_REPLAY_CUDA
void cuda_check(cudaError_t error) { if (error!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(error)); }
template<class T> struct DeviceBuffer {
    T* data=nullptr;
    explicit DeviceBuffer(std::size_t size) { cuda_check(cudaMalloc(reinterpret_cast<void**>(&data),size*sizeof(T))); }
    ~DeviceBuffer() { cudaFree(data); }
    DeviceBuffer(const DeviceBuffer&)=delete;
    DeviceBuffer& operator=(const DeviceBuffer&)=delete;
    void upload(const std::vector<T>& values) { cuda_check(cudaMemcpy(data,values.data(),values.size()*sizeof(T),cudaMemcpyHostToDevice)); }
    std::vector<T> download(std::size_t size) {
        std::vector<T> values(size); cuda_check(cudaMemcpy(values.data(),data,size*sizeof(T),cudaMemcpyDeviceToHost)); return values;
    }
};
__global__ void element_values_kernel(const Geometry* geometry,int count,double* values) {
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<count) element_values(geometry[i],values+138*i);
}
__global__ void face_values_kernel(const Geometry* geometry,const TetBoundaryFace* faces,int count,double* values) {
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<count) face_values(geometry[faces[i].element],faces[i].ordinal,values+18*i);
}
template<int Components>
void reconstruct(TetMeshView<int,double> mesh,const Geometry* geometry,const TetBoundaryFace* faces,int boundaries,
                 const double* field,const double* volumes,double* numerator,double* gradient,int* error) {
    cuda_check(cudaMemset(numerator,0,mesh.node_count*Components*3*sizeof(double)));
    accumulate_interior_gradient<Components><<<(mesh.element_count+127)/128,128>>>(
        mesh,geometry,field,0,mesh.element_count,Components==1,true,numerator);
    cuda_check(cudaGetLastError());
    accumulate_boundary_gradient<Components><<<(boundaries+127)/128,128>>>(
        mesh,geometry,faces,boundaries,field,Components==1,true,numerator,error);
    cuda_check(cudaGetLastError());
    normalize_gradient<Components><<<(mesh.node_count*Components*3+127)/128,128>>>(
        numerator,volumes,mesh.node_count,1.,false,gradient,error);
    cuda_check(cudaGetLastError());
}
void execute(const Fixture& data) {
    const int n=int(data.x.size()), e=int(data.nodes[0].size()), b=int(data.faces.size());
    DeviceBuffer<double> x(n),y(n),z(n),volumes(n),numerator(9*n),gradient(9*n),field(3*n),values(138*e);
    DeviceBuffer<int> n0(e),n1(e),n2(e),n3(e),error(1);
    DeviceBuffer<Geometry> geometry(e); DeviceBuffer<TetBoundaryFace> faces(b);
    x.upload(data.x); y.upload(data.y); z.upload(data.z); faces.upload(data.faces);
    n0.upload(data.nodes[0]); n1.upload(data.nodes[1]); n2.upload(data.nodes[2]); n3.upload(data.nodes[3]);
    TetMeshView<int,double> mesh{{n0.data,n1.data,n2.data,n3.data},x.data,y.data,z.data,n,e};
    cuda_check(cudaMemset(error.data,0,sizeof(int)));
    build_tet_geometry<<<(e+127)/128,128>>>(mesh,geometry.data,error.data); cuda_check(cudaGetLastError());
    require(error.download(1)[0]==0,"invalid native element geometry");
    cuda_check(cudaMemset(volumes.data,0,n*sizeof(double)));
    accumulate_dual_volumes<<<(e+127)/128,128>>>(mesh,geometry.data,0,e,volumes.data); cuda_check(cudaGetLastError());
    element_values_kernel<<<(e+127)/128,128>>>(geometry.data,e,values.data); cuda_check(cudaGetLastError());
    compare(values.download(138*e),data.element_expected,"element geometry");
    face_values_kernel<<<(b+127)/128,128>>>(geometry.data,faces.data,b,values.data); cuda_check(cudaGetLastError());
    compare(values.download(18*b),data.face_expected,"boundary geometry");
    for (const auto& frame:data.frames) {
        field.upload(frame.field);
        if (frame.components==1) reconstruct<1>(mesh,geometry.data,faces.data,b,field.data,volumes.data,numerator.data,gradient.data,error.data);
        else reconstruct<3>(mesh,geometry.data,faces.data,b,field.data,volumes.data,numerator.data,gradient.data,error.data);
        require(error.download(1)[0]==0,"invalid assembled gradient");
        compare(gradient.download(n*frame.components*3),frame.expected,
                std::string(frame.components==1?"pressure":"velocity")+" gradient call="+std::to_string(frame.call));
    }
}
#else
template<int Components>
std::vector<double> reconstruct(const Fixture& data,const std::vector<Geometry>& geometry,
                               const std::vector<double>& volume,const Frame& frame) {
    std::vector<double> sum(data.x.size()*Components*3,0),result(sum.size());
    for (std::size_t e=0;e<geometry.size();++e) {
        double values[4*Components],local[4*Components*3];
        for (int n=0;n<4;++n) for (int c=0;c<Components;++c) values[n*Components+c]=frame.field[data.nodes[n][e]*Components+c];
        tet_gradient_numerator<Components>(geometry[e],values,Components==1,true,local);
        for (int n=0;n<4;++n) for (int c=0;c<Components*3;++c) sum[data.nodes[n][e]*Components*3+c]+=local[n*Components*3+c];
    }
    for (const auto& face:data.faces) {
        double values[3*Components],local[3*Components*3],area[3]; tet_boundary_area(geometry[face.element],face.ordinal,area);
        for (int n=0;n<3;++n) for (int c=0;c<Components;++c)
            values[n*Components+c]=frame.field[data.nodes[tet_face_node(face.ordinal,n)][face.element]*Components+c];
        tri_gradient_numerator<Components>(area,values,Components==1,true,local);
        for (int n=0;n<3;++n) for (int c=0;c<Components*3;++c)
            sum[data.nodes[tet_face_node(face.ordinal,n)][face.element]*Components*3+c]+=local[n*Components*3+c];
    }
    for (std::size_t i=0;i<sum.size();++i) require(finish_gradient(sum[i],volume[i/(Components*3)],0.,1.,false,result[i]),"invalid host gradient");
    return result;
}
void execute(const Fixture& data) {
    const std::size_t e=data.nodes[0].size(); std::vector<Geometry> geometry(e);
    std::vector<double> volume(data.x.size(),0),values(138*e);
    for (std::size_t i=0;i<e;++i) {
        double x[12];
        for (int n=0;n<4;++n) { int node=data.nodes[n][i]; x[3*n]=data.x[node]; x[3*n+1]=data.y[node]; x[3*n+2]=data.z[node]; }
        require(tet_geometry(x,geometry[i]),"invalid host element geometry"); element_values(geometry[i],values.data()+138*i);
        for (int n=0;n<4;++n) volume[data.nodes[n][i]]+=geometry[i].volume/4;
    }
    compare(values,data.element_expected,"element geometry"); values.resize(18*data.faces.size());
    for (std::size_t i=0;i<data.faces.size();++i) face_values(geometry[data.faces[i].element],data.faces[i].ordinal,values.data()+18*i);
    compare(values,data.face_expected,"boundary geometry");
    for (const auto& frame:data.frames) compare(frame.components==1?reconstruct<1>(data,geometry,volume,frame):reconstruct<3>(data,geometry,volume,frame),
        frame.expected,std::string(frame.components==1?"pressure":"velocity")+" gradient call="+std::to_string(frame.call));
}
#endif
int main(int argc,char** argv) {
    try {
        require(argc==2,"usage: mars_segregated_geometry_check PUBLIC_GEOMETRY_FILE"); execute(load(argv[1]));
#ifdef MARS_REPLAY_CUDA
        const char* backend="CUDA";
#else
        const char* backend="host";
#endif
        std::cout << "PASS: " << backend << " native geometry and assembled gradients scalar_checks=" << checks
                  << "; momentum/pressure assembly, MPI and full SIMPLE iteration not tested\n";
        return 0;
    } catch (const std::exception& error) { std::cerr << "FAIL: " << error.what() << '\n'; return 1; }
}
