#ifdef MARS_REPLAY_CUDA
#include "mars_segregated_assembly_device.hpp"
#else
#include "mars_segregated_assembly.hpp"
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace mars::segregated;
using Geometry = TetGeometry<double>;
void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
template<class T> void read(std::istream& s,T& value) {
    require(bool(s>>value) && std::isfinite(double(value)),"invalid/truncated assembly input");
}
template<class T> void read_all(std::istream& s,T& values) { for (auto& v:values) read(s,v); }
struct Frame {
    int components,call;
    std::vector<double> velocity,pressure,alpha,boundary_factor,expected_lhs,expected_rhs,expected_d;
    std::vector<TetInteriorInput> interior;
    std::vector<BoundaryAssemblyInput> boundary;
    std::vector<SteadyMomentumNode> node;
};
struct Fixture {
    std::vector<double> x,y,z;
    std::array<std::vector<int>,4> nodes;
    std::vector<TetBoundaryFace> faces;
    std::vector<int> offsets,columns;
    std::vector<Frame> frames;
};
Fixture load(const char* path) {
    std::ifstream s(path); std::string magic; int n,e,b,f,blocks;
    require(bool(s>>magic>>n>>e>>b>>f>>blocks) && magic=="MARS_PUBLIC_ASSEMBLY_V1"
            && n==425 && e==1536 && b==576 && f==4 && blocks==4921,"expected complete public assembly fixture");
    Fixture data; data.x.resize(n); data.y.resize(n); data.z.resize(n);
    for (int i=0;i<n;++i) { read(s,data.x[i]); read(s,data.y[i]); read(s,data.z[i]); }
    for (auto& nodes:data.nodes) nodes.resize(e);
    for (int i=0;i<e;++i) {
        std::set<int> distinct;
        for (auto& nodes:data.nodes) { read(s,nodes[i]); require(nodes[i]>=0 && nodes[i]<n,"bad connectivity"); distinct.insert(nodes[i]); }
        require(distinct.size()==4,"repeated cell node");
    }
    data.faces.resize(b); std::set<std::pair<int,int>> faces;
    for (auto& face:data.faces) {
        read(s,face.element); read(s,face.ordinal);
        require(face.element>=0 && face.element<e && face.ordinal>=0 && face.ordinal<4,"bad exterior face");
        require(faces.emplace(face.element,face.ordinal).second,"duplicate exterior face");
    }
    data.offsets.resize(n+1); data.columns.resize(blocks); read_all(s,data.offsets); read_all(s,data.columns);
    require(data.offsets.front()==0 && data.offsets.back()==blocks,"bad reference graph offsets");
    for (int row=0;row<n;++row) {
        require(data.offsets[row]>=0 && data.offsets[row]<data.offsets[row+1]
                && data.offsets[row+1]<=blocks,"bad reference graph row");
        int previous=-1;
        for (int k=data.offsets[row];k<data.offsets[row+1];++k) {
            require(data.columns[k]>previous && data.columns[k]<n,"unsorted reference columns"); previous=data.columns[k];
        }
    }
    for (int k=0;k<f;++k) {
        Frame frame; read(s,frame.components); read(s,frame.call);
        require(frame.components==(k%2==0?3:1) && frame.call==k/2+1,"bad assembly state order");
        frame.velocity.resize(3*n); frame.pressure.resize(n);
        for (int i=0;i<n;++i) { for (int j=0;j<3;++j) read(s,frame.velocity[3*i+j]); read(s,frame.pressure[i]); }
        frame.interior.resize(e);
        for (auto& x:frame.interior) {
            x.stage=frame.components==3?1:0;
            read_all(s,x.coordinates); read_all(s,x.velocity); read_all(s,x.density); read_all(s,x.viscosity);
            read_all(s,x.velocity_blend); read_all(s,x.stored_flux); read_all(s,x.pressure);
            read_all(s,x.density_blend); read_all(s,x.density_gradient);
        }
        frame.boundary.resize(b); std::set<std::pair<int,int>> coverage;
        for (auto& input:frame.boundary) {
            auto& x=input.values;
            read(s,x.stage); read(s,input.element); read(s,input.face); read_all(s,input.nodes);
            require(x.stage>=0 && x.stage<=5 && (x.stage>=3)==(frame.components==3),"wrong boundary stage");
            require(faces.count({input.element,input.face}) && coverage.emplace(input.element,input.face).second,"bad boundary coverage");
            read_all(s,x.face_nodes); read_all(s,x.nearest); read_all(s,x.opposing); read_all(s,x.reversal);
            const int count=x.stage==5?3:4;
            std::set<int> face_nodes,nearest,global_nodes;
            for (int j=0;j<count;++j) {
                require(input.nodes[j]>=0 && input.nodes[j]<n,"bad boundary node"); global_nodes.insert(input.nodes[j]);
            }
            require(int(global_nodes.size())==count,"repeated boundary node");
            for (int j=0;j<3;++j) {
                require(x.face_nodes[j]>=0 && x.face_nodes[j]<count && x.nearest[j]>=0 && x.nearest[j]<count,"bad boundary sample map");
                face_nodes.insert(x.face_nodes[j]); nearest.insert(x.nearest[j]);
                require(x.reversal[j]==0 || x.reversal[j]==1,"bad reversal flag");
                if (x.stage!=5) require(x.opposing[j]>=0 && x.opposing[j]<4,"bad opposing node");
            }
            require(face_nodes.size()==3 && face_nodes==nearest,"incomplete boundary samples");
            read_all(s,x.velocity); read_all(s,x.boundary_velocity); read_all(s,x.viscosity);
            read_all(s,x.density); read_all(s,x.pressure); read_all(s,x.stored_flux); read_all(s,x.wall_coefficient);
        }
        if (frame.components==3) {
            frame.node.resize(n); frame.alpha.resize(n); frame.boundary_factor.resize(n);
            for (auto& x:frame.node) {
                read(s,x.density); read(s,x.volume); read(s,x.pseudo_dt); read(s,x.mass_divergence);
                read_all(s,x.velocity); read_all(s,x.pressure_gradient); read_all(s,x.force); read_all(s,x.source);
                require(x.density>0 && x.pseudo_dt>0,"invalid node scales");
            }
            for (int i=0;i<n;++i) { read(s,frame.alpha[i]); read(s,frame.boundary_factor[i]); }
        }
        frame.expected_lhs.resize(blocks*frame.components*frame.components); frame.expected_rhs.resize(n*frame.components);
        read_all(s,frame.expected_lhs); read_all(s,frame.expected_rhs);
        if (frame.components==3) { frame.expected_d.resize(3*n); read_all(s,frame.expected_d); }
        data.frames.push_back(std::move(frame));
    }
    require(!(s>>magic),"trailing assembly data"); return data;
}
std::size_t checks=0;
void compare(const std::vector<double>& actual,const std::vector<double>& expected,const std::string& name) {
    require(actual.size()==expected.size(),"wrong comparison size"); double worst=0; int failures=0;
    for (std::size_t i=0;i<actual.size();++i) {
        ++checks; double error=std::abs(actual[i]-expected[i])/std::max(1.,std::abs(expected[i])); worst=std::max(worst,error);
        if (!std::isfinite(actual[i]) || error>1e-12) {
            if (failures++<3) std::cerr<<name<<" entry="<<i<<" expected="<<expected[i]<<" actual="<<actual[i]<<'\n';
        }
    }
    std::cout<<name<<" scalars="<<actual.size()<<" worst_scaled="<<worst<<'\n';
    require(failures==0,"assembled reference mismatch");
}
#ifndef MARS_REPLAY_CUDA
template<int Components> std::vector<double> reconstruct_host(const Fixture& data,const std::vector<Geometry>& g,
    const std::vector<double>& volume,const std::vector<double>& field) {
    std::vector<double> sum(volume.size()*Components*3,0.);
    for (std::size_t e=0;e<g.size();++e) {
        double values[4*Components],local[12*Components];
        for (int n=0;n<4;++n) for (int c=0;c<Components;++c) values[n*Components+c]=field[data.nodes[n][e]*Components+c];
        tet_gradient_numerator<Components>(g[e],values,Components==1,true,local);
        for (int n=0;n<4;++n) for (int c=0;c<3*Components;++c) sum[data.nodes[n][e]*Components*3+c]+=local[n*Components*3+c];
    }
    for (auto face:data.faces) {
        double values[3*Components],local[9*Components],area[3]; tet_boundary_area(g[face.element],face.ordinal,area);
        for (int n=0;n<3;++n) for (int c=0;c<Components;++c)
            values[n*Components+c]=field[data.nodes[tet_face_node(face.ordinal,n)][face.element]*Components+c];
        tri_gradient_numerator<Components>(area,values,Components==1,true,local);
        for (int n=0;n<3;++n) for (int c=0;c<3*Components;++c)
            sum[data.nodes[tet_face_node(face.ordinal,n)][face.element]*Components*3+c]+=local[n*Components*3+c];
    }
    for (std::size_t i=0;i<sum.size();++i) sum[i]/=volume[i/(Components*3)];
    return sum;
}
template<int Components> void assemble_host(const Fixture& data,const Frame& frame,const std::vector<Geometry>& g,
    const std::vector<double>& volume,const std::vector<double>& vg,const std::vector<double>& pg,std::vector<double>& d) {
    std::vector<double> lhs(data.columns.size()*Components*Components,0),rhs(data.x.size()*Components,0);
    BlockCsrView<Components> matrix{int(data.x.size()),data.offsets.data(),data.columns.data(),lhs.data(),rhs.data()};
    for (std::size_t e=0;e<g.size();++e) {
        int nodes[4]; for (int n=0;n<4;++n) nodes[n]=data.nodes[n][e];
        auto x=frame.interior[e]; native_interior(x,g[e],nodes,vg.data(),pg.data(),d.data());
        TetInteriorOutput y; tet_interior(x,y); require(scatter_block(matrix,nodes,4,y.lhs,y.rhs),"missing interior block");
    }
    for (const auto& input:frame.boundary) {
        int nodes[4]; for (int n=0;n<4;++n) nodes[n]=data.nodes[n][input.element];
        auto x=input.values; require(native_boundary(x,input,g[input.element],nodes,pg.data(),d.data()),"invalid boundary map");
        BoundaryOutput y; boundary_block(x,y);
        require(scatter_block(matrix,input.nodes,x.stage==5?3:4,y.lhs,y.rhs),"missing boundary block");
    }
    if constexpr (Components==3) {
        for (int n=0;n<matrix.nodes;++n) {
            auto x=frame.node[n]; x.volume=volume[n]; for (int j=0;j<3;++j) x.pressure_gradient[j]=pg[3*n+j];
            double a[9],b[3]; steady_momentum_node(x,a,b); require(scatter_block(matrix,&n,1,a,b),"missing node block");
        }
        for (int n=0;n<matrix.nodes;++n) {
            double dt[3]; require(finish_momentum_row(matrix,n,volume[n],frame.alpha[n],frame.boundary_factor[n],false,d.data()+3*n,dt),"bad influence");
        }
    }
    const auto label=std::string(Components==3?"momentum":"pressure")+" call="+std::to_string(frame.call);
    compare(lhs,frame.expected_lhs,label+" matrix"); compare(rhs,frame.expected_rhs,label+" RHS");
    if constexpr (Components==3) compare(d,frame.expected_d,label+" influence");
}
void execute(const Fixture& data) {
    const int n=int(data.x.size()),e=int(data.nodes[0].size()); std::vector<Geometry> g(e);
    std::vector<double> volume(n,0),d(3*n,0);
    for (int i=0;i<e;++i) {
        double coordinates[12]; for (int j=0;j<4;++j) { int node=data.nodes[j][i];
            coordinates[3*j]=data.x[node]; coordinates[3*j+1]=data.y[node]; coordinates[3*j+2]=data.z[node]; }
        require(tet_geometry(coordinates,g[i]),"bad host geometry");
        for (int j=0;j<4;++j) volume[data.nodes[j][i]]+=g[i].volume/4;
    }
    for (const auto& frame:data.frames) {
        auto vg=reconstruct_host<3>(data,g,volume,frame.velocity),pg=reconstruct_host<1>(data,g,volume,frame.pressure);
        if (frame.components==3) assemble_host<3>(data,frame,g,volume,vg,pg,d); else assemble_host<1>(data,frame,g,volume,vg,pg,d);
    }
}
#else
template<class T> struct Buffer {
    T* data=nullptr;
    explicit Buffer(std::size_t n) { assembly_cuda_check(cudaMalloc(reinterpret_cast<void**>(&data),n*sizeof(T))); }
    ~Buffer() { cudaFree(data); }
    Buffer(const Buffer&)=delete; Buffer& operator=(const Buffer&)=delete;
    void upload(const std::vector<T>& x) { assembly_cuda_check(cudaMemcpy(data,x.data(),x.size()*sizeof(T),cudaMemcpyHostToDevice)); }
    std::vector<T> download(std::size_t n) { std::vector<T> x(n); assembly_cuda_check(cudaMemcpy(x.data(),data,n*sizeof(T),cudaMemcpyDeviceToHost)); return x; }
};
template<int Components> void reconstruct_device(TetMeshView<int,double> mesh,const Geometry* g,
    const TetBoundaryFace* faces,int count,const double* volume,const double* field,double* sum,double* result,int* error) {
    assembly_cuda_check(cudaMemset(sum,0,mesh.node_count*Components*3*sizeof(double)));
    accumulate_interior_gradient<Components><<<(mesh.element_count+127)/128,128>>>(mesh,g,field,0,mesh.element_count,Components==1,true,sum);
    assembly_cuda_check(cudaGetLastError());
    accumulate_boundary_gradient<Components><<<(count+127)/128,128>>>(mesh,g,faces,count,field,Components==1,true,sum,error);
    assembly_cuda_check(cudaGetLastError());
    normalize_gradient<Components><<<(mesh.node_count*Components*3+127)/128,128>>>(sum,volume,mesh.node_count,1.,false,result,error);
    assembly_cuda_check(cudaGetLastError());
}
template<int Components> void assemble_device(TetMeshView<int,double> mesh,const Geometry* g,DeviceBlockGraph& graph,
    const TetInteriorInput* interior,const BoundaryAssemblyInput* boundary,int faces,const SteadyMomentumNode* node,
    const double* volume,const double* vg,const double* pg,const double* alpha,const double* factor,
    double* d,double* dt,double* lhs,double* rhs,int* error) {
    auto matrix=graph.view<Components>(lhs,rhs);
    assembly_cuda_check(cudaMemset(lhs,0,graph.columns.size()*Components*Components*sizeof(double)));
    assembly_cuda_check(cudaMemset(rhs,0,mesh.node_count*Components*sizeof(double)));
    assemble_interior_blocks<Components><<<(mesh.element_count+63)/64,64>>>(mesh,g,interior,vg,pg,d,0,mesh.element_count,matrix,error);
    assembly_cuda_check(cudaGetLastError());
    assemble_boundary_blocks<Components><<<(faces+63)/64,64>>>(mesh,g,boundary,faces,pg,d,matrix,error);
    assembly_cuda_check(cudaGetLastError());
    if constexpr (Components==3) {
        assemble_momentum_nodes<<<(mesh.node_count+127)/128,128>>>(node,0,mesh.node_count,volume,pg,matrix,error);
        assembly_cuda_check(cudaGetLastError());
        finish_momentum_rows<<<(mesh.node_count+127)/128,128>>>(matrix,volume,alpha,factor,false,0,mesh.node_count,d,dt,error);
        assembly_cuda_check(cudaGetLastError());
    }
}
void execute(const Fixture& data) {
    const int n=int(data.x.size()),e=int(data.nodes[0].size()),b=int(data.faces.size());
    Buffer<double> x(n),y(n),z(n),volume(n),velocity(3*n),pressure(n),sum(9*n),vg(9*n),pg(3*n),d(3*n),dt(3*n),alpha(n),factor(n);
    Buffer<int> n0(e),n1(e),n2(e),n3(e),error(1); Buffer<Geometry> g(e); Buffer<TetBoundaryFace> faces(b);
    Buffer<TetInteriorInput> interior(e); Buffer<BoundaryAssemblyInput> boundary(b); Buffer<SteadyMomentumNode> node(n);
    x.upload(data.x); y.upload(data.y); z.upload(data.z); faces.upload(data.faces);
    n0.upload(data.nodes[0]); n1.upload(data.nodes[1]); n2.upload(data.nodes[2]); n3.upload(data.nodes[3]);
    TetMeshView<int,double> mesh{{n0.data,n1.data,n2.data,n3.data},x.data,y.data,z.data,n,e};
    assembly_cuda_check(cudaMemset(error.data,0,sizeof(int)));
    build_tet_geometry<<<(e+127)/128,128>>>(mesh,g.data,error.data); assembly_cuda_check(cudaGetLastError());
    require(error.download(1)[0]==0,"invalid CUDA geometry");
    DeviceBlockGraph graph; graph.build(mesh);
    std::vector<int> offsets(graph.offsets.size()),columns(graph.columns.size());
    assembly_cuda_check(cudaMemcpy(offsets.data(),device_data(graph.offsets),offsets.size()*sizeof(int),cudaMemcpyDeviceToHost));
    assembly_cuda_check(cudaMemcpy(columns.data(),device_data(graph.columns),columns.size()*sizeof(int),cudaMemcpyDeviceToHost));
    require(offsets==data.offsets && columns==data.columns,"device CSR pattern differs from topology oracle");
    std::cout<<"device CSR node_blocks="<<columns.size()<<" pattern PASS\n";
    Buffer<double> lhs(columns.size()*9),rhs(3*n);
    assembly_cuda_check(cudaMemset(volume.data,0,n*sizeof(double)));
    accumulate_dual_volumes<<<(e+127)/128,128>>>(mesh,g.data,0,e,volume.data); assembly_cuda_check(cudaGetLastError());
    for (const auto& frame:data.frames) {
        velocity.upload(frame.velocity); pressure.upload(frame.pressure); interior.upload(frame.interior); boundary.upload(frame.boundary);
        reconstruct_device<3>(mesh,g.data,faces.data,b,volume.data,velocity.data,sum.data,vg.data,error.data);
        reconstruct_device<1>(mesh,g.data,faces.data,b,volume.data,pressure.data,sum.data,pg.data,error.data);
        require(error.download(1)[0]==0,"invalid reconstructed gradient");
        if (frame.components==3) {
            node.upload(frame.node); alpha.upload(frame.alpha); factor.upload(frame.boundary_factor);
            assemble_device<3>(mesh,g.data,graph,interior.data,boundary.data,b,node.data,volume.data,vg.data,pg.data,alpha.data,factor.data,d.data,dt.data,lhs.data,rhs.data,error.data);
        } else assemble_device<1>(mesh,g.data,graph,interior.data,boundary.data,b,node.data,volume.data,vg.data,pg.data,alpha.data,factor.data,d.data,dt.data,lhs.data,rhs.data,error.data);
        require(error.download(1)[0]==0,"assembly missing entry, invalid map or invalid influence");
        auto label=std::string(frame.components==3?"momentum":"pressure")+" call="+std::to_string(frame.call);
        compare(lhs.download(frame.expected_lhs.size()),frame.expected_lhs,label+" matrix");
        compare(rhs.download(frame.expected_rhs.size()),frame.expected_rhs,label+" RHS");
        if (frame.components==3) compare(d.download(3*n),frame.expected_d,label+" influence");
    }
}
#endif
int main(int argc,char** argv) {
    try {
        require(argc==2,"usage: mars_segregated_assembly_check PUBLIC_ASSEMBLY_FILE"); execute(load(argv[1]));
#ifdef MARS_REPLAY_CUDA
        const char* backend="CUDA";
#else
        const char* backend="host";
#endif
        std::cout<<"PASS: "<<backend<<" native block assembly scalar_checks="<<checks
                 <<"; frozen stage fields/boundary history; solves, MPI and full SIMPLE iteration not tested\n";
        return 0;
    } catch (const std::exception& error) { std::cerr<<"FAIL: "<<error.what()<<'\n'; return 1; }
}
