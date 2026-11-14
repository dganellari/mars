#include "examples/distributed/unstructured/mars_mixed_nonconformal_runtime.hpp"
#include <limits>

using namespace mixedflow;

static State manufactured(int64_t id,const Vec3<double>& x) {
    State s{};s.p=.7*x.x-.3*x.y+.4*x.z+.013*(id%11);
    for(int j=0;j<3;++j) {
        s.u[j]=.1*(j+1)+.07*x.x-.03*x.y+.05*x.z+.001*(id%7);
        s.d[j]=.2+.03*j+.0001*(id%13);s.gp[j]=.02*(j+1);s.blend[j]=.8;
        for(int k=0;k<3;++k) s.gu[3*j+k]=.01*(j-k+1);
    }
    return s;
}
static Flux poisoned_flux() {
    Flux f{};f.mass=std::numeric_limits<double>::quiet_NaN();
    for(auto& x:f.pressure) x=f.mass;for(auto& x:f.momentum) x=f.mass;for(auto& x:f.coefficient) x=f.mass;
    return f;
}
static bool close(double a,double b) {return std::isfinite(a) && std::abs(a-b)<=2e-13*(1+std::abs(b));}

static void verify_reduced_stencil() {
    Input input;auto& mesh=input.mesh;mesh.globalNodes=8;mesh.partitioned=true;
    Cell cell{};cell.count=8;
    for(int k=0;k<8;++k) {
        cell.node[k]=k;mesh.dofs.push_back(k);
        mesh.x.push_back({.5*(1+sign(k,0)),.5*(1+sign(k,1)),.5*(1+sign(k,2))});
    }
    input.cells.push_back(cell);input.element_ids.push_back(0);mesh.elements.push_back(0);
    for(int f=0;f<6;++f) input.faces.push_back({0,f,2,f});
    for(int r=0;r<rankCount;++r) if(r!=rankId) mesh.candidatePeers.push_back(r);
    Channel full(input,InterfaceMethod::legacy),reduced(input,InterfaceMethod::compatibility);
    Buffer<double> full_rhs(3*full.n_),reduced_rhs(3*reduced.n_);Buffer<int> flags(full.n_);
    for(Channel* channel:{&full,&reduced}) {
        auto states=host<State>(channel->state_);
        for(int k=0;k<8;++k) {double volume=states[k].volume;states[k]=manufactured(k,mesh.x[k]);states[k].volume=volume;}
        channel->state_=Buffer<State>(states.begin(),states.end());channel->view_.state=pointer(channel->state_);
        channel->view_.completed=1;channel->view_.controls.mu=.2;
        auto matrix=channel->pressure_.view();matrix.values=pointer(channel->momentum_.values);
        each<device>(channel->n_,Assemble{channel->view_,matrix,pointer(channel==&full?full_rhs:reduced_rhs),pointer(flags),true});
        ncflow::finish();ncflow::validFlags(flags,"synthetic compatibility assembly");
    }
    const auto a=host<double>(full.momentum_.values),b=host<double>(reduced.momentum_.values);
    const auto fa=host<double>(full_rhs),fb=host<double>(reduced_rhs);
    const auto rows=host<int64_t>(full.pressure_.rowIds),cols=host<int64_t>(full.pressure_.columns);
    const auto offsets=host<int>(full.pressure_.offsets),slots=host<int>(reduced.coefficient_slots_);
    const auto fs=host<State>(full.state_),rs=host<State>(reduced.state_);
    // Cube vertex adjacency, independent of the map built by Channel.
    constexpr int edges[8][3]={{1,3,4},{0,2,5},{1,3,6},{0,2,7},{0,5,7},{1,4,6},{2,5,7},{3,4,6}};
    bool good=true;
    for(int r=0;r<full.n_;++r) {
        int node=int(rows[r]),diagonal=offsets[r]+node;
        for(int p=offsets[r];p<offsets[r+1];++p) {
            int c=int(cols[p]);bool retained=c==node || c==edges[node][0] || c==edges[node][1] || c==edges[node][2];
            good=good && slots[p]==(retained?p:diagonal);
        }
        for(int i=0;i<3;++i) {
            good=good && fa[3*r+i]==fb[3*r+i];
            // Each tangential wall contributes its opposite-corner weight 1/16.
            good=good && close(rs[node].volume/rs[node].d[i]-fs[node].volume/fs[node].d[i],.125*.2);
            for(int j=0;j<3;++j) {
                double sum_a=0,sum_b=0;
                for(int k=0;k<8;++k) {int p=9*offsets[r]+24*i+3*k+j;sum_a+=a[p];sum_b+=b[p];}
                good=good && close(sum_a,sum_b);
            }
        }
    }
    check(good,"reduced block row sums, full residual and analytic wall diagonal");
    if(rankId==0) std::puts("PASS: compatibility coefficient folding preserves block row sums and full residual");
}

static void verify(Channel& channel) {
    const auto& mesh=channel.input_.mesh;
    check(channel.view_.conservative && channel.view_.rank==rankId,"conservative mode and rank initialized");
    auto states=host<State>(channel.state_);
    double local_volume=0,volume=0;
    for(size_t n=0;n<states.size();++n) if(owner(mesh.dofs[n],mesh.globalNodes)==rankId) local_volume+=states[n].volume;
    MPI_Allreduce(&local_volume,&volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    check(std::isfinite(volume) && std::abs(volume-2)<1e-11,"public PI fixture volume equals two");
    for(size_t n=0;n<states.size();++n) {
        double volume=states[n].volume;states[n]=manufactured(mesh.dofs[n],mesh.x[n]);states[n].volume=volume;
        if(owner(mesh.dofs[n],mesh.globalNodes)!=rankId) {
            states[n].p=std::numeric_limits<double>::quiet_NaN();
            for(int j=0;j<3;++j) states[n].u[j]=states[n].d[j]=states[n].gp[j]=states[n].p;
        }
    }
    channel.state_=Buffer<State>(states.begin(),states.end());channel.view_.state=pointer(channel.state_);
    channel.exchange();states=host<State>(channel.state_);bool good=true;
    for(size_t n=0;n<states.size();++n) {
        auto expected=manufactured(mesh.dofs[n],mesh.x[n]);good=good && close(states[n].p,expected.p);
        for(int j=0;j<3;++j) good=good && close(states[n].u[j],expected.u[j]) && close(states[n].d[j],expected.d[j]) && close(states[n].gp[j],expected.gp[j]);
    }
    check(good,"poisoned ghost state replaced by owner data");
    const auto patches=host<Interface>(channel.interfaces_);const auto owners=host<int>(channel.patch_owner_);
    std::vector<double> mass(patches.size()),seed(patches.size());
    int local_ghosts=0,local_owners=0;
    for(size_t q=0;q<patches.size();++q) {
        const auto& p=patches[q];for(int k=0;k<p.count[0]+p.count[1];++k) seed[q]+=double(mesh.dofs[p.node[k]]+1)*1e-6;
        seed[q]*=std::sqrt(dot(p.area,p.area));mass[q]=seed[q];
        if(owners[q]!=rankId) {mass[q]=std::numeric_limits<double>::quiet_NaN();++local_ghosts;}else ++local_owners;
    }
    channel.view_.controls.rho=1.3;channel.view_.controls.mu=.2;channel.view_.controls.alpha_mass=.65;
    channel.interface_mass_=Buffer<double>(mass.begin(),mass.end());channel.view_.interface_mass=pointer(channel.interface_mass_);
    channel.update_interface_mass();mass=host<double>(channel.interface_mass_);good=true;
    for(size_t q=0;q<patches.size();++q) {
        const auto f=interface_flux(patches[q],states.data(),channel.view_.controls,seed[q]);
        good=good && close(mass[q],.65*f.mass+.35*seed[q]);
    }
    check(good,"owned mass update and poisoned ghost mass exchange");
    std::vector<Flux> fluxes(patches.size(),poisoned_flux());
    channel.interface_fluxes_=Buffer<Flux>(fluxes.begin(),fluxes.end());channel.view_.interface_fluxes=pointer(channel.interface_fluxes_);
    channel.prepare_interface_flux();fluxes=host<Flux>(channel.interface_fluxes_);good=true;
    for(size_t q=0;q<patches.size();++q) {
        const auto f=interface_flux(patches[q],states.data(),channel.view_.controls,mass[q]);const auto& got=fluxes[q];
        good=good && close(got.mass,f.mass);
        for(int j=0;j<3;++j) good=good && close(got.momentum[j],f.momentum[j]);
        for(int k=0;k<patches[q].count[0]+patches[q].count[1];++k) {
            good=good && close(got.pressure[k],f.pressure[k]);
            for(int j=0;j<9;++j) good=good && close(got.coefficient[9*k+j],f.coefficient[9*k+j]);
        }
    }
    check(good,"owner flux computed and poisoned ghost flux replaced");
    const auto offsets=host<int>(channel.qo_),indices=host<int>(channel.qi_),signs=host<int>(channel.qs_);
    const auto graph_offsets=host<int>(channel.pressure_.offsets);const auto columns=host<int64_t>(channel.pressure_.columns);
    std::vector<double> matrix(10*mesh.globalNodes),matrix_scale(matrix.size()),global_matrix(matrix.size()),global_matrix_scale(matrix.size());
    double local[12]{},global[12]{};int incidences=0;good=true;
    for(int row=0;row<channel.n_;++row) for(int i=offsets[row];i<offsets[row+1];++i) {
        ++incidences;int q=indices[i],sign=signs[i];const auto& p=patches[q];const auto& f=fluxes[q];
        good=good && (sign==1 || sign==-1);
        double values[5]={f.mass,mass[q],f.momentum[0],f.momentum[1],f.momentum[2]};
        for(int j=0;j<5;++j) {local[j]+=sign*values[j];local[5+j]+=std::abs(values[j]);}
        local[10]+=std::sqrt(dot(p.area,p.area));
        for(int k=0;k<p.count[0]+p.count[1];++k) {
            const int64_t column=mesh.dofs[p.node[k]];
            good=good && std::binary_search(columns.begin()+graph_offsets[row],columns.begin()+graph_offsets[row+1],column);
            matrix[10*column]+=sign*f.pressure[k];matrix_scale[10*column]+=std::abs(f.pressure[k]);
            for(int j=0;j<9;++j) {matrix[10*column+1+j]+=sign*f.coefficient[9*k+j];matrix_scale[10*column+1+j]+=std::abs(f.coefficient[9*k+j]);}
        }
    }
    check(good,"shared incidence signs and complete row graph");
    for(auto query:channel.input_.queries) {
        const auto face=channel.input_.faces[query.face];const auto cell=channel.input_.cells[face.cell];double xyz[24]{};
        for(int k=0;k<cell.count;++k) {
            const auto x=mesh.x[cell.node[k]];xyz[3*k]=x.x;xyz[3*k+1]=x.y;xyz[3*k+2]=x.z;
        }
        Point point;good=good && face_point(cell.count,xyz,face.ordinal,query.sample,true,point,true);
        local[11]+=std::sqrt(dot(point.area,point.area));
    }
    check(good,"host reference interface geometry");
    MPI_Allreduce(local,global,12,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(matrix.data(),global_matrix.data(),int(matrix.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(matrix_scale.data(),global_matrix_scale.data(),int(matrix.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    int counts[3]={local_owners,local_ghosts,incidences},global_counts[3];MPI_Allreduce(counts,global_counts,3,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    good=global_counts[0]>0 && global_counts[2]==2*global_counts[0] && (rankCount==1 || global_counts[1]>0);
    double max_ratio=0;
    for(int j=0;j<5;++j) {good=good && global[5+j]>0;max_ratio=std::max(max_ratio,std::abs(global[j])/global[5+j]);}
    for(size_t j=0;j<matrix.size();++j) {
        const double scale=global_matrix_scale[j];if(scale>0) max_ratio=std::max(max_ratio,std::abs(global_matrix[j])/scale);
        else good=good && global_matrix[j]==0;
    }
    good=good && max_ratio<2e-12 && close(global[10],global[11]);
    check(good,"shared mass momentum and per-column Jacobian cancellation and coverage");
    if(rankId==0) std::printf("PASS: shared conservative runtime ranks=%d patches=%d ghosts=%d incidences=%d area=%.17g mass_scale=%.17g cancellation=%.3e\n",rankCount,global_counts[0],global_counts[1],global_counts[2],global[10],global[5],max_ratio);
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rankId);MPI_Comm_size(MPI_COMM_WORLD,&rankCount);
    try {
        check(argc==2 || (argc==3 && std::string(argv[2])=="--gap"),"partition directory argument and optional --gap");
#if defined(__CUDACC__) || defined(__HIPCC__)
        MPI_Comm local;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rankId,MPI_INFO_NULL,&local);int local_rank,devices;
        MPI_Comm_rank(local,&local_rank);check(cudaGetDeviceCount(&devices)==cudaSuccess && devices>0,"mixed GPU selection");
        const int selected_device=devices==1?0:local_rank%devices;
        check(cudaSetDevice(selected_device)==cudaSuccess,"mixed selected GPU");MPI_Comm_free(&local);
        cudaDeviceProp properties{};check(cudaGetDeviceProperties(&properties,selected_device)==cudaSuccess,"mixed GPU identity");
        char uuid[33]{};for(int i=0;i<16;++i) std::snprintf(uuid+2*i,3,"%02x",static_cast<unsigned char>(properties.uuid.bytes[i]));
        std::printf("GPU rank=%d device=%d uuid=%s\n",rankId,selected_device,uuid);std::fflush(stdout);
#endif
        verify_reduced_stencil();
        auto input=read(argv[1]);
        {
            // The PI mesh has unequal owned-row counts on four ranks.
            Channel compatibility(input,InterfaceMethod::compatibility);
            check(compatibility.view_.coefficient_slots!=nullptr,"compatibility setup on uneven partitions");
        }
        if(argc==3) {
            int local_face=INT_MAX,removed_face=INT_MAX;
            for(auto face:input.faces) if(face.kind==4) local_face=std::min(local_face,face.id);
            MPI_Allreduce(&local_face,&removed_face,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
            check(removed_face!=INT_MAX,"gap control has an interface face");
            std::vector<int> candidates;
            for(auto& query:input.queries) {
                int start=int(candidates.size());
                if(input.faces[query.face].id!=removed_face)
                    for(int i=query.start;i<query.start+query.count;++i)
                        if(input.faces[input.donors[i]].id!=removed_face) candidates.push_back(input.donors[i]);
                query.start=start;query.count=int(candidates.size())-start;
            }
            input.donors=std::move(candidates);
        }
        Channel channel(input,InterfaceMethod::conservative);
        if(argc==3) throw std::runtime_error("FAIL: missing interface coverage was accepted");
        verify(channel);
    } catch(const std::exception& e) {std::fprintf(stderr,"rank %d: %s\n",rankId,e.what());MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
