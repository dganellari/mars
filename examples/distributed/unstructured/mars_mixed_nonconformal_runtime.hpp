#pragma once
#include "mars_nonconformal_flow_runtime.hpp"
#include "backend/distributed/unstructured/fem/mars_mixed_flow.hpp"
#include <filesystem>
#include <fstream>
#include <tuple>

namespace mixedflow {
using namespace mars::mixed;
enum class InterfaceMethod { compatibility, conservative, legacy };
enum class FlowProfile { mars, openaccel };
struct Face { int cell,ordinal,kind,id; };
struct Query { int face,sample,start,count,row; };
struct Input { Mesh mesh;std::vector<Cell> cells;std::vector<int> element_ids;std::vector<Face> faces;std::vector<Query> queries;std::vector<int> donors; };
inline Input read(const std::string& directory) {
    const auto path=std::filesystem::path(directory)/("rank"+std::to_string(rankId)+".mmx");
    std::ifstream file(path,std::ios::binary);
    auto take=[&](void* p,size_t n){if(n && !file.read(static_cast<char*>(p),n)) throw std::runtime_error("truncated mixed partition");};
    char magic[8];uint64_t h[9];take(magic,8);take(h,sizeof(h));
    if(std::memcmp(magic,"MARSMX01",8) || h[0]!=uint64_t(rankCount) || h[1]!=uint64_t(rankId) ||
       !h[2] || h[2]>8000000 || h[3]>h[2] || h[4]>1000000 || h[5]>6*h[4] || h[6]>4*h[5] || h[7]>100*h[6] || h[8]>=h[0])
        throw std::runtime_error("invalid mixed partition header");
    if(std::filesystem::file_size(path)!=80+32*h[3]+40*h[4]+16*h[5]+20*h[6]+4*h[7]+4*h[8]) throw std::runtime_error("mixed partition length mismatch");
    Input out;auto& m=out.mesh;m.globalNodes=int(h[2]);m.partitioned=true;
    m.dofs.resize(h[3]);m.x.resize(h[3]);m.elements.resize(h[4]);out.cells.resize(h[4]);out.element_ids.resize(h[4]);
    out.faces.resize(h[5]);out.queries.resize(h[6]);out.donors.resize(h[7]);m.candidatePeers.resize(h[8]);
    take(m.dofs.data(),8*h[3]);take(m.x.data(),24*h[3]);
    if(!std::is_sorted(m.dofs.begin(),m.dofs.end()) || std::adjacent_find(m.dofs.begin(),m.dofs.end())!=m.dofs.end()) throw std::runtime_error("mixed node IDs are not unique/sorted");
    for(size_t i=0;i<h[3];++i) if(m.dofs[i]<0 || m.dofs[i]>=m.globalNodes || !std::isfinite(m.x[i].x) || !std::isfinite(m.x[i].y) || !std::isfinite(m.x[i].z)) throw std::runtime_error("invalid mixed node");
    for(size_t i=0;i<h[4];++i) {
        int record[10];take(record,sizeof(record));auto& c=out.cells[i];c.count=record[0];out.element_ids[i]=record[1];m.elements[i]=record[1];
        if(c.count!=4 && c.count!=8) throw std::runtime_error("unsupported native element");
        for(int k=0;k<8;++k) {c.node[k]=record[k+2];if(k<c.count && (c.node[k]<0 || c.node[k]>=int(h[3]))) throw std::runtime_error("invalid mixed connectivity");}
    }
    take(out.faces.data(),16*h[5]);take(out.queries.data(),20*h[6]);take(out.donors.data(),4*h[7]);take(m.candidatePeers.data(),4*h[8]);
    std::set<int> face_ids;
    for(auto f:out.faces) if(f.cell<0 || f.cell>=int(h[4]) || f.ordinal<0 || f.ordinal>=face_count(out.cells[f.cell].count) || f.kind<0 || f.kind>4 || f.id<0 || !face_ids.insert(f.id).second) throw std::runtime_error("invalid mixed face or duplicate persistent face ID");
    std::set<std::tuple<int,int,int>> expected_queries;
    for(int f=0;f<int(out.faces.size());++f) if(out.faces[f].kind>=3) {
        const auto face=out.faces[f];const auto& cell=out.cells[face.cell];
        for(int sample=0;sample<face_size(cell.count);++sample) {
            int row=cell.node[face_node(cell.count,face.ordinal,sample)];
            if(owner(m.dofs[row],m.globalNodes)==rankId) expected_queries.emplace(f,sample,row);
        }
    }
    for(auto q:out.queries) {
        if(q.face<0 || q.face>=int(h[5]) || q.sample<0 || q.sample>=face_size(out.cells[out.faces[q.face].cell].count) || q.start<0 || q.count<=0 || uint64_t(q.start)+uint64_t(q.count)>h[7] || q.row<0 || q.row>=int(h[3]) || owner(m.dofs[q.row],m.globalNodes)!=rankId) throw std::runtime_error("invalid mixed interface query");
        if(expected_queries.erase({q.face,q.sample,q.row})!=1) throw std::runtime_error("duplicate or unexpected owned interface query");
    }
    if(!expected_queries.empty()) throw std::runtime_error("missing owned interface query");
    for(int d:out.donors) if(d<0 || d>=int(h[5]) || out.faces[d].kind<3) throw std::runtime_error("invalid mixed donor");
    for(auto id:ownedRows(m)) if(!std::binary_search(m.dofs.begin(),m.dofs.end(),id)) throw std::runtime_error("missing owned mixed node");
    return out;
}
struct View {
    const Cell* cells;const Geometry* geometry;const Face* faces;const Query* queries;const int* donors;
    const double* xyz;const int64_t* ids;const int* owned;
    const int *co,*ci,*fo,*fi,*qo,*qi;
    Interface* interfaces;State* state;double *cell_mass,*face_mass,*interface_mass;int* closed;
    Controls controls;int completed{};
    const int* coefficient_slots=nullptr;
    const int *qs=nullptr,*patch_owner=nullptr;
    Flux* interface_fluxes=nullptr;
    bool conservative=false;int rank=0;
    bool openaccel_profile=false;
    MARS_FS_FN void coordinates(int e,double* x) const {for(int k=0;k<cells[e].count;++k) for(int j=0;j<3;++j) x[3*k+j]=xyz[3*cells[e].node[k]+j];}
    MARS_FS_FN Point facepoint(int f,int q,bool shifted=false) const {
        double x[24];coordinates(faces[f].cell,x);Point p;
        face_point(cells[faces[f].cell].count,x,faces[f].ordinal,q,shifted,p,true);return p;
    }
};
struct BuildGeometry {
    View v;Geometry* geometry;int* valid;
    MARS_FS_FN void operator()(int e) const {double x[24];v.coordinates(e,x);valid[e]=mars::mixed::geometry(v.cells[e].count,x,geometry[e],true);}
};
enum class DonorSelection { persistent_id, openaccel_face_score };
MARS_FS_FN double donor_face_score(int count,int face,const Point& p) {
    // On a coincident face these are Nalu's Tri3/Quad4 parametric distances.
    // The face-sum check below excludes normal offsets; no thickness model is used.
    if(count==4) {
        double minimum=1;
        for(int k=0;k<3;++k) minimum=fmin(minimum,p.shape[face_node(count,face,k)]);
        return fmax(0.,1-3*minimum);
    }
    const double a=p.shape[face_node(count,face,0)],b=p.shape[face_node(count,face,1)];
    const double c=p.shape[face_node(count,face,2)],d=p.shape[face_node(count,face,3)];
    return fmax(fabs(-a+b+c-d),fabs(-a-b+c+d));
}
struct BuildInterface {
    View v;int* valid;DonorSelection selection{DonorSelection::persistent_id};
    MARS_FS_FN void operator()(int index) const {
        const auto query=v.queries[index];const auto f=v.faces[query.face];const auto& c=v.cells[f.cell];
        const Point cp=v.facepoint(query.face,query.sample,true);Interface out;out.row=query.row;out.count[0]=c.count;
        double area=sqrt(mars::mixed::dot(cp.area,cp.area)),normal[3];for(int j=0;j<3;++j) {out.area[j]=cp.area[j];normal[j]=cp.area[j]/area;}
        int best=-1;Point op;double minimum=1e100;
        // First establish the admissible score band, then select by persistent
        // face ID. Neither pass depends on candidate order or owning rank.
        for(int pass=0;pass<2;++pass) for(int i=query.start;i<query.start+query.count;++i) {
            int donor=v.donors[i];auto df=v.faces[donor];auto dc=v.cells[df.cell];double x[24];v.coordinates(df.cell,x);Point p;
            if(df.kind==f.kind || !invert(dc.count,x,cp.position,p)) continue;
            double face_sum=0;for(int k=0;k<face_size(dc.count);++k) face_sum+=p.shape[face_node(dc.count,df.ordinal,k)];
            if(fabs(face_sum-1)>1e-9) continue;
            const auto dp=v.facepoint(donor,0,true);double da=sqrt(mars::mixed::dot(dp.area,dp.area));
            if(mars::mixed::dot(dp.area,normal)/da>-.999999) continue;
            double score=0;for(int k=0;k<dc.count;++k) score=fmax(score,-p.shape[k]);
            if(selection==DonorSelection::openaccel_face_score) score=donor_face_score(dc.count,df.ordinal,p);
            if(pass==0) minimum=fmin(minimum,score);
            else if(prefer_donor(score,minimum,df.id,best<0?-1:v.faces[best].id,1e-14)) {best=donor;op=p;}
        }
        valid[index]=best>=0;if(best<0) return;
        const auto df=v.faces[best];const auto& d=v.cells[df.cell];out.count[1]=d.count;
        if(c.count+d.count>16) {valid[index]=0;return;}
        for(int side=0;side<2;++side) {
            const auto& cell=side==0?c:d;const auto& p=side==0?cp:op;int offset=side==0?0:c.count;
            for(int k=0;k<cell.count;++k) {
                int n=offset+k;out.node[n]=cell.node[k];out.shape[n]=p.shape[k];
                for(int j=0;j<3;++j) out.gradient[3*n+j]=p.gradient[3*k+j];
            }
            int ordinal=side==0?f.ordinal:df.ordinal;
            for(int k=0;k<face_size(cell.count);++k) {
                int n=offset+face_node(cell.count,ordinal,k);out.face_mask[n]=1;
                out.inverse_height[side]+=(side==0?1.:-1.)*mars::mixed::dot(out.gradient+3*n,normal);
            }
            if(!(out.inverse_height[side]>0)) valid[index]=0;
        }
        v.interfaces[index]=out;
    }
};
struct Gradient {
    View v;int field;bool relax;
    MARS_FS_FN double value(int node,int c) const {return field==1?v.state[node].u[c]:(field==2?v.state[node].increment:v.state[node].p);}
    MARS_FS_FN void operator()(int r) const {
        int node=v.owned[r],components=field==1?3:1;double sum[9]{},volume=0;
        for(int t=v.co[r];t<v.co[r+1];++t) {
            int e=v.ci[t];const auto& c=v.cells[e];const auto& g=v.geometry[e];int local=0;while(c.node[local]!=node) ++local;
            volume+=g.volume[local];
            for(int q=0;q<g.samples;++q) {
                const auto& p=g.interior[q];if(p.left!=local && p.right!=local) continue;
                for(int i=0;i<components;++i) {
                    double x=field==1?0:.5*(value(c.node[p.left],i)+value(c.node[p.right],i));
                    if(field==1) for(int k=0;k<c.count;++k) x+=p.shape[k]*value(c.node[k],i);
                    for(int j=0;j<3;++j) sum[3*i+j]+=(p.left==local?1.:-1.)*(x-value(node,i))*p.area[j];
                }
            }
        }
        for(int t=v.fo[r];t<v.fo[r+1];++t) {
            int f=v.fi[t];const auto& c=v.cells[v.faces[f].cell];
            for(int q=0;q<face_size(c.count);++q) if(c.node[face_node(c.count,v.faces[f].ordinal,q)]==node) {
                auto p=v.facepoint(f,q,field!=1);
                for(int i=0;i<components;++i) {double x=-value(node,i);for(int k=0;k<c.count;++k) x+=p.shape[k]*value(c.node[k],i);
                    for(int j=0;j<3;++j) sum[3*i+j]+=x*p.area[j];}
            }
        }
        for(int t=v.qo[r];t<v.qo[r+1];++t) {
            const auto& p=v.interfaces[v.qi[t]];
            for(int i=0;i<components;++i) {double x=-value(node,i);for(int k=0;k<p.count[0]+p.count[1];++k) x+=.5*p.shape[k]*value(p.node[k],i);
                for(int j=0;j<3;++j) sum[3*i+j]+=x*p.area[j]*(v.conservative?v.qs[t]:1);}
        }
        auto& s=v.state[node];s.volume=volume;
        for(int i=0;i<3*components;++i) {
            double* target=field==1?s.gu:(field==2?s.gc:s.gp);double fresh=sum[i]/volume;
            target[i]=relax?(1-v.controls.alpha_gradient)*target[i]+v.controls.alpha_gradient*fresh:fresh;
        }
    }
};
MARS_FS_FN double blend_bound(double value,double lo,double hi,double projection) {
    constexpr double eps=2.2204460492503131e-16;double r=projection+(projection>0?eps:-eps);
    if(fabs(r)<=100*eps) return 1;double y=fmin(hi-value,value-lo)/fabs(r);return y>=2?1:(y*y+2*y)/(y*y+y+2);
}
struct Blend {
    View v;
    MARS_FS_FN void operator()(int r) const {
        int node=v.owned[r];auto& s=v.state[node];double lo[3],hi[3],bound[3]={1,1,1};
        for(int j=0;j<3;++j) lo[j]=hi[j]=s.u[j];
        for(int t=v.co[r];t<v.co[r+1];++t) {
            int e=v.ci[t];const auto& c=v.cells[e];const auto& g=v.geometry[e];
            for(int q=0;q<g.samples;++q) {
                const auto& p=g.interior[q];int other=-1;
                if(c.node[p.left]==node) other=c.node[p.right];
                if(c.node[p.right]==node) other=c.node[p.left];
                if(other>=0) for(int j=0;j<3;++j) {double u=v.state[other].u[j];lo[j]=fmin(lo[j],u);hi[j]=fmax(hi[j],u);}
            }
        }
        for(int t=v.co[r];t<v.co[r+1];++t) {
            int e=v.ci[t];const auto& c=v.cells[e];auto& g=v.geometry[e];
            for(int q=0;q<g.samples;++q) {auto& p=g.interior[q];if(c.node[p.left]!=node && c.node[p.right]!=node) continue;
                for(int i=0;i<3;++i) {double projection=0;for(int j=0;j<3;++j) projection+=s.gu[3*i+j]*(p.position[j]-v.xyz[3*node+j]);bound[i]=fmin(bound[i],blend_bound(s.u[i],lo[i],hi[i],projection));}}
        }
        for(int t=v.fo[r];t<v.fo[r+1];++t) {int f=v.fi[t];auto c=v.cells[v.faces[f].cell];for(int q=0;q<face_size(c.count);++q) if(c.node[face_node(c.count,v.faces[f].ordinal,q)]==node) {
            auto p=v.facepoint(f,q);for(int i=0;i<3;++i) {double projection=0;for(int j=0;j<3;++j) projection+=s.gu[3*i+j]*(p.position[j]-v.xyz[3*node+j]);bound[i]=fmin(bound[i],blend_bound(s.u[i],lo[i],hi[i],projection));}}}
        if(v.openaccel_profile) for(int t=v.qo[r];t<v.qo[r+1];++t) {
            const auto query=v.queries[v.qi[t]];
            const auto p=v.facepoint(query.face,query.sample,true);
            for(int i=0;i<3;++i) {
                double projection=0;
                for(int j=0;j<3;++j) projection+=s.gu[3*i+j]*(p.position[j]-v.xyz[3*node+j]);
                bound[i]=fmin(bound[i],blend_bound(s.u[i],lo[i],hi[i],projection));
            }
        }
        for(int i=0;i<3;++i) s.blend[i]=fmin(bound[i],.75*s.blend[i]+.25*bound[i]);
    }
};
struct Assemble {
    View v;InterfaceCsrView<double> matrix;double* rhs;int* valid;bool momentum;
    MARS_FS_FN void add(int r,int node,int i,int j,double value) const {
        int pos=nonconformal_detail::columnSlot(matrix.columns,matrix.offsets[r],matrix.offsets[r+1],v.ids[node]);
        if(pos<0) {valid[r]=0;return;}
        if(v.coefficient_slots) pos=v.coefficient_slots[pos];
        if(momentum) {int len=matrix.offsets[r+1]-matrix.offsets[r];matrix.values[9*matrix.offsets[r]+3*len*i+3*(pos-matrix.offsets[r])+j]+=value;}
        else matrix.values[pos]+=value;
    }
    MARS_FS_FN void insert(int r,const int* nodes,int count,const Flux& f,double sign) const {
        int width=momentum?3:1;
        for(int i=0;i<width;++i) {rhs[width*r+i]-=sign*(momentum?f.momentum[i]:f.mass);
            for(int k=0;k<count;++k) for(int j=0;j<width;++j) add(r,nodes[k],i,j,sign*(momentum?f.coefficient[9*k+3*i+j]:f.pressure[k]));}
    }
    MARS_FS_FN void operator()(int r) const {
        int node=v.owned[r],width=momentum?3:1;valid[r]=1;
        for(int k=width*width*matrix.offsets[r];k<width*width*matrix.offsets[r+1];++k) matrix.values[k]=0;
        for(int i=0;i<width;++i) rhs[width*r+i]=0;
        for(int t=v.co[r];t<v.co[r+1];++t) {int e=v.ci[t];const auto& c=v.cells[e];auto& g=v.geometry[e];
            for(int q=0;q<g.samples;++q) {auto& p=g.interior[q];if(c.node[p.left]!=node && c.node[p.right]!=node) continue;
                insert(r,c.node,c.count,interior_flux(c,p,v.state,v.xyz,v.controls,v.cell_mass[12*e+q]),c.node[p.left]==node?1.:-1.);}}
        for(int t=v.fo[r];t<v.fo[r+1];++t) {int f=v.fi[t];auto face=v.faces[f];const auto& c=v.cells[face.cell];
            for(int q=0;q<face_size(c.count);++q) if(c.node[face_node(c.count,face.ordinal,q)]==node)
                insert(r,c.node,c.count,boundary_flux(c,face.ordinal,v.facepoint(f,q),v.state,v.xyz,v.controls,face.kind,v.face_mass[4*f+q],v.closed[f],v.completed>0),1.);}
        for(int t=v.qo[r];t<v.qo[r+1];++t) {int q=v.qi[t];auto& p=v.interfaces[q];
            insert(r,p.node,p.count[0]+p.count[1],v.conservative?v.interface_fluxes[q]:interface_flux(p,v.state,v.controls,v.interface_mass[q]),v.conservative?v.qs[t]:1.);}
        if(momentum) {
            auto& s=v.state[node];int diagonal=nonconformal_detail::columnSlot(matrix.columns,matrix.offsets[r],matrix.offsets[r+1],v.ids[node]);
            for(int i=0;i<3;++i) {
                rhs[3*r+i]+=s.divergence*s.u[i]-s.volume*s.gp[i];
                add(r,node,i,i,fmax(-s.divergence,0.)+v.controls.rho*s.volume/v.controls.dt);
                int slot=9*matrix.offsets[r]+3*(matrix.offsets[r+1]-matrix.offsets[r])*i+3*(diagonal-matrix.offsets[r])+i;
                matrix.values[slot]/=v.controls.alpha_u;s.d[i]=s.volume/(matrix.values[slot]+2.2204460492503131e-16);
                if(v.fo[r+1]>v.fo[r]) rhs[3*r+i]*=.75;
                if(!std::isfinite(s.d[i])) valid[r]=0;
            }
        }
    }
};
struct AddIncrement {
    View v;const double* increment;bool momentum;
    MARS_FS_FN void operator()(int r) const {auto& s=v.state[v.owned[r]];if(momentum) for(int j=0;j<3;++j) s.u[j]+=increment[3*r+j];else {s.increment=increment[r];s.p+=v.controls.alpha_p*increment[r];}}
};
struct Correct {
    View v;MARS_FS_FN void operator()(int r) const {auto& s=v.state[v.owned[r]];for(int j=0;j<3;++j) s.u[j]-=s.d[j]*s.gc[j];}
};
struct UpdateCellMass {
    View v;MARS_FS_FN void operator()(int e) const {const auto& c=v.cells[e];const auto& g=v.geometry[e];for(int q=0;q<g.samples;++q) {double& mass=v.cell_mass[12*e+q];auto f=interior_flux(c,g.interior[q],v.state,v.xyz,v.controls,mass);mass=v.controls.alpha_mass*f.mass+(1-v.controls.alpha_mass)*mass;}}
};
struct UpdateInterfaceMass {
    View v;MARS_FS_FN void operator()(int q) const {if(v.conservative && v.patch_owner[q]!=v.rank) return;double& mass=v.interface_mass[q];auto f=interface_flux(v.interfaces[q],v.state,v.controls,mass);mass=v.controls.alpha_mass*f.mass+(1-v.controls.alpha_mass)*mass;}
};
struct PrepareInterfaceFlux {
    View v;MARS_FS_FN void operator()(int q) const {if(v.patch_owner[q]==v.rank) v.interface_fluxes[q]=interface_flux(v.interfaces[q],v.state,v.controls,v.interface_mass[q]);}
};
struct UpdateBoundaryMass {
    View v;
    MARS_FS_FN void operator()(int f) const {
        auto face=v.faces[f];if(face.kind>=3) return;const auto& cell=v.cells[face.cell];const int nf=face_size(cell.count);double net=0,normal[3]{},mean_u[3]{},mean_p=0;
        for(int q=0;q<nf;++q) {auto p=v.facepoint(f,q);double& mass=v.face_mass[4*f+q];auto flux=boundary_flux(cell,face.ordinal,p,v.state,v.xyz,v.controls,face.kind,mass,v.closed[f],v.completed>0);
            mass=v.controls.alpha_mass*flux.mass+(1-v.controls.alpha_mass)*mass;net+=mass;
            const auto& s=v.state[cell.node[p.left]];mean_p+=s.p/nf;
            for(int j=0;j<3;++j) {normal[j]+=p.area[j];mean_u[j]+=s.u[j]/nf;}}
        if(face.kind!=1) return;
        if(!v.closed[f]) {if(net < -2.2204460492503131e-16) v.closed[f]=1;}
        else if(mars::mixed::dot(mean_u,normal)>=0 && 0<=mean_p) v.closed[f]=0;
        // The valve acts on the whole face. An open face can have both signs at
        // its quadrature points; clipping those would change the pressure flux.
        // Only the reference profile deliberately reproduces that mismatch.
        for(int q=0;q<nf;++q) {
            double& mass=v.face_mass[4*f+q];
            if(v.closed[f]) mass=0;
            else if(v.openaccel_profile) mass=fmax(mass,0.);
        }
    }
};
struct Divergence {
    View v;MARS_FS_FN void operator()(int r) const {
        int node=v.owned[r];double sum=0;
        for(int t=v.co[r];t<v.co[r+1];++t) {int e=v.ci[t];auto& c=v.cells[e];auto& g=v.geometry[e];for(int q=0;q<g.samples;++q) {auto p=g.interior[q];if(c.node[p.left]==node) sum+=v.cell_mass[12*e+q];if(c.node[p.right]==node) sum-=v.cell_mass[12*e+q];}}
        for(int t=v.fo[r];t<v.fo[r+1];++t) {int f=v.fi[t];auto face=v.faces[f];auto c=v.cells[face.cell];for(int q=0;q<face_size(c.count);++q) if(c.node[face_node(c.count,face.ordinal,q)]==node) sum+=v.face_mass[4*f+q];}
        for(int t=v.qo[r];t<v.qo[r+1];++t) sum+=(v.conservative?v.qs[t]:1)*v.interface_mass[v.qi[t]];
        v.state[node].divergence=sum;
    }
};
struct Report { double sums[14]{},maxima[4]{}; };
struct ReportSum {
    MARS_FS_FN Report operator()(Report a,const Report& b) const {
        for(int i=0;i<14;++i) a.sums[i]+=b.sums[i];
        for(int i=0;i<4;++i) a.maxima[i]=fmax(a.maxima[i],b.maxima[i]);
        return a;
    }
};
struct Statistics {
    View v;const double* rhs=nullptr;const double* inverse_diagonal=nullptr;int width=0;
    MARS_FS_FN Report operator()(int r) const {
        Report out;int node=v.owned[r];const auto& s=v.state[node];
        out.sums[5]=fabs(s.divergence);out.sums[6]=s.volume;out.maxima[0]=fmax(s.p,0.);out.maxima[3]=fmax(-s.p,0.);double speed=0;
        for(int j=0;j<3;++j) {speed+=s.u[j]*s.u[j];if(!std::isfinite(s.u[j]) || !std::isfinite(s.gp[j])) out.sums[10]+=1;}
        out.maxima[1]=sqrt(speed);if(!std::isfinite(s.p) || !std::isfinite(s.divergence)) out.sums[10]+=1;
        for(int i=0;i<width;++i) {double residual=rhs[width*r+i]*inverse_diagonal[width*r+i];out.sums[7+i]=residual*residual;}
        for(int t=v.fo[r];t<v.fo[r+1];++t) {int f=v.fi[t];auto face=v.faces[f];auto c=v.cells[face.cell];for(int q=0;q<face_size(c.count);++q) if(c.node[face_node(c.count,face.ordinal,q)]==node) {
            out.sums[face.kind]+=v.face_mass[4*f+q];
            if(face.kind==1 && v.closed[f]) {auto p=v.facepoint(f,q);out.sums[13]+=sqrt(mars::mixed::dot(p.area,p.area));}
            if(face.kind==1 && !v.closed[f]) {auto p=v.facepoint(f,q);double area=sqrt(mars::mixed::dot(p.area,p.area));out.sums[3]+=area*s.p;out.sums[4]+=area;}
        }}
        for(int t=v.qo[r];t<v.qo[r+1];++t) {double mass=(v.conservative?v.qs[t]:1)*v.interface_mass[v.qi[t]];out.sums[mass<0?11:12]+=mass;}
        return out;
    }
};
inline Report reduce(int count,Statistics f) {
    Report local,global;
#if defined(__CUDACC__) || defined(__HIPCC__)
    local=thrust::transform_reduce(thrust::device,thrust::counting_iterator<int>(0),thrust::counting_iterator<int>(count),f,Report{},ReportSum{});
#else
    for(int r=0;r<count;++r) local=ReportSum{}(local,f(r));
#endif
    MPI_Allreduce(local.sums,global.sums,14,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(local.maxima,global.maxima,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);return global;
}
class Channel {
public:
    explicit Channel(const Input& input,InterfaceMethod method=InterfaceMethod::compatibility,FlowProfile profile=FlowProfile::mars):input_(input),halo_(MPI_COMM_WORLD),flux_halo_(MPI_COMM_WORLD),mass_halo_(MPI_COMM_WORLD) {
        check(profile!=FlowProfile::openaccel || method==InterfaceMethod::compatibility,
              "OpenAccel profile requires compatibility coupling");
        const auto& m=input.mesh;const auto rows=ownedRows(m);n_=int(rows.size());
        std::vector<int> owned;for(auto id:rows) owned.push_back(int(std::lower_bound(m.dofs.begin(),m.dofs.end(),id)-m.dofs.begin()));
        std::vector<int> local_row(m.x.size(),-1);for(int r=0;r<n_;++r) local_row[owned[r]]=r;
        std::vector<std::vector<int>> by_cell(n_),by_face(n_),by_query(n_);std::vector<std::set<int64_t>> graph(n_);
        for(int e=0;e<int(input.cells.size());++e) {const auto& c=input.cells[e];for(int k=0;k<c.count;++k) {int r=local_row[c.node[k]];if(r<0) continue;by_cell[r].push_back(e);for(int j=0;j<c.count;++j) graph[r].insert(m.dofs[c.node[j]]);}}
        for(int f=0;f<int(input.faces.size());++f) {auto face=input.faces[f];if(face.kind>=3) continue;auto c=input.cells[face.cell];for(int k=0;k<face_size(c.count);++k) {int r=local_row[c.node[face_node(c.count,face.ordinal,k)]];if(r>=0) by_face[r].push_back(f);}}
        for(int q=0;q<int(input.queries.size());++q) {auto query=input.queries[q];int r=local_row[query.row];by_query.at(r).push_back(q);
            for(int k=query.start;k<query.start+query.count;++k) {auto c=input.cells[input.faces[input.donors[k]].cell];for(int j=0;j<c.count;++j) graph[r].insert(m.dofs[c.node[j]]);}}
        // Common patches use the union of face pairs, not each vertex's query
        // alone. Reserve both parents on every owned row that a pair can touch.
        if(method==InterfaceMethod::conservative) for(auto query:input.queries) for(int k=query.start;k<query.start+query.count;++k) {
            const Face pair[2]={input.faces[query.face],input.faces[input.donors[k]]};
            for(auto face:pair) {
                const auto& cell=input.cells[face.cell];
                for(int q=0;q<face_size(cell.count);++q) {
                    int r=local_row[cell.node[face_node(cell.count,face.ordinal,q)]];if(r<0) continue;
                    for(auto parent:pair) {const auto& donor=input.cells[parent.cell];for(int j=0;j<donor.count;++j) graph[r].insert(m.dofs[donor.node[j]]);}
                }
            }
        }
        std::vector<int64_t> cols;std::vector<int> off{0};for(auto& row:graph) {cols.insert(cols.end(),row.begin(),row.end());off.push_back(int(cols.size()));}
        if(method==InterfaceMethod::compatibility) {
            // OpenAccel uses edge neighbours, with full parents at interfaces.
            // Omitted coefficients enter the diagonal; the flux residual stays full.
            std::vector<std::set<int64_t>> implicit(n_);
            for(const auto& c:input.cells) for(int k=0;k<c.count;++k) {
                int r=local_row[c.node[k]];if(r<0) continue;
                for(int j=0;j<c.count;++j) {
                    int differences=0;for(int a=0;a<3;++a) differences+=sign(k,a)!=sign(j,a);
                    if(c.count==4 || differences<=1) implicit[r].insert(m.dofs[c.node[j]]);
                }
            }
            for(auto face:input.faces) if(face.kind>=3) {
                const auto& c=input.cells[face.cell];
                for(int k=0;k<c.count;++k) {
                    int r=local_row[c.node[k]];if(r<0) continue;
                    for(int j=0;j<c.count;++j) implicit[r].insert(m.dofs[c.node[j]]);
                }
            }
            for(auto query:input.queries) {
                int r=local_row[query.row];
                for(int k=query.start;k<query.start+query.count;++k) {
                    const auto& c=input.cells[input.faces[input.donors[k]].cell];
                    for(int j=0;j<c.count;++j) implicit[r].insert(m.dofs[c.node[j]]);
                }
            }
            std::vector<int> slots(cols.size());
            bool valid_diagonals=true;
            for(int r=0;r<n_;++r) {
                int diagonal=int(std::lower_bound(cols.begin()+off[r],cols.begin()+off[r+1],rows[r])-cols.begin());
                valid_diagonals=valid_diagonals && diagonal<off[r+1] && cols[diagonal]==rows[r];
                for(int p=off[r];p<off[r+1];++p) slots[p]=implicit[r].count(cols[p])?p:diagonal;
            }
            check(valid_diagonals,"compatibility diagonal column");
            coefficient_slots_=Buffer<int>(slots.begin(),slots.end());
        }
        pressure_.rowIds=Buffer<int64_t>(rows.begin(),rows.end());pressure_.columns=Buffer<int64_t>(cols.begin(),cols.end());pressure_.offsets=Buffer<int>(off.begin(),off.end());pressure_.values.resize(cols.size());momentum_=ncflow::blockGraph(pressure_);
        cells_=Buffer<Cell>(input.cells.begin(),input.cells.end());faces_=Buffer<Face>(input.faces.begin(),input.faces.end());queries_=Buffer<Query>(input.queries.begin(),input.queries.end());donors_=Buffer<int>(input.donors.begin(),input.donors.end());
        std::vector<double> xyz;for(auto x:m.x) {xyz.push_back(x.x);xyz.push_back(x.y);xyz.push_back(x.z);}xyz_=Buffer<double>(xyz.begin(),xyz.end());ids_=Buffer<int64_t>(m.dofs.begin(),m.dofs.end());owned_=Buffer<int>(owned.begin(),owned.end());
        ncflow::incidence(by_cell,co_,ci_);ncflow::incidence(by_face,fo_,fi_);ncflow::incidence(by_query,qo_,qi_);
        geometry_.resize(cells_.size());interfaces_.resize(queries_.size());state_.resize(m.x.size());cell_mass_.resize(12*cells_.size());face_mass_.resize(4*faces_.size());interface_mass_.resize(queries_.size());closed_.resize(faces_.size());
        view_={pointer(cells_),pointer(geometry_),pointer(faces_),pointer(queries_),pointer(donors_),pointer(xyz_),pointer(ids_),pointer(owned_),pointer(co_),pointer(ci_),pointer(fo_),pointer(fi_),pointer(qo_),pointer(qi_),pointer(interfaces_),pointer(state_),pointer(cell_mass_),pointer(face_mass_),pointer(interface_mass_),pointer(closed_),Controls{},0};
        view_.openaccel_profile=profile==FlowProfile::openaccel;
        if(method==InterfaceMethod::compatibility) view_.coefficient_slots=pointer(coefficient_slots_);
        std::vector<uint64_t> halo_ids(m.dofs.begin(),m.dofs.end());std::vector<int> owners;for(auto id:m.dofs) owners.push_back(owner(id,m.globalNodes));check(halo_.buildFromHost(halo_ids,owners,m.candidatePeers),"mixed state halo");
        Buffer<int> flags(cells_.size());each<device>(int(cells_.size()),BuildGeometry{view_,pointer(geometry_),pointer(flags)});ncflow::finish();ncflow::validFlags(flags,"mixed cell geometry");
        if(method==InterfaceMethod::conservative) build_conservative();
        else {
            flags.resize(queries_.size());each<device>(int(queries_.size()),BuildInterface{view_,pointer(flags),method==InterfaceMethod::legacy?DonorSelection::persistent_id:DonorSelection::openaccel_face_score});ncflow::finish();ncflow::validFlags(flags,"mixed donor search/geometry");
        }
        each<device>(n_,Gradient{view_,0,false});ncflow::finish();auto initial=reduce(n_,Statistics{view_});check(std::isfinite(initial.sums[6]) && initial.sums[6]>0,"mixed domain volume");
        exchange();each<device>(int(faces_.size()),UpdateBoundaryMass{view_});each<device>(n_,Divergence{view_});ncflow::finish();
        if(rankId==0) std::printf("Native mixed geometry and interface setup PASS, volume %.17g\n",initial.sums[6]);
    }
    void exchange() {
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(halo_.forward(pointer(state_),state_.size()),"mixed state exchange");
#else
        check(halo_.forwardHost(state_),"mixed state exchange");
#endif
    }
    void assemble(bool momentum,ChannelOps& ops) {
        prepare_interface_flux();
        Buffer<int> flags(n_);auto matrix=pressure_.view();if(momentum) matrix.values=pointer(momentum_.values);
        each<device>(n_,Assemble{view_,matrix,ops.rhs(),pointer(flags),momentum});ncflow::finish();ncflow::validFlags(flags,"mixed assembly graph/coefficients");
        auto& graph=momentum?momentum_:pressure_;flags.resize(ops.n);each<device>(ops.n,ncflow::Jacobi{graph.view(),pointer(ops.diagonal),pointer(flags)});ncflow::finish();ncflow::validFlags(flags,"mixed diagonal");
    }
    bool solve(int budget,double tolerance,const std::string& output);
    void build_conservative();
    void prepare_interface_flux() {
        if(!view_.conservative) return;
        each<device>(int(interfaces_.size()),PrepareInterfaceFlux{view_});ncflow::finish();
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(flux_halo_.forward(pointer(interface_fluxes_),interface_fluxes_.size()),"shared interface flux exchange");
#else
        check(flux_halo_.forwardHost(interface_fluxes_),"shared interface flux exchange");
#endif
    }
    void update_interface_mass() {
        each<device>(int(interfaces_.size()),UpdateInterfaceMass{view_});ncflow::finish();
        if(!view_.conservative) return;
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(mass_halo_.forward(pointer(interface_mass_),interface_mass_.size()),"shared interface mass exchange");
#else
        check(mass_halo_.forwardHost(interface_mass_),"shared interface mass exchange");
#endif
    }
    void write(const std::string& output) {
        if(output.empty()) return;auto values=host<State>(state_);auto owned=host<int>(owned_);std::string path=output+".rank"+std::to_string(rankId)+".csv";
        if(std::filesystem::exists(path)) throw std::runtime_error("mixed output already exists");std::ofstream f(path);f.precision(17);f<<"node_id,x,y,z,velocity_x,velocity_y,velocity_z,pressure\n";
        for(int i:owned) {auto x=input_.mesh.x[i];auto s=values[i];f<<input_.mesh.dofs[i]+1<<','<<x.x<<','<<x.y<<','<<x.z<<','<<s.u[0]<<','<<s.u[1]<<','<<s.u[2]<<','<<s.p<<'\n';}
        f.close();if(!f) throw std::runtime_error("mixed field write failed");
    }
    const Input& input_;int n_;NonconformalDonorExchange<double,State> halo_;Graph pressure_,momentum_;
    NonconformalDonorExchange<double,Flux> flux_halo_;NonconformalDonorExchange<double,double> mass_halo_;
    Buffer<Flux> interface_fluxes_;Buffer<int> patch_owner_,qs_;
    Buffer<Cell> cells_;Buffer<Geometry> geometry_;Buffer<Face> faces_;Buffer<Query> queries_;Buffer<Interface> interfaces_;Buffer<State> state_;
    Buffer<double> xyz_,cell_mass_,face_mass_,interface_mass_;Buffer<int64_t> ids_;
    Buffer<int> donors_,owned_,co_,ci_,fo_,fi_,qo_,qi_,closed_,coefficient_slots_;View view_;
};
}
#include "mars_mixed_conservative_runtime.hpp"
