#pragma once
#include "mars_mixed_flow.hpp"
#include "mars_nonconformal_overlap.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_MIXED_HD __host__ __device__
#else
#define MARS_MIXED_HD
#endif

namespace mars::mixed {
struct OverlapOptions {
    double distance_tolerance=1e-10; // Physical distance, chosen for the mesh units.
    double normal_tolerance=1e-8;
};
enum class OverlapStatus { valid, invalid_input, invalid_geometry };
struct OverlapPatch : Interface {
    int opposite_row{};
    int topology_index{};
    double position[3]{};
};
inline constexpr int overlap_pair_capacity=64;

namespace overlap_detail {
using Point2=fem::nonconformal_overlap_detail::Point2<double>;
using fem::nonconformal_overlap_detail::cross2;
using fem::nonconformal_overlap_detail::overlap;

struct Face {
    double vertex[12]{},center[3]{},normal[3]{};
    int count{};
};
MARS_MIXED_HD inline bool make_face(const Cell& cell,const double* xyz,int ordinal,
                                    const OverlapOptions& options,Face& out) {
    if((cell.count!=4 && cell.count!=8) || ordinal<0 || ordinal>=face_count(cell.count)) return false;
    Geometry geometry_check;
    if(!geometry(cell.count,xyz,geometry_check,true)) return false;
    out={};out.count=face_size(cell.count);
    for(int k=0;k<cell.count;++k) for(int j=0;j<3;++j)
        if(!segregated::geometry_finite(xyz[3*k+j])) return false;
    for(int k=0;k<out.count;++k) for(int j=0;j<3;++j) {
        out.vertex[3*k+j]=xyz[3*face_node(cell.count,ordinal,k)+j];
        out.center[j]+=out.vertex[3*k+j]/out.count;
    }
    for(int k=1;k+1<out.count;++k) {
        double a[3],b[3],c[3];
        for(int j=0;j<3;++j) {a[j]=out.vertex[3*k+j]-out.vertex[j];b[j]=out.vertex[3*(k+1)+j]-out.vertex[j];}
        cross(a,b,c);for(int j=0;j<3;++j) out.normal[j]+=c[j];
    }
    double norm=sqrt(dot(out.normal,out.normal));
    if(!(norm>0) || !segregated::geometry_finite(norm)) return false;
    for(int j=0;j<3;++j) out.normal[j]/=norm;
    for(int k=0;k<out.count;++k) {
        double delta[3];for(int j=0;j<3;++j) delta[j]=out.vertex[3*k+j]-out.center[j];
        if(fabs(dot(delta,out.normal))>options.distance_tolerance) return false;
    }
    // Reject folded or concave quads: a fan would otherwise double count area.
    for(int k=0;k<out.count;++k) {
        double a[3],b[3],c[3];
        for(int j=0;j<3;++j) {
            a[j]=out.vertex[3*((k+1)%out.count)+j]-out.vertex[3*k+j];
            b[j]=out.vertex[3*((k+2)%out.count)+j]-out.vertex[3*((k+1)%out.count)+j];
        }
        cross(a,b,c);
        if(!(dot(c,out.normal)>0)) return false;
    }
    double outward[3];
    for(int j=0;j<3;++j) {
        outward[j]=out.center[j];
        for(int k=0;k<cell.count;++k) outward[j]-=xyz[3*k+j]/cell.count;
    }
    return dot(outward,out.normal)>0;
}
MARS_MIXED_HD inline void triangle(const Point2* face,int count,int row,int edge,Point2* out) {
    const int neighbor=(row+(edge?count-1:1))%count;
    out[0]=face[row];out[1]={(face[row].x+face[neighbor].x)/2,(face[row].y+face[neighbor].y)/2};
    out[2]={0,0};for(int k=0;k<count;++k) {out[2].x+=face[k].x/count;out[2].y+=face[k].y/count;}
    if(cross2(out[0],out[1],out[2])<0) {auto p=out[1];out[1]=out[2];out[2]=p;}
}
}

// One quadrature point per common median-dual overlap integrates constants and
// physical affine fields exactly. Nonlinear fluxes retain centroid quadrature.
// Emit once and scatter the same flux to row and opposite_row with opposite
// signs. Caller supplies each face pair once and discards emissions on failure.
template<class Emit>
MARS_MIXED_HD inline OverlapStatus overlap_pair(const Cell& a,const double* xa,int fa,
                                               const Cell& b,const double* xb,int fb,
                                               const OverlapOptions& options,Emit emit) {
    using namespace overlap_detail;
    if(!xa || !xb || (a.count!=4 && a.count!=8) || (b.count!=4 && b.count!=8) ||
       fa<0 || fa>=face_count(a.count) || fb<0 || fb>=face_count(b.count) ||
       !segregated::geometry_finite(options.distance_tolerance) || options.distance_tolerance<0 ||
       !segregated::geometry_finite(options.normal_tolerance) || options.normal_tolerance<0 || options.normal_tolerance>=1)
        return OverlapStatus::invalid_input;
    for(int k=0;k<a.count;++k) {
        if(a.node[k]<0) return OverlapStatus::invalid_input;
        for(int j=0;j<k;++j) if(a.node[k]==a.node[j]) return OverlapStatus::invalid_input;
        for(int j=0;j<b.count;++j) if(a.node[k]==b.node[j]) return OverlapStatus::invalid_input;
    }
    for(int k=0;k<b.count;++k) {
        if(b.node[k]<0) return OverlapStatus::invalid_input;
        for(int j=0;j<k;++j) if(b.node[k]==b.node[j]) return OverlapStatus::invalid_input;
    }
    Face ga,gb;
    if(!make_face(a,xa,fa,options,ga) || !make_face(b,xb,fb,options,gb) ||
       1+dot(ga.normal,gb.normal)>options.normal_tolerance) return OverlapStatus::invalid_geometry;
    double u[3],v[3];for(int j=0;j<3;++j) u[j]=ga.vertex[3+j]-ga.vertex[j];
    double length=sqrt(dot(u,u));for(int j=0;j<3;++j) u[j]/=length;
    cross(ga.normal,u,v);
    Point2 pa[4],pb[4];
    for(int side=0;side<2;++side) {
        const auto& g=side?gb:ga;auto* projected=side?pb:pa;
        for(int k=0;k<g.count;++k) {
            double delta[3];for(int j=0;j<3;++j) delta[j]=g.vertex[3*k+j]-ga.vertex[j];
            if(fabs(dot(delta,ga.normal))>options.distance_tolerance) return OverlapStatus::invalid_geometry;
            projected[k]={dot(delta,u),dot(delta,v)};
        }
    }
    for(int ra=0;ra<ga.count;++ra) for(int ea=0;ea<2;++ea)
        for(int rb=0;rb<gb.count;++rb) for(int eb=0;eb<2;++eb) {
            Point2 ta[3],tb[3],center;double area;
            triangle(pa,ga.count,ra,ea,ta);triangle(pb,gb.count,rb,eb,tb);
            if(!overlap(ta,tb,area,center)) return OverlapStatus::invalid_geometry;
            if(area==0) continue;
            OverlapPatch p;
            p.topology_index=(ra*2+ea)*8+rb*2+eb;
            p.row=a.node[face_node(a.count,fa,ra)];p.opposite_row=b.node[face_node(b.count,fb,rb)];
            p.count[0]=a.count;p.count[1]=b.count;
            for(int j=0;j<3;++j) {p.position[j]=ga.vertex[j]+center.x*u[j]+center.y*v[j];p.area[j]=area*ga.normal[j];}
            for(int side=0;side<2;++side) {
                const auto& c=side?b:a;const auto* xyz=side?xb:xa;int face=side?fb:fa,offset=side?a.count:0;
                Point point;if(!invert(c.count,xyz,p.position,point)) return OverlapStatus::invalid_geometry;
                for(int k=0;k<c.count;++k) {
                    p.node[offset+k]=c.node[k];p.shape[offset+k]=point.shape[k];
                    for(int j=0;j<3;++j) p.gradient[3*(offset+k)+j]=point.gradient[3*k+j];
                }
                for(int k=0;k<face_size(c.count);++k) {
                    int node=face_node(c.count,face,k);p.face_mask[offset+node]=1;
                    p.inverse_height[side]+=(side?-1.:1.)*dot(point.gradient+3*node,ga.normal);
                }
                if(!(p.inverse_height[side]>0) || !segregated::geometry_finite(p.inverse_height[side]))
                    return OverlapStatus::invalid_geometry;
            }
            emit(p);
        }
    return OverlapStatus::valid;
}
}
#undef MARS_MIXED_HD
