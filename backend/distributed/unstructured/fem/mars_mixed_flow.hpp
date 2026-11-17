#pragma once
#include "mars_mixed_cvfem.hpp"
#include "segregated/mars_segregated_node.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_MIXED_HD __host__ __device__
#else
#define MARS_MIXED_HD
#endif
namespace mars::mixed {
inline constexpr int interface_capacity=16;
struct Controls {
    double rho=1000,mu=.001,dt=1,speed=.3,alpha_u=1,alpha_p=1,alpha_mass=1,alpha_gradient=.5625;
    bool high_resolution=true;
};
struct State {
    double u[3]{},p{},gp[3]{},gu[9]{},d[3]{},blend[3]{},volume{},divergence{},increment{},gc[3]{};
};
struct Flux {
    double mass{},pressure[interface_capacity]{},momentum[3]{},coefficient[9*interface_capacity]{};
    MARS_MIXED_HD double& a(int k,int i,int j) { return coefficient[9*k+3*i+j]; }
};
struct Interface {
    int count[2]{},node[interface_capacity]{},face_mask[interface_capacity]{},row{};
    double shape[interface_capacity]{},gradient[3*interface_capacity]{},area[3]{},inverse_height[2]{};
};
MARS_MIXED_HD inline double inlet(const double* xyz,const Controls& c) {
    return 16*c.speed*xyz[1]*(1-xyz[1])*xyz[2]*(1-xyz[2]);
}
// Constant-density Rhie--Chow flux, with the same unrelaxed pressure derivative.
MARS_MIXED_HD inline Flux interior_flux(const Cell& cell,const Point& p,const State* state,const double* xyz,
                                      const Controls& c,double advecting_mass) {
    Flux out;double velocity[3]{},influence[3]{},compact[3]{};
    for(int k=0;k<cell.count;++k) {
        const auto& s=state[cell.node[k]];
        for(int j=0;j<3;++j) {
            velocity[j]+=p.shape[k]*s.u[j];influence[j]+=p.shape[k]*s.d[j];compact[j]+=p.gradient[3*k+j]*s.p;
        }
    }
    for(int j=0;j<3;++j) {
        const double reconstructed=.5*(state[cell.node[p.left]].gp[j]+state[cell.node[p.right]].gp[j]);
        out.mass+=c.rho*(velocity[j]-influence[j]*(compact[j]-reconstructed))*p.area[j];
    }
    const int up=advecting_mass>0?p.left:p.right;const auto& su=state[cell.node[up]];
    for(int i=0;i<3;++i) {
        double correction=0;
        if(c.high_resolution)
            for(int j=0;j<3;++j) correction+=(p.position[j]-xyz[3*cell.node[up]+j])*su.gu[3*i+j];
        out.momentum[i]=advecting_mass*su.blend[i]*correction;
    }
    for(int k=0;k<cell.count;++k) {
        double normal_gradient=dot(p.gradient+3*k,p.area);
        for(int j=0;j<3;++j) out.pressure[k]-=c.rho*influence[j]*p.gradient[3*k+j]*p.area[j];
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
            double a=-c.mu*((i==j?normal_gradient:0)+p.gradient[3*k+i]*p.area[j]);
            if(k==up && i==j) a+=advecting_mass;
            out.a(k,i,j)=a;out.momentum[i]+=a*state[cell.node[k]].u[j];
        }
    }
    return out;
}
// The caller chooses independent-side quadrature for OpenAccel compatibility,
// or a common overlap quadrature with equal-and-opposite row contributions.
MARS_MIXED_HD inline Flux interface_flux(const Interface& p,const State* state,const Controls& c,double advecting_mass) {
    Flux out;double velocity[2][3]{},d[2][3]{},compact[2][3]{},g[2][3]{},pressure[2]{};
    const int total=p.count[0]+p.count[1];int face_nodes[2]{};
    for(int k=0;k<total;++k) {
        int side=k<p.count[0]?0:1;const auto& s=state[p.node[k]];
        pressure[side]+=p.shape[k]*s.p;face_nodes[side]+=p.face_mask[k];
        for(int j=0;j<3;++j) {
            velocity[side][j]+=p.shape[k]*s.u[j];d[side][j]+=p.shape[k]*s.d[j];
            compact[side][j]+=p.gradient[3*k+j]*s.p;g[side][j]+=p.face_mask[k]*s.gp[j];
        }
    }
    double penalty=0;
    for(int side=0;side<2;++side) for(int j=0;j<3;++j) penalty+=c.rho*d[side][j]*p.inverse_height[side]/6;
    const double area=sqrt(dot(p.area,p.area));
    out.mass=area*penalty*(pressure[0]-pressure[1]);
    for(int side=0;side<2;++side) for(int j=0;j<3;++j)
        out.mass+=.5*c.rho*(velocity[side][j]-d[side][j]*(compact[side][j]-g[side][j]/face_nodes[side]))*p.area[j];
    const double viscous=.5*c.mu*(p.inverse_height[0]+p.inverse_height[1])*area;
    for(int k=0;k<total;++k) {
        int side=k<p.count[0]?0:1;double sign=side==0?1.:-1.;
        out.pressure[k]=sign*penalty*area*p.shape[k];
        for(int j=0;j<3;++j) out.pressure[k]-=.5*c.rho*d[side][j]*p.gradient[3*k+j]*p.area[j];
        const double advection=side==0?fmax(advecting_mass,0.):fmin(advecting_mass,0.);
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
            double a=-.5*c.mu*((i==j?dot(p.gradient+3*k,p.area):0)+p.gradient[3*k+i]*p.area[j]);
            if(i==j) a+=(sign*viscous+advection)*p.shape[k];
            out.a(k,i,j)=a;out.momentum[i]+=a*state[p.node[k]].u[j];
        }
    }
    return out;
}
MARS_MIXED_HD inline Flux boundary_flux(const Cell& cell,int face,const Point& p,const State* state,const double* xyz,
                                      const Controls& c,int kind,double advecting_mass,bool reversed,bool initialized) {
    Flux out;int mask[8]{};const int nf=face_size(cell.count);
    for(int q=0;q<nf;++q) mask[face_node(cell.count,face,q)]=1;
    if(reversed) return out;
    if(kind!=1) {
        if(kind==0) for(int k=0;k<cell.count;++k) out.mass+=c.rho*p.shape[k]*inlet(xyz+3*cell.node[k],c)*p.area[0];
    } else {
        double u[3]{},d[3]{},compact[3]{};int opposite=0;
        if(cell.count==4) opposite=segregated::tet_opposite_node(face);
        else {
            for(int k=0;k<cell.count;++k) if(!mask[k]) {
                int differences=0;for(int j=0;j<3;++j) differences+=sign(k,j)!=sign(p.left,j);
                if(differences==1) opposite=k;
            }
        }
        for(int k=0;k<cell.count;++k) for(int j=0;j<3;++j) {
            const auto& s=state[cell.node[k]];u[j]+=p.shape[k]*s.u[j];d[j]+=p.shape[k]*s.d[j];
            compact[j]+=p.gradient[3*k+j]*(mask[k]?0.:s.p);
        }
        for(int j=0;j<3;++j) {
            const double g=.5*(state[cell.node[p.left]].gp[j]+state[cell.node[opposite]].gp[j]);
            out.mass+=c.rho*(u[j]-d[j]*(compact[j]-g))*p.area[j];
            for(int k=0;k<cell.count;++k) if(!mask[k]) out.pressure[k]-=c.rho*d[j]*p.gradient[3*k+j]*p.area[j];
        }
    }
    const double area=sqrt(dot(p.area,p.area));double normal[3];for(int j=0;j<3;++j) normal[j]=p.area[j]/area;
    if(kind==2) {
        double inverse_height=0;
        for(int k=0;k<cell.count;++k) if(mask[k]) inverse_height+=dot(p.gradient+3*k,normal);
        const double coefficient=initialized?4*c.mu*area*inverse_height:0;
        for(int k=0;k<cell.count;++k) for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
            double a=coefficient*((i==j?1.:0)-normal[i]*normal[j])*p.shape[k];
            out.a(k,i,j)=a;out.momentum[i]+=a*state[cell.node[k]].u[j];
        }
        return out;
    }
    for(int i=0;i<3;++i) {
        double u=kind==1?state[cell.node[p.left]].u[i]:0;
        if(kind==0 && i==0) for(int k=0;k<cell.count;++k) u+=p.shape[k]*inlet(xyz+3*cell.node[k],c);
        out.momentum[i]+=advecting_mass*u;
        if(kind==1) out.a(p.left,i,i)+=advecting_mass;
        for(int k=0;k<cell.count;++k) for(int j=0;j<3;++j) {
            double a=0;
            for(int l=0;l<3;++l) {
                double projection=(i==l?1.:0)-(kind==1?normal[i]*normal[l]:0);
                a-=c.mu*projection*((l==j?dot(p.gradient+3*k,p.area):0)+p.gradient[3*k+l]*p.area[j]);
            }
            double value=state[cell.node[k]].u[j];
            if(kind==0 && mask[k]) value=j==0?inlet(xyz+3*cell.node[k],c):0;
            out.momentum[i]+=a*value;
            if(kind!=0 || !mask[k]) out.a(k,i,j)+=a;
        }
    }
    return out;
}
}
#undef MARS_MIXED_HD
