#include "backend/distributed/unstructured/fem/mars_mixed_flow.hpp"
#include "backend/distributed/unstructured/fem/segregated/mars_segregated_tet_interior.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>

using namespace mars::mixed;
static int checks=0;
static void require(bool ok,const char* reason) {
    ++checks;if(!ok) {std::fprintf(stderr,"FAIL: %s\n",reason);std::exit(1);}
}
static void near(double a,double b,double tolerance,const char* reason) {
    require(std::isfinite(a) && std::isfinite(b) && fabs(a-b)<=tolerance*std::max({1.,fabs(a),fabs(b)}),reason);
}
static void coordinates(int count,double* x,bool warped=false) {
    if(count==4) {double t[12]={0,0,0,1,0,0,0,1,0,0,0,1};std::copy(t,t+12,x);return;}
    for(int k=0;k<8;++k) for(int j=0;j<3;++j) x[3*k+j]=.5*(1+sign(k,j));
    if(warped) {x[18]+=.13;x[19]-=.07;x[20]+=.11;x[14]+=.06;}
}
static void geometry_checks(int count,bool warped) {
    double x[24]{};coordinates(count,x,warped);Geometry g;
    require(geometry(count,x,g),"positive geometry");double total=0;
    for(int k=0;k<count;++k) {require(g.volume[k]>0,"positive dual volume");total+=g.volume[k];}
    if(!warped) near(total,count==8?1.:1./6,1e-14,"element volume");
    double closure[8][3]{};
    for(int q=0;q<g.samples;++q) {
        const auto& p=g.interior[q];double sum=0;for(int k=0;k<count;++k) sum+=p.shape[k];near(sum,1,1e-14,"shape partition");
        for(int j=0;j<3;++j) {
            closure[p.left][j]+=p.area[j];closure[p.right][j]-=p.area[j];double grad=0;
            for(int k=0;k<count;++k) grad+=p.gradient[3*k+j];near(grad,0,1e-13,"constant derivative");
            for(int a=0;a<3;++a) {double linear=0;for(int k=0;k<count;++k) linear+=p.gradient[3*k+j]*x[3*k+a];near(linear,a==j,1e-13,"linear physical derivative");}
            Point plus,minus;double xp[3],xm[3];std::copy(p.position,p.position+3,xp);std::copy(xp,xp+3,xm);xp[j]+=1e-6;xm[j]-=1e-6;
            require(invert(count,x,xp,plus)&&invert(count,x,xm,minus),"inverse map");
            for(int k=0;k<count;++k) near((plus.shape[k]-minus.shape[k])/2e-6,p.gradient[3*k+j],2e-9,"physical gradient finite difference");
        }
    }
    for(int f=0;f<face_count(count);++f) for(int q=0;q<face_size(count);++q) {
        Point p,shifted;require(face_point(count,x,f,q,false,p)&&face_point(count,x,f,q,true,shifted),"boundary geometry");
        for(int j=0;j<3;++j) {closure[p.left][j]+=p.area[j];near(p.area[j],shifted.area[j],1e-14,"shift preserves area");near(shifted.position[j],x[3*p.left+j],1e-14,"shifted vertex");}
    }
    for(int k=0;k<count;++k) for(int j=0;j<3;++j) near(closure[k][j],0,2e-14,"dual-cell area closure");
    for(int k=0;k<count;++k) x[3*k]=-x[3*k];require(!geometry(count,x,g),"inverted geometry rejected");
}
template<class Function> static void jacobian_checks(State* state,int count,Function flux) {
    const auto base=flux();constexpr double eps=1e-5;
    for(int k=0;k<count;++k) {
        state[k].p+=eps;auto plus=flux();state[k].p-=2*eps;auto minus=flux();state[k].p+=eps;
        near((plus.mass-minus.mass)/(2*eps),base.pressure[k],2e-8,"pressure flux derivative");
        for(int j=0;j<3;++j) {
            state[k].u[j]+=eps;plus=flux();state[k].u[j]-=2*eps;minus=flux();state[k].u[j]+=eps;
            for(int i=0;i<3;++i) near((plus.momentum[i]-minus.momentum[i])/(2*eps),base.coefficient[9*k+3*i+j],2e-8,"momentum flux derivative");
        }
    }
}
static void initialize(State* state,int count) {
    for(int k=0;k<count;++k) {
        state[k].p=.1*k;
        for(int j=0;j<3;++j) {state[k].u[j]=.01*(k+1)*(j+1);state[k].gp[j]=.2+.03*k*j;state[k].d[j]=.003+.0002*k*j;state[k].blend[j]=.2+.03*k;
            for(int a=0;a<3;++a) state[k].gu[3*j+a]=.02*(j+1)*(a+1);}
    }
}
static void flux_checks(int count) {
    double x[24]{};coordinates(count,x);Cell c;c.count=count;for(int k=0;k<count;++k) c.node[k]=k;
    Geometry g;geometry(count,x,g);State state[12];initialize(state,count);Controls controls;
    // With viscosity zero, only the chosen upwind value and reconstruction remain.
    for(double mass:{-2.,0.,3.}) for(int q=0;q<g.samples;++q) {
        const auto& point=g.interior[q];Controls inviscid;inviscid.mu=0;
        const auto high=interior_flux(c,point,state,x,inviscid,mass);
        inviscid.high_resolution=false;
        const auto low=interior_flux(c,point,state,x,inviscid,mass);
        const int up=mass>0?point.left:point.right;
        for(int i=0;i<3;++i) {
            double reconstruction=0;
            for(int j=0;j<3;++j) reconstruction+=(point.position[j]-x[3*up+j])*state[up].gu[3*i+j];
            near(low.momentum[i],mass*state[up].u[i],1e-13,"upwind donor advection");
            near(high.momentum[i]-low.momentum[i],mass*state[up].blend[i]*reconstruction,1e-13,"high-resolution correction");
        }
        near(high.mass,low.mass,1e-13,"advection choice preserves mass flux");
        for(int k=0;k<count;++k) for(int i=0;i<9;++i)
            near(high.coefficient[9*k+i],low.coefficient[9*k+i],1e-13,"advection choice preserves implicit Jacobian");
    }
    for(double mass:{-2.,0.,3.}) {
        for(int q=0;q<g.samples;++q) jacobian_checks(state,count,[&]{return interior_flux(c,g.interior[q],state,x,controls,mass);});
        for(int kind=0;kind<3;++kind) for(int f=0;f<face_count(count);++f) {
            Point p;face_point(count,x,f,0,false,p);
            jacobian_checks(state,count,[&]{return boundary_flux(c,f,p,state,x,controls,kind,mass,false,true);});
        }
    }
    if(count!=4) return;
    // Compare against the previously validated tetrahedral volume algebra.
    mars::segregated::TetInteriorInput in{};
    std::copy(x,x+12,in.coordinates);
    for(int k=0;k<4;++k) {
        in.density[k]=controls.rho;in.viscosity[k]=controls.mu;in.pressure[k]=state[k].p;
        for(int j=0;j<3;++j) {int i=3*k+j;in.velocity[i]=state[k].u[j];in.pressure_gradient[i]=state[k].gp[j];in.influence_lhs[i]=in.influence_rhs[i]=state[k].d[j];in.velocity_blend[i]=state[k].blend[j];}
        std::copy(state[k].gu,state[k].gu+9,in.velocity_gradient+9*k);
    }
    for(int q=0;q<6;++q) {auto p=g.interior[q];in.edges[2*q]=p.left;in.edges[2*q+1]=p.right;in.stored_flux[q]=(q%2?-.3:.4);
        std::copy(p.shape,p.shape+4,in.velocity_shape+4*q);std::copy(p.shape,p.shape+4,in.coordinate_shape+4*q);std::copy(p.gradient,p.gradient+12,in.shape_gradient+12*q);std::copy(p.area,p.area+3,in.area+3*q);}
    for(int stage=0;stage<2;++stage) {
        in.stage=stage;mars::segregated::TetInteriorOutput expected;mars::segregated::tet_interior(in,expected);double rhs[12]{},lhs[144]{};int width=stage?3:1;
        for(int q=0;q<6;++q) {const auto p=g.interior[q];auto f=interior_flux(c,p,state,x,controls,in.stored_flux[q]);
            if(!stage) near(f.mass,expected.flux[q],1e-13,"tet mass oracle");
            for(int side=0;side<2;++side) {int row=side?p.right:p.left;double s=side?-1.:1.;for(int i=0;i<width;++i) {rhs[width*row+i]-=s*(stage?f.momentum[i]:f.mass);for(int k=0;k<4;++k) for(int j=0;j<width;++j) lhs[(width*row+i)*4*width+width*k+j]+=s*(stage?f.coefficient[9*k+3*i+j]:f.pressure[k]);}}
        }
        for(int i=0;i<4*width;++i) {near(rhs[i],expected.rhs[i],2e-13,"tet residual oracle");for(int j=0;j<4*width;++j) near(lhs[i*4*width+j],expected.lhs[i*4*width+j],2e-13,"tet Jacobian oracle");}
    }
}
static void interface_checks() {
    double hex[24],tet[12];coordinates(8,hex);coordinates(4,tet);for(int k=0;k<4;++k) tet[3*k]+=1;
    Point h,t;double point[3]={1,.2,.3};require(invert(8,hex,point,h)&&invert(4,tet,point,t),"interface inverse maps");
    Interface p;p.count[0]=8;p.count[1]=4;p.area[0]=.025;p.inverse_height[0]=p.inverse_height[1]=1;
    for(int k=0;k<12;++k) {p.node[k]=k;p.shape[k]=k<8?h.shape[k]:t.shape[k-8];p.face_mask[k]=k<8?(hex[3*k]==1):(tet[3*(k-8)]==1);for(int j=0;j<3;++j) p.gradient[3*k+j]=k<8?h.gradient[3*k+j]:t.gradient[3*(k-8)+j];}
    State state[12];initialize(state,12);Controls controls;
    for(double mass:{-2.,0.,3.}) jacobian_checks(state,12,[&]{return interface_flux(p,state,controls,mass);});
    for(int k=0;k<12;++k) {state[k]={};state[k].p=2;state[k].u[0]=.2;for(int j=0;j<3;++j) state[k].d[j]=.001;}
    auto f=interface_flux(p,state,controls,5);near(f.mass,5,1e-13,"constant interface mass");near(f.momentum[0],1,1e-13,"constant interface momentum");
    for(int k=0;k<12;++k) {auto x=k<8?hex+3*k:tet+3*(k-8);state[k].p=2+.3*x[0]-.2*x[1]+.1*x[2];state[k].gp[0]=.3;state[k].gp[1]=-.2;state[k].gp[2]=.1;}
    f=interface_flux(p,state,controls,5);near(f.mass,5,1e-13,"linear pressure Rhie-Chow consistency");
    for(int k=0;k<12;++k) {state[k]={};for(int j=0;j<3;++j) state[k].d[j]=.001;state[k].u[0]=k<8?.2:0;}
    f=interface_flux(p,state,controls,0);near(f.momentum[0],5e-6,1e-14,"positive velocity-jump penalty");
    for(int k=0;k<12;++k) {state[k].u[0]=0;state[k].p=k<8?2:0;}
    f=interface_flux(p,state,controls,0);near(f.mass,.05,1e-14,"positive pressure-jump penalty");
    Interface reverse;reverse.count[0]=4;reverse.count[1]=8;reverse.area[0]=-p.area[0];reverse.inverse_height[0]=reverse.inverse_height[1]=1;
    for(int k=0;k<12;++k) {int old=k<4?k+8:k-4;reverse.node[k]=p.node[old];reverse.shape[k]=p.shape[old];reverse.face_mask[k]=p.face_mask[old];for(int j=0;j<3;++j) reverse.gradient[3*k+j]=p.gradient[3*old+j];}
    near(interface_flux(reverse,state,controls,0).mass,-f.mass,1e-14,"reversing both sides reverses the same-point flux");
}
int main() {
    geometry_checks(4,false);geometry_checks(8,false);geometry_checks(8,true);flux_checks(4);flux_checks(8);interface_checks();
    std::printf("PASS: %d mixed geometry, flux, Jacobian and tetrahedral oracle checks\n",checks);
}
