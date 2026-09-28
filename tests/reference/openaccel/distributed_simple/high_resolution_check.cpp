#include "mars_segregated_high_resolution.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
using namespace mars::segregated;
int checks=0;
void check(bool ok) { ++checks; if (!ok) throw std::runtime_error("high-resolution check "+std::to_string(checks)); }
void near(double a,double b) { check(std::isfinite(a) && std::abs(a-b)<3e-13*std::max(1.,std::abs(b))); }

// Literal two-sided reference expression, independent of the simplified production ratio.
double reference_bound(double u,double lo,double hi,double projection) {
    constexpr double eps=std::numeric_limits<double>::epsilon();
    const double r=projection+(projection>0?eps:-eps);
    if (std::abs(r)<=100*eps) return 1;
    const double plus=std::max((hi-u)/r,(lo-u)/r);
    const double minus=std::max(-(hi-u)/r,-(lo-u)/r);
    const double y=std::min(plus,minus);
    return std::min(1.,(y*y+2*y)/(y*y+y+2));
}
int main() {
    near(velocity_blend_bound(0,-1,1,1),.75);
    near(velocity_blend_bound(0,0,2,1),0);
    near(velocity_blend_bound(0,0,2,-1),0);
    near(velocity_blend_bound(3,3,3,0),1);
    near(relax_velocity_blend(.75,0),.1875);
    near(relax_velocity_blend(.2,.9),.2);
    for (double p:{0.,1e-20,98*std::numeric_limits<double>::epsilon(),100*std::numeric_limits<double>::epsilon(),.1,1.,10.,1e20})
        for (double sign:{-1.,1.}) near(velocity_blend_bound(.3,-2,1,sign*p),reference_bound(.3,-2,1,sign*p));

    int n0[]={0,1},n1[]={1,2},n2[]={2,3},n3[]={3,4},error=0;
    double x[]={0,2,0,0,2},y[]={0,0,1,0,2},z[]={0,0,0,3,3};
    TetGeometry<double> geometry[2]; SimpleFace faces[]={{0,0,0},{0,2,2},{0,3,2},{1,0,1},{1,1,2},{1,2,2}};
    SimpleMesh mesh{5,2,6,{n0,n1,n2,n3},x,y,z,faces,geometry};
    int offsets[]={0,4,9,14,19,23},columns[]={0,1,2,3, 0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4, 1,2,3,4};
    double u[15],vg[45],lo[15],hi[15],candidate[15],beta[15]{};
    SimpleState state{}; state.velocity=u; state.velocity_gradient=vg; state.error=&error;
    std::mt19937 rng(19); std::uniform_real_distribution<double> random(-2,2);
    int boundary_dominant=0;
    for (int trial=0;trial<40;++trial) {
        for (double& v:u) v=random(rng);
        for (double& v:vg) v=3*random(rng);
        SimpleBlendBounds bounds{{5,offsets,columns,nullptr,nullptr},u,lo,hi,candidate,&error};
        for (int n=0;n<5;++n) bounds(n);
        double expected[15]; std::fill(expected,expected+15,1.);
        for (int n=0;n<5;++n) for (int c=0;c<3;++c) {
            double lower=u[3*n+c],upper=lower;
            for (int e=0;e<2;++e) {
                bool touches=false; for (int k=0;k<4;++k) touches|=mesh.nodes[k][e]==n;
                if (touches) for (int k=0;k<4;++k) { lower=std::min(lower,u[3*mesh.nodes[k][e]+c]); upper=std::max(upper,u[3*mesh.nodes[k][e]+c]); }
            }
            near(lo[3*n+c],lower); near(hi[3*n+c],upper);
        }
        auto constrain=[&](int n,const double* p) {
            for (int c=0;c<3;++c) {
                const double r=vg[9*n+3*c]*(p[0]-x[n])+vg[9*n+3*c+1]*(p[1]-y[n])+vg[9*n+3*c+2]*(p[2]-z[n]);
                expected[3*n+c]=std::min(expected[3*n+c],reference_bound(u[3*n+c],lo[3*n+c],hi[3*n+c],r));
            }
        };
        // Enumerate all unordered edges independently of MARS's sample numbering.
        for (int e=0;e<2;++e) for (int a=0;a<4;++a) for (int b=a+1;b<4;++b) {
            double p[3]{};
            for (int k=0;k<4;++k) {
                const int n=mesh.nodes[k][e]; const double w=(k==a || k==b)?13./36:5./36;
                p[0]+=w*x[n]; p[1]+=w*y[n]; p[2]+=w*z[n];
            }
            constrain(mesh.nodes[a][e],p); constrain(mesh.nodes[b][e],p);
        }
        const SimpleBlendSamples samples{mesh,state,lo,hi,candidate};
        for (int e=0;e<2;++e) SimpleBlendInterior{samples}(e);
        for (int k=0;k<15;++k) near(candidate[k],expected[k]);
        const std::vector<double> interior(candidate,candidate+15);
        for (const auto& f:faces) for (int s=0;s<3;++s) {
            double p[3]{};
            for (int k=0;k<3;++k) {
                const int n=mesh.nodes[tet_face_node(f.ordinal,k)][f.element]; const double w=k==s?11./18:7./36;
                p[0]+=w*x[n]; p[1]+=w*y[n]; p[2]+=w*z[n];
            }
            constrain(mesh.nodes[tet_face_node(f.ordinal,s)][f.element],p);
        }
        for (int f=0;f<6;++f) SimpleBlendBoundary{samples}(f);
        const std::vector<double> old(beta,beta+15);
        for (int n=0;n<5;++n) SimpleBlendFinish{candidate,beta}(n);
        for (int k=0;k<15;++k) {
            near(candidate[k],expected[k]); boundary_dominant+=candidate[k]<interior[k]-1e-12;
            near(beta[k],std::min(expected[k],.75*old[k]+.25*expected[k]));
            check(beta[k]>=0 && beta[k]<=candidate[k] && candidate[k]<=1);
        }
    }
    check(boundary_dominant>0 && error==0);
    // Reconstruction changes equal-and-opposite residuals, never the upwind Jacobian.
    int nodes[]={0,1,2,3}; double xyz[12]; mesh.cell(0,nodes,xyz);
    check(tet_geometry(xyz,geometry[0])); double p[5]{},pg[15]{},influence[15]{},flux[]={-2,3,-4,5,-6,7};
    auto base=simple_interior(1,nodes,xyz,u,p,flux,{});
    native_interior(base,geometry[0],nodes,vg,pg,influence);
    auto high=base; for (double& v:high.velocity_blend) v=.6;
    TetInteriorOutput a,b; tet_interior(base,a); tet_interior(high,b);
    for (int k=0;k<144;++k) near(a.lhs[k],b.lhs[k]);
    double change=0;
    for (int c=0;c<3;++c) { double sum=0; for (int n=0;n<4;++n) { sum+=b.rhs[3*n+c]-a.rhs[3*n+c]; change+=std::abs(b.rhs[3*n+c]-a.rhs[3*n+c]); } near(sum,0); }
    check(change>1e-3);
    std::cout<<"PASS: "<<checks<<" high-resolution bounds, history, boundary and conservative flux checks\n";
}
