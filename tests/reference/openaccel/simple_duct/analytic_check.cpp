// Independent checks of duct_analytic.hpp. The Dirichlet Poisson problem has one solution, so
// a finite-difference residual of mu lap(u) + G, the wall values and the symmetries together
// verify the series. The flow-rate quadrature then checks the closed-form K, and the tabulated
// Fanning numbers (Shah & London 1978, Table 45) check both against the literature.
// `analytic_check --table` prints K, G and u at sample points for test_duct_analytic.py --cxx,
// which requires the Python module to agree with this header.
#include "duct_analytic.hpp"
#include <cstdio>
#include <string>
#include <vector>

namespace {
int failures=0;
void expect(bool ok,const char* what,double value,double bound) {
    std::printf("%-58s %.3e (bound %.1e) %s\n",what,value,bound,ok?"ok":"FAIL");
    failures+=!ok;
}
void below(const char* what,double value,double bound) { expect(value<=bound,what,value,bound); }

// Gauss-Legendre nodes/weights on [-1,1] (Newton on P_n).
void gauss(int n,std::vector<double>& x,std::vector<double>& w) {
    x.resize(n); w.resize(n);
    for (int i=0;i<n;++i) {
        double t=std::cos(duct::pi*(i+.75)/(n+.5)), dp=0;
        for (int it=0;it<100;++it) {
            double p0=1, p1=t;
            for (int k=2;k<=n;++k) { const double p2=((2*k-1)*t*p1-(k-1)*p0)/k; p0=p1; p1=p2; }
            dp=n*(t*p1-p0)/(t*t-1); const double dt=p1/dp; t-=dt; if (std::abs(dt)<1e-16) break;
        }
        x[i]=t; w[i]=2/((1-t*t)*dp*dp);
    }
}
// Composite Gauss-Legendre over [0,W/2]x[0,H/2] (the profile is even in y and z), panels
// graded towards the walls, where the corner terms r^2 log r live.
double quarter_integral(const duct::Duct& d,double G) {
    std::vector<double> x,w; gauss(10,x,w);
    auto edges=[](double L) {
        std::vector<double> e{0};
        for (int k=1;k<=40;++k) { const double s=double(k)/40; e.push_back(L*(1-(1-s)*(1-s))); }
        return e;
    };
    const auto ey=edges(d.width/2), ez=edges(d.height/2);
    double sum=0;
    for (std::size_t p=0;p+1<ey.size();++p) for (std::size_t q=0;q+1<ez.size();++q) {
        const double hy=(ey[p+1]-ey[p])/2, hz=(ez[q+1]-ez[q])/2, cy=(ey[p+1]+ey[p])/2, cz=(ez[q+1]+ez[q])/2;
        for (std::size_t i=0;i<x.size();++i) for (std::size_t j=0;j<x.size();++j)
            sum+=w[i]*w[j]*hy*hz*d.velocity(cy+hy*x[i],cz+hz*x[j],G);
    }
    return 4*sum;
}
} // namespace

int table() {
    const double mu=.1, U=.1;
    for (auto shape:{std::pair{2.,1.},std::pair{1.,2.},std::pair{1.,1.},std::pair{3.,.75}}) {
        duct::Duct d{shape.first,shape.second,mu}; const double G=d.pressure_gradient(U);
        std::printf("K %.17g %.17g %.17g %.17g\n",d.width,d.height,d.shape_factor(),G);
        for (double fy:{0.,.25,-.5,.75,.97,.999}) for (double fz:{0.,-.3,.6,.95}) {
            const double y=fy*d.width/2, z=fz*d.height/2;
            std::printf("u %.17g %.17g %.17g %.17g %.17g\n",d.width,d.height,y,z,d.velocity(y,z,G));
        }
    }
    return 0;
}

int main(int argc,char** argv) {
    if (argc==2 && std::string(argv[1])=="--table") return table();
    const double mu=.1, U=.1;
    for (auto shape:{std::pair{2.,1.},std::pair{1.,2.},std::pair{1.,1.},std::pair{3.,.75}}) {
        duct::Duct d{shape.first,shape.second,mu};
        const double G=d.pressure_gradient(U), A=d.a(), uc=d.centerline(G);
        std::printf("-- duct %gx%g: a=%g b=%g K=%.12f G=%.12g u_c/U=%.9f fRe=%.6f\n",d.width,d.height,A,d.b(),d.shape_factor(),G,uc/U,d.fanning_re());

        // Poisson residual, 4th-order central differences, relative to G.
        double poisson=0; const double dh=2e-3*A;
        for (double fy:{0.,.3,-.55,.8,.93}) for (double fz:{0.,-.25,.5,.77,.9}) {
            const double y=fy*d.width/2, z=fz*d.height/2;
            auto u=[&](double yy,double zz) { return d.velocity(yy,zz,G); };
            auto d2=[&](double m2,double m1,double c,double p1,double p2) { return (-m2+16*m1-30*c+16*p1-p2)/(12*dh*dh); };
            const double lap=d2(u(y-2*dh,z),u(y-dh,z),u(y,z),u(y+dh,z),u(y+2*dh,z))+d2(u(y,z-2*dh),u(y,z-dh),u(y,z),u(y,z+dh),u(y,z+2*dh));
            poisson=std::max(poisson,std::abs(mu*lap+G)/G);
        }
        below("  |mu lap(u) + G|/G at 25 interior points",poisson,1e-6);

        // No slip: the series cancels the plane profile as |t|->b, and cos(lambda a)=0 at |s|=a.
        double wall=0; const double eps=1e-9*A;
        for (double f:{-.9,-.4,0.,.35,.8}) {
            wall=std::max({wall,std::abs(d.velocity(d.width/2-eps,f*d.height/2,G)),std::abs(d.velocity(-d.width/2+eps,f*d.height/2,G)),
                           std::abs(d.velocity(f*d.width/2,d.height/2-eps,G)),std::abs(d.velocity(f*d.width/2,-d.height/2+eps,G))});
        }
        below("  max |u| at 1e-9 a from the walls, / u_c",wall/uc,1e-7);

        double symmetry=0;
        for (double y:{.1,.37,.6}) for (double z:{.05,.2,.33}) {
            const double v=d.velocity(y*d.width,z*d.height,G);
            symmetry=std::max({symmetry,std::abs(v-d.velocity(-y*d.width,z*d.height,G)),std::abs(v-d.velocity(y*d.width,-z*d.height,G))});
        }
        below("  symmetry y->-y, z->-z, / u_c",symmetry/uc,1e-15);

        // The two expansions (plane profile across a, or across b) agree wherever both converge.
        double expansions=0;
        for (double fs:{0.,.3,.6,.85}) for (double ft:{0.,.4,.7,.9}) {
            const double s=fs*A, t=ft*d.b();
            expansions=std::max(expansions,std::abs(d.expansion(s,t,A,d.b(),G)-d.expansion(t,s,d.b(),A,G)));
        }
        below("  |expansion across a - expansion across b| / u_c",expansions/uc,1e-13);
        duct::Duct r{d.height,d.width,mu}; double rotation=0;
        for (double y:{0.,.21,.43}) for (double z:{0.,.11,.34}) rotation=std::max(rotation,std::abs(d.velocity(y*d.width,z*d.height,G)-r.velocity(z*d.height,y*d.width,G)));
        below("  rotation (W,H)->(H,W), / u_c",rotation/uc,1e-15);

        // Flow rate: quadrature of u against the closed form a^2 G K/(3 mu) * A.
        const double Q=quarter_integral(d,G);
        below("  |quadrature Q - closed-form Q| / Q",std::abs(Q-d.flow_rate(G))/d.flow_rate(G),1e-9);
        below("  |U_mean(G(U)) - U| / U",std::abs(d.mean_velocity(G)-U)/U,1e-15);
    }
    // Literature values: Fanning f Re_Dh (Shah & London, Table 45) and u_max/u_mean.
    auto fre=[&](double w,double h) { return duct::Duct{w,h,mu}.fanning_re(); };
    below("|fRe(1:1) - 14.22708|",std::abs(fre(1,1)-14.22708),5e-5);
    below("|fRe(1:2) - 15.54806|",std::abs(fre(2,1)-15.54806),5e-5);
    below("|fRe(1:4) - 18.23278|",std::abs(fre(4,1)-18.23278),5e-5);
    duct::Duct square{1,1,mu}, plates{1000,1,mu};
    below("|u_max/u_mean(1:1) - 2.09624|",std::abs(square.centerline(square.pressure_gradient(U))/U-2.09624),5e-5);
    // Parallel plates (b/a=1000): u_c/U -> 1.5 and fRe -> 24, up to the O(a/b) side-wall share.
    below("|u_c/U(1:1000) - 1.5/K| (plane profile at the center)",std::abs(plates.centerline(plates.pressure_gradient(U))/U-1.5/plates.shape_factor()),1e-12);
    below("|K(1:1000) - (1 - 0.63 a/b)|",std::abs(plates.shape_factor()-(1-192/std::pow(duct::pi,5)*1.0045237627/1000)),1e-9);
    below("|fRe(1:1000) - 24| (O(a/b) from the side walls)",std::abs(fre(1000,1)-24),.05);
    // Entrance length estimate: creeping limit 0.619 D, and monotone in Re.
    below("|L_Durst(Re->0)/D - 0.619|",std::abs(duct::durst_entrance_length(0,1)-.619),1e-12);
    expect(duct::durst_entrance_length(100,1)>duct::durst_entrance_length(10,1),"L_Durst increases with Re",duct::durst_entrance_length(100,1),0);

    std::printf(failures?"FAIL: %d analytic checks\n":"PASS: analytic duct solution\n",failures);
    return failures?1:0;
}
