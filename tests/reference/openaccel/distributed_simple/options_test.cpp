#include "mars_segregated_simple_options.hpp"
#include "mars_segregated_simple_metrics.hpp"
#include <iostream>
#include <sstream>
using namespace mars::segregated;
int checks=0;
void check(bool pass) { ++checks; if (!pass) throw std::runtime_error("SIMPLE controls check "+std::to_string(checks)); }
SimpleOptions parse(const std::string& args) {
    std::istringstream in("simple --mesh public.exo --output-prefix result "+args);
    std::vector<std::string> words; std::string word;
    while (in>>word) words.push_back(word);
    std::vector<char*> argv; for (auto& w:words) argv.push_back(w.data());
    return simple_options(int(argv.size()),argv.data());
}
int main() {
    const auto defaults=parse("");
    check(defaults.controls.density==1 && defaults.controls.viscosity==.1 && defaults.controls.inlet_speed==.1);
    const auto options=parse("--rho=2 --mu .4 --inlet-velocity=.2 --outlet-pressure=-3 --pseudo-dt .005 "
        "--reference-length=2 --relax-u .4 --relax-p .2 --relax-mass .8 --outlet-beta .1 "
        "--inlet-ss feed --outlet-ss exit --wall-ss casing,cover --iterations 20 --report-every=2");
    const auto c=options.controls;
    check(c.density==2 && c.viscosity==.4 && c.inlet_speed==.2 && c.pressure_reference==-3);
    check(c.pseudo_dt==.005 && c.reference_length==2 && c.alpha_u==.4 && c.alpha_p==.2 && c.alpha_mass==.8 && c.beta==.1);
    check(options.boundaries.kind("feed")==0 && options.boundaries.kind("exit")==1);
    check(options.boundaries.kind("casing")==2 && options.boundaries.kind("cover")==2 && options.boundaries.kind("other")==-1);
    check(options.iterations==20 && options.report==2);
    for (const char* args:{"--rho nan","--rho inf","--rho -1","--mu 0","--inlet-velocity 0","--outlet-pressure inf",
        "--reference-length 0","--pseudo-dt -1","--relax-u 1.1","--relax-p 0","--relax-mass nan","--outlet-beta 0",
        "--residual-tol nan","--mass-tol 0","--change-tol -1","--iterations 2.5","--iterations 2147483648",
        "--report-every 0","--setup-only 2","--rho 2x","--rho 1 --rho 2","--rho=","--mesh-format prepared",
        "--inlet-ss outlet","--wall-ss walls,","--wall-ss walls,walls","--wall-ss walls,inlet","--unknown 1"}) {
        bool rejected=false; try { parse(args); } catch (const std::exception&) { rejected=true; }
        check(rejected);
    }
    // An oblique tetrahedron makes axis-aligned inlet implementations fail.
    double xyz[12]={0,0,0, 2,1,0, 0,2,1, 1,0,3};
    TetGeometry<double> g; check(tet_geometry(xyz,g));
    const int nodes[4]={0,1,2,3}; double velocity[12]{},pressure[4]{},trace[3]{},flux[3]{};
    for (int f=0;f<4;++f) {
        const auto input=simple_boundary(true,{0,f,0},nodes,g,velocity,pressure,trace,flux,c,false);
        double area[3]; tet_boundary_area(g,f,area);
        const double norm=std::sqrt(area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
        for (int k=0;k<3;++k) {
            double dot=0,speed2=0;
            for (int j=0;j<3;++j) {
                const double u=input.values.boundary_velocity[3*k+j];
                dot+=u*area[j]; speed2+=u*u;
                check(std::abs(u+c.inlet_speed*area[j]/norm)<1e-14);
            }
            check(std::abs(dot+c.inlet_speed*norm)<1e-14 && std::abs(speed2-c.inlet_speed*c.inlet_speed)<1e-14);
            check(input.values.density[k]==2 && input.values.viscosity[k]==.4);
        }
        const auto wall=simple_boundary(true,{0,f,2},nodes,g,velocity,pressure,trace,flux,c,false);
        for (double v:wall.values.boundary_velocity) check(v==0);
    }
    SimpleSums sums; sums.volume=2; sums.inlet_area=1; sums.momentum2=3; sums.continuity2=4;
    sums.inlet=-.4; sums.outlet=.4; sums.velocity_change2=.1; sums.pressure_change2=.2;
    const auto scaled=simple_metrics(sums,c),unit=simple_metrics(sums,c,1);
    check(scaled.finite && scaled.momentum==2*unit.momentum && scaled.continuity==2*unit.continuity);
    check(scaled.flux==unit.flux && scaled.velocity_change==unit.velocity_change && scaled.pressure_change==unit.pressure_change);
    std::cout<<"PASS: "<<checks<<" SIMPLE controls, units, oblique inlet normals and metric checks\n";
}
