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
    check(defaults.linear_cache && defaults.halo_overlap && !defaults.profile && defaults.field_output=="gathered");
    check(!defaults.pressure_tolerances && defaults.pressure_rtol==1e-12 && defaults.pressure_atol==0);
    check(!defaults.pressure_refinement);
    check(defaults.pressure_failure_capture.empty());
    check(defaults.pressure_solver_profile.empty());
    check(!defaults.controls.pressure_expansion);
    check(parse("--pressure-expansion 1 --pressure-linear-rtol 1e-10 --pressure-linear-atol 0").controls.pressure_expansion);
    check(parse("--pressure-expansion 1 --pressure-linear-rtol 1e-10 --pressure-linear-atol 0 --pressure-solver-profile settings").controls.pressure_expansion);
    check(parse("--pressure-solver-profile settings --pressure-linear-rtol 1e-4 --pressure-linear-atol 0 "
                "--first-step-audit 1 --iterations 1 --snapshot-iterations 1 --field-output distributed").pressure_solver_profile=="settings");
    check(parse("--pressure-failure-capture private --pressure-linear-rtol 1e-8 --pressure-linear-atol 0").pressure_failure_capture=="private");
    check(parse("--pressure-refinement 1").pressure_refinement);
    check(!parse("--pressure-refinement=0").pressure_refinement);
    check(defaults.snapshot_iterations==0);
    check(!defaults.first_step_audit);
    check(parse("--first-step-audit 1 --iterations 1 --snapshot-iterations 1 --field-output distributed").first_step_audit);
    check(parse("--snapshot-iterations 20 --iterations 20 --field-output distributed").snapshot_iterations==20);
    const auto pressure_options=parse("--pressure-linear-rtol=1e-6 --pressure-linear-atol 0");
    check(pressure_options.pressure_tolerances && pressure_options.pressure_rtol==1e-6 && pressure_options.pressure_atol==0);
    check(pressure_options.residual==defaults.residual && pressure_options.mass==defaults.mass && pressure_options.change==defaults.change);
    const auto performance=parse("--profile 1 --profile-warmup 0 --linear-cache 0 --halo-overlap 0 --field-output distributed");
    check(performance.profile && performance.profile_warmup==0 && !performance.linear_cache && !performance.halo_overlap && performance.field_output=="distributed");
    check(!defaults.controls.high_resolution && !parse("--advection upwind").controls.high_resolution);
    check(parse("--advection=high-resolution").controls.high_resolution);
    check(!defaults.controls.velocity_shifted && !parse("--velocity-interpolation trilinear").controls.velocity_shifted);
    check(parse("--velocity-interpolation=linear-linear").controls.velocity_shifted);
    const auto options=parse("--rho=2 --mu .4 --inlet-velocity=.2 --outlet-pressure=-3 --pseudo-dt .005 "
        "--reference-length=2 --relax-u .4 --relax-p .2 --relax-mass .8 --outlet-beta .1 "
        "--inlet-ss feed --outlet-ss exit --wall-ss casing,cover --iterations 20 --report-every=2");
    const auto c=options.controls;
    check(c.density==2 && c.viscosity==.4 && c.inlet_speed==.2 && c.pressure_reference==-3);
    check(c.pseudo_dt==.005 && c.reference_length==2 && c.alpha_u==.4 && c.alpha_p==.2 && c.alpha_mass==.8 && c.beta==.1);
    check(options.boundaries.kind("feed")==0 && options.boundaries.kind("exit")==1);
    check(options.boundaries.kind("casing")==2 && options.boundaries.kind("cover")==2 && options.boundaries.kind("other")==-1);
    check(options.boundaries.kind("FEED")==0 && options.boundaries.kind("Exit")==1 && options.boundaries.kind("CASING")==2);
    const SimpleBoundaryNames spaced{"Feed Port","Exit Port",{"Outer Wall","Cover"}};
    check(spaced.valid() && spaced.kind("feed_port")==0 && spaced.kind("OUTER_WALL")==2);
    check(!SimpleBoundaryNames{"feed","FEED",{"walls"}}.valid());
    check(!SimpleBoundaryNames{"feed","exit",{"outer wall","OUTER_WALL"}}.valid());
    check(options.iterations==20 && options.report==2);
    for (const char* args:{"--rho nan","--rho inf","--rho -1","--mu 0","--inlet-velocity 0","--outlet-pressure inf",
        "--reference-length 0","--pseudo-dt -1","--relax-u 1.1","--relax-p 0","--relax-mass nan","--outlet-beta 0",
        "--residual-tol nan","--mass-tol 0","--change-tol -1","--iterations 2.5","--iterations 2147483648",
        "--pressure-failure-capture private",
        "--pressure-solver-profile settings",
        "--pressure-solver-profile settings --pressure-linear-rtol 1e-4 --pressure-linear-atol 0",
        "--pressure-solver-profile settings --pressure-linear-rtol 1e-4 --pressure-linear-atol 0 --pressure-refinement 1 "
        "--first-step-audit 1 --iterations 1 --snapshot-iterations 1 --field-output distributed",
        "--pressure-failure-capture private --pressure-linear-rtol 1e-8 --pressure-linear-atol 0 --pressure-refinement 1",
        "--pressure-failure-capture private --pressure-linear-rtol 1e-8 --pressure-linear-atol 0 --setup-only 1",
        "--pressure-failure-capture private --pressure-linear-rtol 1e-8 --pressure-linear-atol 0 --first-step-audit 1 --iterations 1 --snapshot-iterations 1 --field-output distributed",
        "--pressure-linear-rtol 1e-6","--pressure-linear-atol 0",
        "--pressure-linear-rtol 0 --pressure-linear-atol 1e-6",
        "--pressure-linear-rtol 1 --pressure-linear-atol 0",
        "--pressure-linear-rtol -1 --pressure-linear-atol 0",
        "--pressure-linear-rtol 1e-6 --pressure-linear-atol -1",
        "--pressure-linear-rtol nan --pressure-linear-atol 0",
        "--pressure-linear-rtol 1e-6 --pressure-linear-atol inf",
        "--report-every 0","--setup-only 2","--rho 2x","--rho 1 --rho 2","--rho=","--mesh-format prepared",
        "--advection central","--advection 1","--velocity-interpolation other","--velocity-interpolation 1",
        "--inlet-ss outlet","--inlet-ss OUTLET","--wall-ss walls,","--wall-ss walls,walls","--wall-ss walls,inlet","--wall-ss walls,INLET","--unknown 1",
        "--profile 2","--linear-cache -1","--halo-overlap 3","--profile-warmup -1","--profile-warmup 1.5","--field-output vtk",
        "--pressure-expansion 1","--pressure-expansion 2",
        "--pressure-expansion 1 --pressure-linear-rtol 1e-10 --pressure-linear-atol 0 --pressure-refinement 1",
        "--pressure-expansion 1 --pressure-linear-rtol 1e-10 --pressure-linear-atol 0 --snapshot-iterations 1",
        "--pressure-expansion 1 --pressure-linear-rtol 1e-10 --pressure-linear-atol 0 --pressure-failure-capture private",
        "--pressure-refinement 2","--pressure-refinement -1","--pressure-refinement 0.5",
        "--pressure-refinement 0 --pressure-refinement 1",
        "--snapshot-iterations -1","--snapshot-iterations 101","--snapshot-iterations 2.5",
        "--snapshot-iterations 20 --iterations 10","--snapshot-iterations 1 --field-output none",
        "--snapshot-iterations 1 --setup-only 1","--snapshot-iterations 1 --profile 1",
        "--first-step-audit 1","--first-step-audit 2",
        "--first-step-audit 1 --iterations 2 --snapshot-iterations 1 --field-output distributed",
        "--first-step-audit 1 --iterations 1 --snapshot-iterations 1 --field-output gathered"}) {
        bool rejected=false; try { parse(args); } catch (const std::exception&) { rejected=true; }
        check(rejected);
    }
    // An oblique tetrahedron makes axis-aligned inlet implementations fail.
    double xyz[12]={0,0,0, 2,1,0, 0,2,1, 1,0,3};
    TetGeometry<double> g; check(tet_geometry(xyz,g));
    const int nodes[4]={0,1,2,3}; double velocity[12]{},pressure[4]{},trace[3]{},flux[3]{};
    for (int f=0;f<4;++f) {
        double area[3]; tet_boundary_area(g,f,area);
        const double norm=std::sqrt(area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
        double inlet_velocity[12];
        for (int n=0;n<4;++n) for (int j=0;j<3;++j) inlet_velocity[3*n+j]=-c.inlet_speed*area[j]/norm;
        const auto input=simple_boundary(true,{0,f,0},nodes,g,velocity,pressure,trace,flux,c,false,inlet_velocity);
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
        const auto wall=simple_boundary(true,{0,f,2},nodes,g,velocity,pressure,trace,flux,c,false,nullptr);
        for (double v:wall.values.boundary_velocity) check(v==0);
    }
    SimpleSums sums; sums.volume=2; sums.inlet_area=1; sums.momentum2=3; sums.continuity2=4;
    sums.inlet=-.4; sums.outlet=.4; sums.velocity_change2=.1; sums.pressure_change2=.2;
    const auto scaled=simple_metrics(sums,c),unit=simple_metrics(sums,c,1);
    check(scaled.finite && scaled.momentum==2*unit.momentum && scaled.continuity==2*unit.continuity);
    check(scaled.flux==unit.flux && scaled.velocity_change==unit.velocity_change && scaled.pressure_change==unit.pressure_change);
    std::cout<<"PASS: "<<checks<<" SIMPLE controls, units, oblique inlet normals and metric checks\n";
}
