#pragma once
#include "mars_segregated_simple.hpp"
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace mars::segregated {
inline std::string simple_boundary_name(std::string name) {
    // Ioss database names use lowercase ASCII and replace spaces with underscores.
    for (char& c:name) {
        if (c>='A' && c<='Z') c=char(c-'A'+'a');
        else if (c==' ') c='_';
    }
    return name;
}
struct SimpleBoundaryNames {
    std::string inlet="inlet",outlet="outlet";
    std::vector<std::string> walls{"walls"};
    bool valid() const {
        std::set<std::string> names;
        for (const auto& name:{inlet,outlet})
            if (name.empty() || !names.insert(simple_boundary_name(name)).second) return false;
        if (walls.empty()) return false;
        for (const auto& name:walls)
            if (name.empty() || !names.insert(simple_boundary_name(name)).second) return false;
        return true;
    }
    int kind(const std::string& name) const {
        const auto key=simple_boundary_name(name);
        if (key==simple_boundary_name(inlet)) return 0;
        if (key==simple_boundary_name(outlet)) return 1;
        for (const auto& wall:walls) if (key==simple_boundary_name(wall)) return 2;
        return -1;
    }
};
struct SimpleOptions {
    std::string mesh,output,field_output="gathered";
    SimpleControls controls;
    SimpleBoundaryNames boundaries;
    int iterations=2000,report=10,profile_warmup=10;
    double residual=1e-6,mass=1e-6,change=1e-6;
    double pressure_rtol=1e-12,pressure_atol=0;
    bool pressure_tolerances=false;
    bool setup_only=false,help=false,profile=false,linear_cache=true,halo_overlap=true;
};
inline SimpleOptions simple_options(int argc,char** argv) {
    SimpleOptions o;
    std::set<std::string> seen;
    for (int i=1;i<argc;++i) {
        std::string key=argv[i],value;
        if (key=="--help" || key=="-h") { o.help=true; continue; }
        const auto equals=key.find('=');
        if (equals!=std::string::npos) { value=key.substr(equals+1); key.resize(equals); }
        else {
            if (++i==argc) throw std::runtime_error("missing value for "+key);
            value=argv[i];
        }
        if (value.empty() || !seen.insert(key).second) throw std::runtime_error("empty or repeated option: "+key);
        if (key=="--mesh") o.mesh=value;
        else if (key=="--output-prefix") o.output=value;
        else if (key=="--field-output") {
            if (value!="gathered" && value!="distributed" && value!="none")
                throw std::runtime_error("--field-output expects gathered, distributed or none");
            o.field_output=value;
        }
        else if (key=="--mesh-format") {
            if (value!="exodus") throw std::runtime_error("SIMPLE uses native Exodus; prepared input is reserved for the reference gates");
        }
        else if (key=="--inlet-ss") o.boundaries.inlet=value;
        else if (key=="--outlet-ss") o.boundaries.outlet=value;
        else if (key=="--advection") {
            if (value!="upwind" && value!="high-resolution") throw std::runtime_error("--advection expects upwind or high-resolution");
            o.controls.high_resolution=value=="high-resolution";
        }
        else if (key=="--velocity-interpolation") {
            if (value!="trilinear" && value!="linear-linear")
                throw std::runtime_error("--velocity-interpolation expects trilinear or linear-linear");
            o.controls.velocity_shifted=value=="linear-linear";
        }
        else if (key=="--wall-ss") {
            o.boundaries.walls.clear();
            for (size_t start=0;;) {
                const auto end=value.find(',',start);
                o.boundaries.walls.push_back(value.substr(start,end==std::string::npos?end:end-start));
                if (end==std::string::npos) break;
                start=end+1;
            }
        }
        else {
            size_t end=0; const double number=std::stod(value,&end);
            if (end!=value.size() || !std::isfinite(number)) throw std::runtime_error("expected finite number for "+key);
            if (key=="--iterations" || key=="--report-every") {
                if (number<1 || number>std::numeric_limits<int>::max() || number!=std::floor(number))
                    throw std::runtime_error("expected positive integer for "+key);
                (key=="--iterations"?o.iterations:o.report)=int(number);
            }
            else if (key=="--profile-warmup") {
                if (number<0 || number>std::numeric_limits<int>::max() || number!=std::floor(number))
                    throw std::runtime_error("--profile-warmup expects a nonnegative integer");
                o.profile_warmup=int(number);
            }
            else if (key=="--setup-only" || key=="--profile" || key=="--linear-cache" || key=="--halo-overlap") {
                if (number!=0 && number!=1) throw std::runtime_error(key+" expects 0 or 1");
                if (key=="--setup-only") o.setup_only=number!=0;
                else if (key=="--profile") o.profile=number!=0;
                else if (key=="--linear-cache") o.linear_cache=number!=0;
                else o.halo_overlap=number!=0;
            }
            else if (key=="--residual-tol") o.residual=number;
            else if (key=="--mass-tol") o.mass=number;
            else if (key=="--change-tol") o.change=number;
            else if (key=="--pressure-linear-rtol") o.pressure_rtol=number;
            else if (key=="--pressure-linear-atol") o.pressure_atol=number;
            else if (key=="--rho") o.controls.density=number;
            else if (key=="--mu") o.controls.viscosity=number;
            else if (key=="--inlet-velocity") o.controls.inlet_speed=number;
            else if (key=="--outlet-pressure") o.controls.pressure_reference=number;
            else if (key=="--pseudo-dt") o.controls.pseudo_dt=number;
            else if (key=="--relax-u") o.controls.alpha_u=number;
            else if (key=="--relax-p") o.controls.alpha_p=number;
            else if (key=="--relax-mass") o.controls.alpha_mass=number;
            else if (key=="--outlet-beta") o.controls.beta=number;
            else if (key=="--reference-length") o.controls.reference_length=number;
            else throw std::runtime_error("unknown option: "+key);
        }
    }
    o.pressure_tolerances=seen.count("--pressure-linear-rtol")!=0;
    if (o.pressure_tolerances!=(seen.count("--pressure-linear-atol")!=0)
        || !(o.pressure_rtol>0 && o.pressure_rtol<1 && o.pressure_atol>=0))
        throw std::runtime_error("pressure linear tolerances require both 0<rtol<1 and atol>=0");
    if (!valid_simple_controls(o.controls) || !o.boundaries.valid()
        || !(o.residual>0 && o.mass>0 && o.change>0)) throw std::runtime_error("invalid SIMPLE controls, tolerances or boundary names");
    if (!o.help && (o.mesh.empty() || o.output.empty())) throw std::runtime_error("--mesh and --output-prefix are required");
    return o;
}
inline const char* simple_help() {
    return "Native distributed Tet4 SIMPLE (CUDA/Hypre, steady, laminar)\n"
           "  --mesh FILE --output-prefix PREFIX [--mesh-format exodus]\n"
           "  --inlet-ss inlet --outlet-ss outlet --wall-ss walls[,other_wall]\n"
           "  --advection upwind      or high-resolution (OpenAccel limiter, cap 1)\n"
           "  --velocity-interpolation trilinear   or linear-linear (shifted field weights)\n"
           "  --rho 1                 density [kg/m^3]\n"
           "  --mu 0.1                dynamic viscosity [Pa s]; nu=mu/rho\n"
           "  --inlet-velocity 0.1     positive inward-normal speed [m/s]\n"
           "  --outlet-pressure 0     area-mean outlet pressure [Pa]\n"
           "  --pseudo-dt 0.01        steady pseudo-time [s], not physical time\n"
           "  --relax-u 0.3 --relax-p 0.3 --relax-mass 0.75 --outlet-beta 0.05\n"
           "  --reference-length 1    residual normalization length [m]\n"
           "  --iterations 2000 --report-every 10 --setup-only 0\n"
           "  --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6\n"
           "  --pressure-linear-rtol R --pressure-linear-atol A   optional pair; max(A,R*||b||)\n"
           "    Sets pressure Krylov and both true residual targets; momentum stays unchanged.\n"
           "  --linear-cache 1 --halo-overlap 1    set 0 for a performance control\n"
           "  --profile 0 --profile-warmup 10      optional phase timing (adds event fences)\n"
           "  --field-output gathered             distributed writes per-rank CSVs; none skips fields\n"
           "Both --option value and --option=value are accepted. Every exterior\n"
           "face must belong to exactly one selected inlet, outlet or wall set.\n";
}
} // namespace mars::segregated
