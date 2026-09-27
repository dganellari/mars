#pragma once
#include "channel.hpp"
#include "mars_segregated_simple_options.hpp"
namespace dsimple_gate {
inline std::array<double,3> rotate_vector(double x,double y,double z) {
    const double a=std::sqrt(.5);
    return {a*(x-y),.5*(x+y)-a*z,.5*(x+y)+a*z};
}
inline std::array<double,3> unrotate_point(double x,double y,double z) {
    x-=2; y+=1; z-=3;
    const double a=std::sqrt(.5);
    return {a*x+.5*(y+z),-a*x+.5*(y+z),a*(z-y)};
}
inline void rotate_channel(SimpleInput& mesh) {
    for (size_t i=0;i<mesh.x.size();++i) {
        const auto v=rotate_vector(mesh.x[i],mesh.y[i],mesh.z[i]);
        mesh.x[i]=v[0]+2; mesh.y[i]=v[1]-1; mesh.z[i]=v[2]+3;
    }
}
inline mars::segregated::SimpleControls configured_controls() {
    // Use the production parser, so the gate exercises the same controls as the executable.
    std::vector<std::string> words={"simple","--mesh=public.exo","--output-prefix=unused",
        "--rho=2","--mu=.4","--inlet-velocity=.2","--outlet-pressure=.03","--pseudo-dt=.005",
        "--relax-u=.4","--relax-p=.2","--relax-mass=.8","--outlet-beta=.1","--reference-length=1"};
    std::vector<char*> argv; for (auto& word:words) argv.push_back(word.data());
    return mars::segregated::simple_options(int(argv.size()),argv.data()).controls;
}
} // namespace dsimple_gate
