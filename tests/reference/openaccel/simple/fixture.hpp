#pragma once
#include "mars_segregated_simple_input.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
namespace simple_check {
inline void require(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
template<class T> void read(std::istream& s,T& x) { require(bool(s>>x) && std::isfinite(double(x)),"invalid/truncated SIMPLE fixture"); }
struct Expected {
    std::vector<double> momentum,increment,pressure,velocity,influence,trace,interior_flux,boundary_flux,mass_divergence;
};
struct Fixture : mars::segregated::SimpleInput {
    std::vector<Expected> expected;
};
inline Fixture load(const char* path) {
    std::ifstream s(path); std::string magic; int n,e,b,iterations;
    require(bool(s>>magic>>n>>e>>b>>iterations) && magic=="MARS_PUBLIC_SIMPLE_V1"
            && n==425 && e==1536 && b==576 && iterations==2,"expected pinned public SIMPLE fixture");
    Fixture f;
    static_cast<mars::segregated::SimpleInput&>(f)=mars::segregated::read_simple_topology(s,n,e,b);
    for (int i=0;i<iterations;++i) {
        Expected x;
        auto values=[&](std::vector<double>& v,int count) { v.resize(count); for (auto& a:v) read(s,a); };
        values(x.momentum,3*n); values(x.increment,n); values(x.pressure,n); values(x.velocity,3*n); values(x.influence,3*n);
        values(x.trace,3*b); values(x.interior_flux,6*e); values(x.boundary_flux,3*b); values(x.mass_divergence,n);
        f.expected.push_back(std::move(x));
    }
    require(!(s>>magic),"trailing SIMPLE fixture data"); return f;
}
inline std::size_t checks=0;
inline void compare(const std::vector<double>& a,const std::vector<double>& b,const std::string& name) {
    require(a.size()==b.size(),"comparison size"); double scale=0,worst=0; bool finite=true;
    for (double v:b) scale=std::max(scale,std::abs(v));
    for (std::size_t i=0;i<a.size();++i) { ++checks; finite=finite && std::isfinite(a[i]); worst=std::max(worst,std::abs(a[i]-b[i])); }
    std::cout<<name<<" scalars="<<a.size()<<" worst_scaled="<<worst/std::max(1.,scale)<<'\n';
    // The captured pressure solve used rtol=1e-8; direct/tighter solves need not be bit-identical.
    require(finite && worst<=2e-11+1e-7*scale,"SIMPLE state mismatch");
}
} // namespace simple_check
