#pragma once
#include "mars_segregated_simple.hpp"
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
struct Fixture {
    std::vector<double> x,y,z;
    std::array<std::vector<int>,4> nodes;
    std::vector<mars::segregated::SimpleFace> faces;
    std::vector<Expected> expected;
};
inline Fixture load(const char* path) {
    std::ifstream s(path); std::string magic; int n,e,b,iterations;
    require(bool(s>>magic>>n>>e>>b>>iterations) && magic=="MARS_PUBLIC_SIMPLE_V1"
            && n==425 && e==1536 && b==576 && iterations==2,"expected pinned public SIMPLE fixture");
    Fixture f; f.x.resize(n); f.y.resize(n); f.z.resize(n);
    for (int i=0;i<n;++i) { read(s,f.x[i]); read(s,f.y[i]); read(s,f.z[i]); }
    for (auto& a:f.nodes) a.resize(e);
    for (int i=0;i<e;++i) {
        std::set<int> unique;
        for (auto& a:f.nodes) { read(s,a[i]); require(a[i]>=0 && a[i]<n,"bad cell node"); unique.insert(a[i]); }
        require(unique.size()==4,"degenerate connectivity");
    }
    std::map<std::array<int,3>,int> incidence;
    auto face_key=[&](int element,int ordinal) {
        std::array<int,3> key;
        for (int j=0;j<3;++j) key[j]=f.nodes[mars::segregated::tet_face_node(ordinal,j)][element];
        std::sort(key.begin(),key.end()); return key;
    };
    for (int i=0;i<e;++i) for (int j=0;j<4;++j) ++incidence[face_key(i,j)];
    int exterior=0;
    for (const auto& entry:incidence) { require(entry.second<=2,"nonmanifold cell topology"); exterior+=entry.second==1; }
    require(exterior==b,"incomplete exterior topology");
    std::set<std::pair<int,int>> unique;
    f.faces.resize(b); int types[3]{};
    for (auto& face:f.faces) {
        read(s,face.element); read(s,face.ordinal); read(s,face.kind);
        require(face.element>=0 && face.element<e && face.ordinal>=0 && face.ordinal<4 && face.kind>=0 && face.kind<3,"bad boundary tag");
        require(incidence[face_key(face.element,face.ordinal)]==1,"boundary tag on interior face");
        require(unique.emplace(face.element,face.ordinal).second,"duplicate boundary"); ++types[face.kind];
    }
    require(types[0]==32 && types[1]==32 && types[2]==512,"wrong public boundary coverage");
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
