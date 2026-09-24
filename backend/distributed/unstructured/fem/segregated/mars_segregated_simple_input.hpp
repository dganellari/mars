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
namespace mars::segregated {
inline void simple_input_require(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
template<class T> void simple_input_read(std::istream& s,T& x) {
    simple_input_require(bool(s>>x) && std::isfinite(double(x)),"invalid/truncated public SIMPLE mesh");
}
struct SimpleInput {
    std::vector<double> x,y,z;
    std::array<std::vector<int>,4> nodes;
    std::vector<SimpleFace> faces;
};
inline SimpleInput read_simple_topology(std::istream& s,int n,int e,int b) {
    simple_input_require(n==425 && e==1536 && b==576,"expected pinned public channel sizes");
    SimpleInput f; f.x.resize(n); f.y.resize(n); f.z.resize(n);
    for (int i=0;i<n;++i) { simple_input_read(s,f.x[i]); simple_input_read(s,f.y[i]); simple_input_read(s,f.z[i]); }
    for (auto& a:f.nodes) a.resize(e);
    for (int i=0;i<e;++i) {
        std::set<int> unique;
        for (auto& a:f.nodes) { simple_input_read(s,a[i]); simple_input_require(a[i]>=0 && a[i]<n,"bad cell node"); unique.insert(a[i]); }
        simple_input_require(unique.size()==4,"degenerate connectivity");
    }
    std::map<std::array<int,3>,int> incidence;
    auto face_key=[&](int element,int ordinal) {
        std::array<int,3> key;
        for (int j=0;j<3;++j) key[j]=f.nodes[mars::segregated::tet_face_node(ordinal,j)][element];
        std::sort(key.begin(),key.end()); return key;
    };
    for (int i=0;i<e;++i) for (int j=0;j<4;++j) ++incidence[face_key(i,j)];
    int exterior=0;
    for (const auto& entry:incidence) { simple_input_require(entry.second<=2,"nonmanifold cell topology"); exterior+=entry.second==1; }
    simple_input_require(exterior==b,"incomplete exterior topology");
    std::set<std::pair<int,int>> unique;
    f.faces.resize(b); int types[3]{};
    for (auto& face:f.faces) {
        simple_input_read(s,face.element); simple_input_read(s,face.ordinal); simple_input_read(s,face.kind);
        simple_input_require(face.element>=0 && face.element<e && face.ordinal>=0 && face.ordinal<4 && face.kind>=0 && face.kind<3,"bad boundary tag");
        simple_input_require(incidence[face_key(face.element,face.ordinal)]==1,"boundary tag on interior face");
        simple_input_require(unique.emplace(face.element,face.ordinal).second,"duplicate boundary"); ++types[face.kind];
    }
    simple_input_require(types[0]==32 && types[1]==32 && types[2]==512,"wrong public boundary coverage");
    return f;
}
inline SimpleInput load_simple_input(const char* path) {
    std::ifstream s(path); std::string magic; int n,e,b;
    simple_input_require(bool(s>>magic>>n>>e>>b) && magic=="MARS_PUBLIC_SIMPLE_MESH_V1",
                         "expected public channel mesh-only input");
    auto f=read_simple_topology(s,n,e,b);
    simple_input_require(!(s>>magic),"trailing public mesh data"); return f;
}
}
