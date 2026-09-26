#pragma once
// Generated Tet4 channel (inlet x=0, outlet x=Lx, no-slip walls) and its partition into
// rank-local SIMPLE inputs. Test-only: the global mesh is replicated on every rank to build
// explicit ownership, which production gets from ingestion and element-halo completion.
#include "mars_segregated_simple_distributed.hpp"
#include "mars_segregated_simple_input.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace dsimple_gate {
using mars::segregated::SimpleFace;
using mars::segregated::SimpleInput;

inline std::uint64_t mix(std::uint64_t x) {
    x+=0x9e3779b97f4a7c15ull; x=(x^(x>>30))*0xbf58476d1ce4e5b9ull; x=(x^(x>>27))*0x94d049bb133111ebull;
    return x^(x>>31);
}
template<class T> void shuffle(std::vector<T>& v,std::uint64_t seed) {
    for (std::size_t i=v.size();i>1;--i) std::swap(v[i-1],v[mix(seed^(i*0x100000001b3ull))%i]);
}

// nx x ny x nz hexes on [0,nx/ny] x [0,1] x [0,1], six positively oriented tets per hex.
// 16x4x4 gives the public channel's sizes: 425 nodes, 1536 tets, 32/32/512 faces.
inline SimpleInput channel(int nx,int ny,int nz) {
    SimpleInput f;
    const double h=1.0/ny;
    auto id=[&](int i,int j,int k) { return i+(nx+1)*(j+(ny+1)*k); };
    for (int k=0;k<=nz;++k) for (int j=0;j<=ny;++j) for (int i=0;i<=nx;++i) {
        f.x.push_back(i*h); f.y.push_back(j*h); f.z.push_back(double(k)/nz);
    }
    const int axes[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (int k=0;k<nz;++k) for (int j=0;j<ny;++j) for (int i=0;i<nx;++i)
        for (const auto& a:axes) {
            int c[3]={i,j,k}; std::array<int,4> t;
            t[0]=id(c[0],c[1],c[2]);
            for (int s=0;s<2;++s) { ++c[a[s]]; t[s+1]=id(c[0],c[1],c[2]); }
            t[3]=id(i+1,j+1,k+1);
            double m[3][3];
            for (int r=0;r<3;++r) { m[r][0]=f.x[t[r+1]]-f.x[t[0]]; m[r][1]=f.y[t[r+1]]-f.y[t[0]]; m[r][2]=f.z[t[r+1]]-f.z[t[0]]; }
            const double det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                            +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
            if (det<0) std::swap(t[2],t[3]);
            for (int n=0;n<4;++n) f.nodes[n].push_back(t[n]);
        }
    const int e=int(f.nodes[0].size());
    std::map<std::array<int,3>,int> count;
    auto key=[&](int el,int ordinal) {
        std::array<int,3> k3; for (int j=0;j<3;++j) k3[j]=f.nodes[mars::segregated::tet_face_node(ordinal,j)][el];
        std::sort(k3.begin(),k3.end()); return k3;
    };
    for (int el=0;el<e;++el) for (int o=0;o<4;++o) ++count[key(el,o)];
    const double lx=nx*h;
    for (int el=0;el<e;++el) for (int o=0;o<4;++o) {
        if (count[key(el,o)]!=1) continue;
        const auto k3=key(el,o);
        auto all=[&](const std::vector<double>& c,double v) { for (int n:k3) if (std::abs(c[n]-v)>1e-12) return false; return true; };
        f.faces.push_back({el,o,all(f.x,0)?0:all(f.x,lx)?1:2});
    }
    return f;
}

// One rank's SIMPLE input, ownership and the global identity of every local entity.
struct Part {
    SimpleInput input;
    mars::segregated::runtime::SimpleOwnership<long long> ownership;
    std::vector<int> node_global, element_global, face_global;
    std::vector<char> node_owned;
};
struct Partition {
    int ranks;
    std::vector<int> element_owner, node_owner;
    std::vector<std::vector<int>> owned_order;   // solver order per rank
    std::vector<long long> first, solver_node;
    std::vector<int> solver_to_global;
};
// Uneven x slabs (plus a y split at 0.7 for 4 ranks: the outlet and its closing upper band are shared), jagged by moving
// every seventh element to a neighbouring rank; a node belongs to its lowest-rank element owner.
inline Partition partition(const SimpleInput& f,int ranks) {
    Partition p; p.ranks=ranks;
    const int e=int(f.nodes[0].size()), n=int(f.x.size());
    double lx=0; for (double v:f.x) lx=std::max(lx,v);
    p.element_owner.resize(e);
    for (int el=0;el<e;++el) {
        double cx=0, cy=0; for (int k=0;k<4;++k) { cx+=f.x[f.nodes[k][el]]/4; cy+=f.y[f.nodes[k][el]]/4; }
        int r=0;
        if (ranks==2) r=cx>0.43*lx;
        if (ranks==4) r=2*(cx>0.38*lx)+(cy>0.7);
        if (ranks>1 && mix(std::uint64_t(el)*31+7)%7==0) r=(r+1)%ranks;
        p.element_owner[el]=r;
    }
    p.node_owner.assign(n,ranks);
    for (int el=0;el<e;++el) for (int k=0;k<4;++k) p.node_owner[f.nodes[k][el]]=std::min(p.node_owner[f.nodes[k][el]],p.element_owner[el]);
    p.owned_order.resize(ranks); p.first.assign(ranks+1,0); p.solver_node.resize(n); p.solver_to_global.resize(n);
    for (int g=0;g<n;++g) p.owned_order[p.node_owner[g]].push_back(g);
    for (int r=0;r<ranks;++r) {
        shuffle(p.owned_order[r],1000+r); p.first[r+1]=p.first[r]+(long long)p.owned_order[r].size();
        for (std::size_t k=0;k<p.owned_order[r].size();++k) {
            p.solver_node[p.owned_order[r][k]]=p.first[r]+(long long)k; p.solver_to_global[p.first[r]+k]=p.owned_order[r][k];
        }
    }
    return p;
}
inline std::vector<int> present_elements(const SimpleInput& f,const Partition& p,int r) {
    std::vector<int> out;
    for (int el=0;el<int(f.nodes[0].size());++el) {
        bool touch=p.element_owner[el]==r;
        for (int k=0;k<4;++k) touch=touch || p.node_owner[f.nodes[k][el]]==r;
        if (touch) out.push_back(el);
    }
    return out;
}
// Nodes held by rank q: every node of its present elements (without the dropped one on
// drop_rank) plus its owned nodes.
inline std::set<int> holds(const SimpleInput& f,const Partition& p,int q,int drop_element,int drop_rank) {
    std::set<int> s(p.owned_order[q].begin(),p.owned_order[q].end());
    for (int el:present_elements(f,p,q)) if (!(q==drop_rank && el==drop_element)) for (int k=0;k<4;++k) s.insert(f.nodes[k][el]);
    return s;
}
// drop_element on drop_rank: that rank loses the element (fault injection: incomplete star).
inline Part extract(const SimpleInput& f,const Partition& p,int r,int drop_element=-1,int drop_rank=-1) {
    Part part; auto els=present_elements(f,p,r);
    if (r==drop_rank) els.erase(std::remove(els.begin(),els.end(),drop_element),els.end());
    shuffle(els,3000+r); part.element_global=els;
    const std::set<int> nodes=holds(f,p,r,drop_element,drop_rank);
    part.node_global.assign(nodes.begin(),nodes.end()); shuffle(part.node_global,2000+r);
    std::vector<int> local(f.x.size(),-1);
    for (int v=0;v<int(part.node_global.size());++v) local[part.node_global[v]]=v;
    std::vector<int> local_element(f.nodes[0].size(),-1);
    for (int v=0;v<int(els.size());++v) local_element[els[v]]=v;
    auto& in=part.input;
    for (int g:part.node_global) { in.x.push_back(f.x[g]); in.y.push_back(f.y[g]); in.z.push_back(f.z[g]); }
    for (int el:els) for (int k=0;k<4;++k) in.nodes[k].push_back(local[f.nodes[k][el]]);
    std::vector<int> faces;
    for (int i=0;i<int(f.faces.size());++i) if (local_element[f.faces[i].element]>=0) faces.push_back(i);
    shuffle(faces,4000+r); part.face_global=faces;
    auto& o=part.ownership;
    for (int i:faces) {
        in.faces.push_back({local_element[f.faces[i].element],f.faces[i].ordinal,f.faces[i].kind});
        if (p.element_owner[f.faces[i].element]==r) o.owned_faces.push_back(int(in.faces.size())-1);
    }
    for (int v=0;v<int(els.size());++v) if (p.element_owner[els[v]]==r) o.owned_elements.push_back(v);
    for (int g:p.owned_order[r]) o.owned_nodes.push_back(local[g]);
    for (int g:part.node_global) { o.solver_node.push_back(p.solver_node[g]); part.node_owned.push_back(p.node_owner[g]==r); }
    // Halo lists: receive every ghost from its owner; send owned nodes that a peer holds. Both
    // sides order a peer's list by global node, so slots line up.
    for (int q=0;q<p.ranks;++q) {
        if (q==r) continue;
        const auto theirs=holds(f,p,q,drop_element,drop_rank);
        std::vector<int> send, recv;
        for (int g:nodes) {
            if (p.node_owner[g]==q) recv.push_back(local[g]);
            if (p.node_owner[g]==r && theirs.count(g)) send.push_back(local[g]);
        }
        if (send.empty() && recv.empty()) continue;
        o.peers.push_back(q);
        if (o.send_offsets.empty()) { o.send_offsets.push_back(0); o.recv_offsets.push_back(0); }
        o.send_nodes.insert(o.send_nodes.end(),send.begin(),send.end()); o.send_offsets.push_back(int(o.send_nodes.size()));
        o.recv_nodes.insert(o.recv_nodes.end(),recv.begin(),recv.end()); o.recv_offsets.push_back(int(o.recv_nodes.size()));
    }
    if (o.send_offsets.empty()) { o.send_offsets.push_back(0); o.recv_offsets.push_back(0); }
    return part;
}
} // namespace dsimple_gate
